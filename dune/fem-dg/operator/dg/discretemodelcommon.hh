#ifndef DUNE_FEM_DG_DISCRETEMODELCOMMON_HH
#define DUNE_FEM_DG_DISCRETEMODELCOMMON_HH

#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem-dg/operator/fluxes/ldgflux.hh>
#include <dune/fem-dg/operator/limiter/limiter.hh>

#include <dune/fem-dg/operator/dg/passtraits.hh>
#include <dune/fem/gridpart/filter/threadfilter.hh>

namespace Dune {

  template < int idx >
  struct PrintTupleValues
  {
    template <class Tuple>
    static void apply( const Tuple*  )
    {
      std::cout << "Tuple< " << idx << " > = " << std::tuple_element< idx, Tuple >::type::value << std::endl;
    }
  };

  template< class Tuple, int begin, int end>
  struct SelectTupleElements
  {
    typedef typename CutOutTuple< Tuple, begin, end > :: type type;
  };

  template< class Tuple, int begin>
  struct SelectTupleElements< Tuple, begin, begin >
  {
    typedef std::tuple<> type;
  };

  template< class Traits,
            int passUId, int passGradId >
  class FemDGBaseDiscreteModel : public Fem::DGDiscreteModelDefaultWithInsideOutside< Traits, passUId, passGradId >
  {
    typedef Fem::DGDiscreteModelDefaultWithInsideOutside< Traits, passUId, passGradId > BaseType;
    typedef typename BaseType :: Selector BaseSelectorType ;

    typedef typename Traits :: ModelType :: ModelParameter ModelParameter;
    typedef typename Traits :: ExtraParameterTupleType ExtraParameterTupleType;

    static const int modParamSize = std::tuple_size< ModelParameter > :: value ;
    static const int extParamSize = std::tuple_size< ExtraParameterTupleType > :: value ;
    static const int selectedSize = ( extParamSize < modParamSize ) ? extParamSize :  modParamSize ;
    static_assert( selectedSize == 1, " selected size ");
    //static const int modId = typename std::tuple_element< 0, ModelParameter>::type::value;
    static const int newId = 0;//( modParamSize > 0 ) ? modId : 0;
    typedef typename SelectTupleElements< ModelParameter, 0, selectedSize > :: type ExtraParameterType;

  public:
    // overload selector type to add model parameters
    typedef typename Dune::Fem::Selector< newId, passUId, passGradId > ::Type Selector ;
      //Dune::Fem::ElementTuple< ModelParameterpassUId, passGradId, -1, -1, -1, -1, -1, -1, -1, ExtraParameterType > >::Type  Selector;
    //typedef typename Dune::Fem::Selector< passUId, passGradId > :: Type Selector;

    FemDGBaseDiscreteModel()
    {
#ifndef NDEBUG
      if( Dune::Fem::Parameter::verbose() )
      {
        ForLoop< PrintTupleValues, 0, std::tuple_size< Selector >::value-1 >::apply( ( Selector* ) 0 );
      }
#endif
    }
  };

  // AdvectionModel
  //---------------

  template< class OpTraits,
            int passUId, int passGradId,
            bool returnAdvectionPart>
  class AdvectionModel;


  // AdvectionTraits
  //----------------

  template <class Traits,
            int passUId, int passGradId, bool returnAdvectionPart>
  struct AdvectionTraits : public Traits
  {
    typedef AdvectionModel< Traits, passUId, passGradId, returnAdvectionPart >       DGDiscreteModelType;
  };


  // AdvectionModel
  //---------------

  /*  \class AdvectionModel
   *
   *  \tparam Model Mathematical model
   *  \tparam NumFlux Numerical flux
   *  \tparam polOrd Polynomial degree
   *  \tparam passUId The id of a pass whose value is used here
   *  \tparam passGradId The id of a pass whose value is used here
   *  \tparam returnAdvectionPart Switch on/off the advection
   */
  template< class OpTraits,
            int passUId, int passGradId,
            bool returnAdvectionPart>
  class AdvectionModel :
    public FemDGBaseDiscreteModel< AdvectionTraits< OpTraits, passUId, passGradId, returnAdvectionPart >,
                                   passUId, passGradId
                                 >
  {
  public:
    typedef AdvectionTraits< OpTraits, passUId, passGradId, returnAdvectionPart > Traits;

    typedef typename Traits :: ModelType    ModelType ;
    typedef typename Traits :: FluxType     AdvectionFluxType;

    typedef FemDGBaseDiscreteModel< Traits, passUId, passGradId >  BaseType;
    //typedef Fem::DGDiscreteModelDefaultWithInsideOutside
    //    < Traits, passUId, passGradId >                    BaseType;

    // These type definitions allow a convenient access to arguments of pass.
    integral_constant< int, passUId > uVar;

  public:
    enum { dimDomain = Traits :: dimDomain };
    enum { dimRange  = Traits :: dimRange };

    enum { advection = returnAdvectionPart  };
    enum { evaluateJacobian = false };

    // TODO: extract this from the model
    typedef FieldVector< double, dimDomain >               DomainType;
    typedef FieldVector< double, dimDomain-1 >             FaceDomainType;


    typedef typename Traits :: GridPartType                            GridPartType;
    typedef typename Traits :: GridType                                GridType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename BaseType :: EntityType                            EntityType;
    typedef typename Traits :: RangeFieldType                          RangeFieldType;
    typedef typename Traits :: DomainFieldType                         DomainFieldType;
    typedef typename Traits :: RangeType                               RangeType;
    typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

    // discrete function storing the adaptation indicator information
    typedef typename Traits :: AdaptationHandlerType   AdaptationType ;

  public:
    /**
     * @brief constructor
     */
    AdvectionModel(const ModelType& mod,
                   const AdvectionFluxType& numf)
      : model_(mod),
        numflux_( numf )
    {
    }

    //! copy constructor (for thread parallel progs mainly)
    AdvectionModel( const AdvectionModel& other )
      : BaseType( other ),
        model_( other.model_ ),
        numflux_( other.numflux_ )
    {
    }

    void setTime ( double time ) { const_cast< ModelType & >( model_ ).setTime( time ); }

    //! dummy method
    void switchUpwind() const {}

    inline bool hasSource() const
    {
      return model_.hasNonStiffSource();
    }  /*@\label{dm:hasSource}@*/

    inline bool hasFlux() const { return advection; }

    /**
     * @brief Stiff source associated with advection
     */
    template <class LocalEvaluation>
    inline double source( const LocalEvaluation& local,
                          RangeType& s ) const
    {
      return model_.nonStiffSource( local.entity(), local.time(), coordinate( local.point() ), local.values()[ uVar ], s );
    }


    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeIntersection(const Intersection& it,
                                const double time,
                                const QuadratureImp& quadInner,
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec)
    {
    }

    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeBoundary(const Intersection& it,
                            const double time,
                            const QuadratureImp& quadInner,
                            const ArgumentTupleVector& uLeftVec)
    {
    }

  public:
    /**
     * @brief flux function on interfaces between cells for advection and diffusion
     *
     * @param[in] it intersection
     * @param[in] time current time given by TimeProvider
     * @param[in] x coordinate of required evaluation local to \c it
     * @param[in] uLeft DOF evaluation on this side of \c it
     * @param[in] uRight DOF evaluation on the other side of \c it
     * @param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * @param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * @note For dual operators we have \c gDiffLeft = 0 and \c gDiffRight = 0.
     *
     * @return wave speed estimate (multiplied with the integration element of the intersection),
     *              to estimate the time step |T|/wave.
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      gDiffLeft = 0;
      gDiffRight = 0;

      if( advection )
      {
        // returns advection wave speed
        return numflux_.numericalFlux(left.intersection(), left.entity(), right.entity(),
                                      left.time(), left.quadrature(), right.quadrature(), left.index(),
                                      left.values()[ uVar ], right.values()[ uVar ], gLeft, gRight);
      }
      else
      {
        gLeft = 0;
        gRight = 0;
        return 0.0;
      }
    }

    /**
     * @brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
    {
      const FaceDomainType& x = left.localPoint();

      const bool hasBndValue = boundaryValue( left );

      // make sure user sets specific boundary implementation
      gLeft = std::numeric_limits< double >::quiet_NaN();
      gDiffLeft = 0;

      if (advection)
      {
        if( hasBndValue )
        {
          RangeType gRight;
          // returns advection wave speed
          return numflux_.numericalFlux(left.intersection(), left.entity(), left.entity(),
                                        left.time(), left.quadrature(), left.quadrature(), left.index(),
                                        left.values()[ uVar ], uBnd_, gLeft, gRight);
        }
        else
        {
          // returns advection wave speed
          return model_.boundaryFlux( left.intersection(), left.time(), left.localPoint(), left.values()[uVar], gLeft );
        }
      }
      else
      {
        gLeft = 0.;
        return 0.;
      }
    }
                                                  /*@LST0S@*/
    /**
     * @brief analytical flux function for advection only
     */
    template <class LocalEvaluation>
    void analyticalFlux( const LocalEvaluation& local,
                         JacobianRangeType& f ) const
    {
      if( advection )
        model_.advection( local.entity(), local.time(), local.point(), local.values()[ uVar ], f);
      else
        f = 0;
    }


  protected:
    template <class LocalEvaluation>
    bool boundaryValue(const LocalEvaluation& left) const
    {
      const FaceDomainType& x = left.localPoint();
      const bool hasBndValue = model_.hasBoundaryValue( left.intersection(), left.time(), left.localPoint() );
      if( hasBndValue )
      {
        model_.boundaryValue( left.intersection(), left.time(), left.localPoint(),
                              left.values()[ uVar ], uBnd_ );
      }
      else
        // do something bad to uBnd_ as it shouldn't be used
        uBnd_ = std::numeric_limits< double >::quiet_NaN();

      return hasBndValue;
    }

    const ModelType&   model_;
    const AdvectionFluxType& numflux_;
    mutable RangeType uBnd_;
  };                                              /*@LST0E@*/


  //////////////////////////////////////////////////////
  //
  // AdaptiveAdvectionModel
  //
  //////////////////////////////////////////////////////
  template< class OpTraits,
            int passUId, int passGradId,
            bool returnAdvectionPart>
  class AdaptiveAdvectionModel ;

  template <class OpTraits,
            int passUId, int passGradId, bool returnAdvectionPart>
  struct AdaptiveAdvectionTraits
    : public AdvectionTraits< OpTraits, passUId, passGradId,
                              returnAdvectionPart >
  {
    typedef AdaptiveAdvectionModel< OpTraits, passUId, passGradId,
                                    returnAdvectionPart >       DGDiscreteModelType;
  };

  /*  \class AdvectionModel
   *
   *  \tparam Model Mathematical model
   *  \tparam NumFlux Numerical flux
   *  \tparam polOrd Polynomial degree
   *  \tparam passUId The id of a pass whose value is used here
   *  \tparam passGradId The id of a pass whose value is used here
   *  \tparam returnAdvectionPart Switch on/off the advection
   */
  template< class OpTraits,
            int passUId, int passGradId,
            bool returnAdvectionPart>
  class AdaptiveAdvectionModel
    : public AdvectionModel< OpTraits, passUId,  passGradId, returnAdvectionPart >
  {
  public:
    typedef AdaptiveAdvectionTraits< OpTraits, passUId, passGradId, returnAdvectionPart> Traits;

    typedef AdvectionModel< OpTraits, passUId,  passGradId, returnAdvectionPart > BaseType ;

    // These type definitions allow a convenient access to arguments of pass.
    integral_constant< int, passUId > uVar;

  public:
    enum { dimDomain = Traits :: dimDomain };
    enum { dimRange  = Traits :: dimRange };

    typedef typename BaseType :: ModelType          ModelType;
    typedef typename BaseType :: AdvectionFluxType  AdvectionFluxType;

    typedef typename BaseType :: DomainType      DomainType ;
    typedef typename BaseType :: FaceDomainType  FaceDomainType;

    typedef typename Traits :: GridPartType                            GridPartType;
    typedef typename Traits :: GridType                                GridType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename BaseType :: EntityType                            EntityType;
    typedef typename Traits :: RangeFieldType                          RangeFieldType;
    typedef typename Traits :: DomainFieldType                         DomainFieldType;
    typedef typename Traits :: RangeType                               RangeType;
    typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

    // discrete function storing the adaptation indicator information
    typedef typename Traits :: AdaptationHandlerType    AdaptationType ;

    typedef typename AdaptationType :: LocalIndicatorType  LocalIndicatorType;

    // type of thread filter in thread parallel runs
    typedef Fem :: ThreadFilter<GridPartType> ThreadFilterType;

  public:
    /**
     * @brief constructor
     */
    AdaptiveAdvectionModel(const ModelType& mod,
                           const AdvectionFluxType& numf)
      : BaseType( mod, numf ),
        adaptation_( 0 ),
        threadFilter_( 0 ),
        enIndicator_(),
        nbIndicator_(),
        weight_( 1 )
    {
    }

    //! copy constructor (for thread parallel progs mainly)
    AdaptiveAdvectionModel( const AdaptiveAdvectionModel& other )
      : BaseType( other ),
        adaptation_( other.adaptation_ ),
        threadFilter_( other.threadFilter_ ),
        enIndicator_( other.enIndicator_ ),
        nbIndicator_( other.nbIndicator_ ),
        weight_( other.weight_ )
    {
    }

    void setEntity( const EntityType& entity )
    {
      BaseType :: setEntity( entity );

      if( adaptation_ )
        enIndicator_ = adaptation_->localIndicator( entity );
    }

    void setNeighbor( const EntityType& neighbor )
    {
      BaseType :: setNeighbor( neighbor );

      if( adaptation_ )
      {
#ifdef USE_SMP_PARALLEL
        // if the neighbor does not belong to our
        // thread domain reset the pointer
        // to avoid update of indicators
        if( threadFilter_ &&
            ! threadFilter_->contains( neighbor ) )
        {
          // remove neighbor indicator
          nbIndicator_.reset();
        }
        else
#endif
        {
          nbIndicator_ = adaptation_->localIndicator( neighbor );
        }
      }
    }

    //! set pointer to adaptation indicator
    void setAdaptation( AdaptationType& adaptation, double weight )
    {
      adaptation_ = & adaptation;
      threadFilter_ = 0;
      weight_ = weight ;
    }

    //! set pointer to adaptation indicator
    void setAdaptation( AdaptationType& adaptation,
                        const ThreadFilterType& threadFilter,
                        double weight )
    {
      adaptation_   = & adaptation;
      threadFilter_ = & threadFilter ;
      weight_ = weight ;
    }

    //! remove pointer to adaptation indicator
    void removeAdaptation()
    {
      adaptation_   = 0 ;
      threadFilter_ = 0 ;
    }

  public:
    /**
     * @brief flux function on interfaces between cells for advection and diffusion
     *
     * @param[in] it intersection
     * @param[in] time current time given by TimeProvider
     * @param[in] x coordinate of required evaluation local to \c it
     * @param[in] uLeft DOF evaluation on this side of \c it
     * @param[in] uRight DOF evaluation on the other side of \c it
     * @param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * @param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * @note For dual operators we have \c gDiffLeft = 0 and \c gDiffRight = 0.
     *
     * @return wave speed estimate (multiplied with the integration element of the intersection),
     *              to estimate the time step |T|/wave.
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      if( ! model_.allowsRefinement( left.intersection(), left.time(), left.localPoint() ) )
        return 0.;

      double ldt = BaseType :: numericalFlux( left, gLeft, gRight, gDiffLeft, gDiffRight );

      if( BaseType :: advection && adaptation_ )
      {
        RangeType error ;
        RangeType v ;
        // v = g( ul, ul ) = f( ul )
        numflux_.numericalFlux(left.intersection(), left.entity(), right.entity(),
                               left.time(), left.quadrature(), right.quadrature(), left.index(),
                               left.values()[ uVar ], left.values(), v, error);

        RangeType w ;
        // v = g( ur, ur ) = f( ur )
        numflux_.numericalFlux(left.intersection(), left.entity(), right.entity(),
                               left.time(), left.quadrature(), right.quadrature(), left.index(),
                               right.values()[ uVar ], right.values()[ uVar ], w, error);

        // err = 2 * g(u,v) - g(u,u) - g(v,v)
        // 2 * g(u,v) = gLeft + gRight
        error  = gLeft ;
        error += gRight ;  // gRight +  gLeft  = 2*g(v,w)

        error -= v;
        error -= w;

        for( int i=0; i<dimRange; ++i )
        {
          error[ i ] = std::abs( error[ i ] );
          /*
          // otherwise ul = ur
          if( std::abs( error[ i ] > 1e-12 ) )
            error[ i ] /= (uLeft[ uVar ][i] - uRight[ uVar ][i] );
            */
        }

        // get face volume
        const double faceVol = left.intersection().geometry().integrationElement( left.localPoint() );

        // calculate grid width
        double weight = weight_ *
          (0.5 * ( this->enVolume() + this->nbVolume() ) / faceVol );

        // add error to indicator
        enIndicator_.add( error, weight );
        nbIndicator_.addChecked( error, weight );
      }

      return ldt ;
    }

    /**
     * @brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      return BaseType :: boundaryFlux( left, gLeft, gDiffLeft );
    }

  protected:
    using BaseType :: inside ;
    using BaseType :: outside ;
    using BaseType :: model_ ;
    using BaseType :: numflux_ ;

    AdaptationType*  adaptation_;
    const ThreadFilterType*  threadFilter_;

    mutable LocalIndicatorType enIndicator_;
    mutable LocalIndicatorType nbIndicator_;

    double weight_ ;
  };

} // end namespace Dune

#endif
