#ifndef DUNE_FEM_DG_DISCRETEMODELCOMMON_HH
#define DUNE_FEM_DG_DISCRETEMODELCOMMON_HH

#include <type_traits>

#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem-dg/pass/discretemodel.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem-dg/operator/limiter/limiter.hh>

#include <dune/fem/gridpart/filter/domainfilter.hh>

namespace Dune
{
namespace Fem
{

  template< class >
  struct Sequence2Tuple;

  template< int... i >
  struct Sequence2Tuple< std::integer_sequence< int, i... > >
  {
    typedef std::tuple< std::integral_constant< int, i >... > type;
  };

  /**
   * \brief A discrete model that can take additional model parameters
   *
   * This discrete model allows us to "insert" additional model parameters.
   *
   * \ingroups DiscreteModels
   * \ingroups Pass
   */
  template< class Traits, int... passIds >
  class FemDGBaseDiscreteModel
    : public Fem::DGDiscreteModelDefaultWithInsideOutside< Traits, passIds... >
  {
    static const int modParamSize = Traits::ModelType::modelParameterSize;
    static const int extParamSize = std::tuple_size< typename Traits::ExtraParameterTupleType >::value ;

    static_assert( modParamSize <= extParamSize, "Specify enough ExtraParameterTuple elements!" );

    static_assert( modParamSize >= 0, "invalid integer_sequence: Throw this additional assertion here because \
                                       gcc-6 won't stop compiling for a very long time... :(" );
    typedef typename Sequence2Tuple< std::make_integer_sequence< int, modParamSize > >::type ExtraParameterIds;
  public:

    // overload selector type to add model parameters
    typedef typename Dune::Fem::VariadicSelectorBase< ExtraParameterIds, passIds... >::Type  Selector;
  };

  // AdvectionModel
  //---------------

  template< class OpTraits, bool enableAdvection, int passUId, int... passIds >
  class AdvectionModel;


  // AdvectionTraits
  //----------------

  template <class Traits, bool enableAdvection, int... passIds >
  struct AdvectionTraits : public Traits
  {
    typedef AdvectionModel< Traits, enableAdvection, passIds... >       DGDiscreteModelType;
  };


  /*  \brief Discrete model for advection terms
   *
   *  \ingroup DiscreteModels
   *  \ingroup PassBased
   *
   *  \tparam OpTraits Operator traits describing the operator
   *  \tparam passUId The id of a pass whose value is used here
   *  \tparam passGradId The id of a pass whose value is used here
   *  \tparam enableAdvection Switch on/off the advection
   */
  template< class OpTraits, bool enableAdvection, int passUId, int... passIds >
  class AdvectionModel :
    public FemDGBaseDiscreteModel< AdvectionTraits< OpTraits, enableAdvection, passUId, passIds... >,
                                   passUId, passIds... >
  {
  public:
    typedef AdvectionTraits< OpTraits, enableAdvection, passUId, passIds... > Traits;

    typedef typename Traits::ModelType                                        ModelType;
    typedef typename Traits::AdvectionFluxType                                AdvectionFluxType;

    typedef FemDGBaseDiscreteModel< Traits, passUId, passIds... >             BaseType;

    // These type definitions allow a convenient access to arguments of pass.
    std::integral_constant< int, passUId > uVar;

  public:
    enum { advection = enableAdvection };
    enum { evaluateJacobian = false };

    typedef typename Traits::GridPartType                            GridPartType;
    typedef typename Traits::GridType                                GridType;
    typedef typename Traits::DiscreteFunctionSpaceType               DiscreteFunctionSpaceType;

    typedef typename GridPartType::IntersectionIteratorType          IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType                  IntersectionType;
    typedef typename BaseType::EntityType                            EntityType;
    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType    FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeFieldType               RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType              DomainFieldType;
    typedef typename FunctionSpaceType::DomainType                   DomainType;
    typedef typename FunctionSpaceType::RangeType                    RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType            JacobianRangeType;


    // discrete function storing the adaptation indicator information
    typedef typename Traits::AdaptationHandlerType                   AdaptationType;

  public:
    /**
     * \brief constructor
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

    void setEntity( const EntityType& entity )
    {
      BaseType::setEntity( entity );
      model_.setEntity( entity );
    }

    void setTime ( double time ) { const_cast< ModelType & >( model_ ).setTime( time ); }

    //! dummy method
    void switchUpwind() const {}

    //! return true if source term is present
    inline bool hasSource() const
    {
      return model_.hasNonStiffSource();
    }

    //! return true if flux term is present
    inline bool hasFlux() const { return advection; }

    //! return true if the mass term is not the identity
    inline bool hasMass() const { return model_.hasMass(); }

    /**
     * \brief Stiff source associated with advection
     */
    template <class LocalEvaluation>
    inline double source( const LocalEvaluation& local,
                          RangeType& s ) const
    {
      return model_.nonStiffSource( local, local.values()[uVar], s );
    }

    template <class LocalEvaluationVec >
    void initializeIntersection(const LocalEvaluationVec& left,
                                const LocalEvaluationVec& right )
    {}

    template <class LocalEvaluationVec >
    void initializeBoundary(const LocalEvaluationVec& local )
    {}

  public:
    /**
     * \brief flux function on interfaces between cells for advection and diffusion
     *
     * \param[in]  left local evaluation context of inside cell
     * \param[in]  right local evaluation context of outside cell
     * \param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * \param[out] gRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \phi \f$
     * \param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * \param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * \note For dual operators we have \c gDiffLeft = 0 and \c gDiffRight = 0.
     *
     * \return wave speed estimate (multiplied with the integration element of the intersection),
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
        return numflux_.numericalFlux(left, right,
                                      left.values()[uVar], right.values()[uVar],
                                      left.jacobians()[uVar], right.jacobians()[uVar],
                                      gLeft, gRight);
      }
      else
      {
        gLeft = 0;
        gRight = 0;
        return 0.0;
      }
    }

    /**
     * \brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
    {
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
          return numflux_.numericalFlux(left, left,
                                        left.values()[uVar], uBnd_, left.jacobians()[uVar], left.jacobians()[uVar],
                                        gLeft, gRight);
        }
        else
        {
          // returns advection wave speed
          return model_.boundaryFlux( left, left.values()[uVar], left.jacobians()[uVar], gLeft );
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
     * \brief analytical flux function for advection only
     */
    template <class LocalEvaluation>
    void analyticalFlux( const LocalEvaluation& local,
                         JacobianRangeType& f ) const
    {
      if( advection )
        model_.advection( local, local.values()[uVar],local.jacobians()[uVar], f);
      else
        f = 0;
    }

    template <class LocalEvaluation>
    void mass (const LocalEvaluation& local,
               RangeType& m ) const
    {
      assert( hasMass() );
      model_.mass( local, local.values()[uVar], m );
    }

    const ModelType& model() const
    {
      return model_;
    }

  protected:
    template <class LocalEvaluation>
    bool boundaryValue(const LocalEvaluation& left) const
    {
      const bool hasBndValue = model_.hasBoundaryValue( left ) || model_.hasRobinBoundaryValue( left );
      if( hasBndValue )
      {
        model_.boundaryValue( left, left.values()[uVar], uBnd_ );
      }
      else
        // do something bad to uBnd_ as it shouldn't be used
        uBnd_ = std::numeric_limits< double >::quiet_NaN();

      return hasBndValue;
    }

    // store an instance here so that for thread parallel
    // runs class variables are thread private since discrete models
    // are thread private
    ModelType  model_;

    AdvectionFluxType numflux_;
    mutable RangeType uBnd_;
  };                                              /*@LST0E@*/


  //////////////////////////////////////////////////////
  //
  // AdaptiveAdvectionModel
  //
  //////////////////////////////////////////////////////
  template< class OpTraits, bool enableAdvection, int passUId, int... passIds >
  class AdaptiveAdvectionModel;

  template <class OpTraits, bool enableAdvection, int... passIds >
  struct AdaptiveAdvectionTraits
    : public AdvectionTraits< OpTraits, enableAdvection, passIds... >
  {
    typedef AdaptiveAdvectionModel< OpTraits, enableAdvection, passIds... > DGDiscreteModelType;
  };

  /*  \brief discrete model for adaptive operator
   *
   *  \ingroup DiscreteModels
   *  \ingroup PassBased
   *
   *  \tparam OpTraits Operator traits describing the operator
   *  \tparam passUId The id of a pass whose value is used here
   *  \tparam passGradId The id of a pass whose value is used here
   *  \tparam enableAdvection Switch on/off the advection
   */
  template< class OpTraits, bool enableAdvection, int passUId, int... passIds >
  class AdaptiveAdvectionModel
    : public AdvectionModel< OpTraits, enableAdvection, passUId, passIds... >
  {
  public:
    typedef AdaptiveAdvectionTraits< OpTraits, enableAdvection, passUId, passIds... > Traits;

    typedef AdvectionModel< OpTraits, enableAdvection, passUId, passIds... > BaseType ;

    // These type definitions allow a convenient access to arguments of pass.
    std::integral_constant< int, passUId > uVar;

  public:
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::AdvectionFluxType  AdvectionFluxType;


    typedef typename BaseType::GridPartType                          GridPartType;
    typedef typename BaseType::GridType                              GridType;
    typedef typename BaseType::IntersectionIteratorType              IntersectionIteratorType;
    typedef typename BaseType::IntersectionType                      IntersectionType;
    typedef typename BaseType::FunctionSpaceType                     FunctionSpaceType;
    typedef typename BaseType::EntityType                            EntityType;
    typedef typename BaseType::DomainType                            DomainType ;
    typedef typename BaseType::RangeFieldType                        RangeFieldType;
    typedef typename BaseType::DomainFieldType                       DomainFieldType;
    typedef typename BaseType::RangeType                             RangeType;
    typedef typename BaseType::JacobianRangeType                     JacobianRangeType;

    // discrete function storing the adaptation indicator information
    typedef typename Traits::AdaptationHandlerType    AdaptationType ;

    typedef typename AdaptationType::LocalIndicatorType  LocalIndicatorType;

    // type of thread filter in thread parallel runs
    typedef Fem::DomainFilter<GridPartType> ThreadDomainFilterType;

  public:
    /**
     * \brief constructor
     */
    AdaptiveAdvectionModel(const ModelType& mod,
                           const AdvectionFluxType& numf)
      : BaseType( mod, numf ),
        adaptation_( 0 ),
        threadDomainFilter_( 0 ),
        enIndicator_(),
        nbIndicator_(),
        weight_( 1 )
    {
    }

    //! copy constructor (for thread parallel progs mainly)
    AdaptiveAdvectionModel( const AdaptiveAdvectionModel& other )
      : BaseType( other ),
        adaptation_( other.adaptation_ ),
        threadDomainFilter_( other.threadDomainFilter_ ),
        enIndicator_( other.enIndicator_ ),
        nbIndicator_( other.nbIndicator_ ),
        weight_( other.weight_ )
    {
    }

    void setEntity( const EntityType& entity )
    {
      BaseType::setEntity( entity );
      model_.setEntity( entity );

      if( adaptation_ )
        enIndicator_ = adaptation_->localIndicator( entity );
    }

    void setNeighbor( const EntityType& neighbor )
    {
      BaseType::setNeighbor( neighbor );

      if( adaptation_ )
      {
        // if the neighbor does not belong to our
        // thread domain reset the pointer
        // to avoid update of indicators
        if( threadDomainFilter_ &&
            ! threadDomainFilter_->contains( neighbor ) )
        {
          // remove neighbor indicator
          nbIndicator_.reset();
        }
        else
        {
          nbIndicator_ = adaptation_->localIndicator( neighbor );
        }
      }
    }

    //! set pointer to adaptation indicator
    void setAdaptation( AdaptationType& adaptation, const double weight,
                        const ThreadDomainFilterType* threadDomainFilter = nullptr )
    {
      adaptation_ = & adaptation;
      weight_ = weight ;
      threadDomainFilter_ = threadDomainFilter;
    }

    //! remove pointer to adaptation indicator
    void removeAdaptation()
    {
      adaptation_   = 0 ;
      threadDomainFilter_ = 0 ;
    }

  public:
    /**
     * \brief flux function on interfaces between cells for advection and diffusion
     *
     * \param[in]  left local evaluation context of inside cell
     * \param[in]  right local evaluation context of outside cell
     * \param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * \param[out] gRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \phi \f$
     * \param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * \param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * \note For dual operators we have \c gDiffLeft = 0 and \c gDiffRight = 0.
     *
     * \return wave speed estimate (multiplied with the integration element of the intersection),
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
      if( ! model_.allowsRefinement( left ) )
        return 0.;

      double ldt = BaseType::numericalFlux( left, gLeft, gRight, gDiffLeft, gDiffRight );

      if( BaseType::advection && adaptation_ )
      {
        RangeType error ;
        RangeType v ;
        // v = g( ul, ul ) = f( ul )
        numflux_.numericalFlux(left, left,
                               left.values()[uVar], left.values()[uVar], v, error);

        RangeType w ;
        // v = g( ur, ur ) = f( ur )
        numflux_.numericalFlux(right, right,
                               right.values()[uVar], right.values()[uVar], w, error);

        // err = 2 * g(u,v) - g(u,u) - g(v,v)
        // 2 * g(u,v) = gLeft + gRight
        error  = gLeft ;
        error += gRight ;  // gRight +  gLeft  = 2*g(v,w)

        error -= v;
        error -= w;

        const int dimRange = FunctionSpaceType::dimRange;
        for( int i=0; i<dimRange; ++i )
        {
          error[ i ] = std::abs( error[ i ] );
          /*
          // otherwise ul = ur
          if( std::abs( error[ i ] > 1e-12 ) )
            error[ i ] /= (uLeft[uVar][i] - uRight[uVar][i] );
            */
        }

        // get face volume
        const double faceVol = left.intersection().geometry().integrationElement( left.localPosition() );

        // calculate grid width
        double weight = weight_ *
          (0.5 * ( this->enVolume() + this->nbVolume() ) / faceVol );

        // add error to indicator
        enIndicator_.add( error, weight );
        nbIndicator_.addChecked( error, weight );
      }

      return ldt;
    }

    /**
     * \brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      return BaseType::boundaryFlux( left, gLeft, gDiffLeft );
    }

    const ModelType& model() const
    {
      return model_;
    }

  protected:
    using BaseType::inside ;
    using BaseType::outside ;
    using BaseType::model_ ;
    using BaseType::numflux_ ;

    AdaptationType*  adaptation_;
    const ThreadDomainFilterType*  threadDomainFilter_;

    mutable LocalIndicatorType enIndicator_;
    mutable LocalIndicatorType nbIndicator_;

    double weight_ ;
  };

} // end namespace
} // end namespace Dune

#endif
