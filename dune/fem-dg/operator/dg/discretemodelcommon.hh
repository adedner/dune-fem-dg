#ifndef DUNE_FEM_DG_DISCRETEMODELCOMMON_HH
#define DUNE_FEM_DG_DISCRETEMODELCOMMON_HH

// Dune includes
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/fvspace.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/operator/fluxes/ldgflux.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>

namespace Dune {

  //PassTraits
  //----------

  template <class Model,int dimRange,int polOrd>
  class PassTraits
  {
  public:
    typedef typename Model :: Traits                                 ModelTraits;
    typedef typename ModelTraits :: GridPartType                     GridPartType;
    typedef typename GridPartType :: GridType                        GridType;
    typedef typename GridType :: ctype                               ctype;
    static const int dimDomain = Model :: Traits :: dimDomain;

    //typedef ElementQuadrature< GridPartType, 0 >                     VolumeQuadratureType;
    typedef CachingQuadrature< GridPartType, 0 >                     VolumeQuadratureType;
    typedef CachingQuadrature< GridPartType, 1 >                     FaceQuadratureType;
    //typedef ElementQuadrature< GridPartType, 1 >                     FaceQuadratureType;

    // Allow generalization to systems
    typedef FunctionSpace< ctype, double, dimDomain, dimRange >      FunctionSpaceType;
    typedef DiscontinuousGalerkinSpace< FunctionSpaceType,
                                        GridPartType, polOrd,
                                        CachingStorage >             DiscreteFunctionSpaceType;
    //typedef LegendreDiscontinuousGalerkinSpace< FunctionSpaceType,
    //                                    GridPartType, polOrd,
    //                                    CachingStorage >             DiscreteFunctionSpaceType;
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >    DestinationType;

    // Indicator for Limiter
    typedef FunctionSpace< ctype, double, dimDomain, 3> FVFunctionSpaceType;
    typedef FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, SimpleStorage> IndicatorSpaceType;
    typedef AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;

    typedef AdaptationHandler< GridType, FunctionSpaceType >  AdaptationHandlerType ;
  };


  // AdvectionModel
  //---------------

  template< class Model, 
            class NumFlux, 
            int polOrd, int passUId, int passGradId,
            bool returnAdvectionPart> 
  class AdvectionModel;


  // AdvectionTraits
  //----------------

  template <class Model, class NumFlux,
            int polOrd, int passUId, int passGradId, bool returnAdvectionPart>
  struct AdvectionTraits
  {
    typedef typename Model :: Traits                                 ModelTraits;
    typedef typename ModelTraits :: GridType                         GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits< Model, dimRange, polOrd >                    Traits;
    typedef typename Traits :: FunctionSpaceType                     FunctionSpaceType;

    typedef typename Traits :: VolumeQuadratureType                  VolumeQuadratureType;
    typedef typename Traits :: FaceQuadratureType                    FaceQuadratureType;
    typedef typename Traits :: GridPartType                          GridPartType;
    typedef typename Traits :: DiscreteFunctionSpaceType             DiscreteFunctionSpaceType;
    typedef typename Traits :: DestinationType                       DestinationType;
    typedef DestinationType                                          DiscreteFunctionType;
    typedef typename Traits :: IndicatorType                         IndicatorType;

    typedef typename DestinationType :: DomainType                   DomainType;
    typedef typename DestinationType :: RangeType                    RangeType;
    typedef typename DestinationType :: RangeFieldType               RangeFieldType;
    typedef typename DestinationType :: DomainFieldType              DomainFieldType;
    typedef typename DestinationType :: JacobianRangeType            JacobianRangeType;

    typedef typename Traits :: AdaptationHandlerType  AdaptationHandlerType ;

    typedef AdvectionModel
      < Model, NumFlux, polOrd, passUId, passGradId, returnAdvectionPart >       DGDiscreteModelType;
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
  template< class Model, 
            class NumFlux, 
            int polOrd, int passUId, int passGradId,
            bool returnAdvectionPart> 
  class AdvectionModel :
    public DGDiscreteModelDefaultWithInsideOutside
      <AdvectionTraits<Model, NumFlux, polOrd, passUId, passGradId, returnAdvectionPart>,
       passUId, passGradId>
  {
  public:
    typedef AdvectionTraits 
      <Model, NumFlux, polOrd, passUId, passGradId, returnAdvectionPart> Traits;

    typedef DGDiscreteModelDefaultWithInsideOutside
              < Traits, passUId, passGradId >                          BaseType;

    // These type definitions allow a convenient access to arguments of pass.
#if DUNE_VERSION_NEWER_REV(DUNE_COMMON,2,1,0)
    integral_constant< int, passUId > uVar;
#else
    const Int2Type<passUId> uVar; /*@\label{dm:int2type0}@*/
#endif

  public:
    enum { dimDomain = Traits :: dimDomain };
    enum { dimRange  = Traits :: dimRange };

    enum { advection = returnAdvectionPart  };
    enum { evaluateJacobian = false };

    typedef FieldVector< double, dimDomain >               DomainType;
    typedef FieldVector< double, dimDomain-1 >             FaceDomainType;


    typedef typename Traits :: GridPartType                            GridPartType;
    typedef typename Traits :: GridType                                GridType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename GridType :: template Codim< 0 > :: Entity         EntityType;
    typedef typename GridType :: template Codim< 0 > :: EntityPointer  EntityPointerType;
    typedef typename Traits :: RangeFieldType                          RangeFieldType;
    typedef typename Traits :: DomainFieldType                         DomainFieldType;
    typedef typename Traits :: RangeType                               RangeType;
    typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

    // discrete function storing the adaptation indicator information 
    typedef typename Traits :: IndicatorType          IndicatorTpye;

    // discrete function storing the adaptation indicator information 
    typedef typename Traits :: AdaptationHandlerType  AdaptationHandlerType ;

  public:
    /**
     * @brief constructor
     */
    AdvectionModel(const Model& mod,
                   const NumFlux& numf)
      : model_(mod),
        numflux_( numf ),
        adaptation_( 0 ),
        weight_( 1 ),
        maxAdvTimeStep_( 0.0 ),
        maxDiffTimeStep_( 0.0 )
    {
    }

    //! copy constructor (for thread parallel progs mainly)
    AdvectionModel( const AdvectionModel& other )
      : BaseType( other ),
        model_( other.model_ ),
        numflux_( other.numflux_ ),
        weight_( other.weight_ ),
        adaptation_( other.adaptation_ ),
        maxAdvTimeStep_( other.maxAdvTimeStep_ ),
        maxDiffTimeStep_( other.maxDiffTimeStep_ )
    {
    }

#ifdef APOST_ERROR_INDICATOR
    void setEntity( const EntityType& entity ) 
    {
      BaseType :: setEntity( entity );

      if( adaptation_ ) 
        adaptation_->setEntity( entity );
    }

    void setNeighbor( const EntityType& neighbor ) 
    {
      BaseType :: setNeighbor( neighbor );

      if( adaptation_ ) 
        adaptation_->setNeighbor( neighbor );
    }
#endif

    //! dummy method 
    void switchUpwind() const {
      maxAdvTimeStep_ = 0;
      maxDiffTimeStep_ = 0;
    } 

    //! set pointer to adaptation indicator 
    void setAdaptationHandler( AdaptationHandlerType& adaptation, double weight ) 
    {
      adaptation_ = & adaptation;       
      weight_ = weight ;
    }

    //! remove pointer to adaptation indicator 
    void removeAdaptationHandler() 
    {
      adaptation_ = 0 ;
    }

    inline bool hasSource() const 
    { 
      return model_.hasNonStiffSource(); 
    }  /*@\label{dm:hasSource}@*/

    inline bool hasFlux() const { return advection; }

    //! return true if diffusion time step is defining the time step size (default is false)
    double maxAdvectionTimeStep() const { return maxAdvTimeStep_; }
    double maxDiffusionTimeStep() const { return maxDiffTimeStep_; }

    //! this method is needed for thread pass 
    void setMaxTimeSteps( const double advStep, const double diffStep ) 
    {
      maxAdvTimeStep_  = advStep ;
      maxDiffTimeStep_ = diffStep ;
    }

    /**
     * @brief Stiff source associated with advection
     */
    template <class ArgumentTuple, class JacobianTuple >
    inline double source( const EntityType& en,
                 const double time, 
                 const DomainType& x,
                 const ArgumentTuple& u, 
                 const JacobianTuple& jac, 
                 RangeType& s ) const
    {
      return model_.nonStiffSource( en, time, x, u[ uVar ], s );
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
    template <class QuadratureImp,
              class ArgumentTuple, 
              class JacobianTuple >          /*@LST0S@*/
    double numericalFlux(const Intersection& it,
                         const double time,
                         const QuadratureImp& faceQuadInner,
                         const QuadratureImp& faceQuadOuter,
                         const int quadPoint, 
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         const JacobianTuple& jacLeft,
                         const JacobianTuple& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      gDiffLeft = 0;
      gDiffRight = 0;

      if( advection ) 
      {
        double ldt = numflux_.numericalFlux(it, this->inside(), this->outside(),
                                           time, faceQuadInner, faceQuadOuter, quadPoint, 
                                           uLeft[ uVar ], uRight[ uVar ], gLeft, gRight);

#ifdef APOST_ERROR_INDICATOR
        if( adaptation_ ) 
        {
          RangeType v, w, error ;
          // v = g( ul, ul )
          numflux_.numericalFlux(it, this->inside(), this->inside(),
                                 time, faceQuadInner, faceQuadInner, quadPoint, 
                                 uLeft[ uVar ], uLeft[ uVar ], v, error);

          // w = g( ur, ur )
          numflux_.numericalFlux(it, this->outside(), this->outside(),
                                 time, faceQuadOuter, faceQuadInner, quadPoint, 
                                 uRight[ uVar ], uRight[ uVar ], w, error);

          // err = 2 * g(u,v) - g(u,u) - g(v,v) 
          // 2 * g(u,v) = gLeft + gRight
          error  = gLeft;
          error += gRight;

          error -= v;
          error -= w;

          const DomainType normal = it.integrationOuterNormal( faceQuadInner.localPoint( quadPoint ) );

          // calculate grid width 
          double weight = weight_ *
            (0.5 * ( this->enVolume() + this->nbVolume() ) / normal.two_norm() );

          error *= weight ;

          // add error to indicator 
          adaptation_->addToLocalIndicator(error);
          adaptation_->addToNeighborIndicator(error);
        }
#endif
        return ldt ;
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
    template <class QuadratureImp, 
              class ArgumentTuple, class JacobianTuple>
    double boundaryFlux(const Intersection& it,
                        const double time, 
                        const QuadratureImp& faceQuadInner,
                        const int quadPoint,
                        const ArgumentTuple& uLeft,
                        const JacobianTuple& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
    {
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );

      const bool hasBndValue = boundaryValue( it, time, 
                                              faceQuadInner, quadPoint,
                                              uLeft );

      // make sure user sets specific boundary implementation
      gLeft = std::numeric_limits< double >::quiet_NaN();
      gDiffLeft = 0;

      if (advection)
      {
        if( hasBndValue )
        {
          RangeType gRight;
          return numflux_.numericalFlux(it, this->inside(), this->inside(),
                                        time, faceQuadInner, faceQuadInner, quadPoint, 
                                        uLeft[ uVar ], uBnd_, gLeft, gRight);
        }
        else 
        {
          return model_.boundaryFlux( it, time, x, uLeft[uVar], gLeft );
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
    template <class ArgumentTuple, class JacobianTuple >
    void analyticalFlux( const EntityType& en,
                         const double time, 
                         const DomainType& x,
                         const ArgumentTuple& u, 
                         const JacobianTuple& jac, 
                         JacobianRangeType& f ) const
    {
      if( advection ) 
        model_.advection(en, time, x, u[ uVar ], f);
      else 
        f = 0;
    }


  protected:
    template <class QuadratureImp, 
              class ArgumentTuple>
    bool boundaryValue(const Intersection& it,
                       const double time, 
                       const QuadratureImp& faceQuadInner,
                       const int quadPoint,
                       const ArgumentTuple& uLeft ) const 
    {
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
      const bool hasBndValue = model_.hasBoundaryValue(it, time, x);
      if( hasBndValue ) 
      {
        model_.boundaryValue(it, time, x, uLeft[ uVar ], uBnd_ );
      }
      else 
        // do something bad to uBnd_ as it shouldn't be used
        uBnd_ = std::numeric_limits< double >::quiet_NaN();

      return hasBndValue;
    }

    const Model&   model_;
    const NumFlux& numflux_;
    AdaptationHandlerType* adaptation_;
    mutable RangeType uBnd_;
    double weight_ ;
    mutable double maxAdvTimeStep_;
    mutable double maxDiffTimeStep_;
  };                                              /*@LST0E@*/



  // LimiterDiscreteModel
  //---------------------

  template <class GlobalTraitsImp, class Model, int passId = -1 >
  class LimiterDiscreteModel;
  

  // LimiterTraits
  //--------------

  template <class GlobalTraitsImp, class Model, int passId >
  struct LimiterTraits 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef GlobalTraitsImp Traits; 
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef typename Traits::IndicatorType IndicatorType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainFieldType DomainFieldType;
    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;
    typedef FieldVector<DomainFieldType, dimDomain - 1> FaceDomainType;

    typedef LimiterDiscreteModel <GlobalTraitsImp,Model,passId> DGDiscreteModelType;

    typedef MinModLimiter< DomainFieldType > LimiterFunctionType;
    //typedef SuperBeeLimiter< DomainFieldType > LimiterFunctionType;
    //typedef VanLeerLimiter< DomainFieldType > LimiterFunctionType;
  };
  
  
  // LimiterDiscreteModel
  //---------------------

  template <class GlobalTraitsImp, class Model, int passId >
  class LimiterDiscreteModel :
    public LimiterDefaultDiscreteModel
      < LimiterTraits< GlobalTraitsImp,Model,passId >, Model, passId >  
  {
    typedef LimiterDefaultDiscreteModel
          < LimiterTraits< GlobalTraitsImp,Model,passId >, Model, passId > BaseType;
  public:
    typedef LimiterTraits<GlobalTraitsImp,Model,passId> Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::FaceDomainType FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename Traits::IndicatorType IndicatorType;

    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    typedef typename DomainType :: field_type DomainFieldType;

    typedef typename Traits :: LimiterFunctionType  LimiterFunctionType;

    enum { dimRange = RangeType :: dimension };
    
  public:
    /** \brief default limiter discrete model */
    LimiterDiscreteModel( const Model& mod, 
                          const DomainFieldType veloEps = 1e-8 )
      : BaseType( mod, veloEps )
      , refTol_(1), crsTol_(0.1), finLevel_(0), crsLevel_(0)
    {}

    /** \brief copy constructor */
    LimiterDiscreteModel(const LimiterDiscreteModel& org) 
      : BaseType( org )
    {}
    
    void setIndicator( IndicatorType* ind )
    {
      indicator_ = ind;
    }


    enum { Coarsen = -1 , None = 0, Refine = 1 };


    void adaptation(GridType& grid, 
                    const EntityType& en,
                    const RangeType& shockIndicator, 
                    const RangeType& adaptIndicator) const 
    {
      const double val = adaptIndicator[0] / max_[1];
      // set indicator 
      if( indicator_ )
      {
        typedef typename IndicatorType :: LocalFunctionType  LocalFunctionType;
        LocalFunctionType lf = indicator_->localFunction( en );
        assert( lf.numDofs() == 3 );
        lf[0] = shockIndicator[0];
        lf[1] = adaptIndicator[0];
        lf[2] = val; // store real adapt indicator 
      }

      // get refinement marker 
      int refinement = grid.getMark( en );
      
      // if element already is marked for refinement then do nothing 
      if( refinement >= Refine ) return ;
      
      // use component max 1 (max of adapt indicator)
      {
        const int level = en.level();
        if( (( val > refTol_ ) || (val < 0)) && level < finLevel_) 
        {
          refinement = Refine;
        }
        else if ( (val >= 0) && (val < crsTol_)  && level > crsLevel_ ) 
        {
          if( refinement == None ) 
          {
            refinement = Coarsen;
          }
        }
      }

      // set new refinement marker 
      grid.mark( refinement, en );
    }

    /** \brief returns true if model provides physical check */
    bool hasPhysical() const { return this->model_.hasPhysical(); } 

    /** \brief check physical values */
    bool physical( const EntityType& en,
                   const DomainType& x,
                   const RangeType& u ) const
    {
      return this->model_.physical( en, x, u );
    }

  public:
    using BaseType :: model_;
    IndicatorType* indicator_;
    mutable RangeType max_;

    double refTol_;
    double crsTol_;
    int finLevel_;
    int crsLevel_;

  protected:
    //const Model& model_;
    using BaseType :: velocity_;
    using BaseType :: veloEps_;
  };



} // end namespace Dune

#endif
