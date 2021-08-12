#ifndef DUNE_FEM_LIMITERDISCRETEMODEL_HH
#define DUNE_FEM_LIMITERDISCRETEMODEL_HH

#include <type_traits>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>

namespace Dune
{
namespace Fem
{

  template <class GlobalPassTraitsImp, class Model, int passId = -1>
  class LimiterDefaultDiscreteModel;

  template <class GlobalTraitsImp, class Model, int passId >
  struct LimiterDefaultTraits
  {
    typedef Model ModelType;
    typedef typename ModelType::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };
    enum { dimGrid = GridType::dimension };

    typedef GlobalTraitsImp Traits;
    typedef typename ModelTraits::FunctionSpaceType FunctionSpaceType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainFieldType DomainFieldType;
    typedef typename DestinationType::DomainType DomainType;
    typedef FieldVector<DomainFieldType, dimGrid >  LocalDomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;
    typedef FieldVector<DomainFieldType, dimGrid - 1> FaceLocalDomainType;

    // Indicator function type for Limiter (for output mainly)
    typedef Fem::FunctionSpace< DomainFieldType, double, dimDomain, 3> FVFunctionSpaceType;
    typedef Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Fem::SimpleStorage> IndicatorSpaceType;
    typedef Fem::AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;

    typedef LimiterDefaultDiscreteModel <GlobalTraitsImp,Model,passId> DGDiscreteModelType;

    //typedef typename ExistsLimiterFunction< Model, DomainFieldType > :: LimiterFunctionType  LimiterFunctionType;
    //typedef MinModLimiter< DomainFieldType > LimiterFunctionType;
    typedef typename Model :: Traits :: LimiterFunctionType  LimiterFunctionType;
  };

  /**
   * \brief Standard implementation of limiter discret model.
   *
   * \ingroup Pass
   * \ingroup DiscreteModel
   */
  template <class GlobalTraitsImp, class Model, int passId >
  class LimiterDefaultDiscreteModel :
    public Fem::DGDiscreteModelDefaultWithInsideOutside< LimiterDefaultTraits<GlobalTraitsImp,Model,passId >, passId >
  {
    typedef Fem::DGDiscreteModelDefaultWithInsideOutside< LimiterDefaultTraits<GlobalTraitsImp,Model,passId >, passId >  BaseType;

    // These type definitions allow a convenient access to arguments of pass.
    std::integral_constant< int, passId > uVar;
  public:
    typedef Model                                                       ModelType;
    typedef LimiterDefaultTraits<GlobalTraitsImp,Model, passId >        Traits;

    typedef typename Traits::RangeType                                  RangeType;
    typedef typename Traits::DomainType                                 DomainType;
    typedef typename Traits::LocalDomainType                            LocalDomainType;
    typedef typename Traits::DomainFieldType                            DomainFieldType;
    typedef typename Traits::GridType                                   GridType;
    typedef typename Traits::GridPartType                               GridPartType;
    typedef typename Traits::JacobianRangeType                          JacobianRangeType;
    typedef typename GridPartType::IntersectionType                     IntersectionType;
    typedef typename GridPartType::template Codim<0>::EntityType        EntityType;
    typedef typename GridType::template Codim<0>::Entity                GridEntityType;

    enum { dimGrid = GridType :: dimension };

    // type of surface domain type
    typedef FieldVector< DomainFieldType, dimGrid - 1 >                 FaceLocalDomainType;

    typedef typename Traits::DestinationType                            DestinationType;

    // Indicator for Limiter
    typedef typename Traits::IndicatorType                              IndicatorType ;

    // type of limiter function, No, MinMod, VanLeer, SuperBee
    typedef typename Traits::LimiterFunctionType                        LimiterFunctionType;

    enum { dimRange = Traits :: dimRange };
    enum { evaluateJacobian = false };

  protected:
    // store an instance here so that for thread parallel
    // runs class variables are thread private since discrete models
    // are thread private
    Model model_;

    const LimiterFunctionType limiterFunction_;
    mutable DomainType velocity_;
    const DomainFieldType veloEps_;

    IndicatorType* indicator_;

    double refTol_;
    double crsTol_;
    int finLevel_;
    int crsLevel_;
    const AdaptationParameters param;
    const bool shockIndicatorAdaptivty_;

  public:
    using BaseType::inside;
    using BaseType::outside;

    //! constructor
    LimiterDefaultDiscreteModel(const Model& mod,
                                const int polOrd,
                                const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container(),
                                const DomainFieldType veloEps = 1e-8 )
      : model_( mod ),
        limiterFunction_(parameter),
        velocity_(0),
        veloEps_( veloEps ),
        indicator_( 0 ),
        refTol_( -1 ),
        crsTol_( -1 ),
        finLevel_( 0 ),
        crsLevel_( 0 ),
        param(parameter),
        shockIndicatorAdaptivty_( param.shockIndicator() )
    {
      if( shockIndicatorAdaptivty_ )
      {
        refTol_   = param.refinementTolerance();
        crsTol_   = param.coarsenTolerance();
        finLevel_ = param.finestLevel( DGFGridInfo<GridType>::refineStepsForHalf() );
        crsLevel_ = param.coarsestLevel( DGFGridInfo<GridType>::refineStepsForHalf() );
      }
    }

    bool calculateIndicator () const { return model_.calculateIndicator(); }

    void setIndicator(IndicatorType* indicator)
    {
      indicator_ = indicator;
    }

    void setEntity(const EntityType& en)
    {
      BaseType::setEntity( en );
      model_.setEntity ( en );
    }

    //! \brief returns false
    bool hasSource() const { return false; }
    //! \brief returns true
    bool hasFlux() const   { return true;  }

    template <class LocalEvaluationVec >
    void initializeIntersection( const LocalEvaluationVec& left,
                                 const LocalEvaluationVec& right )
    {}

    template <class LocalEvaluationVec >
    void initializeBoundary(const LocalEvaluationVec& local )
    {}

    // old version
    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeIntersection(const IntersectionType& it,
                                const QuadratureImp& quadInner,
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec)
    {
    }

    // old version
    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeBoundary(const IntersectionType& it,
                            const QuadratureImp& quadInner,
                            const ArgumentTupleVector& uLeftVec)
    {
    }

    void indicatorMax() const
    {
    }

    bool isConstant( const RangeType& minVal, const RangeType& maxVal ) const
    {
      return model_.isConstant( minVal, maxVal );
    }

    void obtainBounds( RangeType& globalMin, RangeType& globalMax ) const
    {
      model_.obtainBounds( globalMin, globalMax );
    }

    enum { Coarsen = -1 , None = 0, Refine = 1 };

    void clearIndicator() { if( indicator_ ) indicator_->clear(); }

    void markIndicator()
    {
      // set indicator
      if( indicator_ )
      {
        typedef typename IndicatorType :: RangeType IndicatorRangeType;
        // set all components to 1
        IndicatorRangeType localIndicator( 1 );
        indicator_->setLocalDofs( inside(), localIndicator );
      }
    }

    //! mark element for refinement or coarsening
    void adaptation(GridType& grid,
                    const EntityType& entity,
                    const RangeType& shockIndicator,
                    const RangeType& adaptIndicator) const
    {
      const double val = adaptIndicator[0];

      // set indicator
      if( indicator_ )
      {
        typedef typename IndicatorType :: RangeType IndicatorRangeType;
        IndicatorRangeType localIndicator( 0 );

        localIndicator[0] = shockIndicator[0];
        localIndicator[1] = adaptIndicator[0];
        localIndicator[2] = val; // store real adapt indicator
        indicator_->setLocalDofs( entity, localIndicator );
      }

      // only mark for adaptation if this type of adaptation is enabled
      if( ! shockIndicatorAdaptivty_ ) return ;

      // get real grid entity from possibly wrapped entity
      const GridEntityType& gridEntity = Fem :: gridEntity( entity );

      // get refinement marker
      int refinement = grid.getMark( gridEntity );

      // if element already is marked for refinement then do nothing
      if( refinement >= Refine ) return ;

      // use component 1 (max of adapt indicator)
      {
        const int level = gridEntity.level();
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
      grid.mark( refinement, gridEntity );
    }

    /** \brief numericalFlux of for limiter evaluateing the difference
         of the solution in the current integration point if we are a an
         inflow intersection.
         This is needed for the shock detection.
    */
    template <class FaceQuadratureImp,
              class ArgumentTuple,
              class JacobianTuple>
    double numericalFlux(const IntersectionType& it,
                         const double time,
                         const FaceQuadratureImp& innerQuad,
                         const FaceQuadratureImp& outerQuad,
                         const int quadPoint,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         const JacobianTuple& jacLeft,
                         const JacobianTuple& jacRight,
                         RangeType& shockIndicator,
                         RangeType& adaptIndicator,
                         JacobianRangeType&,
                         JacobianRangeType& ) const
    {
      std::abort();
      const FaceLocalDomainType& x = innerQuad.localPoint( quadPoint );

      if( checkDirection(it,time, x, uLeft) )
      {
        shockIndicator  = uLeft[ uVar ];
        shockIndicator -= uRight[ uVar ];
        adaptIndicator = shockIndicator;
        return it.integrationOuterNormal( x ).two_norm();
      }
      else
      {
        adaptIndicator = shockIndicator = 0;
        return 0.0;
      }
    }

    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& shockIndicator,
                         RangeType& adaptIndicator,
                         JacobianRangeType&,
                         JacobianRangeType& ) const
    {
      const FaceLocalDomainType& x = left.localPosition();

      if (! physical(left.entity(),  left.position(),  left.values()[ uVar ] ) ||
          ! physical(right.entity(), right.position(), right.values()[ uVar ] ) )
      {
        adaptIndicator = shockIndicator = 1e10;
        return -1.;
      }
      else
      {
        // evaluate adaptation indicator
        model_.adaptationIndicator( left.intersection(), x,
                                    left.values()[ uVar ], right.values()[ uVar ],
                                    adaptIndicator );

        model_.jump( left.intersection(), x,
                     left.values()[ uVar ], right.values()[ uVar ],
                     shockIndicator );
        return 1.;
      }
    }

    //! returns difference between internal value and boundary
    //! value
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& jump,
                        JacobianRangeType& dummy) const
    {

      if( model_.hasBoundaryValue( left ) )
      {
        RangeType uRight;
        // evaluate boundary value
        model_.boundaryValue( left, left.values()[ uVar ], uRight );

        // use boundaryValue to check physical and jump
        if (! physical(left.entity(), left.position(), left.values()[ uVar ] ) ||
            ! physical(left.entity(), left.position(), uRight ) )
        {
          jump = 1e10;
          return -1.;
        }
        else
        {
          model_.jump( left.intersection(), left.localPosition(), left.values()[ uVar ], uRight, jump);
          return 1.;
        }
      }
      else
      {
        // otherwise evaluate boundary flux
        //model_.boundaryFlux( left, left.values()[ uVar ], left.jacobians()[ uVar ], jump );
        //std::cout << jump << " boundary flux" << std::endl;
        //return ( std::abs( jump[ 0 ] ) > 1e-10 ) ? -1. : 1.;
        jump = 0;
        return 0;
      }

    }

    //! return true if method boundaryValue returns something meaningful.
    template <class LocalEvaluation>
    bool hasBoundaryValue( const LocalEvaluation& local ) const
    {
      return model_.hasBoundaryValue( local );
    }

    //! returns difference between internal value and boundary
    //! value
    template <class LocalEvaluation>
    inline void boundaryValue(const LocalEvaluation& local,
                              const RangeType& uLeft,
                              RangeType& uRight) const
    {
      model_.boundaryValue( local, uLeft, uRight);
    }

    /** \brief returns true if model provides physical check */
    bool hasPhysical() const { return model_.hasPhysical(); }

    /** \brief check physical values */
    template <class ArgumentTuple>
    bool checkPhysical( const EntityType& entity,
                        const LocalDomainType& xLocal,
                        const ArgumentTuple& u ) const
    {
      return physical( entity, xLocal, u[ uVar ] );
    }

    /** \brief check physical values */
    bool physical( const EntityType& entity,
                   const LocalDomainType& xLocal,
                   const RangeType& u ) const
    {
      return model_.physical( entity, xLocal, u );
    }

    /** \brief adjust average values, e.g. transform to primitive or something similar */
    void adjustAverageValue( const EntityType& entity,
                             const LocalDomainType& xLocal,
                             RangeType& u ) const
    {
      model_.adjustAverageValue( entity, xLocal, u );
    }

    //! return reference to model
    const Model& model() const { return model_; }

    // g = grad L ( w_E,i - w_E ) ,  d = u_E,i - u_E
    // default limiter function is minmod
    const LimiterFunctionType& limiterFunction() const
    {
      return limiterFunction_;
    }

    //! returns true, if we have an inflow boundary
    template <class ArgumentTuple>
    bool checkDirection(const IntersectionType& it,
                        const FaceLocalDomainType& xFace,
                        const LocalDomainType& xLocal,
                        const ArgumentTuple& uLeft) const
    {
      // evaluate velocity
      model_.velocity(this->inside(), xLocal, uLeft[ uVar ], velocity_);
      return checkDirection(it, xFace, velocity_);
    }

  protected:
    //! returns true, if we have an inflow boundary
    bool checkDirection(const IntersectionType& it,
                        const FaceLocalDomainType& x,
                        const DomainType& velocity) const
    {
      // calculate scalar product of normal and velocity
      const double scalarProduct = it.outerNormal(x) * velocity;

      // if scalar product is larger than veloEps it's considered to be
      // an outflow intersection, otherwise an inflow intersection
      return (scalarProduct < veloEps_);
    }
  };

  template <class GlobalPassTraitsImp, class Model, int passId = -1 >
  using StandardLimiterDiscreteModel = LimiterDefaultDiscreteModel< GlobalPassTraitsImp, Model, passId >;

} // end namespace Fem
} // end namespace Dune
#endif
