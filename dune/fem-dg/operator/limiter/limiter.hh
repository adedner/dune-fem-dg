#ifndef DUNE_FEM_LIMITER_HH
#define DUNE_FEM_LIMITER_HH

#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include <dune/fem-dg/operator/limiter/scalinglimitpass.hh>
#include <dune/fem-dg/operator/limiter/limiterdiscretemodel.hh>
#include <dune/fem-dg/operator/limiter/limitermodel.hh>
#include <dune/fem/function/adaptivefunction.hh>

namespace Dune
{
namespace Fem
{

  namespace detail
  {
    template <class DiscreteFunction>
    struct LimiterTraits
    {
      typedef DiscreteFunction     DiscreteFunctionType;
      typedef DiscreteFunctionType DestinationType;

      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType :: GridPartType  GridPartType ;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridType :: ctype ctype;
      enum { dimDomain = FunctionSpaceType :: dimDomain };

      typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
      typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;

      typedef typename FunctionSpaceType :: RangeFieldType FieldType;

      // Indicator for Limiter
      typedef FunctionSpace< ctype, FieldType, dimDomain, 3> FVFunctionSpaceType;
      typedef FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, SimpleStorage> IndicatorSpaceType;
      typedef AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;
    };
  }

  /**
   * \brief Limited reconstruction.
   *
   * \ingroup PassBased
   */
  template <class DomainFunction, class RangeFunction = DomainFunction,
            class Model = LimiterDefaultModel< typename DomainFunction::GridType,
                                               typename DomainFunction::FunctionSpaceType >
           >
  class Limiter
  {
    typedef DomainFunction DiscreteFunction ;

  public:
    typedef DomainFunction DomainFunctionType;
    typedef RangeFunction  RangeFunctionType;

    typedef typename DomainFunctionType :: DiscreteFunctionSpaceType  DomainSpaceType;
    typedef typename RangeFunctionType  :: DiscreteFunctionSpaceType  RangeSpaceType;

    typedef DiscreteFunction DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    typedef typename GridPartType :: GridType   GridType;

    typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;

    // pass ids
    enum PassIdType{ u = 0 , limitPass = 1 };
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

    typedef typename DiscreteFunctionSpaceType :: RangeFieldType ftype;

    typedef detail::LimiterTraits< DomainFunctionType >  PassTraitsType;
    //typedef LimiterDefaultModel< GridType, FunctionSpaceType > Model;
    typedef StandardLimiterDiscreteModel<PassTraitsType, Model, u > LimiterDiscreteModelType;

    // same type as DiscreteFunction
    typedef typename LimiterDiscreteModelType :: DestinationType DestinationType ;

    typedef StartPass< DiscreteFunctionType , u > StartPassType;
    typedef LimitDGPass< LimiterDiscreteModelType, StartPassType, limitPass > LimitPassType;

    typedef typename PassTraitsType :: IndicatorSpaceType  IndicatorSpaceType;
    typedef typename PassTraitsType :: IndicatorType       IndicatorType;

  public:
    Limiter( const DomainSpaceType& domainSpace,
             const double lowerBound, const double upperBound )
      : Limiter( domainSpace, domainSpace, lowerBound, upperBound ) {}

    Limiter( const DomainSpaceType& domainSpace,
             const RangeSpaceType&  rangeSpace,
             const double lowerBound, const double upperBound )
      : domainSpace_( domainSpace )
      , rangeSpace_( rangeSpace )
      , indicatorSpace_( domainSpace_.gridPart() )
      , indicator_("indicator", indicatorSpace_ )
      , model_( lowerBound, upperBound )
      , discreteModel_( model_, domainSpace_.order() )
      , startPass_()
      , limitPass_( discreteModel_ , startPass_, domainSpace_ )
    {
      discreteModel_.setIndicator( &indicator_ );
    }

    //! calculate internal reconstruction
    void operator () ( const DomainFunctionType& arg, RangeFunctionType& dest )
    {
      /*
      values_.resize( size );
      gradients_.resize( size );

      // helper class for evaluation of average value of discrete function
      EvalAverage average( *this, arg, discreteModel_);

      // obtain average values for all cells
      const auto endit = spc_.end();
      for( auto it = spc_.begin(); (it != endit); ++it )
      {
        const auto& en = *it;
        average.evaluate( en, values_[ gridPart_.indexSet().index( en ) ] );
      }
      */

      // apply limit pass in any case
      limitPass_.enable();

      // calculate reconstruction
      limitPass_( arg, dest );
    }

    bool activated( const EntityType& entity ) const
    {
      return std::abs( indicator_.localFunction( entity )[ 0 ] ) > 0.0;
    }

    const IndicatorType& indicator() const { return indicator_; }

  private:
    Limiter( const Limiter& );

  protected:
    const DomainSpaceType& domainSpace_;
    const RangeSpaceType&  rangeSpace_;

    IndicatorSpaceType indicatorSpace_;
    IndicatorType      indicator_;

    Model model_;
    LimiterDiscreteModelType discreteModel_;

    StartPassType startPass_;
    LimitPassType limitPass_;
  };




  /**
   * \brief Limited reconstruction.
   *
   * \ingroup PassBased
   */
  template <class DomainFunction, class RangeFunction = DomainFunction>
  class ScalingLimiter
  {
    typedef DomainFunction DiscreteFunction ;

    // MethodOrderTraits
    struct PassTraits
    {
      typedef DiscreteFunction     DiscreteFunctionType;
      typedef DiscreteFunctionType DestinationType;

      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;

      typedef typename DiscreteFunctionSpaceType :: GridPartType  GridPartType ;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridType :: ctype ctype;
      enum { dimDomain = FunctionSpaceType :: dimDomain };

      typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
      typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;

      typedef typename FunctionSpaceType :: RangeFieldType FieldType;

      // Indicator for Limiter
      typedef FunctionSpace< ctype, FieldType, dimDomain, 3> FVFunctionSpaceType;
      typedef FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, SimpleStorage> IndicatorSpaceType;
      typedef AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;
    };

  public:
    typedef DomainFunction DomainFunctionType;
    typedef RangeFunction  RangeFunctionType;

    typedef typename DomainFunctionType :: DiscreteFunctionSpaceType  DomainSpaceType;
    typedef typename RangeFunctionType  :: DiscreteFunctionSpaceType  RangeSpaceType;

    typedef DiscreteFunction DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    typedef typename GridPartType :: GridType   GridType;

    typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;

    // pass ids
    enum PassIdType{ u = 0 , limitPass = 1 };
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

    typedef typename DiscreteFunctionSpaceType :: RangeType      RangeType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType ftype;

    typedef PassTraits PassTraitsType;
    typedef LimiterDefaultModel< GridType, FunctionSpaceType > Model;
    typedef StandardLimiterDiscreteModel<PassTraitsType, Model, u > LimiterDiscreteModelType;

    // same type as DiscreteFunction
    typedef typename LimiterDiscreteModelType :: DestinationType DestinationType ;

    typedef StartPass< DiscreteFunctionType , u > StartPassType;
    typedef ScalingLimitDGPass< LimiterDiscreteModelType, StartPassType, limitPass > LimitPassType;

    typedef typename PassTraitsType :: IndicatorSpaceType  IndicatorSpaceType;
    typedef typename PassTraitsType :: IndicatorType       IndicatorType;

  public:
    ScalingLimiter( const DomainSpaceType& domainSpace,
                    const double lowerBound, const double upperBound )
      : ScalingLimiter( domainSpace, domainSpace, lowerBound, upperBound ) {}

    ScalingLimiter( const DomainSpaceType& domainSpace,
                    const RangeSpaceType& rangeSpace,
                    const double lowerBound, const double upperBound )
      : domainSpace_( domainSpace )
      , rangeSpace_( rangeSpace )
      , indicatorSpace_( domainSpace_.gridPart() )
      , indicator_("indicator", indicatorSpace_ )
      , model_( lowerBound, upperBound )
      , discreteModel_( model_, domainSpace_.order() )
      , startPass_()
      , limitPass_( discreteModel_ , startPass_, domainSpace_ )
    {
      discreteModel_.setIndicator( &indicator_ );
    }

    //! calculate internal reconstruction
    void operator () ( const DomainFunctionType& arg, RangeFunctionType& dest )
    {
      // apply limit pass in any case
      limitPass_.enable();

      // calculate reconstruction
      limitPass_( arg, dest );
    }

    bool activated( const EntityType& entity ) const
    {
      return std::abs( indicator_.localFunction( entity )[ 0 ] ) > 0.0;
    }

    const IndicatorType& indicator() const { return indicator_; }

  private:
    ScalingLimiter( const ScalingLimiter& );

  protected:
    const DomainSpaceType& domainSpace_;
    const RangeSpaceType&  rangeSpace_;

    IndicatorSpaceType indicatorSpace_;
    IndicatorType      indicator_;

    Model model_;
    LimiterDiscreteModelType discreteModel_;

    StartPassType startPass_;
    LimitPassType limitPass_;
  };

  /**
   * \brief Limited reconstruction.
   *
   * \ingroup PassBased
   */
  template <class Model, class DiscreteFunction>
  class LimitedReconstruction
  {
    // MethodOrderTraits
    template <class GridPartImp, class ftype, int dimRange,int polOrd>
    class PassTraits
    {
    public:
      typedef ftype        FieldType;
      typedef GridPartImp  GridPartType ;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridType :: ctype ctype;
      enum { dimDomain = GridType :: dimensionworld };

      typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
      typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;

      typedef FunctionSpace< ctype, FieldType, dimDomain, dimRange> FunctionSpaceType;
      typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType,
                 polOrd> DiscreteFunctionSpaceType;

      typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
      typedef DestinationType DiscreteFunctionType;

      // Indicator for Limiter
      typedef FunctionSpace< ctype, FieldType, dimDomain, 3> FVFunctionSpaceType;
      typedef FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, SimpleStorage> IndicatorSpaceType;
      typedef AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;
    };

  public:
    typedef DiscreteFunction  DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    typedef typename GridPartType :: GridType   GridType;

    typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;

    // pass ids
    enum PassIdType{ u = 0 , limitPass = 1 };
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

    typedef typename DiscreteFunctionSpaceType :: RangeFieldType ftype;

    typedef PassTraits<GridPartType, ftype, dimRange, 1> PassTraitsType;
    typedef StandardLimiterDiscreteModel<PassTraitsType, Model, u > LimiterDiscreteModelType;

    typedef typename LimiterDiscreteModelType :: DestinationType ReconstructionType;

    typedef typename ReconstructionType :: DiscreteFunctionSpaceType
      ReconstructionSpaceType;

    //! define local function, this could also be a temporary local function
    typedef typename ReconstructionType :: LocalFunctionType  LocalFunctionType;

    typedef StartPass< DiscreteFunctionType , u > StartPassType;
    typedef LimitDGPass< LimiterDiscreteModelType, StartPassType, limitPass > LimitPassType;

  public:
    LimitedReconstruction( const Model& model, const DiscreteFunctionSpaceType& space )
      : space_( space.gridPart() )
      , reconstruction_( "LimitedReconstruction" , space_ )
      , problem_( model, 1 )
      , startPass_()
      , limitPass_( problem_ , startPass_, space_ )
    {}

    //! calculate internal reconstruction
    void update( const DiscreteFunctionType& arg )
    {
      // apply limit pass in any case
      limitPass_.enable();
      // calculate reconstruction
      limitPass_( arg, reconstruction_ );
    }

    //! return local reconstruction
    LocalFunctionType localFunction( const EntityType& entity )
    {
      return reconstruction_.localFunction( entity );
    }

    //! return local reconstruction
    const LocalFunctionType localFunction( const EntityType& entity ) const
    {
      return reconstruction_.localFunction( entity );
    }

    ReconstructionType& function()
    {
      return reconstruction_ ;
    }

    const ReconstructionType& function() const
    {
      return reconstruction_ ;
    }

  private:
    LimitedReconstruction( const LimitedReconstruction& );

  protected:
    ReconstructionSpaceType space_;
    ReconstructionType reconstruction_;

    LimiterDiscreteModelType problem_;

    StartPassType startPass_;
    LimitPassType limitPass_;

  };

} // end namespace Fem

} // end namespace Dune

#endif // DUNE_FEM_LIMITER_HH
