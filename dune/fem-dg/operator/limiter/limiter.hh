#ifndef DUNE_FEM_LIMITER_HH
#define DUNE_FEM_LIMITER_HH

#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include <dune/fem-dg/operator/limiter/limiterdiscretemodel.hh>
#include <dune/fem-dg/operator/limiter/limitermodel.hh>
#include <dune/fem/function/adaptivefunction.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Limited reconstruction.
   *
   * \ingroup PassBased
   */
  template <class DiscreteFunction>
  class Limiter
  {
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

    typedef PassTraits PassTraitsType;
    typedef LimiterDefaultModel< GridType, FunctionSpaceType > Model;
    typedef StandardLimiterDiscreteModel<PassTraitsType, Model, u > LimiterDiscreteModelType;

    // same type as DiscreteFunction
    typedef typename LimiterDiscreteModelType :: DestinationType DestinationType ;

    typedef StartPass< DiscreteFunctionType , u > StartPassType;
    typedef LimitDGPass< LimiterDiscreteModelType, StartPassType, limitPass > LimitPassType;

  public:
    Limiter( const DiscreteFunctionSpaceType& space )
      : space_( space )
      , model_( )
      , discreteModel_( model_, space_.order() )
      , startPass_()
      , limitPass_( discreteModel_ , startPass_, space_ )
    {}

    //! calculate internal reconstruction
    void operator () ( const DiscreteFunctionType& arg, DiscreteFunctionType& dest )
    {
      // apply limit pass in any case
      limitPass_.enable();
      // calculate reconstruction
      limitPass_( arg, dest );
    }

  private:
    Limiter( const Limiter& );

  protected:
    const DiscreteFunctionSpaceType& space_;
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
