#ifndef DUNE_FEM_DG_OPERATORTRAITS_HH
#define DUNE_FEM_DG_OPERATORTRAITS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/fluxes.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem/space/finitevolume/space.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>

#include <dune/fem-dg/operator/fluxes/diffusion/parameters.hh>
namespace Dune
{
namespace Fem
{

  // traits for the operator passes
  template< class GridPart,
            int polOrd,
            class AnalyticalTraits,
            class DiscreteFunctionImp,
            class AdvectionFluxImp,
            class DiffusionFluxImp,
            class LimiterIndicatorFunctionImp,
            class AdaptationHandlerImp,
            class ExtraParameterTupleImp = std::tuple<>
          >
  struct OperatorTraits
  {
    typedef GridPart                                                     GridPartType;
    typedef typename GridPartType::GridType                              GridType;
    typedef AnalyticalTraits                                             AnalyticalTraitsType;

    typedef typename AnalyticalTraitsType::InitialDataType               InitialDataType;
    typedef typename AnalyticalTraitsType::ModelType                     ModelType ;
    typedef AdvectionFluxImp                                             AdvectionFluxType;
    typedef DiffusionFluxImp                                             DiffusionFluxType;

    static const int polynomialOrder = polOrd == -1 ? 0 : polOrd;

    typedef DiscreteFunctionImp                                          DestinationType ;
    typedef typename DestinationType::DiscreteFunctionSpaceType          DiscreteFunctionSpaceType;

    typedef Fem::CachingQuadrature< GridPartType, 0 >                    VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 >                    FaceQuadratureType;

    typedef LimiterIndicatorFunctionImp                                  LimiterIndicatorType;

    typedef AdaptationHandlerImp                                         AdaptationHandlerType ;

    static const int limiterPolynomialOrder = polOrd == -1 ? 1 : polOrd;
    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };

  // traits for the operator passes
  template< int polOrd,
            class AnalyticalTraits,
            class DiscreteFunctionImp,
            class AdvectionFluxImp,
            class DiffusionFluxImp,
            class ExtraParameterTupleImp = std::tuple<>,
            class AdaptationHandlerFunctionSpaceImp = typename DiscreteFunctionImp::DiscreteFunctionSpaceType::FunctionSpaceType
          >
  struct DefaultOperatorTraits
  {
    typedef typename DiscreteFunctionImp::GridPartType                   GridPartType;
    typedef typename GridPartType::GridType                              GridType;
    typedef AnalyticalTraits                                             AnalyticalTraitsType;

    typedef typename AnalyticalTraitsType::InitialDataType               InitialDataType;
    typedef typename AnalyticalTraitsType::ModelType                     ModelType ;
    typedef AdvectionFluxImp                                             AdvectionFluxType;
    typedef DiffusionFluxImp                                             DiffusionFluxType;

    static const int polynomialOrder = polOrd == -1 ? 0 : polOrd;

    typedef DiscreteFunctionImp                                          DestinationType ;
    typedef typename DestinationType::DiscreteFunctionSpaceType          DiscreteFunctionSpaceType;

    typedef Fem::CachingQuadrature< GridPartType, 0 >                    VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 >                    FaceQuadratureType;

  private:
    typedef Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraits::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Fem::SimpleStorage> IndicatorSpaceType;
  public:
    typedef Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>            LimiterIndicatorType;

    typedef AdaptationHandler< GridType, AdaptationHandlerFunctionSpaceImp >
                                                                         AdaptationHandlerType;

    static const int limiterPolynomialOrder = polOrd == -1 ? 1 : polOrd;
    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };

} // end namespace
} // end namespace
#endif
