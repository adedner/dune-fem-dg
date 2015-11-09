#ifndef DUNE_FEM_DG_OPERATORTRAITS_HH
#define DUNE_FEM_DG_OPERATORTRAITS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem/space/finitevolume/space.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>

namespace Dune {

  // traits for the operator passes
  template< class GridPart,
            int polOrd,
            class AnalyticalTraits,
            class DiscreteFunctionImp,
            class AdvectionFluxImp,
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
    typedef AdvectionFluxImp                                             FluxType;

    static const Dune::DGDiffusionFluxIdentifier PrimalDiffusionFluxId = Dune::method_general;

    static const int polynomialOrder = polOrd == -1 ? 0 : polOrd;

    typedef DiscreteFunctionImp DestinationType ;
    typedef typename DestinationType::DiscreteFunctionSpaceType          DiscreteFunctionSpaceType;

    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >              VolumeQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 1 >              FaceQuadratureType;

    typedef LimiterIndicatorFunctionImp                                  LimiterIndicatorType;

    typedef AdaptationHandlerImp                                         AdaptationHandlerType ;

    static const int limiterPolynomialOrder = polOrd == -1 ? 1 : polOrd;
    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };

  // traits for the operator passes
  template< class GridPart,
            int polOrd,
            class AnalyticalTraits,
            class DiscreteFunctionImp,
            class AdvectionFluxImp,
            class ExtraParameterTupleImp = std::tuple<>,
            class AdaptationHandlerFunctionSpaceImp = typename DiscreteFunctionImp::DiscreteFunctionSpaceType::FunctionSpaceType
          >
  struct DefaultOperatorTraits
  {
    typedef GridPart                                                     GridPartType;
    typedef typename GridPartType::GridType                              GridType;
    typedef AnalyticalTraits                                             AnalyticalTraitsType;

    typedef typename AnalyticalTraitsType::InitialDataType               InitialDataType;
    typedef typename AnalyticalTraitsType::ModelType                     ModelType ;
    typedef AdvectionFluxImp                                             FluxType;

    static const Dune::DGDiffusionFluxIdentifier PrimalDiffusionFluxId = Dune::method_general;

    static const int polynomialOrder = polOrd == -1 ? 0 : polOrd;

    typedef DiscreteFunctionImp DestinationType ;
    typedef typename DestinationType::DiscreteFunctionSpaceType          DiscreteFunctionSpaceType;

    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >              VolumeQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 1 >              FaceQuadratureType;

  private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraits::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
  public:
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>      LimiterIndicatorType;

    typedef Dune::AdaptationHandler< GridType, AdaptationHandlerFunctionSpaceImp >
                                                                         AdaptationHandlerType;

    static const int limiterPolynomialOrder = polOrd == -1 ? 1 : polOrd;
    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };
} // end namespace Dune
#endif
