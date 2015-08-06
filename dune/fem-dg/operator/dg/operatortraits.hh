#ifndef DUNE_FEM_DG_OPERATORTRAITS_HH
#define DUNE_FEM_DG_OPERATORTRAITS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>

namespace Dune {

  // traits for the operator passes
  template< class GridPart,
            int polOrd,
            class AnalyticalTraits,
            class DiscreteFunctionImp,
            class IndicatorFunctionImp,
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
    typedef UpwindFlux< ModelType >                                      FluxType;
    //typedef LLFFlux< ModelType >                    FluxType;
    static const Dune::DGDiffusionFluxIdentifier PrimalDiffusionFluxId = Dune::method_general;

    static const int polynomialOrder = polOrd == -1 ? 0 : polOrd;

    typedef DiscreteFunctionImp DestinationType ;
    typedef typename DestinationType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >              VolumeQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 1 >              FaceQuadratureType;

    typedef IndicatorFunctionImp                                         IndicatorType;

    typedef AdaptationHandlerImp                                         AdaptationHandlerType ;

    static const int limiterPolynomialOrder = polOrd == -1 ? 1 : polOrd;
    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };

} // end namespace Dune
#endif