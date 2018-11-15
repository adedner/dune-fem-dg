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
            class ModelImp,
            class DiscreteFunctionImp,
            class AdvectionFluxImp,
            class DiffusionFluxImp,
            class LimiterIndicatorFunctionImp,
            class AdaptationHandlerImp,
            class ExtraParameterTupleImp = std::tuple<>,
            template <class F, int d> class QuadratureTraits = Dune::Fem::DefaultQuadratureTraits
          >
  struct OperatorTraits
  {
    typedef GridPart                                                     GridPartType;
    typedef typename GridPartType::GridType                              GridType;

    typedef ModelImp                                                     ModelType;
    typedef AdvectionFluxImp                                             AdvectionFluxType;
    typedef DiffusionFluxImp                                             DiffusionFluxType;

    // polynomial order of ansatz space
    static const int polynomialOrder = polOrd;

    // if polynomialOrder is 0 then limiting will be linear (higher order FV case)
    static const int limiterPolynomialOrder = (polynomialOrder == 0) ? 1 : polynomialOrder;

    typedef DiscreteFunctionImp                                          DestinationType ;
    //static_assert( std::is_same<typename  ModelType::RangeType, typename DiscreteFunctionType::RangeType>::value, "range type does not fit.");
    typedef typename DestinationType::DiscreteFunctionSpaceType          DiscreteFunctionSpaceType;

    typedef Fem::CachingQuadrature< GridPartType, 0, QuadratureTraits >  VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1, QuadratureTraits >  FaceQuadratureType;

    typedef LimiterIndicatorFunctionImp                                  LimiterIndicatorType;

    typedef AdaptationHandlerImp                                         AdaptationHandlerType ;

    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };

  // traits for the operator passes
  template< class ModelImp,
            class DiscreteFunctionImp,
            class AdvectionFluxImp,
            class DiffusionFluxImp,
            class ExtraParameterTupleImp = std::tuple<>,
            class AdaptationHandlerFunctionSpaceImp = typename DiscreteFunctionImp::DiscreteFunctionSpaceType::FunctionSpaceType,
            template <class F, int d> class QuadratureTraits = Dune::Fem::DefaultQuadratureTraits,
            bool enableThreaded =
    // static cmake variables provided by dune-fem
#ifdef USE_SMP_PARALLEL
              true
#else
              false
#endif
          >
  struct DefaultOperatorTraits
  {
    typedef typename DiscreteFunctionImp::GridPartType                   GridPartType;
    typedef typename GridPartType::GridType                              GridType;

    typedef ModelImp                                                     ModelType;
    typedef AdvectionFluxImp                                             AdvectionFluxType;
    typedef DiffusionFluxImp                                             DiffusionFluxType;

    typedef DiscreteFunctionImp                                          DestinationType;
    typedef typename DestinationType::DiscreteFunctionSpaceType          DiscreteFunctionSpaceType;

    // polynomial order of ansatz space
    static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;

    // if polynomialOrder is 0 then limiting will be linear (higher order FV case)
    static const int limiterPolynomialOrder = (polynomialOrder == 0) ? 1 : polynomialOrder;

    // enables the possibility to run in threaded mode
    static const bool threading = enableThreaded ;

    static_assert( std::is_same<typename ModelType::RangeType, typename DestinationType::RangeType>::value, "range type does not fit.");

    typedef Fem::CachingQuadrature< GridPartType, 0, QuadratureTraits >  VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1, QuadratureTraits >  FaceQuadratureType;

  private:
    typedef Fem::FunctionSpace< typename GridType::ctype, double, ModelImp::dimDomain, 3> FVFunctionSpaceType;
    typedef Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Fem::SimpleStorage> IndicatorSpaceType;
  public:
    typedef Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>            LimiterIndicatorType;

    typedef AdaptationHandler< GridType, AdaptationHandlerFunctionSpaceImp >
                                                                         AdaptationHandlerType;

    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };

} // end namespace
} // end namespace
#endif
