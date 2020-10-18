#ifndef DUNE_FEM_DG_OPERATORTRAITS_HH
#define DUNE_FEM_DG_OPERATORTRAITS_HH

#include <dune/fem/quadrature/dunequadratures.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/interpolationquadrature.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/fluxes.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem/space/finitevolume/space.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/space/common/capabilities.hh>

#include <dune/fem-dg/operator/fluxes/diffusion/parameters.hh>

#if HAVE_DUNE_POLYGONGRID
#include <dune/polygongrid/declaration.hh>
#endif

#if HAVE_OPM_GRID
#include <opm/grid/polyhedralgrid/declaration.hh>
#endif // #if HAVE_OPM_GRID

namespace Dune
{
namespace Fem
{
  // traits for the operator passes
  template< class ModelImp,
            class DiscreteFunctionImp,
            class AdvectionFluxImp,
            class DiffusionFluxImp,
            class ExtraParameterTupleImp = std::tuple<>,
            class AdaptationHandlerFunctionSpaceImp = typename DiscreteFunctionImp::DiscreteFunctionSpaceType::FunctionSpaceType,
            template <class F, int d> class QuadratureTraits = Capabilities::DefaultQuadrature< typename DiscreteFunctionImp::DiscreteFunctionSpaceType >::template DefaultQuadratureTraits,
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

    // enables the possibility to run in threaded mode
    static const bool threading = enableThreaded ;

    static const bool scalingLimiter = ModelType::scalingLimiter;

    static_assert( std::is_same<typename ModelType::RangeType, typename DestinationType::RangeType>::value, "range type does not fit.");

    // default quadrature selection should be CachingQuadrature
    template <class Grid, bool caching>
    struct SelectQuadrature
    {
      typedef Fem::CachingQuadrature< GridPartType, 0, QuadratureTraits >  VolumeQuadratureType;
      typedef Fem::CachingQuadrature< GridPartType, 1, QuadratureTraits >  FaceQuadratureType;
    };

    // if selected discrete function space has no caching storage,
    // then select ElementQuadrature, and not  CachingQuadrature
    // also, if the grid is polygonal or polyhedral
    // CachingQuadarature cannot be used
    template <class Grid>
    struct SelectQuadrature< Grid, false >
    {
      typedef Fem::ElementQuadrature< GridPartType, 0, QuadratureTraits >  VolumeQuadratureType;
      typedef Fem::ElementQuadrature< GridPartType, 1, QuadratureTraits >  FaceQuadratureType;
    };

    template <class Grid>
    struct CheckGrid { static const bool value = true; };

    // also, if the grid is polygonal or polyhedral
    // CachingQuadarature cannot be used
#if HAVE_DUNE_POLYGONGRID
    template <class ct>
    struct CheckGrid< Dune::PolygonGrid< ct > > { static const bool value = false; };
#endif

#if HAVE_OPM_GRID
    template <int dim, int dimworld, class ct>
    struct CheckGrid< Dune::PolyhedralGrid< dim, dimworld, ct > > { static const bool value = false; };
#endif

    template <class DFS>
    struct CheckSpace { static const bool value = true; };

    // if selected discrete function space has no caching storage,
    // then select ElementQuadrature, and not  CachingQuadrature
    template <class FS, class GP, int polOrd,
              template <class, class, int, class> class DFS>
    struct CheckSpace< DFS< FS, GP, polOrd, Dune::Fem::SimpleStorage > >
    {
      static const bool value = false;
    };

    // define quad selector which contains the correct types for quadratures
    typedef SelectQuadrature< GridType,
                              CheckSpace< DiscreteFunctionSpaceType > :: value  &&
                              CheckGrid < GridType > :: value >  SelectQuadratureType;

    typedef typename SelectQuadratureType::VolumeQuadratureType   VolumeQuadratureType;
    typedef typename SelectQuadratureType::FaceQuadratureType     FaceQuadratureType;

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
