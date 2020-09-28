#ifndef FEMDG_POISSONSTEPPER_HH
#define FEMDG_POISSONSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
#include "gridinitializer.hh"
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
#include <dune/fem-dg/operator/adaptation/poissonestimator.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/elliptic.hh>
#include <dune/fem-dg/algorithm/steadystate.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/configurator.hh>

namespace Dune
{
namespace Fem
{

  template< class GridSelectorGridType >
  struct PoissonAlgorithmCreator
  {
    typedef AlgorithmConfigurator< GridSelectorGridType,
                                   Galerkin::Enum::dg,
                                   Adaptivity::Enum::yes,
                                   DiscreteFunctionSpaces::Enum::orthonormal,
                                   Solver::Enum::istl,
                                   AdvectionLimiter::Enum::unlimited,
                                   Matrix::Enum::assembled,
                                   AdvectionFlux::Enum::upwind,
                                   DiffusionFlux::Enum::primal > ACPoisson;

    template< class AC >
    struct SubPoissonAlgorithmCreator
    {


      typedef typename AC::GridType                         GridType;
      typedef typename AC::GridParts                        HostGridPartType;
      typedef HostGridPartType                              GridPartType;

      // define problem type here if interface should be avoided
      typedef ProblemInterface< typename AC::template FunctionSpaces<DIMRANGE> >
                                                                    ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType      FunctionSpaceType;

      typedef Poisson::Model< GridType, ProblemInterfaceType >      ModelType;

      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem()
      {
        int probNr = Parameter::getValue< int > ( "problem" );
        return new Poisson::Problem< GridType, DIMRANGE > ( probNr );
      }

      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctions< FunctionSpaceType, polOrd > DiscreteFunctionType;

        class Operator
        {
          typedef typename AC::template DefaultAssembTraits< ModelType, FunctionSpaceType, polOrd >
                                                                                     OpTraits;
        public:
          typedef typename AC::template Operators< OpTraits >                        AssemblerType;
          typedef typename AssemblerType::LinearOperatorType                         type;
        };

        struct Solver
        {
          typedef typename AC::template LinearSolvers< DiscreteFunctionType, false > type;
        };

      private:
        typedef typename AC::template WrappedDiscreteFunctions< SigmaDiscreteFunctionSelector, DiscreteFunctionType >
                                                                                     SigmaDiscreteFunctionType;

        typedef ErrorEstimator< DiscreteFunctionType, SigmaDiscreteFunctionType, typename Operator::AssemblerType >
                                                                                     ErrorEstimatorType;
        typedef PoissonSigmaEstimator< ErrorEstimatorType >                          SigmaEstimatorType;

        typedef PAdaptivity<typename DiscreteFunctionType::DiscreteFunctionSpaceType, polOrd, SigmaEstimatorType >
                                                                                     PAdaptivityType;
      public:

        typedef SubSolverMonitor< SolverMonitor >                                    SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                        DiagnosticsType;
        typedef PAdaptIndicator< PAdaptivityType, ModelType >                        AdaptIndicatorType;
        typedef SubDataWriter< SolutionOutput<DiscreteFunctionType>, ExactSolutionOutput<DiscreteFunctionType> >
                                                                                     DataWriterType;
      };

      template <int polOrd>
      using Algorithm = SubEllipticAlgorithm< GridType, SubPoissonAlgorithmCreator<AC>, polOrd >;

    };

    template <int polOrd>
    using Algorithm = SteadyStateAlgorithm< polOrd, SubPoissonAlgorithmCreator<ACPoisson> >;

    static inline std::string moduleName() { return ""; }


    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef typename SubPoissonAlgorithmCreator<ACPoisson>::GridType  GridType;
      //Discrete Functions
      typedef typename SubPoissonAlgorithmCreator<ACPoisson>::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType;

      //Item1
      typedef _t< SubSteadyStateContainerItem >                         Steady;
      typedef std::tuple< Steady >                                      Item1TupleType;

      //Item2
      typedef _t<SubEllipticContainerItem, ACPoisson::template Containers > Def;
      //typedef _t<SubEllipticContainerItem, SparseRowLinearOperator >  Def;
      //typedef _t< SubEllipticContainerItem, SparseRowLinearOperator >   Sp;

      typedef std::tuple< std::tuple< Def > >                           Item2TupleType;

      //Sub (discrete function argument ordering)
      typedef std::tuple<__0 >                                          PoissonOrder;

      typedef std::tuple< PoissonOrder >                                SubOrderRowType;
      typedef SubOrderRowType                                           SubOrderColType;

      //external params lists
      typedef ExtraArg<>                                                ExtraType;

      //Global container
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, ExtraType, DFType > GlobalContainerType;

      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( moduleName(), gridptr );
    }

  };

}
}

#endif // FEMHOWTO_HEATSTEPPER_HH
