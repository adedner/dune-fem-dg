#ifndef FEMDG_EULERPROBLEMCREATOR_HH
#define FEMDG_EULERPROBLEMCREATOR_HH
#include <config.h>

#ifndef GRIDDIM
#define GRIDDIM 2
#endif

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM + 2
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/l1eocerror.hh>
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
  struct EulerAlgorithmCreator
  {
    typedef AlgorithmConfigurator< GridSelectorGridType,
                                   Galerkin::Enum::dg,
                                   Adaptivity::Enum::yes,
                                   DiscreteFunctionSpaces::Enum::orthonormal,
                                   Solver::Enum::fem,
                                   AdvectionLimiter::Enum::limited,
                                   Matrix::Enum::matrixfree,
                                   AdvectionFlux::Enum::euler_hllc,
                                   PrimalDiffusionFlux::Enum::none > ACEuler;

    template< class AC >
    struct SubEulerAlgorithmCreator
    {
      // define problem type here if interface should be avoided
      typedef typename AC::GridType                                GridType;
      typedef typename AC::GridParts                               HostGridPartType;
      typedef typename AC::GridParts                               GridPartType;

      typedef ProblemBase< GridType >                              ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType     FunctionSpaceType;


      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                               ProblemType;
        typedef ProblemInterfaceType                               InitialDataType;
        typedef EulerModel< GridType, InitialDataType >            ModelType;

        template< class Solution, class Model, class ExactFunction, class TimeProvider >
        static void addEOCErrors ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
        {
          static L1EOCError l1EocError( "$L^1$-Error");
          l1EocError.add( tp, u, model, f );
          static L2EOCError l2EocError( "$L^2$-Error");
          l2EocError.add( tp, u, model, f );
        }
      };


      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem()
      {
        return AnalyticalEulerProblemCreator<GridType>::apply();
      }

      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctions< FunctionSpaceType, polOrd >       DiscreteFunctionType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                 IOTupleType;

        class Operator
        {
          typedef typename AC::template DefaultOpTraits< AnalyticalTraits, FunctionSpaceType, polOrd >
                                                                                           OpTraits;
        public:
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::full >    type;
        };

        struct Solver
        {
          // type of linear solver for implicit ode
          typedef typename AC::template LinearSolvers< DiscreteFunctionType >               LinearSolverType;
          typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                       type;
          //typedef RungeKuttaSolver< typename Operator::type, typename Operator::ExplicitType, typename Operator::ImplicitType,
          //                          LinearSolverType >                                    type;
        };

      private:
        typedef typename AC::template DefaultOpTraits< AnalyticalTraits, FunctionSpaceType, polOrd >
                                                                                            OpTraits;
        typedef DGAdaptationIndicatorOperator< OpTraits >                                   IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >   GradientIndicatorType ;
      public:

        typedef AdaptIndicator< IndicatorType, GradientIndicatorType >           AdaptIndicatorType;
        typedef SubSolverMonitor< SolverMonitor >                                SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                    DiagnosticsType;
        typedef ExactSolutionOutput< DiscreteFunctionType >                      AdditionalOutputType;
      };

      template <int polOrd>
      using Algorithm = SubAdvectionAlgorithm< GridType, SubEulerAlgorithmCreator<AC>, polOrd >;
    };

    template <int polOrd>
    using Algorithm = EvolutionAlgorithm< polOrd, SubEulerAlgorithmCreator<ACEuler> >;

    static inline std::string moduleName() { return ""; }

    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef typename SubEulerAlgorithmCreator<ACEuler>::GridType    GridType;
      //Discrete Functions
      typedef typename SubEulerAlgorithmCreator<ACEuler>::template DiscreteTraits<polOrd>::DiscreteFunctionType
                                                                      DFType;

      //Item1
      typedef _t< SubEvolutionContainerItem >                         Steady;
      typedef std::tuple< Steady >                                    Item1TupleType;

      //Item2
      typedef _t< EmptyContainerItem >                                Empty;
      typedef std::tuple< std::tuple< Empty > >                       Item2TupleType;


      //Sub (discrete function argument ordering)
      typedef std::tuple<__0 >                                        AdvDiffOrder;

      typedef std::tuple< AdvDiffOrder >                              SubOrderRowType;
      typedef SubOrderRowType                                         SubOrderColType;

      //external params lists
      typedef ExtraArg<>                                              ExtraType;

      //Global container
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, ExtraType, DFType >
                                                                      GlobalContainerType;

      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( moduleName(), gridptr );
    }
  };

}
}
#endif
