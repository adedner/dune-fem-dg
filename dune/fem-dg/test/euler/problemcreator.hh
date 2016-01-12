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

//--------- HANDLER --------------------------------
#include <dune/fem-dg/algorithm/handler/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/sub/adapt.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
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

  template< class GridImp >
  struct EulerProblemCreator
  {

    struct SubEulerProblemCreator
    {
      typedef AlgorithmConfigurator< GridImp,
                                     Galerkin::Enum::dg,
                                     Adaptivity::Enum::yes,
                                     DiscreteFunctionSpaces::Enum::legendre,
                                     Solver::Enum::fem,
                                     AdvectionLimiter::Enum::limited,
                                     Matrix::Enum::matrixfree,
                                     AdvectionFlux::Identifier<AdvectionFlux::Enum::upwind>,
                                     PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::general > > AC;

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
        typedef Fem::EulerModel< GridPartType, InitialDataType >   ModelType;

        template< class Solution, class Model, class ExactFunction, class TimeProvider >
        static void addEOCErrors ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
        {
          static L2EOCError l2EocError( "$L^2$-Error");
          l2EocError.add( tp, u, model, f );
          static L1EOCError l1EocError( "$L^1$-Error");
          l1EocError.add( tp, u, model, f );
        }
      };


      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem()
      {
        return AnalyticalEulerProblemCreator<GridType>::apply();
      }

      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType >
                                                                                           DFSpaceType;
        typedef std::tuple<>                                                               ExtraParameterTuple;
      public:
        typedef typename AC::template DiscreteFunctions< DFSpaceType >                     DiscreteFunctionType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                 IOTupleType;

        class Operator
        {
          typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                           OpTraits;
        public:
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::full >    type;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::expl >    ExplicitType;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::impl >    ImplicitType;
        };

        struct Solver
        {
          // type of linear solver for implicit ode
          typedef typename AC::template LinearSolvers< DFSpaceType >                        BasicLinearSolverType;
          typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                       type;
          //typedef RungeKuttaSolver< typename Operator::type, typename Operator::ExplicitType, typename Operator::ImplicitType,
          //                          BasicLinearSolverType >                                 type;
        };

      private:
        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                            OpTraits;
        typedef DGAdaptationIndicatorOperator< OpTraits >                                   IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >   GradientIndicatorType ;
      public:

        typedef Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                 AdaptIndicatorType;
        typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                          SolverMonitorHandlerType;
        typedef Fem::SubDiagnosticsHandler< Diagnostics >                                   DiagnosticsHandlerType;
        typedef Fem::ExactSolutionOutputHandler< DiscreteFunctionType >                     AdditionalOutputHandlerType;
      };

      template <int polOrd>
      struct Stepper
      {
       // this should be ok but could lead to a henn-egg problem
        typedef Fem::AdvectionAlgorithm< GridType, SubEulerProblemCreator, polOrd > Type;
      };

    };

    template <int polOrd>
    struct Stepper
    {
      typedef Fem::EvolutionAlgorithm< polOrd, DefaultEvolutionCreator, SubEulerProblemCreator > Type;
    };

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return Fem::DefaultGridInitializer< GridType >::initialize(); }


  };

}
}
#endif
