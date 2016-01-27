#ifndef FEMDG_NAVIERSTOKESSTEPPER_HH
#define FEMDG_NAVIERSTOKESSTEPPER_HH
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
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
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
  struct NavierStokesProblemCreator
  {


    struct SubNavierStokesProblemCreator
    {
      typedef AlgorithmConfigurator< GridImp,
                                     Galerkin::Enum::dg,
                                     Adaptivity::Enum::yes,
                                     DiscreteFunctionSpaces::Enum::legendre,
                                     Solver::Enum::fem,
                                     AdvectionLimiter::Enum::unlimited,
                                     Matrix::Enum::matrixfree,
                                     AdvectionFlux::Enum::euler_llf2,
                                     PrimalDiffusionFlux::Enum::general > AC;

      typedef typename AC::GridType                                 GridType;
      typedef typename AC::GridParts                                HostGridPartType;
      typedef HostGridPartType                                      GridPartType;

      // define problem type here if interface should be avoided
      typedef NSWaves< GridType >                                   ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType      FunctionSpaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                ProblemType;
        typedef ProblemInterfaceType                                InitialDataType;
        typedef NSModel< GridPartType, InitialDataType >            ModelType;

        template< class Solution, class Model, class ExactFunction, class TimeProvider >
        static void addEOCErrors ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
        {
          static L2EOCError l2EocError( "$L^2$-Error");
          l2EocError.add( tp, u, model, f );
        }
      };

      static inline std::string moduleName() { return ""; }

      static inline GridPtr<GridType>
      initializeGrid() { return Fem::DefaultGridInitializer< GridType >::initialize(); }

      static ProblemInterfaceType* problem() { return new ProblemInterfaceType(); }


      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
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
          typedef typename AC::template LinearSolvers< DFSpaceType >                       BasicLinearSolverType;
          typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                      type;
        };

      private:
        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                           OpTraits;
        typedef DGAdaptationIndicatorOperator< OpTraits >                                  IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >  GradientIndicatorType ;
      public:

        typedef Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                AdaptIndicatorType;
        typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                         SolverMonitorHandlerType;
        typedef Fem::SubDiagnosticsHandler< Diagnostics >                                  DiagnosticsHandlerType;
        typedef Fem::ExactSolutionOutputHandler< DiscreteFunctionType >                    AdditionalOutputHandlerType;
      };

      template <int polOrd>
      struct Stepper
      {
        // this should be ok but could lead to a henn-egg problem
        typedef Fem::SubAdvectionDiffusionAlgorithm< GridType, SubNavierStokesProblemCreator, polOrd > Type;
      };


    };

    template <int polOrd>
    struct Stepper
    {
      typedef Fem::EvolutionAlgorithm< polOrd, DefaultEvolutionCreator, SubNavierStokesProblemCreator > Type;
    };

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return Fem::DefaultGridInitializer< GridType >::initialize(); }


  };

}
}

#endif // FEMHOWTO_HEATSTEPPER_HH
