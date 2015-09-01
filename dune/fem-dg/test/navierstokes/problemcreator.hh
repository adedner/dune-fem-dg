#ifndef FEMDG_NAVIERSTOKESSTEPPER_HH
#define FEMDG_NAVIERSTOKESSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM + 2
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

//--------- HANDLER --------------------------------
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/advectiondiffusionstepper.hh>
#include <dune/fem-dg/algorithm/advectionstepper.hh>
#include <dune/fem-dg/algorithm/combinedevolution.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/l1eocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>

template< class GridImp >
struct NavierStokesProblemCreator
{


  struct SubNavierStokesProblemCreator
  {

    typedef GridImp                                         GridType;
    typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
    typedef HostGridPartType                                GridPartType;

    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

    // define problem type here if interface should be avoided
    typedef Dune::NSWaves< GridType >                             ProblemInterfaceType;

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

    static inline Dune::GridPtr<GridType>
    initializeGrid() { return Dune::Fem::DefaultGridInitializer< GridType >::initialize(); }

    static ProblemInterfaceType* problem() { return new ProblemInterfaceType(); }


    //Stepper Traits
    template< int polOrd >
    struct DiscreteTraits
    {
    private:
      static const SolverType solverType = fem;
    public:
      typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type    DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

      typedef typename InitialProjectors< typename AnalyticalTraits::ProblemType::TimeDependentFunctionType, DiscreteFunctionType, dg >::type   InitialProjectorType;

      typedef std::tuple<> ExtraParameterTuple;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

      typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 OdeSolverType;
      // type of linear solver for implicit ode
      typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >            BasicLinearSolverType;

      class Operator
      {
        friend DiscreteTraits;
        typedef LLFFlux< typename AnalyticalTraits::ModelType >                                    FluxType;

        typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                             DiscreteFunctionType, FluxType, ExtraParameterTuple > OperatorTraitsType;

        static const int hasAdvection = AnalyticalTraits::ModelType::hasAdvection;
        static const int hasDiffusion = AnalyticalTraits::ModelType::hasDiffusion;
        typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;
      public:
        typedef typename AdvectionDiffusionOperatorType::FullOperatorType                           type;
        typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                       ImplicitType;
        typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                       ExplicitType;
      };

      struct Solver
      {
        // type of linear solver for implicit ode
        typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >            BasicLinearSolverType;

        typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 type;
      };

    private:
      typedef Dune::DGAdaptationIndicatorOperator< typename Operator::OperatorTraitsType, Operator::hasAdvection, Operator::hasDiffusion >
                                                                                                    IndicatorType;
      typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >             GradientIndicatorType ;
    public:

      typedef Dune::Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
      typedef Dune::Fem::SubSolverMonitorHandler< Dune::Fem::SolverMonitor >                        SolverMonitorHandlerType;
      typedef Dune::Fem::SubDiagnosticsHandler< Dune::Diagnostics >                                 DiagnosticsHandlerType;
      typedef Dune::Fem::ExactSolutionOutputHandler< DiscreteFunctionType >                         AdditionalOutputHandlerType;
    };

    template <int polOrd>
    struct Stepper
    {
      // this should be ok but could lead to a henn-egg problem
      typedef Dune::Fem::AdvectionDiffusionStepper< GridType, SubNavierStokesProblemCreator, polOrd > Type;
    };


  };

  template <int polOrd>
  struct Stepper
  {
    typedef Dune::Fem::EvolutionAlgorithm< polOrd, SubNavierStokesProblemCreator > Type;
  };

  typedef GridImp                                         GridType;

  static inline std::string moduleName() { return ""; }

  static inline Dune::GridPtr<GridType>
  initializeGrid() { return Dune::Fem::DefaultGridInitializer< GridType >::initialize(); }


};

#endif // FEMHOWTO_HEATSTEPPER_HH
