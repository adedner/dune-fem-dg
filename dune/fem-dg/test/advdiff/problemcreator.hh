#ifndef FEMHOWTO_HEATSTEPPER_HH
#define FEMHOWTO_HEATSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
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
#include <dune/fem-dg/algorithm/monitor.hh>

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
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>

template< class GridImp >
struct AdvectionDiffusionProblemCreator
{

  struct SubAdvectionDiffusionProblemCreator
  {

    typedef GridImp                                         GridType;
    typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
    typedef HostGridPartType                                GridPartType;

    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

    // define problem type here if interface should be avoided
    typedef Dune::EvolutionProblemInterface< FunctionSpaceType,false >      ProblemInterfaceType;

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                    ProblemType;
      typedef ProblemInterfaceType                                    InitialDataType;
      typedef HeatEqnModel< GridPartType, InitialDataType >           ModelType;

      template< class Solution, class Model, class ExactFunction, class TimeProvider >
      static void addEOCErrors ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
      {
        static L2EOCError l2EocError( "$L^2$-Error");
        l2EocError.add( tp, u, model, f );
      }
    };

    static inline std::string moduleName() { return ""; }

    static ProblemInterfaceType* problem()
    {
      // choice of explicit or implicit ode solver
      static const std::string probString[]  = { "heat" ,"quasi", "pulse", "sin" };
      const int probNr = Dune::Fem::Parameter::getEnum( "problem", probString, 0 );
      if( probNr == 0 )
        return new Dune :: U0< GridType, DIMRANGE > ();
      else if ( probNr == 1 )
        return new Dune :: QuasiHeatEqnSolution< GridType, DIMRANGE > ();
      else if ( probNr == 2 )
        return new Dune :: Pulse< GridType, DIMRANGE > ();
      else if ( probNr == 3 )
        return new Dune :: U0Sin< GridType, DIMRANGE > ();
      else
      {
        abort();
        return 0;
      }
    }


    //Stepper Traits
    template< int polOrd >
    struct DiscreteTraits
    {
    private:
      static const SolverType solverType = fem ;
    public:

      //static const int polynomialOrder = polOrd;

      typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type             DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

      typedef typename InitialProjectors< typename AnalyticalTraits::ProblemType::TimeDependentFunctionType, DiscreteFunctionType, dg >::type   InitialProjectorType;

      typedef std::tuple<> ExtraParameterTuple;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                         IOTupleType;

      typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 OdeSolverType;
      // type of linear solver for implicit ode
      typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >            BasicLinearSolverType;

      class OperatorType
      {
        friend DiscreteTraits;
        typedef Dune::UpwindFlux< typename AnalyticalTraits::ModelType >                            FluxType;

        typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits, DiscreteFunctionType, FluxType, ExtraParameterTuple >
                                                                                                  OperatorTraitsType;

        // TODO: advection/diffusion should not be precribed by model
        static const int hasAdvection = AnalyticalTraits::ModelType::hasAdvection;
        static const int hasDiffusion = AnalyticalTraits::ModelType::hasDiffusion;
        typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;
      public:
        typedef typename AdvectionDiffusionOperatorType::FullOperatorType                         FullType;
        typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                     ImplicitType;
        typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                     ExplicitType;
      };

    private:
      typedef Dune::DGAdaptationIndicatorOperator< typename OperatorType::OperatorTraitsType, OperatorType::hasAdvection, OperatorType::hasDiffusion >
                                                                                                    IndicatorType;
      typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >             GradientIndicatorType ;
    public:

      typedef Dune::Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
      typedef Dune::Fem::SubSolverMonitorHandler< Dune::Fem::SolverMonitor< 1 > >                   SolverMonitorHandlerType;
      typedef Dune::Fem::SubDiagnosticsHandler< Dune::Diagnostics >                                 DiagnosticsHandlerType;
      typedef Dune::Fem::ExactSolutionOutputHandler< DiscreteFunctionType >                         AdditionalOutputHandlerType;
    };

    template <int polOrd>
    struct Stepper
    {
     // this should be ok but could lead to a henn-egg problem
      typedef Dune::Fem::AdvectionDiffusionStepper< GridType, SubAdvectionDiffusionProblemCreator, polOrd > Type;
    };

  };

  template <int polOrd>
  struct Stepper
  {
    typedef Dune::Fem::EvolutionAlgorithm< polOrd, SubAdvectionDiffusionProblemCreator > Type;
  };

  typedef GridImp                                         GridType;

  static inline std::string moduleName() { return ""; }

  static inline Dune::GridPtr<GridType>
  initializeGrid() { return Dune::Fem::DefaultGridInitializer< GridType >::initialize(); }


};

#endif // FEMHOWTO_HEATSTEPPER_HH
