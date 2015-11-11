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
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>

namespace Dune
{

/**
 *  \brief problem creator for an advection diffusion problem
 */
template< class GridImp >
struct AdvectionDiffusionProblemCreator
{

  struct SubAdvectionDiffusionProblemCreator
  {

    typedef GridImp                                       GridType;
    typedef Fem::DGAdaptiveLeafGridPart< GridType >       HostGridPartType;
    typedef HostGridPartType                              GridPartType;

    // define problem type here if interface should be avoided
    typedef EvolutionProblemInterface< Fem::FunctionSpace< double, double, GridImp::dimension, DIMRANGE > >
                                                                       ProblemInterfaceType;

    typedef typename ProblemInterfaceType::FunctionSpaceType           FunctionSpaceType;

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
      const int probNr = Fem::Parameter::getEnum( "problem", probString, 0 );
      if( probNr == 0 )
        return new U0< GridType, DIMRANGE > ();
      else if ( probNr == 1 )
        return new QuasiHeatEqnSolution< GridType, DIMRANGE > ();
      else if ( probNr == 2 )
        return new Pulse< GridType, DIMRANGE > ();
      else if ( probNr == 3 )
        return new U0Sin< GridType, DIMRANGE > ();
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
      typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type             DiscreteFunctionSpaceType;
    public:
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

      typedef std::tuple<> ExtraParameterTuple;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                         IOTupleType;

      class Operator
      {
        friend DiscreteTraits;
        typedef UpwindFlux< typename AnalyticalTraits::ModelType >                                FluxType;

        typedef DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits, DiscreteFunctionType, FluxType, ExtraParameterTuple >
                                                                                                  OperatorTraitsType;

        // TODO: advection/diffusion should not be precribed by model
        static const int hasAdvection = AnalyticalTraits::ModelType::hasAdvection;
        static const int hasDiffusion = AnalyticalTraits::ModelType::hasDiffusion;
        typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;
      public:
        typedef typename AdvectionDiffusionOperatorType::FullOperatorType                         type;
        typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                     ImplicitType;
        typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                     ExplicitType;
      };

      struct Solver
      {
        // type of linear solver for implicit ode
        typedef Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >                  BasicLinearSolverType;

        typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 type;
      };

    private:
      typedef DGAdaptationIndicatorOperator< typename Operator::OperatorTraitsType, Operator::hasAdvection, Operator::hasDiffusion >
                                                                                                    IndicatorType;
      typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >             GradientIndicatorType ;
    public:

      typedef Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
      typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                              SolverMonitorHandlerType;
      typedef Fem::SubDiagnosticsHandler< Diagnostics >                                       DiagnosticsHandlerType;
      typedef Fem::ExactSolutionOutputHandler< DiscreteFunctionType >                         AdditionalOutputHandlerType;
    };


    template <int polOrd>
    struct Stepper
    {
     // this should be ok but could lead to a henn-egg problem
      typedef Fem::AdvectionDiffusionAlgorithm< GridType, SubAdvectionDiffusionProblemCreator, polOrd > Type;
    };

  };

  template <int polOrd>
  struct Stepper
  {
    typedef Fem::EvolutionAlgorithm< polOrd, SubAdvectionDiffusionProblemCreator > Type;
  };

  typedef GridImp                                         GridType;

  static inline std::string moduleName() { return ""; }

  static inline GridPtr<GridType>
  initializeGrid() { return Fem::DefaultGridInitializer< GridType >::initialize(); }


};

}

#endif // FEMHOWTO_HEATSTEPPER_HH
