#ifndef FEMDG_NAVIERSTOKESSTEPPER_HH
#define FEMDG_NAVIERSTOKESSTEPPER_HH
#include <config.h>

#ifndef POLORDER
#define POLORDER 1
#endif

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
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/l1eocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>
//--------- HANDLER -------------------------
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>


template< class GridImp >
struct NavierStokesProblemCreator
{
  typedef GridImp                                         GridType;
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
  typedef HostGridPartType                                GridPartType;

  // define problem type here if interface should be avoided
  typedef Dune::NSWaves< GridType >                             ProblemInterfaceType;

  typedef typename ProblemInterfaceType :: FunctionSpaceType    FunctionSpaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef NSModel< GridPart, InitialDataType >                ModelType;

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
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
  private:
    static const SolverType solverType = fem;
  public:
    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static const int quadOrder = polynomialOrder * 3 + 1;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef typename InitialProjectors< typename AnalyticalTraitsType::ProblemType::TimeDependentFunctionType, DiscreteFunctionType, dg >::type   InitialProjectorType;

    typedef std::tuple<> ExtraParameterTuple;

  private:
    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType                       ExactSolutionType;
  public:
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >                   GridExactSolutionType;
    typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

    typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 OdeSolverType;
    // type of linear solver for implicit ode
    typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >            BasicLinearSolverType;

    class HandlerTraits;

    class OperatorType
    {
      friend HandlerTraits;
      typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;

      typedef LLFFlux< typename AnalyticalTraitsType::ModelType >                                 FluxType;

      typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
      typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
      typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

      typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                    DiscreteFunctionType, FluxType, IndicatorType,
                                    AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

      static const int hasAdvection = AnalyticalTraitsType::ModelType::hasAdvection;
      static const int hasDiffusion = AnalyticalTraitsType::ModelType::hasDiffusion;
      typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;
    public:
      typedef typename AdvectionDiffusionOperatorType::FullOperatorType                           FullType;
      typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                       ImplicitType;
      typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                       ExplicitType;
    };

    //------HANDLER-----------------------------------------------------
    class HandlerTraits
    {
      //adaptivity
      typedef Dune::DGAdaptationIndicatorOperator< typename OperatorType::OperatorTraitsType,
                                                   OperatorType::hasAdvection,
                                                   OperatorType::hasDiffusion >                   IndicatorType;
      typedef Estimator< DiscreteFunctionType, typename  AnalyticalTraitsType::ProblemType >      GradientIndicatorType ;
      typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                           RestrictionProlongationType;
      //limiting
      typedef typename OperatorType::FullType                                                     LimiterOperatorType;
      public:
      typedef Dune::Fem::DefaultDiagnosticsHandler                                                DiagnosticsHandlerType;
      typedef Dune::Fem::DefaultSolverMonitorHandler                                              SolverMonitorHandlerType;
      typedef Dune::Fem::DefaultCheckPointHandler< GridType >                                     CheckPointHandlerType;
      typedef Dune::Fem::DefaultDataWriterHandler< GridType, IOTupleType >                        DataWriterHandlerType;
      typedef Dune::Fem::NoAdditionalOutputHandler                                                AdditionalOutputHandlerType;
      typedef Dune::Fem::DefaultSolutionLimiterHandler< LimiterOperatorType >                     SolutionLimiterHandlerType;
      typedef Dune::Fem::DefaultAdaptHandler< IndicatorType,
                                              GradientIndicatorType,
                                              SolutionLimiterHandlerType >                        AdaptHandlerType;
    };


  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::AdvectionDiffusionStepper< GridType, NavierStokesProblemCreator<GridType>, polOrd > Type;
  };


};

#endif // FEMHOWTO_HEATSTEPPER_HH
