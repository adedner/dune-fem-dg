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

  template< class GridImp2 >
  struct SubAdvectionDiffusionProblemCreator
  {

    typedef GridImp2                                        GridType;
    typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
    typedef HostGridPartType                                GridPartType;

    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

    // define problem type here if interface should be avoided
    typedef Dune::EvolutionProblemInterface< FunctionSpaceType,false >      ProblemInterfaceType;

    template< class GridPart > // TODO: is this template parameter needed?
    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                ProblemType;
      typedef ProblemInterfaceType                                InitialDataType;
      typedef HeatEqnModel< GridPart, InitialDataType >           ModelType;

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
    template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
    struct DiscreteTraits
    {
    private:
      static const SolverType solverType = fem ;
      typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;
    public:

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
        public: //TODO private later
        friend HandlerTraits;
        typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;

        typedef Dune::UpwindFlux< typename AnalyticalTraitsType::ModelType >                        FluxType;

        typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
        typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
        typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             LimiterIndicatorType;

        typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                      DiscreteFunctionType, FluxType, LimiterIndicatorType,
                                      AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

        // TODO: advection/diffusion should not be precribed by model
        static const int hasAdvection = AnalyticalTraitsType::ModelType::hasAdvection;
        static const int hasDiffusion = AnalyticalTraitsType::ModelType::hasDiffusion;
        typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;
      public:
        typedef typename AdvectionDiffusionOperatorType::FullOperatorType                         FullType;
        typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                     ImplicitType;
        typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                     ExplicitType;
      };

    private:
      typedef Dune::DGAdaptationIndicatorOperator< typename OperatorType::OperatorTraitsType, OperatorType::hasAdvection, OperatorType::hasDiffusion >
                                                                                                    IndicatorType;
      typedef Estimator< DiscreteFunctionType, typename AnalyticalTraitsType::ProblemType >         GradientIndicatorType ;
    public:

      typedef Dune::Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
      typedef Dune::Fem::DefaultSolverMonitorHandler                                                SolverMonitorHandlerType;
      typedef Dune::Fem::DefaultDiagnosticsHandler                                                  DiagnosticsHandlerType;
    };

    template <int polOrd>
    struct Stepper
    {
     // this should be ok but could lead to a henn-egg problem
      typedef Dune::Fem::AdvectionDiffusionStepper< GridType, SubAdvectionDiffusionProblemCreator<GridType>, polOrd > Type;
    };

  };




  typedef GridImp                                         GridType;
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
  typedef HostGridPartType                                GridPartType;

  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

  // define problem type here if interface should be avoided
  typedef Dune::EvolutionProblemInterface< FunctionSpaceType,false >      ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef HeatEqnModel< GridPart, InitialDataType >           ModelType;

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

  template< class GridPart, int polOrd >
  struct DiscreteTraits;


  template <int polOrd>
  struct Stepper
  {
    typedef Dune::Fem::CombinedEvolutionAlgorithm< polOrd, typename DiscreteTraits< GridPartType, polOrd>::HandlerTraits, SubAdvectionDiffusionProblemCreator<GridType> > Type;
  };


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
  private:
    static const SolverType solverType = fem ;
    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;
  public:

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
      public: //private later
      friend HandlerTraits;
      typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;

      typedef Dune::UpwindFlux< typename AnalyticalTraitsType::ModelType >                        FluxType;

      typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
      typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
      typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             LimiterIndicatorType;

      typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                    DiscreteFunctionType, FluxType, LimiterIndicatorType,
                                    AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

      // TODO: advection/diffusion should not be precribed by model
      static const int hasAdvection = AnalyticalTraitsType::ModelType::hasAdvection;
      static const int hasDiffusion = AnalyticalTraitsType::ModelType::hasDiffusion;
      typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;
    public:
      typedef typename AdvectionDiffusionOperatorType::FullOperatorType                         FullType;
      typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                     ImplicitType;
      typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                     ExplicitType;
    };

  private:
    typedef Dune::DGAdaptationIndicatorOperator< typename OperatorType::OperatorTraitsType, OperatorType::hasAdvection, OperatorType::hasDiffusion >
                                                                                                  IndicatorType;
    typedef Estimator< DiscreteFunctionType, typename AnalyticalTraitsType::ProblemType >         GradientIndicatorType ;
  public:
    typedef Dune::Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                                AdaptIndicatorType;

    //------HANDLER-----------------------------------------------------
    class HandlerTraits
    {
            //limiting
      typedef typename OperatorType::FullType                                                        LimiterOperatorType;
      typedef typename SubAdvectionDiffusionProblemCreator<GridType>::template Stepper<polOrd>::Type SubStepperType;
    public:
      typedef Dune::Fem::CombinedDefaultDiagnosticsHandler< SubStepperType >                         DiagnosticsHandlerType;
      typedef Dune::Fem::CombinedDefaultSolverMonitorHandler< SubStepperType >                       SolverMonitorHandlerType;
      typedef Dune::Fem::CombinedDefaultCheckPointHandler< SubStepperType >                          CheckPointHandlerType;
      typedef Dune::Fem::CombinedDefaultDataWriterHandler< SubStepperType >                          DataWriterHandlerType;
      typedef Dune::Fem::NoAdditionalOutputHandler                                                   AdditionalOutputHandlerType;
      typedef Dune::Fem::CombinedDefaultSolutionLimiterHandler< SubStepperType >                     SolutionLimiterHandlerType;
      typedef Dune::Fem::CombinedDefaultAdaptHandler< SubStepperType >                               AdaptHandlerType;
    };

  };




};

#endif // FEMHOWTO_HEATSTEPPER_HH
