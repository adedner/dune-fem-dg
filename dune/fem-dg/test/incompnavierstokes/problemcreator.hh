#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#include <config.h>

#ifndef GRIDDIM
#define GRIDDIM 2
#endif

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem/misc/l2norm.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/noflux.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/test/stokes/stokesalgorithm.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
#include "../stokes/problems.hh"
//--------- MODELS --------------------------
#include "stokesmodel.hh"
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
struct IncompressibleNavierStokesProblemCreator
{

//=======================================
//     STOKES PROBLEM CREATOR
//=======================================

  struct SubStokesProblemCreator
  {

    struct SubPoissonProblemCreator
    {
      typedef GridImp                                         GridType;
      typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
      typedef HostGridPartType                                GridPartType;

      // define problem type here if interface should be avoided
      typedef Dune::NavierStokesProblemDefault< GridType >           ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType       FunctionSpaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                   ProblemType;
        typedef ProblemInterfaceType                                   InitialDataType;
        typedef StokesModel< GridPartType, InitialDataType, false >    ModelType;

        template< class Solution, class Model, class ExactFunction, class SigmaFunction>
        static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
        {}
      };

      static inline std::string moduleName() { return "";}

      static ProblemInterfaceType* problem()
      {
        return new Dune::NavierStokesProblemDefault< GridType > ();
      }


      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
        static const bool symmetricSolver = true ;
        static const SolverType solverType = istl ;
        typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type    DiscreteFunctionSpaceType;
      public:
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                          DiscreteFunctionType;
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                      JacobianOperatorType;

        typedef std::tuple<>                                                                        ExtraParameterTuple;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

        class Operator
        {
          friend DiscreteTraits;
          friend SolverType;
          typedef Dune::NoFlux< typename AnalyticalTraits::ModelType >                                FluxType;

          typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                               DiscreteFunctionType, FluxType, ExtraParameterTuple >  OperatorTraitsType;

          typedef Dune::DGAdvectionDiffusionOperator< OperatorTraitsType >                            AssemblyOperatorType;
          typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>                     SolversType;



          struct RhsAnalyticalTraits
          {
            typedef ProblemInterfaceType                                   ProblemType;
            typedef ProblemInterfaceType                                   InitialDataType;
            typedef StokesModel< GridPartType, InitialDataType, true >     ModelType;
          };

          typedef Dune::NoFlux< typename RhsAnalyticalTraits::ModelType >                                RhsFluxType;
          typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, RhsAnalyticalTraits,
                                        DiscreteFunctionType, RhsFluxType, ExtraParameterTuple,
                                        FunctionSpaceType >                                           RhsOperatorTraitsType;


        public:
          typedef Dune::DGPrimalMatrixAssembly< AssemblyOperatorType >                                AssemblerType;
          typedef typename SolversType::LinearOperatorType                                            type;
          typedef Dune::DGAdvectionDiffusionOperator< RhsOperatorTraitsType >                         RhsType;
        };

        struct Solver
        {
          typedef typename Operator::SolversType::LinearInverseOperatorType                           type;
        };

        typedef Dune::Fem::SubSolverMonitorHandler< Dune::Fem::SolverMonitor >                    SolverMonitorHandlerType;
        typedef Dune::Fem::SubDiagnosticsHandler< Dune::Diagnostics >                             DiagnosticsHandlerType;
      };

      template <int polOrd>
      struct Stepper
      {
        // this should be ok but could lead to a henn-egg problem
        typedef Dune::Fem::EllipticAlgorithm< GridType, SubPoissonProblemCreator, polOrd > Type;
      };
    };

    typedef typename SubPoissonProblemCreator::GridType             GridType;
    typedef typename SubPoissonProblemCreator::HostGridPartType     HostGridPartType;
    typedef typename SubPoissonProblemCreator::GridPartType         GridPartType;

    typedef Dune::NavierStokesProblemDefault< GridType >             ProblemInterfaceType;

    typedef typename ProblemInterfaceType::PressureFunctionSpaceType FunctionSpaceType;

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                      ProblemType;
      typedef ProblemInterfaceType                                      InitialDataType;
      typedef StokesModel< GridPartType, InitialDataType, false >       ModelType;

      template< class Solution, class Model, class ExactFunction >
      static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f )
      {}
    };

    static inline std::string moduleName() { return "";}

    static ProblemInterfaceType* problem() { return new Dune::NavierStokesProblemDefault< GridType > (); }

    //Stepper Traits
    template< int polOrd >
    struct DiscreteTraits
    {
    private:
      typedef typename SubPoissonProblemCreator::template DiscreteTraits< polOrd >          PoissonDiscreteTraits;
      static const SolverType solverType = PoissonDiscreteTraits::solverType;
      static const bool symmetricSolver = true ;
      typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type    DiscreteFunctionSpaceType;
    public:
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

      typedef std::tuple<>                                                                        ExtraParameterTuple;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

      class Operator
      {
        friend DiscreteTraits;
        friend SolverType;
        typedef Dune::NoFlux< typename AnalyticalTraits::ModelType >                                FluxType;

        typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                      DiscreteFunctionType, FluxType, ExtraParameterTuple,
                                      typename SubPoissonProblemCreator::FunctionSpaceType >        OperatorTraitsType;

        struct RhsAnalyticalTraits
        {
          typedef ProblemInterfaceType                                   ProblemType;
          typedef ProblemInterfaceType                                   InitialDataType;
          typedef StokesModel< GridPartType, InitialDataType, true >     ModelType;
        };

        typedef Dune::NoFlux< typename RhsAnalyticalTraits::ModelType >                                RhsFluxType;
        typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, RhsAnalyticalTraits,
                                      DiscreteFunctionType, RhsFluxType, ExtraParameterTuple,
                                      typename SubPoissonProblemCreator::FunctionSpaceType >        RhsOperatorTraitsType;

      public:

        typedef typename Dune::StokesAssembler< typename PoissonDiscreteTraits::DiscreteFunctionType,
                                                DiscreteFunctionType,
                                                OperatorTraitsType>                                 AssemblerType;
        typedef typename PoissonDiscreteTraits::Operator::type                                      type;
        //typedef Dune::DGAdvectionDiffusionOperator< RhsOperatorTraitsType >                         RhsType;
        typedef void                         RhsType;
      };

      struct Solver
      {
        typedef Dune::UzawaSolver< typename PoissonDiscreteTraits::DiscreteFunctionType, DiscreteFunctionType, typename Operator::AssemblerType,
                                   typename PoissonDiscreteTraits::Solver::type >                   type;
      };

      static_assert( (int)DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

      typedef Dune::Fem::SubSolverMonitorHandler< Dune::Fem::SolverMonitor >                    SolverMonitorHandlerType;
      typedef Dune::Fem::SubDiagnosticsHandler< Dune::Diagnostics >                             DiagnosticsHandlerType;
    };

    template <int polOrd>
    struct Stepper
    {
      //typedef Dune::Fem::SubSteadyState2EvolutionAlgorithm< GridType, Dune::Fem::StokesAlgorithm< GridType, SubStokesProblemCreator, SubPoissonProblemCreator, polOrd >, polOrd >Type;
      typedef Dune::Fem::StokesAlgorithm< GridType, SubStokesProblemCreator, SubPoissonProblemCreator, polOrd > Type;
    };
  };


//=======================================
//     NAVIER-STOKES PROBLEM CREATOR
//=======================================

  struct SubNavierStokesProblemCreator
  {

    typedef GridImp                                         GridType;
    typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
    typedef HostGridPartType                                GridPartType;

    // define problem type here if interface should be avoided
    typedef Dune:: NavierStokesProblemDefault< GridType >                ProblemInterfaceType;

    typedef typename ProblemInterfaceType::FunctionSpaceType             FunctionSpaceType;

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                       ProblemType;
      typedef ProblemInterfaceType                                       InitialDataType;
      typedef NavierStokesModel< GridPartType, InitialDataType, false >  ModelType;

      template< class Solution, class Model, class ExactFunction, class TimeProvider >
      static void addEOCErrors ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
      {
        static L2EOCError l2EocError( "$L^2$-Error");
        l2EocError.add( tp, u, model, f );
      }
    };

    static inline std::string moduleName() { return ""; }

    static ProblemInterfaceType* problem() { return new ProblemInterfaceType(); }

    //Stepper Traits
    template< int polOrd >
    struct DiscreteTraits
    {
    private:
      static const SolverType solverType = fem;
      typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type             DiscreteFunctionSpaceType;
    public:
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

      typedef typename InitialProjectors< typename AnalyticalTraits::ProblemType::TimeDependentFunctionType, DiscreteFunctionType, dg >::type   InitialProjectorType;

      typedef std::tuple<> ExtraParameterTuple;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

      class Operator
      {
        friend DiscreteTraits;
        typedef LLFFlux< typename AnalyticalTraits::ModelType >                                    FluxType;

        typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits, DiscreteFunctionType, FluxType, ExtraParameterTuple >
                                                                                                    OperatorTraitsType;

        static const int hasAdvection = AnalyticalTraits::ModelType::hasAdvection;
        static const int hasDiffusion = AnalyticalTraits::ModelType::hasDiffusion;
        typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;

        struct RhsAnalyticalTraits
        {
          typedef ProblemInterfaceType                                       ProblemType;
          typedef ProblemInterfaceType                                       InitialDataType;
          typedef NavierStokesModel< GridPartType, InitialDataType, true >   ModelType;
        };

        typedef LLFFlux< typename RhsAnalyticalTraits::ModelType >                                  RhsFluxType;

        typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, RhsAnalyticalTraits,
                                      DiscreteFunctionType, RhsFluxType, ExtraParameterTuple >         RhsOperatorTraitsType;


      public:
        typedef typename AdvectionDiffusionOperatorType::FullOperatorType                           type;
        typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                       ExplicitType;
        typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                       ImplicitType;
        typedef Dune::DGDiffusionOperator< RhsOperatorTraitsType >                         RhsType;
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
     // typedef Dune::Fem::ExactSolutionOutputHandler< DiscreteFunctionType >                         AdditionalOutputHandlerType;
    };

    template <int polOrd>
    struct Stepper
    {
      typedef Dune::Fem::AdvectionDiffusionAlgorithm< GridType, SubNavierStokesProblemCreator, polOrd > Type;
    };

  };


  template <int polOrd>
  struct Stepper
  {
    typedef Dune::Fem::EvolutionAlgorithm< polOrd, SubStokesProblemCreator, SubNavierStokesProblemCreator, SubStokesProblemCreator > Type;
  };

  typedef GridImp                                         GridType;

  static inline std::string moduleName() { return ""; }

  static inline Dune::GridPtr<GridType>
  initializeGrid() { return Dune::Fem::DefaultGridInitializer< GridType >::initialize(); }

};

#endif
