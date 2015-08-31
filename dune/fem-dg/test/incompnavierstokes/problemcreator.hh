#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#include <config.h>

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
#include <dune/fem-dg/algorithm/advectiondiffusionstepper.hh>
#include <dune/fem-dg/algorithm/advectionstepper.hh>
#include <dune/fem-dg/test/stokes/stokesalgorithm.hh>
#include <dune/fem-dg/algorithm/combinedevolution.hh>
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

  template< bool rhsModel = false >
  struct SubStokesProblemCreator
  {

    struct SubPoissonProblemCreator
    {
      typedef GridImp                                         GridType;
      typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
      typedef HostGridPartType                                GridPartType;

      // define problem type here if interface should be avoided
      typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE >
                                                                    FunctionSpaceType;
      typedef Dune::ProblemInterface< FunctionSpaceType >           ProblemInterfaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                   ProblemType;
        typedef ProblemInterfaceType                                   InitialDataType;
        typedef StokesModel< GridPartType, InitialDataType, rhsModel > ModelType;

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
        static const SolverType solverType = istl ;
      public:
        typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type    DiscreteFunctionSpaceType;
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                          DiscreteFunctionType;
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                      JacobianOperatorType;

        typedef std::tuple<>                                                                        ExtraParameterTuple;

      private:
        typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;
        //typedef Dune::UpwindFlux< typename AnalyticalTraits::ModelType >                          FluxType;
        typedef Dune::NoFlux< typename AnalyticalTraits::ModelType >                                FluxType;

        typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraits::ModelType::dimDomain, 3> FVFunctionSpaceType;
        typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
        typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

        typedef Dune::OperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                      DiscreteFunctionType, FluxType,  IndicatorType,
                                      AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

        typedef Dune::DGAdvectionDiffusionOperator< OperatorTraitsType >                            AssemblyOperatorType;
      public:
        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

        static const bool symmetricSolver = true ;
        typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>                     SolversType;

        typedef Dune::DGPrimalMatrixAssembly< AssemblyOperatorType >                                AssemblerType;
        typedef typename SolversType::LinearOperatorType                                            OperatorType;

        typedef typename SolversType::LinearInverseOperatorType                                     BasicLinearSolverType;


        typedef Dune::Fem::SubSolverMonitorHandler< Dune::Fem::SolverMonitor< 1 > >               SolverMonitorHandlerType;
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


  private:
    typedef typename SubPoissonProblemCreator::FunctionSpaceType  VelocityFunctionSpaceType;
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, 1 >
                                                                  PressureFunctionSpaceType;
  public:
    // define problem type here if interface should be avoided
    typedef typename SubPoissonProblemCreator::FunctionSpaceType FunctionSpaceType;

    typedef Dune::StokesProblemInterface< VelocityFunctionSpaceType, PressureFunctionSpaceType >       ProblemInterfaceType;

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                      ProblemType;
      typedef ProblemInterfaceType                                      InitialDataType;
      typedef StokesModel< GridPartType, InitialDataType, rhsModel >    ModelType;

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
    public:
      typedef typename DiscreteFunctionSpaces< PressureFunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type    DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

      typedef std::tuple<>                                                                        ExtraParameterTuple;

    private:
      typedef Dune::AdaptationHandler< GridType, VelocityFunctionSpaceType >                      AdaptationHandlerType;

      typedef Dune::NoFlux< typename AnalyticalTraits::ModelType >                                FluxType;

      typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraits::ModelType::dimDomain, 3> FVFunctionSpaceType;
      typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
      typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

      typedef Dune::OperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                    DiscreteFunctionType, FluxType, IndicatorType,
                                    AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;


    public:
      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;
      static const bool symmetricSolver = true ;

      typedef typename Dune::StokesAssembler< typename PoissonDiscreteTraits::DiscreteFunctionType,
                                              DiscreteFunctionType,
                                              OperatorTraitsType>                                 AssemblerType;
      typedef typename PoissonDiscreteTraits::OperatorType                                        OperatorType;
      typedef Dune::UzawaSolver< typename PoissonDiscreteTraits::DiscreteFunctionType, DiscreteFunctionType, AssemblerType,
                                 typename PoissonDiscreteTraits::BasicLinearSolverType >          BasicLinearSolverType;

      static_assert( (int)DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

      typedef Dune::Fem::SubSolverMonitorHandler< Dune::Fem::SolverMonitor< 1 > >               SolverMonitorHandlerType;
      typedef Dune::Fem::SubDiagnosticsHandler< Dune::Diagnostics >                             DiagnosticsHandlerType;
    };

    template <int polOrd>
    struct Stepper
    {
      typedef Dune::Fem::StokesAlgorithm< GridType, SubStokesProblemCreator, SubPoissonProblemCreator, polOrd > Type;
    };
  };


//=======================================
//     NAVIER-STOKES PROBLEM CREATOR
//=======================================

  template< bool rhsModel = false >
  struct SubNavierStokesProblemCreator
  {

    typedef GridImp                                         GridType;
    typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
    typedef HostGridPartType                                GridPartType;

    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

    // define problem type here if interface should be avoided
    typedef Dune:: NavierStokesProblemDefault< GridType >                    ProblemInterfaceType;

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                ProblemType;
      typedef ProblemInterfaceType                                InitialDataType;
      typedef  NavierStokesModel< GridPartType, InitialDataType, rhsModel >  ModelType;

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

      class HandlerTraits;

      class OperatorType
      {
        friend DiscreteTraits;
        typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;

        typedef LLFFlux< typename AnalyticalTraits::ModelType >                                 FluxType;

        typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraits::ModelType::dimDomain, 3> FVFunctionSpaceType;
        typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
        typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

        typedef Dune::OperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                      DiscreteFunctionType, FluxType, IndicatorType,
                                      AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

        static const int hasAdvection = AnalyticalTraits::ModelType::hasAdvection;
        static const int hasDiffusion = AnalyticalTraits::ModelType::hasDiffusion;
        typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;
      public:
        typedef typename AdvectionDiffusionOperatorType::FullOperatorType                           FullType;
        typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                       ImplicitType;
        typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                       ExplicitType;
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
      typedef Dune::Fem::AdvectionDiffusionStepper< GridType, SubNavierStokesProblemCreator, polOrd > Type;
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
