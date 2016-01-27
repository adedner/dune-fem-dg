#ifndef FEMDG_STOKESSTEPPER_HH
#define FEMDG_STOKESSTEPPER_HH
#include <config.h>

#include <dune/fem-dg/misc/static_warning.hh>


#ifndef GRIDDIM
#define GRIDDIM 2
#endif

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM
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
#include <dune/fem-dg/test/poisson/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/test/stokes/stokesalgorithm.hh>
#include <dune/fem-dg/algorithm/steadystate.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
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

  static const Galerkin::Enum               galerkinEnum   = Galerkin::Enum::dg;
  static const Adaptivity::Enum             adaptivityEnum = Adaptivity::Enum::yes;
  static const DiscreteFunctionSpaces::Enum dfSpaceEnum    = DiscreteFunctionSpaces::Enum::hierarchic_legendre;
  static const Solver::Enum                 solverEnum     = Solver::Enum::istl;
  static const AdvectionLimiter::Enum       advLimiterEnum = AdvectionLimiter::Enum::unlimited;
  static const Matrix::Enum                 matrixEnum     = Matrix::Enum::assembled;
  static const AdvectionFlux::Enum          advFluxEnum    = AdvectionFlux::Enum::none;
  static const PrimalDiffusionFlux::Enum    diffFluxEnum   = PrimalDiffusionFlux::Enum::general;

  //produce some static compiler warnings in case we are using an uninstalled solver
  static const AvailableSolvers< solverEnum > checkSolverInstalled;

  template< class GridImp >
  struct StokesProblemCreator
  {

    struct SubStokesProblemCreator
    {
      typedef AlgorithmConfigurator< GridImp, galerkinEnum, adaptivityEnum, dfSpaceEnum, solverEnum, advLimiterEnum,
                                     matrixEnum, advFluxEnum, diffFluxEnum > AC;

      struct SubPoissonProblemCreator
      {

        typedef typename AC::GridType                         GridType;
        typedef typename AC::GridParts                        HostGridPartType;
        typedef HostGridPartType                              GridPartType;

        // define problem type here if interface should be avoided
        typedef typename Stokes::ProblemInterface< GridType >::PoissonProblemType  ProblemInterfaceType;

        typedef typename ProblemInterfaceType::FunctionSpaceType      FunctionSpaceType;


        struct AnalyticalTraits
        {
          typedef ProblemInterfaceType                                  ProblemType;
          typedef ProblemInterfaceType                                  InitialDataType;
          typedef Stokes::PoissonModel< GridPartType, InitialDataType > ModelType;

          template< class Solution, class Model, class ExactFunction, class SigmaFunction>
          static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
          {
            static L2EOCError l2EocError( "$L^2$-Error" );
            l2EocError.add( u, f );
            static DGEOCError dgEocError( "DG-Error" );
            dgEocError.add( u, f );
            static H1EOCError sigmaEocError( "sigma-norm" );
            sigmaEocError.add( sigma, f );
          }
        };

        static inline std::string moduleName() { return "";}

        static ProblemInterfaceType* problem()
        {
          return new typename Stokes::Problem< GridType, Stokes::GeneralizedProblem >::PoissonProblemType();
          //return new typename Stokes::GeneralizedProblem< GridType >::PoissonProblemType();
        }


        template< int polOrd >
        struct DiscreteTraits
        {
          typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                     DFSpaceType;
        public:
          typedef typename AC::template DiscreteFunctions< DFSpaceType >             DiscreteFunctionType;

          typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >         IOTupleType;

          class Operator
          {
            typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd, AnalyticalTraits >
                                                                                     OpTraits;
          public:
            typedef typename AC::template Operators< OpTraits >                      AssemblerType;
            typedef typename AssemblerType::MatrixType                               type;
          };

          struct Solver
          {
            typedef typename AC::template LinearSolvers< DFSpaceType, true> type;
          };

          typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                 SolverMonitorHandlerType;
          typedef Fem::SubDiagnosticsHandler< Diagnostics >                          DiagnosticsHandlerType;
        };

        template <int polOrd>
        struct Stepper
        {
          // this should be ok but could lead to a henn-egg problem
          typedef Fem::SubEllipticAlgorithm< GridType, SubPoissonProblemCreator, polOrd > Type;
        };
      };

      typedef typename AC::GridType                         GridType;
      typedef typename AC::GridParts                        HostGridPartType;
      typedef HostGridPartType                              GridPartType;

      // define problem type here if interface should be avoided
      typedef Stokes::ProblemInterface< GridType >             ProblemInterfaceType;

      typedef typename ProblemInterfaceType::StokesProblemType::FunctionSpaceType FunctionSpaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                      ProblemType;
        typedef ProblemInterfaceType                                      InitialDataType;
        typedef Stokes::StokesModel< GridPartType, InitialDataType >      ModelType;

        template< class Solution, class Model, class ExactFunction >
        static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f )
        {
          static L2EOCError l2EocError( "$L^2$-p-Error" );
          l2EocError.add( u, f );
        }
      };

      static inline std::string moduleName() { return "";}

      static ProblemInterfaceType* problem()
      {
        return new Stokes::Problem< GridType, Stokes::GeneralizedProblem > ();
      }

      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        typedef typename SubPoissonProblemCreator::template DiscreteTraits< polOrd > PoissonDiscreteTraits;
        typedef typename PoissonDiscreteTraits::DiscreteFunctionType                 VelDiscreteFunctionType;
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                     DFSpaceType;
      public:
        typedef typename AC::template DiscreteFunctions< DFSpaceType >               DiscreteFunctionType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >           IOTupleType;
        typedef std::tuple<>                                                         ExtraParameterTuple;

        class Operator
        {
          typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd, AnalyticalTraits >
                                                                                     OpTraits;

          typedef AssemblerTraitsList< std::tuple< VelDiscreteFunctionType, DiscreteFunctionType >, AC::template Containers > AssTraits;
        public:
          typedef StokesAssembler< AssTraits, OpTraits >                             AssemblerType;
          //the following typedef is not needed by stokes algorithm atm
          //typedef typename AssemblerType::MatrixType                               type;
        };

        struct Solver
        {
          typedef UzawaSolver< typename Operator::AssemblerType, typename PoissonDiscreteTraits::Solver::type >
                                                                                     type;
        };

        static_assert( (int)DFSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

        typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                   SolverMonitorHandlerType;
        typedef Fem::SubDiagnosticsHandler< Diagnostics >                            DiagnosticsHandlerType;
      };

      template <int polOrd>
      struct Stepper
      {
        typedef Fem::SubStokesAlgorithm< GridType, SubStokesProblemCreator, SubPoissonProblemCreator, polOrd > Type;
      };
    };


    template <int polOrd>
    struct Stepper
    {
      typedef Fem::SteadyStateAlgorithm< polOrd, DefaultSteadyStateCreator, SubStokesProblemCreator > Type;
    };

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return Fem::PoissonGridInitializer< GridType >::initializeGrid(); }


  };
}
}
#endif // FEMHOWTO_HEATSTEPPER_HH
