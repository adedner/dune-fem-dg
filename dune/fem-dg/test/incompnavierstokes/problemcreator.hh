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
//--------- HANDLER --------------------------------
#include <dune/fem-dg/algorithm/handler/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/sub/adapt.hh>
//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
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
//#include "../stokes/problems.hh"
//--------- MODELS --------------------------
#include "stokesmodel.hh"
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/configurator.hh>


namespace Dune
{
namespace Fem
{

  template< class GridImp >
  struct IncompressibleNavierStokesProblemCreator
  {

  //=======================================
  //     STOKES PROBLEM CREATOR
  //=======================================

    struct SubStokesProblemCreator
    {
      typedef AlgorithmConfigurator< GridImp,
                                     Galerkin::Enum::dg,
                                     Adaptivity::Enum::yes,
                                     DiscreteFunctionSpaces::Enum::legendre,
                                     Solver::Enum::istl,
                                     AdvectionLimiter::Enum::unlimited,
                                     Matrix::Enum::assembled,
                                     AdvectionFlux::Identifier<AdvectionFlux::Enum::none>,
                                     PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::general > > AC;


      struct SubPoissonProblemCreator
      {

        typedef typename AC::GridType                         GridType;
        typedef typename AC::GridParts                        HostGridPartType;
        typedef HostGridPartType                              GridPartType;

        // define problem type here if interface should be avoided
        //typedef ThetaProblemInterface< typename AC::template FunctionSpaces<DIMRANGE> >  ProblemInterfaceType;
        typedef typename NavierStokesProblemInterface< GridType >::PoissonProblemType  ProblemInterfaceType;

        typedef typename ProblemInterfaceType::FunctionSpaceType       FunctionSpaceType;

        struct AnalyticalTraits
        {
          typedef ProblemInterfaceType                                 ProblemType;
          typedef ProblemInterfaceType                                 InitialDataType;
          typedef PoissonModel< GridPartType, InitialDataType, false > ModelType;

          template< class Solution, class Model, class ExactFunction, class SigmaFunction>
          static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
          {}
        };

        static inline std::string moduleName() { return "";}

        static ProblemInterfaceType* problem()
        {
          return new typename NavierStokesProblem< GridType, NavierStokesProblemDefault >::PoissonProblemType();
        }

        template< int polOrd >
        struct DiscreteTraits
        {
          typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                  DFSpaceType;
        public:
          typedef typename AC::template DiscreteFunctions< DFSpaceType >          DiscreteFunctionType;

          typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >      IOTupleType;
          typedef std::tuple<>                                                    ExtraParameterTuple;

          class Operator
          {
            typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd, AnalyticalTraits >
                                                                                  OpTraits;
            typedef typename AC::template RhsAnalyticalTraits< AnalyticalTraits, PoissonModel< GridPartType, typename AnalyticalTraits::InitialDataType, true > >
                                                                                  RhsAnalyticalTraitsType;

            typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, RhsAnalyticalTraitsType, ExtraParameterTuple >
                                                                                  RhsOpTraits;
          public:
            typedef typename AC::template Operators< OpTraits >                   AssemblerType;
            typedef typename AssemblerType::MatrixType                            type;

            typedef DGAdvectionDiffusionOperator< RhsOpTraits >                   RhsType;
          };

          struct Solver
          {
            typedef typename AC::template LinearSolvers< DFSpaceType, true > type;
          };

          typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >               SolverMonitorHandlerType;
          typedef Fem::SubDiagnosticsHandler< Diagnostics >                        DiagnosticsHandlerType;
        };

        template <int polOrd>
        struct Stepper
        {
          // this should be ok but could lead to a henn-egg problem
          typedef Fem::EllipticAlgorithm< GridType, SubPoissonProblemCreator, polOrd > Type;
        };
      };

      typedef typename AC::GridType                                    GridType;
      typedef typename AC::GridParts                                   HostGridPartType;
      typedef HostGridPartType                                         GridPartType;

      typedef NavierStokesProblemInterface< GridType >                 ProblemInterfaceType;
      //typedef StokesProblemInterface< typename SubPoissonProblemCreator::ProblemInterfaceType /*velocity*/,
      //                                ThetaProblemInterface< typename AC::template FunctionSpaces<1> > >
      //                                                                    ProblemInterfaceType;

      typedef typename ProblemInterfaceType::StokesProblemType::FunctionSpaceType
                                                                          FunctionSpaceType;

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

      static ProblemInterfaceType* problem()
      {
        return new NavierStokesProblem< GridType, NavierStokesProblemDefault > ();
        //return new NavierStokesProblemDefault< GridType > ();
      }

      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        typedef typename SubPoissonProblemCreator::template DiscreteTraits< polOrd > PoissonDiscreteTraits;
        typedef typename PoissonDiscreteTraits::DiscreteFunctionType                 VelDiscreteFunctionType;
        typedef typename SubPoissonProblemCreator::FunctionSpaceType                 VelFunctionSpaceType;
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                     DFSpaceType;
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, VelFunctionSpaceType>
                                                                                     VelDFSpaceType;
      public:
        typedef typename AC::template DiscreteFunctions< DFSpaceType >               DiscreteFunctionType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >           IOTupleType;
        typedef std::tuple<>                                                         ExtraParameterTuple;

        class Operator
        {
          typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd, AnalyticalTraits >
                                                                                     OpTraits;

          typedef AssemblerTraitsList< std::tuple< VelDiscreteFunctionType, DiscreteFunctionType >, AC::template Containers > AssTraits;

          typedef typename AC::template RhsAnalyticalTraits< AnalyticalTraits, StokesModel< GridPartType, typename AnalyticalTraits::InitialDataType, true > >
                                                                                RhsAnalyticalTraitsType;
          typedef typename AC::template DefaultOpTraits< VelDFSpaceType, polOrd, RhsAnalyticalTraitsType, ExtraParameterTuple >
                                                                                RhsOpTraits;
        public:
          typedef StokesAssembler< AssTraits, OpTraits >                             AssemblerType;
          //the following typedef is not needed by stokes algorithm atm
          //typedef typename AssemblerType::MatrixType                               type;

          //typedef DGAdvectionDiffusionOperator< RhsOpTraits >                        RhsType;
        };

        struct Solver
        {
          typedef UzawaSolver< VelDiscreteFunctionType, DiscreteFunctionType, typename Operator::AssemblerType,
                               typename PoissonDiscreteTraits::Solver::type >        type;
        };

        static_assert( (int)DFSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

        typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                    SolverMonitorHandlerType;
        typedef Fem::SubDiagnosticsHandler< Diagnostics >                             DiagnosticsHandlerType;
      };

      template <int polOrd>
      struct Stepper
      {
        typedef Fem::StokesAlgorithm< GridType, SubStokesProblemCreator, SubPoissonProblemCreator, polOrd > Type;
      };
    };


  //=======================================
  //     NAVIER-STOKES PROBLEM CREATOR
  //=======================================

    struct SubNavierStokesProblemCreator
    {
      typedef AlgorithmConfigurator< GridImp,
                                     Galerkin::Enum::dg,
                                     Adaptivity::Enum::yes,
                                     DiscreteFunctionSpaces::Enum::legendre,
                                     Solver::Enum::fem,
                                     AdvectionLimiter::Enum::unlimited,
                                     Matrix::Enum::matrixfree,
                                     AdvectionFlux::Identifier<AdvectionFlux::Enum::llf>,
                                     PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::general > > AC;

      typedef typename AC::GridType                         GridType;
      typedef typename AC::GridParts                        HostGridPartType;
      typedef HostGridPartType                              GridPartType;

      // define problem type here if interface should be avoided
      //typedef NavierStokesProblemDefault< GridType >                       ProblemInterfaceType;
      typedef typename NavierStokesProblemInterface< GridType >::NavierStokesProblemType ProblemInterfaceType;

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

      static ProblemInterfaceType* problem()
      {
        return new typename NavierStokesProblem< GridType, NavierStokesProblemDefault >::NavierStokesProblemType();
        //return new ProblemInterfaceType();
      }

      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                           DFSpaceType;
      public:
        typedef typename AC::template DiscreteFunctions< DFSpaceType >                     DiscreteFunctionType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                 IOTupleType;
        typedef std::tuple<>                                                               ExtraParameterTuple;

        class Operator
        {
          typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                           OpTraits;

          typedef typename AC::template RhsAnalyticalTraits< AnalyticalTraits, NavierStokesModel< GridPartType, typename AnalyticalTraits::InitialDataType, true > >
                                                                                           RhsAnalyticalTraitsType;
          typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, RhsAnalyticalTraitsType, ExtraParameterTuple >
                                                                                           RhsOpTraits;

        public:
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::full >    type;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::expl >    ExplicitType;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::impl >    ImplicitType;

          typedef DGAdvectionDiffusionOperator< RhsOpTraits >                              RhsType;
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
        typedef Fem::AdvectionDiffusionAlgorithm< GridType, SubNavierStokesProblemCreator, polOrd > Type;
      };

    };


    template <int polOrd>
    struct Stepper
    {
      typedef Fem::EvolutionAlgorithm< polOrd, SubStokesProblemCreator, SubNavierStokesProblemCreator, SubStokesProblemCreator > Type;
    };

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return Fem::DefaultGridInitializer< GridType >::initialize(); }

  };

}
}
#endif
