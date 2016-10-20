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
//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>
//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/examples/stokes/stokesalgorithm.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
#include "incompnavierstokesalgorithm.hh"
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

  struct IncompressibleNavierStokesAlgorithmCreator
  {

  //=======================================
  //     STOKES PROBLEM CREATOR
  //=======================================

    struct SubStokesAlgorithmCreator
    {
      typedef AlgorithmConfigurator< Dune::GridSelector::GridType,
                                     Galerkin::Enum::dg,
                                     Adaptivity::Enum::yes,
                                     DiscreteFunctionSpaces::Enum::legendre,
                                     Solver::Enum::istl,
                                     AdvectionLimiter::Enum::unlimited,
                                     Matrix::Enum::assembled,
                                     AdvectionFlux::Enum::none,
                                     DiffusionFlux::Enum::primal > AC;


      struct SubPoissonAlgorithmCreator
      {

        typedef typename AC::GridType                         GridType;
        typedef typename AC::GridParts                        HostGridPartType;
        typedef HostGridPartType                              GridPartType;

        // define problem type here if interface should be avoided
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
            typedef typename AC::template LinearSolvers< DFSpaceType, true >      type;
          };

        private:
          //small helper class
          template< class SigmaDFSpaceType > struct SigmaFunctionChooser
          { typedef typename AC::template DiscreteFunctions< SigmaDFSpaceType > type; };
        public:
          typedef typename SigmaDiscreteFunctionSelector< DiscreteFunctionType, SigmaFunctionChooser >::type SigmaDiscreteFunctionType;

          typedef ErrorEstimator< DiscreteFunctionType, SigmaDiscreteFunctionType, typename Operator::AssemblerType >
                                                                                  ErrorEstimatorType;
          typedef PoissonSigmaEstimator< ErrorEstimatorType >                     SigmaEstimatorType;

          typedef PAdaptivity<DFSpaceType, polOrd, SigmaEstimatorType >           PAdaptivityType;

          typedef PAdaptIndicator< PAdaptivityType, ProblemInterfaceType >        AdaptIndicatorType;
          // typedef NoPAdaptIndicator                                               AdaptIndicatorType;

          typedef SubSolverMonitor< SolverMonitor >                               SolverMonitorType;
          typedef SubDiagnostics< Diagnostics >                                   DiagnosticsType;
        };

        template <int polOrd>
        using Algorithm = SubEllipticAlgorithm< GridType, SubPoissonAlgorithmCreator, polOrd >;
      };

      typedef typename AC::GridType                                       GridType;
      typedef typename AC::GridParts                                      HostGridPartType;
      typedef HostGridPartType                                            GridPartType;

      typedef NavierStokesProblemInterface< GridType >                    ProblemInterfaceType;

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
      }

      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        typedef typename SubPoissonAlgorithmCreator::template DiscreteTraits< polOrd > PoissonDiscreteTraits;
        typedef typename PoissonDiscreteTraits::DiscreteFunctionType                 VelDiscreteFunctionType;
        typedef typename SubPoissonAlgorithmCreator::FunctionSpaceType                 VelFunctionSpaceType;
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
          typedef UzawaSolver< typename Operator::AssemblerType, typename PoissonDiscreteTraits::Solver::type >
                                                                                      type;
        };

        static_assert( (int)DFSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

      private:
        typedef typename SubPoissonAlgorithmCreator::template DiscreteTraits<polOrd>::ErrorEstimatorType  PoissonErrorEstimatorType;
        typedef typename SubPoissonAlgorithmCreator::template DiscreteTraits<polOrd>::SigmaEstimatorType  PoissonSigmaEstimatorType;
        typedef typename SubPoissonAlgorithmCreator::template DiscreteTraits<polOrd>::PAdaptivityType     PoissonPAdaptivityType;

        typedef StokesSigmaEstimator< PoissonErrorEstimatorType, typename Operator::AssemblerType > StokesSigmaEstimatorType;
        typedef typename StokesSigmaEstimatorType::StokesErrorEstimatorType             StokesErrorEstimatorType;

        typedef StokesPAdaptivity< PoissonPAdaptivityType,
                                   DFSpaceType, polOrd, StokesSigmaEstimatorType >       StokesPAdaptivityType;
      public:

        typedef SubSolverMonitor< SolverMonitor >                                     SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                         DiagnosticsType;
        typedef StokesPAdaptIndicator< StokesPAdaptivityType, ProblemInterfaceType >    AdaptIndicatorType;
      };

      template <int polOrd>
      using Algorithm = SubStokesAlgorithm< GridType, SubStokesAlgorithmCreator, SubPoissonAlgorithmCreator, polOrd >;
    };


  //=======================================
  //     NAVIER-STOKES PROBLEM CREATOR
  //=======================================

    struct SubNavierStokesAlgorithmCreator
    {
      typedef AlgorithmConfigurator< Dune::GridSelector::GridType,
                                     Galerkin::Enum::dg,
                                     Adaptivity::Enum::yes,
                                     DiscreteFunctionSpaces::Enum::legendre,
                                     Solver::Enum::fem,
                                     AdvectionLimiter::Enum::unlimited,
                                     Matrix::Enum::matrixfree,
                                     AdvectionFlux::Enum::llf,
                                     DiffusionFlux::Enum::primal > AC;

      typedef typename AC::GridType                         GridType;
      typedef typename AC::GridParts                        HostGridPartType;
      typedef HostGridPartType                              GridPartType;

      // define problem type here if interface should be avoided
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
      }

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
          typedef typename AC::template LinearSolvers< DFSpaceType >                       LinearSolverType;
          typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                      type;
        };

      private:
        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                           OpTraits;
        typedef DGAdaptationIndicatorOperator< OpTraits >                                  IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >  GradientIndicatorType ;
      public:

        typedef AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
        typedef SubSolverMonitor< SolverMonitor >                                          SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                              DiagnosticsType;
        typedef ExactSolutionOutput< DiscreteFunctionType >                                AdditionalOutputType;
      };

      template <int polOrd>
      using Algorithm = SubAdvectionDiffusionAlgorithm< GridType, SubNavierStokesAlgorithmCreator, polOrd >;

    };

    template< int polOrd >
    using GlobalDiscreteFunction = std::tuple< typename SubStokesAlgorithmCreator::SubPoissonAlgorithmCreator::template DiscreteTraits< polOrd >::DiscreteFunctionType,
                                               typename SubStokesAlgorithmCreator::template DiscreteTraits< polOrd >::DiscreteFunctionType >;

    struct GlobalContainer
    {

      template< int polOrd, class Grid >
      static GlobalContainer& init( Grid& grid )
      {
        static const int numObjects = 2;


        typedef typename std::tuple_element< 0, GlobalDiscreteFunction< polOrd > >::type DFType1;
        typedef typename std::tuple_element< 1, GlobalDiscreteFunction< polOrd > >::type DFType2;


        typedef ContainerItem< DFType1 > ContainerItemType1;
        typedef ContainerItem< DFType2 > ContainerItemType2;

        typedef typename DFType1::DiscreteFunctionSpaceType SpaceType1;
        typedef typename DFType2::DiscreteFunctionSpaceType SpaceType2;

        typedef typename SpaceType1::GridPartType GPType1;
        typedef typename SpaceType2::GridPartType GPType2;

        typedef typename GPType1::GridType GType1;
        typedef typename GPType2::GridType GType2;

        static_assert( std::is_same< GType1, GridType >::value, "type of grids does not match" );
        static_assert( std::is_same< GType2, GridType >::value, "type of grids does not match" );

        GPType1 gridPart1( grid );
        GPType2 gridPart2( grid );

        SpaceType1 space1( gridPart1 );
        SpaceType2 space2( gridPart2 );

        DFType1 df1( "space1", space1 );
        DFType2 df2( "space2", space2 );


        ContainerItemType1 c1( "c", grid );

        ContainerItemType2 c7( "c7" ); c7 = c1.gridPart();

        ContainerItemType2 sol2( "sol2", c7 );
        ContainerItemType2 exactSol2( "exact", c7 );
        ContainerItemType2 rhs2( "rhs", c7 );
      }

    };


    template< class Grid, class... Args >
    using MyThetaSchemeCoupling=ThetaSchemeCoupling< Grid, GlobalContainer, Args... >;


    template <int polOrd>
    using Algorithm = IncompNavierStokesAlgorithm< polOrd, MyThetaSchemeCoupling, SubStokesAlgorithmCreator, SubNavierStokesAlgorithmCreator, SubStokesAlgorithmCreator >;

    typedef typename SubStokesAlgorithmCreator::GridType                  GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return DefaultGridInitializer< GridType >::initialize(); }

  };

}
}



#endif
