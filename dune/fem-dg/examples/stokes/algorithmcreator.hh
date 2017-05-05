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

//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
#include <dune/fem-dg/examples/poisson/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/examples/stokes/stokesalgorithm.hh>
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

//#define ROBERT_WANTS_USABILITY


namespace Dune
{
namespace Fem
{

  // for Taylorhood (P2,P1) equal one
  static const int pressureOrderReduction = 1;

  ////produce some static compiler warnings in case we are using an uninstalled solver
  //static const AvailableSolvers< solverEnum > checkSolverInstalled;

  template< class GridSelectorGridType >
  struct StokesAlgorithmCreator
  {
    typedef AlgorithmConfigurator< GridSelectorGridType,
                                   Galerkin::Enum::dg,
                                   Adaptivity::Enum::yes,
                                   DiscreteFunctionSpaces::Enum::hierarchic_legendre,
                                   Solver::Enum::istl,
                                   AdvectionLimiter::Enum::unlimited,
                                   Matrix::Enum::assembled,
                                   AdvectionFlux::Enum::none,
                                   DiffusionFlux::Enum::primal > ACStokes;
    template< class AC >
    struct SubPoissonAlgorithmCreator
    {

      typedef typename AC::GridType                         GridType;
      typedef typename AC::GridParts                        HostGridPartType;
      typedef HostGridPartType                              GridPartType;

      // define problem type here if interface should be avoided
      typedef typename Stokes::ProblemInterface< GridType >::PoissonProblemType  ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType               FunctionSpaceType;

      typedef Stokes::PoissonModel< GridType, ProblemInterfaceType >         ModelType;

      static inline std::string moduleName() { return "";}

      static ProblemInterfaceType* problem()
      {
        int problemNr = Parameter::getValue< int > ( "problem" );
        switch( problemNr )
        {
          case 1:
            return new typename Stokes::Problem< GridType, Stokes::DrivenCavityProblem >::PoissonProblemType ();
          default:
            return new typename Stokes::Problem< GridType, Stokes::GeneralizedProblem >::PoissonProblemType ();
        }
      }


      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctions< FunctionSpaceType, polOrd > DiscreteFunctionType;

        class Operator
        {
          typedef typename AC::template DefaultAssembTraits< ModelType, FunctionSpaceType, polOrd >
                                                                                     OpTraits;
        public:
          typedef typename AC::template Operators< OpTraits >                        AssemblerType;
          typedef typename AssemblerType::MatrixType                                 type;
        };

        struct Solver
        {
          typedef typename AC::template LinearSolvers< DiscreteFunctionType, true>   type;
        };

        typedef typename AC::template WrappedDiscreteFunctions< SigmaDiscreteFunctionSelector, DiscreteFunctionType >
                                                                                     SigmaDiscreteFunctionType;

        typedef ErrorEstimator< DiscreteFunctionType, SigmaDiscreteFunctionType, typename Operator::AssemblerType >
                                                                                     ErrorEstimatorType;
        typedef PoissonSigmaEstimator< ErrorEstimatorType >                          SigmaEstimatorType;

        typedef PAdaptivity<typename DiscreteFunctionType::DiscreteFunctionSpaceType, polOrd, SigmaEstimatorType >
                                                                                     PAdaptivityType;

        typedef PAdaptIndicator< PAdaptivityType, ModelType >                        AdaptIndicatorType;
        // typedef NoPAdaptIndicator                                                 AdaptIndicatorType;

        typedef SubSolverMonitor< SolverMonitor >                                    SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                        DiagnosticsType;
        typedef SubDataWriter< SolutionOutput<DiscreteFunctionType>, ExactSolutionOutput<DiscreteFunctionType> >
                                                                                     DataWriterType;

      };

      template <int polOrd>
      using Algorithm = SubEllipticAlgorithm< GridType, SubPoissonAlgorithmCreator<AC>, polOrd >;
    };


    template< class AC >
    struct SubStokesAlgorithmCreator
    {
      typedef SubPoissonAlgorithmCreator<AC>                   SubCreatorType;

      template< int polOrd >
      using SubAlgorithms = typename SubCreatorType::template Algorithm<polOrd>;


      typedef typename AC::GridType                            GridType;
      typedef typename AC::GridParts                           HostGridPartType;
      typedef HostGridPartType                                 GridPartType;

      // define problem type here if interface should be avoided
      typedef Stokes::ProblemInterface< GridType >             ProblemInterfaceType;

      typedef typename ProblemInterfaceType::StokesProblemType::FunctionSpaceType FunctionSpaceType;

      typedef Stokes::StokesModel< GridType, ProblemInterfaceType >              ModelType;

      static inline std::string moduleName() { return "";}

      static ProblemInterfaceType* problem()
      {
        int problemNr = Parameter::getValue< int > ( "problem" );
        switch( problemNr )
        {
          case 1:
            return new Stokes::Problem< GridType, Stokes::DrivenCavityProblem > ();
          default:
            return new Stokes::Problem< GridType, Stokes::GeneralizedProblem > ();
        }
      }

      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        static const int redPolOrd = polOrd-pressureOrderReduction;
      public:
        typedef typename AC::template DiscreteFunctions< FunctionSpaceType, redPolOrd > DiscreteFunctionType;

        class Operator
        {
          typedef typename AC::template DefaultAssembTraits< ModelType, FunctionSpaceType, redPolOrd >
                                                                                        OpTraits;
        public:
          typedef StokesAssembler< OpTraits, AC::template Containers,
                                   typename SubCreatorType::template DiscreteTraits< polOrd >::DiscreteFunctionType,
                                   DiscreteFunctionType >                               AssemblerType;

          //template< template<class,class>class MatrixImp >
          //using AssemblerType = StokesAssembler< OpTraits, OpTraits, MatrixImp >
          //the following typedef is not needed by stokes algorithm atm
          //typedef typename AssemblerType::MatrixType                                  type;
        };

        struct Solver
        {
          typedef UzawaSolver< typename Operator::AssemblerType,typename SubCreatorType::template Algorithm<polOrd> >
                                                                                        type;
        };

        static_assert( (int)FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

      private:
        typedef typename SubCreatorType::template DiscreteTraits<polOrd>::ErrorEstimatorType  PoissonErrorEstimatorType;
        typedef typename SubCreatorType::template DiscreteTraits<polOrd>::PAdaptivityType     PoissonPAdaptivityType;

        typedef StokesSigmaEstimator< PoissonErrorEstimatorType, typename Operator::AssemblerType > StokesSigmaEstimatorType;

        typedef StokesPAdaptivity< PoissonPAdaptivityType,typename DiscreteFunctionType::DiscreteFunctionSpaceType,
                                   polOrd, StokesSigmaEstimatorType >                   StokesPAdaptivityType;

      public:
        typedef SubSolverMonitor< SolverMonitor >                                       SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                           DiagnosticsType;
        typedef StokesPAdaptIndicator< StokesPAdaptivityType, ModelType >               AdaptIndicatorType;
        typedef SubDataWriter< SolutionOutput<DiscreteFunctionType>, ExactSolutionOutput<DiscreteFunctionType> >
                                                                                        DataWriterType;
      };

      template <int polOrd>
      using Algorithm = SubStokesAlgorithm< GridType, SubStokesAlgorithmCreator<AC>, SubCreatorType, polOrd >;
    };



    template <int polOrd>
    using Algorithm = SteadyStateAlgorithm< polOrd, SubStokesAlgorithmCreator<ACStokes> >;


    static inline std::string moduleName() { return ""; }

  public:
    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef typename SubStokesAlgorithmCreator<ACStokes>::GridType        GridType;

      //Discrete Functions
      typedef typename SubStokesAlgorithmCreator<ACStokes>::SubCreatorType::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType1;
      typedef typename SubStokesAlgorithmCreator<ACStokes>::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType2;

#ifdef ROBERT_WANTS_USABILITY
      typedef std::tuple< typename SubStokesAlgorithmCreator<ACStokes>::template Algorithm<polOrd>::ContainerType >
                                                                           ContainerTupleType;

      typedef std::tuple< std::tuple<__0,__1> >                            SubOrderType;

      //external params lists
      typedef std::tuple< ExtraArg<> >                                     ExtraType;

      //Global container
      typedef NewCombinedGlobalContainer< ContainerTupleType, SubOrderType, SubOrderType, ExtraType, DFType1, DFType2 >
                                                                           GlobalContainerType;
#else
      //Item1
      typedef _t< SubSteadyStateContainerItem >                             Steady;
      typedef std::tuple< Steady, Steady >                                  Item1TupleType;

      //Item2
      typedef _t< SubEllipticContainerItem, ACStokes::template Containers > Def;
      typedef _t< SubEllipticContainerItem, SparseRowLinearOperator >       Sp;
      typedef _t< SubEllipticContainerItem, NoPreconditioner, ACStokes::template Containers > DefNo;

      typedef std::tuple< std::tuple< Def, Sp >,
                          std::tuple< Sp, Sp > >                            Item2TupleType;

      //Sub (discrete function argument ordering)
      typedef std::tuple<__0,__1 >                                          StokesOrder;
      typedef std::tuple< StokesOrder >                                     SubOrderRowType;
      typedef SubOrderRowType                                               SubOrderColType;

      //external params lists
      typedef ExtraArg<>                                                    ExtraType;


      //Global container
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, ExtraType, DFType1, DFType2 >
                                                                            GlobalContainerType;
#endif
      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( moduleName(), gridptr );

    }


  };
}
}
#endif // FEMHOWTO_HEATSTEPPER_HH
