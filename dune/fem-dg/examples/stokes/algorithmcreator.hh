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
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>
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
    struct SubStokesAlgorithmCreator
    {


      struct SubPoissonAlgorithmCreator
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
          typedef Stokes::PoissonModel< GridType, InitialDataType >     ModelType;

          template< class Solution, class Model, class ExactFunction, class SigmaFunction>
          static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
          {
            static L2EOCError l2EocError( "$L^2$-Error" );
            l2EocError.add( u, f );
            static DGEOCError dgEocError( "DG-Error" );
            dgEocError.add( u, f );
            //static H1EOCError sigmaEocError( "sigma-norm" );
            //sigmaEocError.add( sigma, f );
          }
        };

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

        private:
          //small helper class
          template< class SigmaDFSpaceType > struct SigmaFunctionChooser
          { typedef typename AC::template DiscreteFunctions< SigmaDFSpaceType > type; };
        public:
          typedef typename SigmaDiscreteFunctionSelector< DiscreteFunctionType, SigmaFunctionChooser >::type SigmaDiscreteFunctionType;

          typedef ErrorEstimator< DiscreteFunctionType, SigmaDiscreteFunctionType, typename Operator::AssemblerType >
                                                                                     ErrorEstimatorType;
          typedef PoissonSigmaEstimator< ErrorEstimatorType >                        SigmaEstimatorType;

          typedef PAdaptivity<DFSpaceType, polOrd, SigmaEstimatorType >              PAdaptivityType;

          typedef PAdaptIndicator< PAdaptivityType, ProblemInterfaceType >           AdaptIndicatorType;
          // typedef NoPAdaptIndicator                                               AdaptIndicatorType;

          typedef SubSolverMonitor< SolverMonitor >                                  SolverMonitorType;
          typedef SubDiagnostics< Diagnostics >                                      DiagnosticsType;
        };

        template <int polOrd>
        using Algorithm = SubEllipticAlgorithm< GridType, SubPoissonAlgorithmCreator, polOrd >;
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
        typedef Stokes::StokesModel< GridType, InitialDataType >          ModelType;

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
        typedef typename SubPoissonAlgorithmCreator::template DiscreteTraits< polOrd > PoissonDiscreteTraits;
        typedef typename PoissonDiscreteTraits::DiscreteFunctionType                 VelDiscreteFunctionType;
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd-pressureOrderReduction, FunctionSpaceType>
                                                                                     DFSpaceType;
      public:
        typedef typename AC::template DiscreteFunctions< DFSpaceType >               DiscreteFunctionType;

        typedef std::tuple< VelDiscreteFunctionType, DiscreteFunctionType >          DiscreteFunctionsType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >           IOTupleType;
        typedef std::tuple<>                                                         ExtraParameterTuple;

        class Operator
        {
          typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd-pressureOrderReduction, AnalyticalTraits >
                                                                                     OpTraits;

          //typedef AssemblerTraitsList< std::tuple< typename std::tuple_element< 0, DiscreteFunctionsType>::type,
          //                                         typename std::tuple_element< 1, DiscreteFunctionsType>::type >,
          //                                         AC::template Containers >         AssTraits;
        public:
          typedef StokesAssembler< OpTraits, AC::template Containers,
                                   typename std::tuple_element< 0, DiscreteFunctionsType>::type,
                                   typename std::tuple_element< 1, DiscreteFunctionsType>::type>
                                                                                     AssemblerType;

          //template< template<class,class>class MatrixImp >
          //using AssemblerType = StokesAssembler< OpTraits, OpTraits, MatrixImp >
          //the following typedef is not needed by stokes algorithm atm
          //typedef typename AssemblerType::MatrixType                               type;
        };

        struct Solver
        {
          typedef UzawaSolver< typename Operator::AssemblerType,typename SubPoissonAlgorithmCreator::template Algorithm<polOrd> >
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
        typedef SubSolverMonitor< SolverMonitor >                                       SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                           DiagnosticsType;
        typedef StokesPAdaptIndicator< StokesPAdaptivityType, ProblemInterfaceType >    AdaptIndicatorType;
      };

      template <int polOrd>
      using Algorithm = SubStokesAlgorithm< GridType, SubStokesAlgorithmCreator<ACStokes>, SubPoissonAlgorithmCreator, polOrd >;
    };



    template <int polOrd>
    using Algorithm = SteadyStateAlgorithm< polOrd, SubStokesAlgorithmCreator<ACStokes> >;

    typedef typename SubStokesAlgorithmCreator<ACStokes>::GridType       GridType;

    static inline std::string moduleName() { return ""; }

  private:
    //helper struct
    template< class DDF, class RDF >
    using DefaultContainer = typename ACStokes::template Containers< typename DDF::DiscreteFunctionSpaceType, typename RDF::DiscreteFunctionSpaceType >;

  public:
    template< int polOrd >
    static decltype(auto) initContainer()
    {

      //Discrete Functions
      typedef typename SubStokesAlgorithmCreator<ACStokes>::SubPoissonAlgorithmCreator::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType1;
      typedef typename SubStokesAlgorithmCreator<ACStokes>::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType2;

      //Item1
      typedef _t< SubSteadyStateContainerItem >                         Steady;
      typedef std::tuple< Steady, Steady >                              Item1TupleType;

      //Item2
      typedef _t< SubEllipticContainerItem, DefaultContainer >          Def;
      typedef _t< SubEllipticContainerItem, NoPreconditioner, DefaultContainer > DefNo;

      typedef std::tuple< std::tuple< Def,   DefNo >,
                          std::tuple< DefNo, Def > >                        Item2TupleType;

      //Sub (discrete function argument ordering)
      typedef std::tuple<__0,__1 >                                      StokesOrder;

      typedef std::tuple< StokesOrder >                                 SubOrderRowType;
      typedef SubOrderRowType                                           SubOrderColType;



      //Global container
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, DFType1, DFType2 > GlobalContainerType;

      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( gridptr );

    }


  };
}
}
#endif // FEMHOWTO_HEATSTEPPER_HH
