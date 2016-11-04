#ifndef FEMDG_POISSONSTEPPER_HH
#define FEMDG_POISSONSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
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
#include "gridinitializer.hh"
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
#include <dune/fem-dg/operator/adaptation/poissonestimator.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/elliptic.hh>
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

  // only to ensure unique types in sub algorithm creators
  //namespace Global
  //{
  //  typedef std::tuple< Dune::GridSelector::GridType > GridTypes;

  //  typedef AlgorithmConfigurator< typename std::tuple_element< 0, GridTypes >::type,
  //                                 Galerkin::Enum::dg,
  //                                 Adaptivity::Enum::yes,
  //                                 DiscreteFunctionSpaces::Enum::legendre,
  //                                 Solver::Enum::istl,
  //                                 AdvectionLimiter::Enum::unlimited,
  //                                 Matrix::Enum::assembled,
  //                                 AdvectionFlux::Enum::upwind,
  //                                 DiffusionFlux::Enum::primal > AC;

  //  typedef std::tuple< AC >       AlgorithmConfiguratorTypes;

  //  typedef std::tuple< typename AC::GridParts >                  GridPartTypes;

  //  //typedef std::tuple< AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType> > DiscreteFunctionSpaceTypes;

  //  typedef std::tuple< AC::template DiscreteFunctionSpaceTemplates< GridPartType > > DiscreteFunctionSpaceTypes;

  //  typedef typename AC::template DiscreteFunctions< DiscreteFunctionSpaceTypes >             DiscreteFunctionType;


  //}

  struct PoissonAlgorithmCreator
  {
    struct SubPoissonAlgorithmCreator
    {
      typedef AlgorithmConfigurator< Dune::GridSelector::GridType,
                                     Galerkin::Enum::dg,
                                     Adaptivity::Enum::yes,
                                     DiscreteFunctionSpaces::Enum::legendre,
                                     Solver::Enum::istl,
                                     AdvectionLimiter::Enum::unlimited,
                                     Matrix::Enum::assembled,
                                     AdvectionFlux::Enum::upwind,
                                     DiffusionFlux::Enum::primal > AC;

      typedef typename AC::GridType                         GridType;
      typedef typename AC::GridParts                        HostGridPartType;
      typedef HostGridPartType                              GridPartType;

      // define problem type here if interface should be avoided
      typedef ProblemInterface< typename AC::template FunctionSpaces<DIMRANGE> >
                                                                    ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType      FunctionSpaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                ProblemType;
        typedef ProblemInterfaceType                                InitialDataType;
        typedef Poisson::Model< GridType, InitialDataType >         ModelType;

        template< class Solution, class Model, class ExactFunction, class SigmaFunction>
        static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
        {
          static L2EOCError l2EocError( "$L^2$-Error" );
          l2EocError.add( u, f );
          static DGEOCError dgEocError( "DG-Error" );
          dgEocError.add( u, f );

          //TODO what dimRange is sigma? This does not fit, yet
          //static H1EOCError sigmaEocError( "sigma-norm" );
          //sigmaEocError.add( sigma, f );
        }
      };

      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem()
      {
        int probNr = Parameter::getValue< int > ( "problem" );
        return new Poisson::Problem< GridType, DIMRANGE > ( probNr );
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
          typedef typename AssemblerType::LinearOperatorType                       type;
        };

        struct Solver
        {
          typedef typename AC::template LinearSolvers< DFSpaceType, false/*true*/> type;
        };

      private:
        //small helper class
        template< class SigmaDFSpaceType > struct SigmaFunctionChooser
        { typedef typename AC::template DiscreteFunctions< SigmaDFSpaceType > type; };

        //TODO improve -> algorithm configurator


        typedef typename SigmaDiscreteFunctionSelector< DiscreteFunctionType, SigmaFunctionChooser >::type SigmaDiscreteFunctionType;

        typedef ErrorEstimator< DiscreteFunctionType, SigmaDiscreteFunctionType, typename Operator::AssemblerType >
                                                                                   ErrorEstimatorType;
        typedef PoissonSigmaEstimator< ErrorEstimatorType >                        SigmaEstimatorType;

        typedef PAdaptivity<DFSpaceType, polOrd, SigmaEstimatorType >              PAdaptivityType;
      public:

        typedef SubSolverMonitor< SolverMonitor >                                  SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                      DiagnosticsType;
        typedef PAdaptIndicator< PAdaptivityType, ProblemInterfaceType >           AdaptIndicatorType;
      };

      template <int polOrd>
      using Algorithm = SubEllipticAlgorithm< GridType, SubPoissonAlgorithmCreator, polOrd >;

    };

    // structure of sub containers

    template <int polOrd>
    using Algorithm = SteadyStateAlgorithm< polOrd, UncoupledSubAlgorithms, SubPoissonAlgorithmCreator >;

    typedef typename SubPoissonAlgorithmCreator::GridType         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return DefaultGridInitializer< GridType >::initialize(); }



    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef std::tuple< std::tuple<_index<0> > >                                          SubOrderType;

      typedef std::tuple< _t< SubSteadyStateContainerItem > >                               Item1TupleType;
      typedef std::tuple< std::tuple< _t<SubEllipticContainerItem, ISTLLinearOperator > > > Item2TupleType;
      typedef typename SubPoissonAlgorithmCreator::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType;

      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderType, DFType > GlobalContainerType;

      //create grid
      std::unique_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( *gridptr );
    }

  };

}
}

#endif // FEMHOWTO_HEATSTEPPER_HH
