#ifndef FEMDG_POISSONSTEPPER_HH
#define FEMDG_POISSONSTEPPER_HH
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
#include "gridinitializer.hh"
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/noflux.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/ellipticalgorithm.hh>
#include <dune/fem-dg/algorithm/combinedsteadystate.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>

template< class GridImp >
struct PoissonProblemCreator
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
      typedef ProblemInterfaceType                                ProblemType;
      typedef ProblemInterfaceType                                InitialDataType;
      typedef PoissonModel< GridPartType, InitialDataType >       ModelType;

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

    static inline std::string moduleName() { return ""; }

    static ProblemInterfaceType* problem()
    {
      int probNr = Dune::Fem::Parameter::getValue< int > ( "problem" );
      return new Dune :: PoissonProblem< GridType, DIMRANGE > ( probNr );
    }


    //Stepper Traits
    template< int polOrd >
    struct DiscreteTraits
    {
    private:
      static const SolverType solverType = istl ;
      static const bool symmetricSolver = true ;
    public:
      typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type    DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                          DiscreteFunctionType;
      typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                      JacobianOperatorType;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

      typedef std::tuple<>                                                                        ExtraParameterTuple;

      class Operator
      {
        friend DiscreteTraits;
        friend SolverType;
        typedef Dune::UpwindFlux< typename AnalyticalTraits::ModelType >                            FluxType;

        typedef Dune::DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                             DiscreteFunctionType, FluxType, ExtraParameterTuple >  OperatorTraitsType;

        typedef Dune::DGAdvectionDiffusionOperator< OperatorTraitsType >                            AssemblyOperatorType;
        typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>                     SolversType;
      public:
        typedef Dune::DGPrimalMatrixAssembly< AssemblyOperatorType >                                AssemblerType;
        typedef typename SolversType::LinearOperatorType                                            type;
      };

      struct Solver
      {
        typedef typename Operator::SolversType::LinearInverseOperatorType                           type;
      };

    private:
      //typedef Dune::DGAdaptationIndicatorOperator< OperatorTraitsType, true, true >                 IndicatorType;
      //typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >             GradientIndicatorType ;
    public:

      //typedef Dune::Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
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

  template <int polOrd>
  struct Stepper
  {
    typedef Dune::Fem::SteadyStateAlgorithm< polOrd, SubPoissonProblemCreator > Type;
  };

  typedef GridImp                                         GridType;

  static inline std::string moduleName() { return ""; }

  static inline Dune::GridPtr<GridType>
  initializeGrid() { return Dune::Fem::DefaultGridInitializer< GridType >::initialize(); }



};

#endif // FEMHOWTO_HEATSTEPPER_HH
