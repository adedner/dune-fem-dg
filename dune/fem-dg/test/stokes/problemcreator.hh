#ifndef FEMDG_STOKESSTEPPER_HH
#define FEMDG_STOKESSTEPPER_HH
#include <config.h>

#ifndef POLORDER
#define POLORDER 1
#endif

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
#include <dune/fem-dg/test/poisson/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/noflux.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/test/stokes/stokesalgorithm.hh>
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
  typedef GridImp                                         GridType;
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
  typedef HostGridPartType                                GridPartType;

  // define problem type here if interface should be avoided
  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, GridType::dimensionworld >
                                                                FunctionSpaceType;
  typedef Dune::ProblemInterface< FunctionSpaceType >           ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef StokesModel< GridPart, InitialDataType >           ModelType;

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

  static inline Dune::GridPtr<GridType>
  initializeGrid() { return Dune::Fem::PoissonGridInitializer< GridType >::initializeGrid(); }

  static ProblemInterfaceType* problem()
  {
    //int probNr = Dune::Fem::Parameter::getValue< int > ( "problem" );
    //return new Dune :: GeneralizedStokesProblem< GridType, DIMRANGE > ( probNr );
    return new Dune::GeneralizedStokesProblem< GridType > ();
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
    static const SolverType solverType = istl ;
  public:
    typedef AnalyticalTraits< GridPartType >                                                   AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static const int quadOrder = polynomialOrder * 3 + 1;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef std::tuple<>                                                                        ExtraParameterTuple;

  private:
    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType               ExactSolutionType;

    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                      AdaptationHandlerType;
    //typedef Dune::UpwindFlux< typename AnalyticalTraitsType::ModelType >              FluxType;
    typedef Dune::NoFlux< typename AnalyticalTraitsType::ModelType >                    FluxType;
    static const Dune::DGDiffusionFluxIdentifier PrimalDiffusionFluxId =                Dune::method_general;


    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType,
                                  PrimalDiffusionFluxId,
                                  IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

    typedef Dune::DGAdvectionDiffusionOperator< OperatorTraitsType >                            AssemblyOperatorType;
  public:
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >                   GridExactSolutionType;

    typedef std::tuple< DiscreteFunctionType* >                                                 IOTupleType;

    static const bool symmetricSolver = true ;
    typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>                     SolversType;

    typedef Dune::DGPrimalMatrixAssembly< AssemblyOperatorType >                                AssemblerType;
    typedef typename SolversType::LinearOperatorType                                            OperatorType;

    typedef typename SolversType::LinearInverseOperatorType                                     BasicLinearSolverType;

    //------HANDLER-----------------------------------------------------
    struct HandlerTraits
    {
      typedef Dune::Fem::DefaultSteadyStateSolverMonitorHandler                                 SolverMonitorHandlerType;
      typedef Dune::Fem::DefaultDataWriterHandler< GridType, IOTupleType >                      DataWriterHandlerType;
    };

  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::EllipticAlgorithm< GridType, PoissonProblemCreator<GridType>, polOrd > Type;
  };


};




template< class GridImp >
struct StokesProblemCreator
{
  typedef PoissonProblemCreator<GridImp>                           PoissonProblemCreatorType;

  typedef typename PoissonProblemCreatorType::GridType             GridType;
  typedef typename PoissonProblemCreatorType::HostGridPartType     HostGridPartType;
  typedef typename PoissonProblemCreatorType::GridPartType         GridPartType;


private:
  typedef typename PoissonProblemCreatorType::FunctionSpaceType  VelocityFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, 1 >
                                                                PressureFunctionSpaceType;
public:
  // define problem type here if interface should be avoided
  typedef typename PoissonProblemCreatorType::FunctionSpaceType FunctionSpaceType;


  typedef Dune::StokesProblemInterface< VelocityFunctionSpaceType, PressureFunctionSpaceType >       ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                      ProblemType;
    typedef ProblemInterfaceType                                      InitialDataType;
    typedef StokesModel< GridPartType, InitialDataType >              ModelType;

    template< class Solution, class Model, class ExactFunction >
    static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f )
    {
      static L2EOCError l2EocError( "$L^2$-p-Error" );
      l2EocError.add( u, f );
    }
  };

  static inline std::string moduleName() { return "";}

  static inline Dune::GridPtr<GridType>
  initializeGrid() { return PoissonProblemCreatorType::initializeGrid(); }

  static ProblemInterfaceType* problem() { return new Dune::GeneralizedStokesProblem< GridType > (); }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
  private:
    typedef typename PoissonProblemCreatorType::template DiscreteTraits< GridPart, polOrd >       PoissonDiscreteTraits;
    static const SolverType solverType = PoissonDiscreteTraits::solverType;
  public:
    typedef AnalyticalTraits< GridPart >                                                        AnalyticalTraitsType;

    static const int polynomialOrder = PoissonDiscreteTraits::polynomialOrder;

    static const int quadOrder = PoissonDiscreteTraits::quadOrder;

    typedef typename DiscreteFunctionSpaces< PressureFunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef std::tuple<>                                                                        ExtraParameterTuple;

  private:
    typedef typename AnalyticalTraitsType::ProblemType::ExactPressureType                       ExactSolutionType;

    typedef Dune::AdaptationHandler< GridType, VelocityFunctionSpaceType >                      AdaptationHandlerType;

    typedef Dune::NoFlux< typename AnalyticalTraitsType::ModelType >                            FluxType;

    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    static const Dune::DGDiffusionFluxIdentifier PrimalDiffusionFluxId =                        Dune::method_general;

    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType, PrimalDiffusionFluxId, IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;


  public:
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >                   GridExactSolutionType;
    typedef std::tuple< typename std::tuple_element<0,typename PoissonDiscreteTraits::IOTupleType>::type, DiscreteFunctionType* >
                                                                                                IOTupleType;
    static const bool symmetricSolver = true ;

    typedef typename Dune::StokesAssembler< typename PoissonDiscreteTraits::DiscreteFunctionType,
                                            DiscreteFunctionType,
                                            OperatorTraitsType>                                 AssemblerType;
    typedef typename PoissonDiscreteTraits::OperatorType                                        OperatorType;
    typedef Dune::UzawaSolver< typename PoissonDiscreteTraits::DiscreteFunctionType, DiscreteFunctionType, AssemblerType,
                               typename PoissonDiscreteTraits::BasicLinearSolverType >          BasicLinearSolverType;

    static_assert( (int)DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");
    //------HANDLER-----------------------------------------------------
    struct HandlerTraits
    {
      typedef Dune::Fem::DefaultSteadyStateSolverMonitorHandler                                 SolverMonitorHandlerType;
      typedef Dune::Fem::DefaultDataWriterHandler< GridType, IOTupleType >                      DataWriterHandlerType;
    };

  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::StokesAlgorithm< GridType, StokesProblemCreator<GridType>, PoissonProblemCreatorType, polOrd > Type;
  };


};

#endif // FEMHOWTO_HEATSTEPPER_HH
