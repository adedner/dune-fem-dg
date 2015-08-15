#ifndef FEMDG_POISSONSTEPPER_HH
#define FEMDG_POISSONSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem/misc/l2norm.hh>

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
  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE >
                                                                FunctionSpaceType;
  typedef Dune::ProblemInterface< FunctionSpaceType >           ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef PoissonModel< GridPart, InitialDataType >           ModelType;

    typedef std::vector< int >                                  EOCErrorIDs;

    static EOCErrorIDs initEoc ()
    {
      EOCErrorIDs ids;
      ids.resize( 3 );
      ids[ 0 ] = Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-error" ) );
      ids[ 1 ] = Dune::Fem::FemEoc::addEntry( std::string( "DG-error" ) );
      ids[ 2 ] = Dune::Fem::FemEoc::addEntry( std::string( "sigma-norm" ) );
      return ids;
    }

    template< class SolutionImp, class Model, class ExactFunction, class SigmaEstimatorImp >
    static void addEOCErrors ( const EOCErrorIDs &ids, SolutionImp &u, const Model &model, const ExactFunction &f, SigmaEstimatorImp& sigma  )
    {
      // calculate L2 - Norm
      Dune::Fem::L2Norm< GridPartType > l2norm( u.space().gridPart() );
      const double l2error = l2norm.distance( f, u );

      Dune::Fem::DGNorm< GridPartType > dgnorm( u.space().gridPart() );
      const double dgerror = dgnorm.distance( f, u );

      Dune::Fem::H1Norm< GridPartType > sigmanorm( u.space().gridPart() );
      const double sigmaerror = sigmanorm.distance( f, sigma );

      Dune::Fem::FemEoc::setErrors( ids[ 0 ], l2error );
      Dune::Fem::FemEoc::setErrors( ids[ 1 ], dgerror );
      Dune::Fem::FemEoc::setErrors( ids[ 2 ], sigmaerror );
    }
  };

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    return Dune::Fem::PoissonGridInitializer< GridType >::initializeGrid();
  }

  static ProblemInterfaceType* problem()
  {
    int probNr = Dune::Fem::Parameter::getValue< int > ( "problem" );
    return new Dune :: PoissonProblem< GridType, DIMRANGE > ( probNr );
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
  private:
    static const SolverType solverType = istl ;
  public:
    typedef AnalyticalTraits< GridPartType >                                                    AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static const int quadOrder = polynomialOrder * 3 + 1;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef std::tuple<>                                                                        ExtraParameterTuple;

  private:
    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType                       ExactSolutionType;

    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;

    typedef Dune::UpwindFlux< typename AnalyticalTraitsType::ModelType >                        FluxType;

    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType,  IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

    typedef Dune::DGAdvectionDiffusionOperator< OperatorTraitsType >                            AssemblyOperatorType;
  public:
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >                   GridExactSolutionType;

    typedef std::tuple< DiscreteFunctionType*, GridExactSolutionType* >                         IOTupleType;

    static const bool symmetricSolver = true ;
    typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>                     SolversType;

    typedef Dune::DGPrimalMatrixAssembly< AssemblyOperatorType >                                AssemblerType;
    typedef typename SolversType::LinearOperatorType                                            OperatorType;

    typedef typename SolversType::LinearInverseOperatorType                                     BasicLinearSolverType;

    //------HANDLER-----------------------------------------------------
    struct HandlerTraits
    {
    private:
      // type of restriction/prolongation projection for adaptive simulations
      typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                         RestrictionProlongationType;
    public:
      typedef Dune::Fem::DefaultSteadyStateSolverMonitorHandler                                 SolverMonitorHandlerType;
    };

  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::EllipticAlgorithm< GridType, PoissonProblemCreator<GridType>, polOrd > Type;
  };


};

#endif // FEMHOWTO_HEATSTEPPER_HH
