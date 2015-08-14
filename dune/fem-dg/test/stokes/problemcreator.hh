#ifndef FEMDG_STOKESSTEPPER_HH
#define FEMDG_STOKESSTEPPER_HH
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
    typedef StokesModel< GridPart, InitialDataType >           ModelType;

    typedef std::vector< int >                                  EOCErrorIDs;

    static EOCErrorIDs initEoc()
    {
      EOCErrorIDs ids;
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-error" ) ) );
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "DG-error" ) ) );
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "sigma-norm" ) ) );
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
    //int probNr = Dune::Fem::Parameter::getValue< int > ( "problem" );
    //return new Dune :: GeneralizedStokesProblem< GridType, DIMRANGE > ( probNr );
    return new Dune::GeneralizedStokesProblem< GridType > ();
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
public:
    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static const int quadOrder = polynomialOrder * 3 + 1;
    static const SolverType solverType = istl ;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef std::tuple<> ExtraParameterTuple;

    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;
    typedef std::tuple< DiscreteFunctionType*, GridExactSolutionType* > IOTupleType;

private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
public:
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    // type of restriction/prolongation projection for adaptive simulations
    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                           RestrictionProlongationType;
    // type of linear solver for implicit ode
    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;


    //############################ poisson issues ############################
    static const bool symmetricSolver = true ;
    typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver> SolversType;

    //typedef Dune::UpwindFlux< typename AnalyticalTraitsType::ModelType >              FluxType;
    typedef Dune::NoFlux< typename AnalyticalTraitsType::ModelType >              FluxType;

    typedef typename SolversType::LinearOperatorType         FullOperatorType;
    typedef typename SolversType::LinearInverseOperatorType  BasicLinearSolverType;

    //----------- passes! ------------------------
    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType,  IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;
    typedef Dune::DGAdvectionDiffusionOperator< OperatorTraitsType >  AssemblyOperatorType;
    typedef Dune::DGPrimalMatrixAssembly< AssemblyOperatorType >            AssemblerType;
    //#########################################################################

    //------------- Limiter ---------------------------------------------
    typedef FullOperatorType                                                   LimiterOperatorType;
    //------------ Limiter ---------------------------------------------


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
    typedef StokesModel< GridPart, InitialDataType >                  ModelType;

    typedef typename PoissonProblemCreatorType::template AnalyticalTraits< GridPart >::EOCErrorIDs     EOCErrorIDs;

    static EOCErrorIDs initEoc ()
    {
      EOCErrorIDs ids;
      ids.push_back(Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-p-error" ) ) );
      return ids;
    }

    template< class SolutionImp, class Model, class ExactFunction >
    static void addEOCErrors ( const EOCErrorIDs &ids, SolutionImp &p, const Model &model, const ExactFunction &g  )
    {
      // calculate L2 - p-Norm
      Dune::Fem::L2Norm< GridPartType > l2pnorm( p.space().gridPart() );
      const double l2perror = l2pnorm.distance( g, p );
      Dune::Fem::FemEoc::setErrors( ids[ 0 ], l2perror );
    }
  };

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    return PoissonProblemCreatorType::initializeGrid();
  }

  static ProblemInterfaceType* problem()
  {
    return new Dune::GeneralizedStokesProblem< GridType > ();
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
    typedef typename PoissonProblemCreatorType::template DiscreteTraits< GridPart, polOrd >       PoissonDiscreteTraits;
public:
    typedef AnalyticalTraits<GridPart>                                                            AnalyticalTraitsType;

    static const int polynomialOrder = PoissonDiscreteTraits::polynomialOrder;

    static const int quadOrder = PoissonDiscreteTraits::quadOrder;
    static const SolverType solverType = PoissonDiscreteTraits::solverType;

    typedef typename DiscreteFunctionSpaces< PressureFunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef typename  std::tuple<> ExtraParameterTuple;


private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
public:
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                           RestrictionProlongationType;

    typedef Dune::AdaptationHandler< GridType, VelocityFunctionSpaceType >                              AdaptationHandlerType;


    //############################ stokes ############################
    static const bool symmetricSolver = true ;
    //typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver> SolversType;
    typedef Dune::NoFlux< typename AnalyticalTraitsType::ModelType >              FluxType;
    typedef typename PoissonDiscreteTraits::FullOperatorType                      FullOperatorType;

    //----------- passes! ------------------------
    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType, IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;
    typedef typename Dune::StokesAssembler< typename PoissonDiscreteTraits::DiscreteFunctionType,
                                            DiscreteFunctionType,
                                            OperatorTraitsType>                   AssemblerType;
    //------------------------------------------------------------

    static_assert( (int)DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

    typedef Dune::UzawaSolver< typename PoissonDiscreteTraits::DiscreteFunctionType, DiscreteFunctionType, AssemblerType,
                               typename PoissonDiscreteTraits::BasicLinearSolverType >    BasicLinearSolverType;

    //#########################################################################

    //------------- Limiter ---------------------------------------------
    typedef typename PoissonDiscreteTraits::LimiterOperatorType                                         LimiterOperatorType;
    //------------ Limiter ---------------------------------------------

    typedef typename AnalyticalTraitsType::ProblemType::ExactPressureType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;

    typedef std::tuple< typename std::tuple_element<0,typename PoissonDiscreteTraits::IOTupleType>::type, typename std::tuple_element<1,typename PoissonDiscreteTraits::IOTupleType>::type,
                        GridExactSolutionType*, DiscreteFunctionType* > IOTupleType;
  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::StokesAlgorithm< GridType, StokesProblemCreator<GridType>, PoissonProblemCreatorType, polOrd > Type;
  };


};

#endif // FEMHOWTO_HEATSTEPPER_HH
