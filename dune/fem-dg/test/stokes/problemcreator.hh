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
#include "corner.hh"

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/noflux.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/test/stokes/stokes.hh>
//--------- PROBLEMS ------------------------
#include <dune/fem-dg/models/stokesprobleminterfaces.hh>
#include "problem.hh"
//--------- MODELS --------------------------
#include "models.hh"


//--------- PROBLEMCREATOR ------------------
#include "../poisson/problemcreator.hh"



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
      //todo: correction
      EOCErrorIDs ids = PoissonProblemCreatorType::template AnalyticalTraits< GridPart >::initEoc();
      ids[ 3 ] = Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-p-error" ) );
      return ids;
    }

    template< class SolutionImp, class Model, class ExactFunction, class SigmaEstimatorImp >
    static void addEOCErrors ( const EOCErrorIDs &ids, SolutionImp &u, const Model &model, const ExactFunction &f, SigmaEstimatorImp& sigma  )
    {
      typename PoissonProblemCreatorType::template AnalyticalTraits< GridPart > addEOCErrors( ids, u, model, f, sigma );
    }

    template< class SolutionPressureImp, class Model, class ExactPressureFunction >
    static void addEOCErrors ( const EOCErrorIDs &ids, SolutionPressureImp &p, const Model &model, const ExactPressureFunction &g  )
    {
      // calculate L2 - p-Norm
      Dune::Fem::L2Norm< GridPartType > l2pnorm( p.space().gridPart() );
      const double l2perror = l2pnorm.distance( g, p );
      Dune::Fem::FemEoc::setErrors( ids[ 1 ], l2perror );
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

    static inline std::string advectionFluxName()
    {
      return PoissonDiscreteTraits::advectionFluxName();
    }

    static inline std::string diffusionFluxName()
    {
      return PoissonDiscreteTraits::diffusionFluxName();
    }

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

    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;


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

    typedef typename AnalyticalTraitsType::ProblemType::ExactPressureType ExactPressureType;
    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactPressureType, GridPartType >  GridExactPressureType;

    typedef std::tuple< typename std::tuple_element<0,typename PoissonDiscreteTraits::IOTupleType>::type, GridExactSolutionType*,
                        DiscreteFunctionType*, GridExactPressureType* > IOTupleType;
  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::StokesAlgorithm< GridType, StokesProblemCreator<GridType>, PoissonProblemCreatorType, polOrd > Type;
  };


};

#define NEW_STEPPER_SELECTOR_USED
#endif // FEMHOWTO_HEATSTEPPER_HH
