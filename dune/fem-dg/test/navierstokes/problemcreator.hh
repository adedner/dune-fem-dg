#ifndef FEMDG_NAVIERSTOKESSTEPPER_HH
#define FEMDG_NAVIERSTOKESSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM + 2
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

// misc
#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem/misc/l2norm.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/advectiondiffusionstepper.hh>
#include <dune/fem-dg/algorithm/advectionstepper.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>



template< class GridImp >
struct NavierStokesProblemCreator
{
  typedef GridImp                                         GridType;
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
  typedef HostGridPartType                                GridPartType;

  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

  // define problem type here if interface should be avoided
    typedef Dune::NSWaves< GridType >                           ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef NSModel< GridPart, InitialDataType >                ModelType;

    typedef std::vector< int >                                  EOCErrorIDs;

    static EOCErrorIDs initEoc ()
    {
      EOCErrorIDs ids;
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-error" ) ) );
      return ids;
    }

    template< class Solution, class Model, class ExactFunction, class TimeProvider >
    static void addEOCErrors ( const EOCErrorIDs &ids, TimeProvider& tp,
                               Solution &u, Model &model, ExactFunction &f )
    {
      if( model.problem().calculateEOC( tp, u, 0 ) )
      {
        const int order = 2*u.space().order()+4;
        Dune::Fem::L2Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > l2norm( u.space().gridPart(), order );
        const double error = l2norm.distance( model.problem().fixedTimeFunction( tp.time() ), u );
        Dune::Fem::FemEoc::setErrors( ids[ 0 ], error );
      }
    }
  };

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    // use default implementation
    return Dune::Fem::DefaultGridInitializer< GridType >::initialize();
  }

  static ProblemInterfaceType* problem()
  {
    return new ProblemInterfaceType();
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
public:
    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static const int quadOrder = polynomialOrder * 3 + 1;
    static const SolverType solverType = fem ;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef typename InitialProjectors< typename AnalyticalTraitsType::ProblemType::TimeDependentFunctionType, DiscreteFunctionType, dg >::type   InitialProjectorType;

    typedef std::tuple<> ExtraParameterTuple;

    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;
    typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* > IOTupleType;

    typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 OdeSolverType;
    // type of restriction/prolongation projection for adaptive simulations
    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                           RestrictionProlongationType;
    // type of linear solver for implicit ode
    typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >            BasicLinearSolverType;

    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;


    // --------- Operators using PASSES --------------------------
    //============================================================
    typedef LLFFlux< typename AnalyticalTraitsType::ModelType >                                 FluxType;

private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;
public:


    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType, IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

    // TODO: advection/diffusion should not be precribed by model
    static const int hasAdvection = AnalyticalTraitsType::ModelType::hasAdvection;
    static const int hasDiffusion = AnalyticalTraitsType::ModelType::hasDiffusion;

    typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;

    typedef typename AdvectionDiffusionOperatorType::FullOperatorType                           FullOperatorType;
    typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                       ImplicitOperatorType;
    typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                       ExplicitOperatorType;
    // --------- Operators using PASSES --------------------------
    //============================================================



    //---------Adaptivity ----------------------------------------------
    // advection = true , diffusion = true
    typedef Dune :: DGAdaptationIndicatorOperator< OperatorTraitsType, hasAdvection, hasDiffusion >  DGIndicatorType;
    // gradient estimator
    typedef Estimator< DiscreteFunctionType, typename  AnalyticalTraitsType::ProblemType >                   GradientIndicatorType ;
    typedef std::tuple< DGIndicatorType*, GradientIndicatorType* >             IndicatorTupleType;
    // --------Adaptivity ----------------------------------------------

    //------------- Limiter ---------------------------------------------
    typedef FullOperatorType                                                   LimiterOperatorType;
    //------------ Limiter ---------------------------------------------


  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::AdvectionDiffusionStepper< GridType, NavierStokesProblemCreator<GridType>, polOrd > Type;
  };


};

#endif // FEMHOWTO_HEATSTEPPER_HH
