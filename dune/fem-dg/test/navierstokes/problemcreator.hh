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
//--------- HANDLER -------------------------
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>


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
  private:
    static const SolverType solverType = fem;
  public:
    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static const int quadOrder = polynomialOrder * 3 + 1;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef typename InitialProjectors< typename AnalyticalTraitsType::ProblemType::TimeDependentFunctionType, DiscreteFunctionType, dg >::type   InitialProjectorType;

    typedef std::tuple<> ExtraParameterTuple;

  private:
    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType                       ExactSolutionType;
  public:
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >                   GridExactSolutionType;
    typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

    typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 OdeSolverType;
    // type of linear solver for implicit ode
    typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >            BasicLinearSolverType;

  private:
    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;

    typedef LLFFlux< typename AnalyticalTraitsType::ModelType >                                 FluxType;

    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType, IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

    static const int hasAdvection = AnalyticalTraitsType::ModelType::hasAdvection;
    static const int hasDiffusion = AnalyticalTraitsType::ModelType::hasDiffusion;
  public:
    typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, _unlimited > AdvectionDiffusionOperatorType;

    typedef typename AdvectionDiffusionOperatorType::FullOperatorType                           FullOperatorType;
    typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                       ImplicitOperatorType;
    typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                       ExplicitOperatorType;

    //------HANDLER-----------------------------------------------------
    struct HandlerTraits
    {
      private:
      //adaptivity
      typedef Dune::DGAdaptationIndicatorOperator< OperatorTraitsType, hasAdvection, hasDiffusion >  IndicatorType;
      typedef Estimator< DiscreteFunctionType, typename  AnalyticalTraitsType::ProblemType >         GradientIndicatorType ;
      typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                              RestrictionProlongationType;
      //limiting
      typedef FullOperatorType                                                                       LimiterOperatorType;
      public:
      typedef Dune::Fem::DefaultDiagnosticsHandler                                                   DiagnosticsHandlerType;
      typedef Dune::Fem::DefaultSolverMonitorHandler                                                 SolverMonitorHandlerType;
      typedef Dune::Fem::DefaultCheckPointHandler< GridType >                                        CheckPointHandlerType;
      typedef Dune::Fem::DefaultDataWriterHandler< GridType, IOTupleType >                           DataWriterHandlerType;
      typedef Dune::Fem::NoAdditionalOutputHandler                                                   AdditionalOutputHandlerType;
      typedef Dune::Fem::DefaultSolutionLimiterHandler< LimiterOperatorType >                        SolutionLimiterHandlerType;
      typedef Dune::Fem::DefaultAdaptHandler< GridPartType, DiscreteFunctionType,
                                              RestrictionProlongationType, IndicatorType,
                                              GradientIndicatorType, SolutionLimiterHandlerType >    AdaptHandlerType;
    };


  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::AdvectionDiffusionStepper< GridType, NavierStokesProblemCreator<GridType>, polOrd > Type;
  };


};

#endif // FEMHOWTO_HEATSTEPPER_HH
