#ifndef FEMHOWTO_HEATSTEPPER_HH
#define FEMHOWTO_HEATSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

// dune-fem includes
#include <dune/fem/io/parameter.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>

// local includes
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

// overload default stepper traits
#include <dune/fem-dg/algorithm/advectiondiffusionstepper.hh>

#include <dune/fem-dg/algorithm/gridinitializer.hh>

#include <dune/fem/operator/linear/spoperator.hh>


#include "problems/problem.hh"
#include "problems/problemQuasiHeatEqn.hh"
#include "problems/pulse.hh"
#include "problems/sin.hh"
#include "problems/deformationalflow.hh"
#include "models.hh"


// traits for the operator passes
template< class GridPart,
          int polOrd,
          class AnalyticalTraits,
          class DiscreteFunctionImp,
          class IndicatorFunctionImp,
          class AdaptationHandlerImp,
          class ExtraParameterTupleImp = std::tuple<> >
struct OperatorTraits
{
  typedef GridPart                                                     GridPartType;
  typedef typename GridPartType::GridType                              GridType;
  typedef AnalyticalTraits                                             AnalyticalTraitsType;


  typedef typename AnalyticalTraitsType::InitialDataType               InitialDataType;
  typedef typename AnalyticalTraitsType::ModelType                     ModelType ;
  typedef UpwindFlux< ModelType >                                      FluxType;
  //typedef LLFFlux< ModelType >                    FluxType;
  static const Dune::DGDiffusionFluxIdentifier PrimalDiffusionFluxId = Dune::method_general;


  static const int polynomialOrder = polOrd == -1 ? 0 : polOrd;

  typedef typename ModelType::Traits::FaceDomainType                   FaceDomainType;

  typedef DiscreteFunctionImp                                          DiscreteFunctionType;
  typedef DiscreteFunctionType                                         DestinationType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::FunctionSpaceType             FunctionSpaceType;

  static const int dimRange  = FunctionSpaceType::dimRange;
  static const int dimDomain = FunctionSpaceType::dimDomain ;

  typedef typename FunctionSpaceType :: DomainType                     DomainType;
  typedef typename FunctionSpaceType :: RangeType                      RangeType;
  typedef typename FunctionSpaceType :: JacobianRangeType              JacobianRangeType;
  typedef typename FunctionSpaceType :: RangeFieldType                 RangeFieldType ;
  typedef typename FunctionSpaceType :: DomainFieldType                DomainFieldType ;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >              VolumeQuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 1 >              FaceQuadratureType;

  typedef IndicatorFunctionImp                                         IndicatorType;

  typedef AdaptationHandlerImp                                         AdaptationHandlerType ;

  static const int limiterPolynomialOrder = polOrd == -1 ? 1 : polOrd;
  typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
};

template< class GridTypeImp >
struct AdvectionDiffusionProblemCreator
{

  typedef Dune::Fem::FunctionSpace< typename GridTypeImp::ctype, double, GridTypeImp::dimension, DIMRANGE> FunctionSpaceType;

  // define problem type here if interface should be avoided
    typedef Dune::EvolutionProblemInterface< FunctionSpaceType,false >      ProblemInterfaceType;

  template< class GridPart >
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef HeatEqnModel< GridPart, InitialDataType >           ModelType;
    typedef typename InitialDataType::TimeDependentFunctionType TimeDependentFunctionType;
  };

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridTypeImp>
  initializeGrid()
  {
    // use default implementation
    return Dune::Fem::DefaultGridInitializer< GridTypeImp >::initialize();
  }

  static ProblemInterfaceType* problem()
  {
    // choice of explicit or implicit ode solver
    static const std::string probString[]  = { "heat" ,"quasi", "pulse", "sin" };
    const int probNr = Dune::Fem::Parameter::getEnum( "problem", probString, 0 );
    if( probNr == 0 )
      return new Dune :: U0< GridTypeImp, DIMRANGE > ();
    else if ( probNr == 1 )
      return new Dune :: QuasiHeatEqnSolution< GridTypeImp, DIMRANGE > ();
    else if ( probNr == 2 )
      return new Dune :: Pulse< GridTypeImp, DIMRANGE > ();
    else if ( probNr == 3 )
      return new Dune :: U0Sin< GridTypeImp, DIMRANGE > ();
    else
    {
      abort();
      return 0;
    }
  }


  //Stepper Traits
  template< class GridPart, int polOrd >
  struct DiscreteTraits
  {
public:

    typedef GridPart                                                      GridPartType;
    typedef typename GridPartType::GridType                               GridType;
    typedef GridPartType                                                  HostGridPartType;

    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static inline std::string advectionFluxName()
    {
      return "LLF";
    }

    static inline std::string diffusionFluxName()
    {
#ifdef EULER
      return "";
#elif (defined PRIMALDG)
      return Dune::Fem::Parameter::getValue< std::string >("dgdiffusionflux.method");
#else
      return "LDG";
#endif
    }

    static const int quadOrder = polynomialOrder * 3 + 1;

#ifdef CONTINOUSFUNCTIONS
    typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polynomialOrder, Dune::Fem::CachingStorage > DiscreteFunctionSpaceType;
#elif defined USE_MINIELEMENT
    typedef Dune::Fem::MiniElementDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polynomialOrder, Dune::Fem::CachingStorage > DiscreteFunctionSpaceType;
#else
    typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polynomialOrder, Dune::Fem::CachingStorage > DiscreteFunctionSpaceType;
#endif

#ifdef ISTLVEC
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
    typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > JacobianOperatorType;
#elif defined PETSCVEC
    typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
    typedef Dune::Fem::PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType > JacobianOperatorType;
#else
    typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > JacobianOperatorType;
#endif


#ifdef CONTINOUSFUNCTIONS
    typedef Dune::Fem::LagrangeInterpolation< typename AnalyticalTraitsType::TimeDependentFunctionType, DiscreteFunctionType > InitialProjectorType;
#else
    typedef Dune::Fem::L2Projection< typename AnalyticalTraitsType::TimeDependentFunctionType, DiscreteFunctionType > InitialProjectorType;
#endif

    typedef std::tuple<> ExtraParameterTuple;


private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
public:
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;



    typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                           OdeSolverType;
    // type of restriction/prolongation projection for adaptive simulations
    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                     RestrictionProlongationType;
    // type of linear solver for implicit ode
    typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >      BasicLinearSolverType;

    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                        AdaptationHandlerType;


    // --------- Operators using PASSES --------------------------
    //============================================================
    typedef OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                            DiscreteFunctionType, IndicatorType,
                            AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;

private:
#ifdef LIMITER
    #if (not defined EULER) and (defined FLUXDG)
    #warning "DGAdvectionDiffusionOperator: using LIMITER."
      typedef Dune :: DGLimitedAdvectionDiffusionOperator< OperatorTraitsType >  DgType;
    #else
    #warning "DGAdvectionDiffusionOperator: LIMITER can NOT be used. Not supported -> LIMITER, no EULER, no FLUXDG."
      typedef Dune :: DGAdvectionDiffusionOperator< OperatorTraitsType >  DgType;
    #endif
    #ifndef HIGHER_ORDER_FV
    #warning "DGAdvectionOperator: using LIMITER."
      typedef Dune :: DGLimitedAdvectionOperator< OperatorTraitsType >    DgAdvectionType;
    #else
    #warning "DGAdvectionOperator: using HIGHER ORDER FV."
      typedef AnalyticalTraitsType< GridPartType, -1 >                    FVOperatorTraitsType;
      typedef Dune :: DGLimitedAdvectionOperator< FVOperatorTraitsType >  DgAdvectionType;
    #endif
#else // no LIMITER
#warning "No limiter is applied to the numerical solution !!"
    typedef Dune :: DGAdvectionDiffusionOperator< OperatorTraitsType >    DgType;
    typedef Dune :: DGAdvectionOperator< OperatorTraitsType >             DgAdvectionType;
#endif
    typedef Dune :: DGDiffusionOperator< OperatorTraitsType >             DgDiffusionType;
    // default is that both are enabled
    template < bool advection, bool diffusion >
    struct OperatorChooser
    {
      typedef DgType          FullOperatorType;
      typedef DgDiffusionType ImplicitOperatorType;
      typedef DgAdvectionType ExplicitOperatorType;
    };
    // advection operator only, i.e. linear advection equation
    template < bool advection >
    struct OperatorChooser< advection, false >
    {
      typedef DgAdvectionType  FullOperatorType;
      typedef FullOperatorType ImplicitOperatorType;
      typedef FullOperatorType ExplicitOperatorType;
    };
    // diffusion operator only, i.e. for the heat equation
    template < bool diffusion >
    struct OperatorChooser< false, diffusion >
    {
      typedef DgDiffusionType  FullOperatorType;
      typedef FullOperatorType ImplicitOperatorType;
      typedef FullOperatorType ExplicitOperatorType;
    };

    static const bool advection = AnalyticalTraitsType::ModelType::hasAdvection;
    static const bool diffusion = AnalyticalTraitsType::ModelType::hasDiffusion;
    typedef OperatorChooser< advection, diffusion > OperatorChooserType;

public:
    typedef typename OperatorChooserType :: FullOperatorType      FullOperatorType;
    typedef typename OperatorChooserType :: ImplicitOperatorType  ImplicitOperatorType;
    typedef typename OperatorChooserType :: ExplicitOperatorType  ExplicitOperatorType;
    // --------- Operators using PASSES --------------------------
    //============================================================
  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::AdvectionDiffusionStepper< GridTypeImp, AdvectionDiffusionProblemCreator<GridTypeImp>, polOrd > Type;
  };


};

#ifndef COMBINED_PROBLEM_CREATOR
#define ProblemCreator AdvectionDiffusionProblemCreator
#endif

#define NEW_STEPPER_SELECTOR_USED
#endif // FEMHOWTO_HEATSTEPPER_HH
