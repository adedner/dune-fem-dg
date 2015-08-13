#ifndef FEMDG_EULERPROBLEMCREATOR_HH
#define FEMDG_EULERPROBLEMCREATOR_HH
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
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/advectionstepper.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"

enum AdvectionDiffusionOperatorType
{
  _unlimited = 0,
  _limited = 1,
};


template< class Op, class DiffusionOp, class AdvectionOp, bool advection, bool diffusion >
struct OperatorChooser
{
  typedef Op                   FullOperatorType;
  typedef DiffusionOp          ImplicitOperatorType;
  typedef AdvectionOp          ExplicitOperatorType;
};
template< class Op, class DiffusionOp, class AdvectionOp, bool advection >
struct OperatorChooser< Op, DiffusionOp, AdvectionOp, advection, false >
{
  typedef AdvectionOp          FullOperatorType;
  typedef FullOperatorType     ImplicitOperatorType;
  typedef FullOperatorType     ExplicitOperatorType;
};
template<class Op, class DiffusionOp, class AdvectionOp, bool diffusion >
struct OperatorChooser< Op, DiffusionOp, AdvectionOp, false, diffusion >
{
  typedef DiffusionOp          FullOperatorType;
  typedef FullOperatorType     ImplicitOperatorType;
  typedef FullOperatorType     ExplicitOperatorType;
};



template< class OperatorTraits, bool advection, bool diffusion, AdvectionDiffusionOperatorType op >
class AdvectionDiffusionOperators;


template< class OperatorTraits, bool advection, bool diffusion >
class AdvectionDiffusionOperators< OperatorTraits, advection, diffusion, _unlimited >
{
  typedef Dune::DGAdvectionDiffusionOperator< OperatorTraits >             DgType;
  typedef Dune::DGAdvectionOperator< OperatorTraits >                      DgAdvectionType;
  typedef Dune::DGDiffusionOperator< OperatorTraits >                      DgDiffusionType;
  typedef OperatorChooser< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                           OperatorChooserType;
public:
  typedef typename OperatorChooserType::FullOperatorType                            FullOperatorType;
  typedef typename OperatorChooserType::ImplicitOperatorType                        ImplicitOperatorType;
  typedef typename OperatorChooserType::ExplicitOperatorType                        ExplicitOperatorType;
};

template< class OperatorTraits, bool advection, bool diffusion >
class AdvectionDiffusionOperators< OperatorTraits, advection, diffusion, _limited >
{
  typedef Dune::DGLimitedAdvectionDiffusionOperator< OperatorTraits >      DgType;
  typedef Dune::DGLimitedAdvectionOperator< OperatorTraits >               DgAdvectionType;
  typedef Dune::DGDiffusionOperator< OperatorTraits >                      DgDiffusionType;
  typedef OperatorChooser< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                           OperatorChooserType;
public:
  typedef typename OperatorChooserType::FullOperatorType                            FullOperatorType;
  typedef typename OperatorChooserType::ImplicitOperatorType                        ImplicitOperatorType;
  typedef typename OperatorChooserType::ExplicitOperatorType                        ExplicitOperatorType;
};


enum DiscreteFunctionSpacesType
{
  _lagrange = 0,
  _legendre = 1,
};


template< class FunctionSpaceImp, class GridPartImp, int polOrder, DiscreteFunctionSpacesType dfType, GalerkinType opType >
struct DiscreteFunctionSpaces;

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
struct DiscreteFunctionSpaces< FunctionSpaceImp, GridPartImp, polOrder, _lagrange, cg >
{
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, Dune::Fem::CachingStorage > type;
};

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
struct DiscreteFunctionSpaces< FunctionSpaceImp, GridPartImp, polOrder, _legendre, dg >
{
  typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, Dune::Fem::CachingStorage > type;
};



template< class TimeDependentFunctionImp, class DiscreteFunctionImp, GalerkinType gt >
struct InitialProjectors;

template< class TimeDependentFunctionImp, class DiscreteFunctionImp >
struct InitialProjectors< TimeDependentFunctionImp, DiscreteFunctionImp, cg >
{
  typedef Dune::Fem::LagrangeInterpolation< TimeDependentFunctionImp, DiscreteFunctionImp >                           type;
};

template< class TimeDependentFunctionImp, class DiscreteFunctionImp >
struct InitialProjectors< TimeDependentFunctionImp, DiscreteFunctionImp, dg >
{
  typedef Dune::Fem::L2Projection< TimeDependentFunctionImp, DiscreteFunctionImp >                                    type;
};


template< class DiscreteFunctionSpaceImp, SolverType solverType >
struct DiscreteFunctions;


template< class DiscreteFunctionSpaceImp >
struct DiscreteFunctions< DiscreteFunctionSpaceImp, fem >
{
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > type;
  typedef Dune::Fem::SparseRowLinearOperator< type, type >                jacobian;
};

#if HAVE_DUNE_ISTL
template< class DiscreteFunctionSpaceImp >
struct DiscreteFunctions< DiscreteFunctionSpaceImp, istl >
{
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >  type;
  typedef Dune::Fem::ISTLLinearOperator< type, type >                             jacobian;
};
#endif

#if HAVE_PETSC
template< class DiscreteFunctionSpaceImp >
struct DiscreteFunctions< DiscreteFunctionSpaceImp, petsc >
{
  typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceImp > type;
  typedef Dune::Fem::PetscLinearOperator< type, type >                 jacobian;
};
#endif



template< class GridImp >
struct EulerProblemCreator
{
  typedef GridImp                                         GridType;
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
  typedef HostGridPartType                                GridPartType;

  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

  // define problem type here if interface should be avoided
  typedef ProblemBase< GridType >  ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef Dune::Fem::EulerModel< GridPart, InitialDataType >  ModelType;
    typedef typename InitialDataType::TimeDependentFunctionType TimeDependentFunctionType;

    typedef std::vector< int >                                  EOCErrorIDs;

    static EOCErrorIDs initEoc ()
    {
      EOCErrorIDs ids;
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-error" ) ) );
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "$L^1$-error" ) ) );
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
        const double l2error = l2norm.distance( model.problem().fixedTimeFunction( tp.time() ), u );
        Dune::Fem::FemEoc::setErrors( ids[0], l2error );

        Dune::Fem::L1Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > l1norm( u.space().gridPart() );
        const double l1error = l1norm.distance( model.problem().fixedTimeFunction( tp.time() ), u );
        Dune::Fem::FemEoc::setErrors( ids[1], l1error );

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
    const std::string problemNames []
        = { "sod" , "withman", "withmansmooth", "smooth1d" , "ffs" , "diffraction" , "shockbubble", "p123" };

    const int problemNumber = Dune :: Fem :: Parameter :: getEnum ( "problem" , problemNames );

    if( problemNames[ problemNumber ] == "sod" )
    {
      return new U0Sod< GridType > ( );
    }
    else if( problemNames[ problemNumber ] == "smooth1d" )
    {
      return new U0Smooth1D< GridType > ();
    }
    else if( problemNames[ problemNumber ] == "ffs" )
    {
      return new U0FFS< GridType > ();
    }
    else if( problemNames[ problemNumber ] == "p123" )
    {
      return new U0P123< GridType >();
    }
    std::cerr << "Error: Problem " << problemNames[ problemNumber ]
              << " not implemented." << std::endl;

    // choice of explicit or implicit ode solver
    return new U0Smooth1D< GridType > ();
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
public:
    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static inline std::string advectionFluxName()
    {
      return "LLF";
    }

    static inline std::string diffusionFluxName()
    {
      return Dune::Fem::Parameter::getValue< std::string >("dgdiffusionflux.method");
    }

    static const int quadOrder = polynomialOrder * 3 + 1;
    static const SolverType solverType = fem ;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef typename InitialProjectors< typename AnalyticalTraitsType::TimeDependentFunctionType, DiscreteFunctionType, dg >::type   InitialProjectorType;

    typedef std::tuple<> ExtraParameterTuple;


    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;
    typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* > IOTupleType;

private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
public:
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;



    typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 OdeSolverType;
    // type of restriction/prolongation projection for adaptive simulations
    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                           RestrictionProlongationType;
    // type of linear solver for implicit ode
    typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >            BasicLinearSolverType;

    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;


    // --------- Operators using PASSES --------------------------
    //============================================================
    typedef LLFFlux< typename AnalyticalTraitsType::ModelType >                                 FluxType;
    //typedef HLLNumFlux< typename AnalyticalTraitsType::ModelType >                            FluxType;
    //typedef HLLCNumFlux< typename AnalyticalTraitsType::ModelType >                           FluxType;

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
    typedef Dune::Fem::AdvectionStepper< GridType, EulerProblemCreator<GridType>, polOrd > Type;
  };


};

#endif
