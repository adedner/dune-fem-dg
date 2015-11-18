#ifndef FEMDG_PROBLEMCREATOR_SELECTOR_HH
#define FEMDG_PROBLEMCREATOR_SELECTOR_HH

#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>
#include <dune/fem/operator/projection/l2projection.hh>

namespace Dune
{
namespace Fem
{

  struct AdvectionLimiterIdentifier
  {
    typedef enum
    {
      unlimited = 0,
      limited = 1,
    } id;
  };

  struct DiscreteFunctionSpaceIdentifier
  {
    typedef enum
    {
      lagrange = 0,
      legendre = 1,
    } id;
  };


  struct GalerkinIdentifier
  {
    typedef enum
    {
      cg = 0,
      dg = 1
    } id;
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



  template< class OperatorTraits, bool advection, bool diffusion, AdvectionLimiterIdentifier::id op >
  class AdvectionDiffusionOperators;


  template< class OperatorTraits, bool advection, bool diffusion >
  class AdvectionDiffusionOperators< OperatorTraits, advection, diffusion, AdvectionLimiterIdentifier::unlimited >
  {
    typedef DGAdvectionDiffusionOperator< OperatorTraits >                   DgType;
    typedef DGAdvectionOperator< OperatorTraits >                            DgAdvectionType;
    typedef DGDiffusionOperator< OperatorTraits >                            DgDiffusionType;
    typedef OperatorChooser< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                             OperatorChooserType;
  public:
    typedef typename OperatorChooserType::FullOperatorType                            FullOperatorType;
    typedef typename OperatorChooserType::ImplicitOperatorType                        ImplicitOperatorType;
    typedef typename OperatorChooserType::ExplicitOperatorType                        ExplicitOperatorType;
  };

  template< class OperatorTraits, bool advection, bool diffusion >
  class AdvectionDiffusionOperators< OperatorTraits, advection, diffusion, AdvectionLimiterIdentifier::limited >
  {
    typedef DGLimitedAdvectionDiffusionOperator< OperatorTraits >            DgType;
    typedef DGLimitedAdvectionOperator< OperatorTraits >                     DgAdvectionType;
    typedef DGDiffusionOperator< OperatorTraits >                            DgDiffusionType;
    typedef OperatorChooser< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                             OperatorChooserType;
  public:
    typedef typename OperatorChooserType::FullOperatorType                            FullOperatorType;
    typedef typename OperatorChooserType::ImplicitOperatorType                        ImplicitOperatorType;
    typedef typename OperatorChooserType::ExplicitOperatorType                        ExplicitOperatorType;
  };



  template< class FunctionSpaceImp, class GridPartImp, int polOrder, DiscreteFunctionSpaceIdentifier::id dfType, GalerkinIdentifier::id opType >
  struct DiscreteFunctionSpaces;

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaces< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaceIdentifier::lagrange, GalerkinIdentifier::cg >
  {
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, CachingStorage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaces< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaceIdentifier::legendre, GalerkinIdentifier::dg >
  {
    typedef DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, CachingStorage > type;
  };


  struct NoImplementation { static const int implemented = false; };
  struct Implementation { static const int implemented = true; };

  template< class DiscreteFunctionSpaceImp, SolverType solverType >
  struct DiscreteFunctions : NoImplementation
  {};


  template< class DiscreteFunctionSpaceImp >
  struct DiscreteFunctions< DiscreteFunctionSpaceImp, fem > : Implementation
  {
    typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > type;
    typedef SparseRowLinearOperator< type, type >                jacobian;
  };

  #if HAVE_DUNE_ISTL
  template< class DiscreteFunctionSpaceImp >
  struct DiscreteFunctions< DiscreteFunctionSpaceImp, istl > : Implementation
  {
    typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >  type;
    typedef ISTLLinearOperator< type, type >                             jacobian;
  };
  #endif

  #if HAVE_PETSC
  template< class DiscreteFunctionSpaceImp >
  struct DiscreteFunctions< DiscreteFunctionSpaceImp, petsc > : Implementation
  {
    typedef PetscDiscreteFunction< DiscreteFunctionSpaceImp > type;
    typedef PetscLinearOperator< type, type >                 jacobian;
  };
  #endif

  template< class DiscreteFunctionSpaceImp, SolverType solverType >
  struct DiscreteFunctionSelector
  {
    typedef DiscreteFunctions< DiscreteFunctionSpaceImp, solverType > BaseType;
    static_assert( BaseType::implemented,
                   "\n==============================================================================="
                   "\nRequested implementation of DiscreteFunctions not available."
                   "\nPlease try one of the the following SolverTypes in your ProblemCreator: "
                   "\n"
                   "\n + matrixFree   use the matrix free version of the dune-fem solvers"
                   "\n + fem          use the matrix based version of the dune-fem solvers"
                   "\n + femoem       use the matrix based version of the dune-fem solvers with blas"
  #if HAVE_DUNE_ISTL
                   "\n + istl         use the dune-istl solvers"
  #else
                   "\n + istl         use the dune-istl solvers (WARNING: dune-ist is not installed!)"
  #endif
                   "\n + umfpack      use the direct solver umfpack"
  #if HAVE_PETSC
                   "\n + petsc        use the petsc package"
  #else
                   "\n + petsc        use the petsc package (WARNING: petsc is not installed!)"
  #endif
                   "\n===============================================================================" );

    typedef typename BaseType::type type;
    typedef typename BaseType::jacobian jacobian;
  };


//template< bool passBased >
//class AlgorithmConfigurator
//{
//
//
//
//};
//
//
//template< DiscreteFunctionSpacesType dtType,
//          GalerkinType opType
//          SolverType sType,
//          GalerkinType galerkin,
//          DGDiffusionFluxIdentifier diffId,
//          bool adaptivity,
//          int dimRange,
//
//  >
//class AlgorithmConfigurator
//{
//
//
//
//};
////DiffusionFlux 	BR, CDG, CDG2...
////AdvectionFlux 	LLF, HLLC...
////GalerkinType 	CG, DG
////Solver 	istl, fem...
////polOrder 	0, 1, 2...
////GridDim 	1,2,3...
////Adaptation 	no, yes (Estimator...)
////DimRange 	1,2,3..
////Assembly 	matrix-free/matrix-based
////Codegen 	yes/no
////Parallel 	yes/no

//  template< class GridImp,
//  >
//  struct AlgorithmConfigurator
//  {
//
//
//    typedef GridImp                                         GridType;
//    typedef typename GridPartChooser<GridType,dg>::Type     HostGridPartType;
//    typedef HostGridPartType                                GridPartType;
//
//    // define problem type here if interface should be avoided
//    typedef ProblemInterface< Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE > >
//                                                                  ProblemInterfaceType;
//
//    typedef typename ProblemInterfaceType::FunctionSpaceType      FunctionSpaceType;
//
//    struct AnalyticalTraits
//    {
//      typedef ProblemInterfaceType                                ProblemType;
//      typedef ProblemInterfaceType                                InitialDataType;
//      typedef PoissonModel< GridPartType, InitialDataType >       ModelType;
//
//      template< class Solution, class Model, class ExactFunction, class SigmaFunction>
//      static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
//      {
//        static L2EOCError l2EocError( "$L^2$-Error" );
//        l2EocError.add( u, f );
//        static DGEOCError dgEocError( "DG-Error" );
//        dgEocError.add( u, f );
//        static H1EOCError sigmaEocError( "sigma-norm" );
//        sigmaEocError.add( sigma, f );
//      }
//
//    };
//
//    static inline std::string moduleName() { return ""; }
//
//    static ProblemInterfaceType* problem()
//    {
//      int probNr = Fem::Parameter::getValue< int > ( "problem" );
//      return new PoissonProblem< GridType, DIMRANGE > ( probNr );
//    }
//
//
//    //Stepper Traits
//    template< int polOrd >
//    struct DiscreteTraits
//    {
//    private:
//      static const SolverType solverType = istl;
//      static const bool symmetricSolver = true;
//      typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, _legendre, dg >::type    DiscreteFunctionSpaceType;
//    public:
//      typedef typename DiscreteFunctionSelector< DiscreteFunctionSpaceType, solverType >::type                          DiscreteFunctionType;
//      typedef typename DiscreteFunctionSelector< DiscreteFunctionSpaceType, solverType >::jacobian                      JacobianOperatorType;
//
//      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;
//
//      typedef std::tuple<>                                                                        ExtraParameterTuple;
//
//      class Operator
//      {
//        friend DiscreteTraits;
//        friend SolverType;
//        typedef UpwindFlux< typename AnalyticalTraits::ModelType >                                 FluxType;
//
//        typedef DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
//                                       DiscreteFunctionType, FluxType, ExtraParameterTuple >        OperatorTraitsType;
//
//        typedef DGAdvectionDiffusionOperator< OperatorTraitsType >                                  AssemblyOperatorType;
//        typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>                     SolversType;
//      public:
//        typedef DGPrimalMatrixAssembly< AssemblyOperatorType >                                      AssemblerType;
//        typedef typename SolversType::LinearOperatorType                                            type;
//      };
//
//      struct Solver
//      {
//        typedef typename Operator::SolversType::LinearInverseOperatorType                           type;
//      };
//
//    private:
//      //typedef DGAdaptationIndicatorOperator< OperatorTraitsType, true, true >                       IndicatorType;
//      //typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >             GradientIndicatorType ;
//    public:
//
//      //typedef Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
//      typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                                SolverMonitorHandlerType;
//      typedef Fem::SubDiagnosticsHandler< Diagnostics >                                         DiagnosticsHandlerType;
//    };
//
//    template <int polOrd>
//    struct Stepper
//    {
//      // this should be ok but could lead to a henn-egg problem
//      typedef Fem::EllipticAlgorithm< GridType, SubPoissonProblemCreator, polOrd > Type;
//    };
//
//  };
//
}
}
#endif
