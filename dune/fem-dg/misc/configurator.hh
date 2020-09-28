#ifndef FEMDG_PROBLEMCREATOR_CONFIGURATOR_HH
#define FEMDG_PROBLEMCREATOR_CONFIGURATOR_HH

#include "algorithmcreatorselector.hh"

//default values
#define GALERKINENUMDEFAULT 2
#define ADAPTIVITYENUMDEFAULT 2
#define DISCRETEFUNCTIONSPACESENUMDEFAULT 1
#define SOLVERENUMDEFAULT 1
#define ADVECTIONLIMITERENUMDEFAULT 1
#define MATRIXENUMDEFAULT 1
#define ADVECTIONFLUXENUMDEFAULT 1
#define DIFFUSIONFLUXENUMDEFAULT 1

//Default Preprocessor flags (if not specified)
#ifndef GALERKINENUM
#define GALERKINENUM GALERKINENUMDEFAULT
#endif

#ifndef ADAPTIVITYENUM
#define ADAPTIVITYENUM ADAPTIVITYENUMDEFAULT
#endif


#ifndef DISCRETEFUNCTIONSPACESENUM
#define DISCRETEFUNCTIONSPACESENUM DISCRETEFUNCTIONSPACESENUMDEFAULT
#endif

#ifndef SOLVERENUM
#define SOLVERENUM SOLVERENUMDEFAULT
#endif

#ifndef ADVECTIONLIMITERENUM
#define ADVECTIONLIMITERENUM ADVECTIONLIMITERENUMDEFAULT
#endif

#ifndef MATRIXENUM
#define MATRIXENUM MATRIXENUMDEFAULT
#endif

#ifndef ADVECTIONFLUXENUM
#define ADVECTIONFLUXENUM ADVECTIONFLUXENUMDEFAULT
#endif

#ifndef DIFFUSIONFLUXENUM
#define DIFFUSIONFLUXENUM DIFFUEIONSFLUXENUMDEFAULT
#endif

namespace Dune
{
namespace Fem
{

  static const Galerkin::Enum galerkinEnum =
#if GALERKINENUM == 1
  Galerkin::Enum::cg;
#else //GALERKINENUM == 1
  Galerkin::Enum::dg;
#endif

  static const Adaptivity::Enum adaptivityEnum =
#if ADAPTIVITYENUM == 1
  Adaptivity::Enum::no;
#else //ADAPTIVITYENUM == 1
  Adaptivity::Enum::yes;
#endif

  static const DiscreteFunctionSpaces::Enum dfSpaceEnum =
#if DISCRETEFUNCTIONSPACESENUM == 2
  DiscreteFunctionSpaces::Enum::legendre;
#elif DISCRETEFUNCTIONSPACESENUM == 3
  DiscreteFunctionSpaces::Enum::hierarchic_legendre;
#elif DISCRETEFUNCTIONSPACESENUM == 4
  DiscreteFunctionSpaces::Enum::lagrange;
#elif DISCRETEFUNCTIONSPACESENUM == 5
  DiscreteFunctionSpaces::Enum::padaptive;
#elif DISCRETEFUNCTIONSPACESENUM == 6
  DiscreteFunctionSpaces::Enum::gausslobatto;
#elif DISCRETEFUNCTIONSPACESENUM == 7
  DiscreteFunctionSpaces::Enum::gausslegendre;
#else // DISCRETEFUNCTIONSPACESENUM == 1
  DiscreteFunctionSpaces::Enum::orthonormal;
#endif

  static const Solver::Enum solverEnum =
#if SOLVERENUM == 2
  Solver::Enum::femoem;
#elif SOLVERENUM == 3
  Solver::Enum::istl;
#elif SOLVERENUM == 4
  Solver::Enum::umfpack;
#elif SOLVERENUM == 5
  Solver::Enum::petsc;
#elif SOLVERENUM == 6
  Solver::Enum::eigen;
#else //SOLVERENUM == 1
  Solver::Enum::fem;
#endif

  static const AdvectionLimiter::Enum advLimiterEnum =
#if ADVECTIONLIMITERENUM == 2
  AdvectionLimiter::Enum::limited;
#else //ADVECTIONLIMITERENUM == 1
  AdvectionLimiter::Enum::unlimited;
#endif

  static const Matrix::Enum matrixEnum =
#if MATRIXENUM == 2
  Matrix::Enum::assembled;
#else //MATRIXENUM == 1
  Matrix::Enum::matrixfree;
#endif

  //TODO insert all the available fluxes...
  static const AdvectionFlux::Enum advFluxEnum  = AdvectionFlux::Enum::general;
  static const DiffusionFlux::Enum diffFluxEnum = DiffusionFlux::Enum::primal;



  /**
   * \brief Convenience class for AlgorithmCreator classes.
   *
   * The module dune-fem-dg needs a lot of templates and typedefs which might confuse some users.
   * Furthermore, redundancy of templates expressions which are not updated correctly
   * may lead to nasty compiler errors.
   *
   * This class protects the user from this situation and enables all available implementations.
   *
   * \tparam GridImp type of the grid
   * \tparam Galerkin::Enum enum of the galerkin type
   * \tparam Adaptivity::Enum enum indicating whether we have adaptivity or not
   * \tparam DiscreteFunctionSpaces::Enum enum defining a discrete function space
   * \tparam Solver::Enum enum defining a solver
   * \tparam AdvectionLimiter::Enum enum defining the limiting of the advection operator
   * \tparam Matrix::Enum enum describing whether to assemble or not
   * \tparam AdvectionFlux::Enum enum describing the chosen numerical advection flux,
     \tparam DiffusionFlux::Enum enum describing the chosen primal numerical diffusion flux
   */
  template< class GridImp,
            Galerkin::Enum dgEnumId,
            Adaptivity::Enum adapEnumId,
            DiscreteFunctionSpaces::Enum spaceEnumId,
            Solver::Enum solverEnumId,
            AdvectionLimiter::Enum advLimitEnumId,
            Matrix::Enum matrixEnumId,
            AdvectionFlux::Enum advFluxEnumId,
            DiffusionFlux::Enum diffFluxEnumId>
  class AlgorithmConfigurator
  {
    static const Galerkin::Enum dgId                 =(dgEnumId == Galerkin::Enum::default_)?                   galerkinEnum   : dgEnumId;
    static const Adaptivity::Enum adapId             =(adapEnumId == Adaptivity::Enum::default_)?               adaptivityEnum : adapEnumId;
    static const DiscreteFunctionSpaces::Enum spaceId=(spaceEnumId == DiscreteFunctionSpaces::Enum::default_ )? dfSpaceEnum    : spaceEnumId;
    static const Solver::Enum solverId               =(solverEnumId == Solver::Enum::default_)?                 solverEnum     : solverEnumId;
    static const AdvectionLimiter::Enum advLimitId   =(advLimitEnumId == AdvectionLimiter::Enum::default_)?     advLimiterEnum : advLimitEnumId;
    static const Matrix::Enum matrixId               =(matrixEnumId == Matrix::Enum::default_)?                 matrixEnum     : matrixEnumId;
    static const AdvectionFlux::Enum advFluxId       =(advFluxEnumId == AdvectionFlux::Enum::default_)?         advFluxEnum    : advFluxEnumId;
    static const DiffusionFlux::Enum diffFluxId      =(diffFluxEnumId == DiffusionFlux::Enum::default_)?        diffFluxEnum   : diffFluxEnumId;


    template< bool symmetric=false>
    using SolverSelect = typename SolverSelector< solverId, symmetric, Matrix::Enum::matrixfree == matrixId ? true : false >::type;


  public:
    typedef GridImp                                           GridType;

    // select formulation depending on chosen diffusion flux
    static const Formulation::Enum formId = FormulationSelector< diffFluxId > :: formId ;

    using GridParts = typename GridPartSelector< GridType, dgId, adapId >::type;

    template< int dimRange >
    using FunctionSpaces = Fem::FunctionSpace< typename GridType::ctype, typename GridType::ctype, GridType::dimension, dimRange >;

    template< class FunctionSpaceImp, int polOrd, class GridPartImp = GridParts >
    using DiscreteFunctionSpaces = typename DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrd, spaceId, dgId >::type;

    template< class FunctionSpaceImp, int polOrd, class GridPartImp = GridParts >
    using DiscreteFunctions = typename SolverSelect<>::template DiscreteFunctionType<DiscreteFunctionSpaces<FunctionSpaceImp,polOrd,GridPartImp> >;

  public:
    template<class DF>
    struct FunctionChooser
    {
      typedef DiscreteFunctionSpaces< typename DF::DiscreteFunctionSpaceType::FunctionSpaceType, DF::DiscreteFunctionSpaceType::polynomialOrder > type;
    };
  public:

    template< template<class,template<class> class > class DFSelector, class DF >
    using WrappedDiscreteFunctions = typename SolverSelect<>::template DiscreteFunctionType< typename DFSelector<DF,FunctionChooser >::type >;

    template< class ModelImp, AdvectionFlux::Enum id = advFluxId >
    using AdvectionFluxes = DGAdvectionFlux< ModelImp, id >;

    template< class ModelImp, class DF, DiffusionFlux::Enum id = diffFluxId >
    using DiffusionFluxes = typename DiffusionFluxSelector< ModelImp, typename DF::DiscreteFunctionSpaceType, id, formId >::type;

    template< class ModelImp,
              class DFunctionSpaceImp,
              int dPolOrd,
              class RFunctionSpaceImp = DFunctionSpaceImp,
              int rPolOrd = dPolOrd,
              class ExtraParameterTupleImp = std::tuple<>,
              AdvectionFlux::Enum advId = advFluxId,
              DiffusionFlux::Enum diffId = diffFluxId >
    using DefaultAssembTraits = DefaultAssemblerTraits< ModelImp,
                                                        typename SolverSelect<>::template LinearOperatorType<DiscreteFunctionSpaces< DFunctionSpaceImp, dPolOrd >, DiscreteFunctionSpaces< RFunctionSpaceImp, rPolOrd > >,
                                                        AdvectionFluxes< ModelImp, advId >,
                                                        DiffusionFluxes< ModelImp, DiscreteFunctions< DFunctionSpaceImp, dPolOrd >, diffId >,
                                                        DiscreteFunctions< DFunctionSpaceImp, dPolOrd >,
                                                        DiscreteFunctions< RFunctionSpaceImp, rPolOrd >,
                                                        ExtraParameterTupleImp >;
    template< class ModelImp,
              class FunctionSpaceImp,
              int polOrd,
              class ExtraParameterTupleImp = std::tuple<>,
              class AdaptationIndicatorFunctionSpaceImp = typename DiscreteFunctionSpaces< FunctionSpaceImp, polOrd >::FunctionSpaceType,
              AdvectionFlux::Enum advId = advFluxId,
              DiffusionFlux::Enum diffId = diffFluxId>
    using DefaultOpTraits = DefaultOperatorTraits< ModelImp,
                                                   DiscreteFunctions< FunctionSpaceImp, polOrd >,
                                                   AdvectionFluxes< ModelImp, advId >,
                                                   DiffusionFluxes< ModelImp, DiscreteFunctions< FunctionSpaceImp, polOrd >, diffId >,
                                                   ExtraParameterTupleImp,
                                                   AdaptationIndicatorFunctionSpaceImp >;

    //Operator/Assembler
    template< class OpTraits, OperatorSplit::Enum opSplit = OperatorSplit::Enum::full, Matrix::Enum matrix = matrixId >
    using Operators = typename OperatorSelector< OpTraits, formId, advLimitId, opSplit, matrix >::type;

    //Matrix Containers
    template< class DDF, class RDF  >
    using Containers = typename SolverSelect<>::template LinearOperatorType<typename DDF::DiscreteFunctionSpaceType,typename RDF::DiscreteFunctionSpaceType>;

    //Solver
    template< class DF, bool symmetric = false >
    using LinearSolvers = typename SolverSelect<symmetric>::template LinearInverseOperatorType<typename DF::DiscreteFunctionSpaceType>;


    //Extra Parameters
    template< class ModelImp, class AC, int... polOrds >
    using ExtraParameters = typename ExtraParameterSelector< typename ModelImp::ParameterSpacesType, AC, polOrds... >::type;

    ////Real Jacobian for Runge-Kutta
    //template< class OperatorImp >
    //using JacobianOperators = DGHelmholtzOperator< OperatorImp >; //i.e. JacobianImplicitOperators

    //template< class OperatorImp >
    //using JacobianOperators = DGHelmholtzOperator< DGPrimalMatrixAssembly >; //i.e. JacobianImplicitOperators

    //general, rungekutta, rungekutta_pardg, rungekutta_assembled
    //template< class DiscreteFunctionImp >
    //using OdeSolvers = DuneODE::OdeSolverInterface< DiscreteFunctionImp >;

    //using OdeSolvers = RungeKuttaSolver< OperatorImp, ExplicitOperatorImp, ImplicitOperatorImp, BasicLinearSolverImp >;
    //using OdeSolvers = RungeKuttaSolver< OperatorImp, ExplicitOperatorImp, ImplicitOperatorImp, JacobianImplicitOperatorImp, BasicLinearSolverImp >;

  };



}
}
#endif
