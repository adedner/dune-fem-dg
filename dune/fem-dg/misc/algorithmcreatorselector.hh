#ifndef FEMDG_ALGORITHMCREATOR_SELECTOR_HH
#define FEMDG_ALGORITHMCREATOR_SELECTOR_HH

// iostream includes
#include <iostream>
#include <type_traits>

#include <dune/common/dynvector.hh>

#include <dune/fem-dg/misc/static_warning.hh>

// include gridpart
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include discrete function
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/vectorfunction/managedvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>

#if HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlinverseoperators.hh>
#endif

#if HAVE_UMFPACK || HAVE_SUITESPARSE_UMFPACK
//#define USE_UMFPACKSOLVER 1
#include <dune/fem/solver/umfpacksolver.hh>
#endif

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscinverseoperators.hh>
#endif

//include operators
#include <dune/fem-dg/operator/dg/primaloperator.hh>
#include <dune/fem-dg/operator/dg/fluxoperator.hh>

#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem-dg/operator/dg/operatortraits.hh>

#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/fluxes.hh>
#include <dune/fem-dg/operator/limiter/limiterutility.hh>

#include <dune/fem-dg/misc/algorithmcreatorenums.hh>


namespace Dune
{
namespace Fem
{
  /* Enums have been moved to dune/fem-dg/misc/algorithmcreatorenums.hh */

///////////////////////////////////////////////////////////////////////////
// GridPartSelector
///////////////////////////////////////////////////////////////////////////

  template <class Grid, Galerkin::Enum op, Adaptivity::Enum adap >
  struct GridPartSelector;

  template <class Grid>
  struct GridPartSelector<Grid, Galerkin::Enum::cg, Adaptivity::Enum::no >
  {
    //typedef Dune::Fem::LeafGridPart< Grid > type;
    typedef Dune::Fem::AdaptiveLeafGridPart< Grid, Dune::InteriorBorder_Partition > type;
  };

  template <class Grid>
  struct GridPartSelector<Grid, Galerkin::Enum::dg, Adaptivity::Enum::no>
  {
    //typedef Dune::Fem::LeafGridPart< Grid > type;
    typedef Dune::Fem::DGAdaptiveLeafGridPart< Grid > type;
  };

  template <class Grid>
  struct GridPartSelector<Grid, Galerkin::Enum::cg, Adaptivity::Enum::yes >
  {
    typedef Dune::Fem::AdaptiveLeafGridPart< Grid, Dune::InteriorBorder_Partition > type;
  };

  template <class Grid>
  struct GridPartSelector<Grid, Galerkin::Enum::dg, Adaptivity::Enum::yes>
  {
    typedef Dune::Fem::DGAdaptiveLeafGridPart< Grid > type;
  };

///////////////////////////////////////////////////////////////////////////
// ImplExplOperatorSelector
///////////////////////////////////////////////////////////////////////////

  // advection and diffusion or none
  template< class Op, class AdvectionOp, class DiffusionOp, bool advection, bool diffusion >
  struct ImplExplOperatorSelector
  {
    typedef Op                   FullOperatorType;
    typedef DiffusionOp          ImplicitOperatorType;
    typedef AdvectionOp          ExplicitOperatorType;
  };
  // advection only
  template< class Op, class AdvectionOp, class DiffusionOp >
  struct ImplExplOperatorSelector< Op, AdvectionOp, DiffusionOp, true, false >
  {
    typedef AdvectionOp          FullOperatorType;
    typedef FullOperatorType     ImplicitOperatorType;
    typedef FullOperatorType     ExplicitOperatorType;
  };
  // diffusion only
  template<class Op, class AdvectionOp, class DiffusionOp>
  struct ImplExplOperatorSelector< Op, AdvectionOp, DiffusionOp, false, true >
  {
    typedef DiffusionOp          FullOperatorType;
    typedef FullOperatorType     ImplicitOperatorType;
    typedef FullOperatorType     ExplicitOperatorType;
  };

///////////////////////////////////////////////////////////////////////////
// AdvectionLimiterFunctionSelector
///////////////////////////////////////////////////////////////////////////
  // default is minmod
  template< class DomainFieldType, AdvectionLimiterFunction::Enum limiterFct >
  struct AdvectionLimiterFunctionSelector
  {
    typedef MinModLimiter< DomainFieldType > type;
  };

  template< class DomainFieldType >
  struct AdvectionLimiterFunctionSelector< DomainFieldType, AdvectionLimiterFunction::Enum::superbee >
  {
    typedef SuperBeeLimiter< DomainFieldType > type;
  };

  template< class DomainFieldType >
  struct AdvectionLimiterFunctionSelector< DomainFieldType, AdvectionLimiterFunction::Enum::vanleer >
  {
    typedef VanLeerLimiter< DomainFieldType > type;
  };

  template< class DomainFieldType >
  struct AdvectionLimiterFunctionSelector< DomainFieldType, AdvectionLimiterFunction::Enum::none >
  {
    typedef NoLimiter< DomainFieldType > type;
  };


///////////////////////////////////////////////////////////////////////////
// AdvectionDiffusionOperatorSelector
///////////////////////////////////////////////////////////////////////////

  template< class OperatorTraits, Formulation::Enum form, AdvectionLimiter::Enum op >
  class AdvectionDiffusionOperatorSelector;

  template< class OperatorTraits >
  class AdvectionDiffusionOperatorSelector< OperatorTraits, Formulation::Enum::primal, AdvectionLimiter::Enum::unlimited >
  {
    static const int advection = OperatorTraits::ModelType::hasAdvection;
    static const int diffusion = OperatorTraits::ModelType::hasDiffusion;
    typedef DGAdvectionDiffusionOperator< OperatorTraits >                   DgType;
    typedef DGAdvectionOperator< OperatorTraits >                            DgAdvectionType;
    typedef DGDiffusionOperator< OperatorTraits >                            DgDiffusionType;
    typedef ImplExplOperatorSelector< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                             ImplExplOperatorSelectorType;
  public:
    typedef typename ImplExplOperatorSelectorType::FullOperatorType          FullOperatorType;
    typedef typename ImplExplOperatorSelectorType::ImplicitOperatorType      ImplicitOperatorType;
    typedef typename ImplExplOperatorSelectorType::ExplicitOperatorType      ExplicitOperatorType;

  };

  template< class OperatorTraits >
  class AdvectionDiffusionOperatorSelector< OperatorTraits, Formulation::Enum::primal, AdvectionLimiter::Enum::limited >
  {
    static const int advection = OperatorTraits::ModelType::hasAdvection;
    static const int diffusion = OperatorTraits::ModelType::hasDiffusion;
    typedef DGLimitedAdvectionDiffusionOperator< OperatorTraits >            DgType;
    typedef DGLimitedAdvectionOperator< OperatorTraits >                     DgAdvectionType;
    typedef DGDiffusionOperator< OperatorTraits >                            DgDiffusionType;
    typedef ImplExplOperatorSelector< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                             ImplExplOperatorSelectorType;
  public:
    typedef typename ImplExplOperatorSelectorType::FullOperatorType          FullOperatorType;
    typedef typename ImplExplOperatorSelectorType::ImplicitOperatorType      ImplicitOperatorType;
    typedef typename ImplExplOperatorSelectorType::ExplicitOperatorType      ExplicitOperatorType;
  };

  template< class OperatorTraits >
  class AdvectionDiffusionOperatorSelector< OperatorTraits, Formulation::Enum::primal, AdvectionLimiter::Enum::scalinglimited >
   : public AdvectionDiffusionOperatorSelector< OperatorTraits, Formulation::Enum::primal, AdvectionLimiter::Enum::limited >
  {
  };

  template< class OperatorTraits >
  class AdvectionDiffusionOperatorSelector< OperatorTraits, Formulation::Enum::local, AdvectionLimiter::Enum::unlimited >
  {
    static const int advection = OperatorTraits::ModelType::hasAdvection;
    static const int diffusion = OperatorTraits::ModelType::hasDiffusion;
    typedef LDGAdvectionDiffusionOperator< OperatorTraits >                  DgType;
    typedef LDGAdvectionDiffusionOperator< OperatorTraits >                  DgAdvectionType;
    typedef LDGAdvectionDiffusionOperator< OperatorTraits >                  DgDiffusionType;
    typedef ImplExplOperatorSelector< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                             ImplExplOperatorSelectorType;
  public:
    typedef typename ImplExplOperatorSelectorType::FullOperatorType          FullOperatorType;
    typedef typename ImplExplOperatorSelectorType::ImplicitOperatorType      ImplicitOperatorType;
    typedef typename ImplExplOperatorSelectorType::ExplicitOperatorType      ExplicitOperatorType;

  };

  template< class OperatorTraits >
  class AdvectionDiffusionOperatorSelector< OperatorTraits, Formulation::Enum::local, AdvectionLimiter::Enum::limited >
  {
    static const int advection = OperatorTraits::ModelType::hasAdvection;
    static const int diffusion = OperatorTraits::ModelType::hasDiffusion;
    typedef LDGLimitedAdvectionDiffusionOperator< OperatorTraits >           DgType;
    typedef LDGLimitedAdvectionDiffusionOperator< OperatorTraits >           DgAdvectionType;
    typedef LDGLimitedAdvectionDiffusionOperator< OperatorTraits >           DgDiffusionType;
    typedef ImplExplOperatorSelector< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                             ImplExplOperatorSelectorType;
  public:
    typedef typename ImplExplOperatorSelectorType::FullOperatorType          FullOperatorType;
    typedef typename ImplExplOperatorSelectorType::ImplicitOperatorType      ImplicitOperatorType;
    typedef typename ImplExplOperatorSelectorType::ExplicitOperatorType      ExplicitOperatorType;
  };


  template< class OperatorTraits, Formulation::Enum form, AdvectionLimiter::Enum op, OperatorSplit::Enum split, Matrix::Enum ass >
  struct OperatorSelector;

  //matrixfree
  template< class OperatorTraits, Formulation::Enum form, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, form, op, OperatorSplit::Enum::expl, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, form, op >::ExplicitOperatorType type;
  };

  template< class OperatorTraits, Formulation::Enum form, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, form, op, OperatorSplit::Enum::impl, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, form, op >::ImplicitOperatorType type;
  };

  template< class OperatorTraits, Formulation::Enum form, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, form, op, OperatorSplit::Enum::full, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, form, op >::FullOperatorType type;
  };

  template< class OperatorTraits, Formulation::Enum form, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, form, op, OperatorSplit::Enum::rhs, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, form, op >::FullOperatorType type;
  };

  //assembled
  //TODO improve DGPrimalMatrixAssembly for correct splitting
  template< class AssemblerTraitsImp, AdvectionLimiter::Enum op >
  struct OperatorSelector< AssemblerTraitsImp, Formulation::Enum::primal, op, OperatorSplit::Enum::full, Matrix::Enum::assembled >
  {
    typedef DGPrimalMatrixAssembly< AssemblerTraitsImp > type;
  };

  template< class AssemblerTraitsImp, AdvectionLimiter::Enum op >
  struct OperatorSelector< AssemblerTraitsImp, Formulation::Enum::primal, op, OperatorSplit::Enum::rhs, Matrix::Enum::assembled >
  {
    typedef DGPrimalMatrixAssembly< AssemblerTraitsImp > type;
  };

///////////////////////////////////////////////////////////////////////////
// MatrixSolverSelector
///////////////////////////////////////////////////////////////////////////

  template< Solver::Enum solver, bool dummy = false >
  struct AvailableSolvers
  {};

  template <Solver::Enum solver, bool symmetric>
  struct MatrixSolverSelector
  {
    static const bool solverConfigured = false; // this implementation is used for not installed packages
    // choose type of discrete function, Matrix implementation and solver implementation
    // this should work with any discrete function implementation
    typedef Dune::DynamicVector<double>                              DofVectorType;
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction<Space,DofVectorType> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType = Dune::Fem::CgInverseOperator< DiscreteFunctionType<DSpace > >;
  };

  template <bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::fem,symmetric>
  {
    static const bool solverConfigured = true;
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction<Space>;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType
      = typename std::conditional<symmetric,
                                  Dune::Fem::CgInverseOperator< DiscreteFunctionType<DSpace> >,
                                  Dune::Fem::KrylovInverseOperator< DiscreteFunctionType<DSpace> > >::type;
  };

#if HAVE_DUNE_ISTL
  template <bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::istl,symmetric>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    // here we need the special ISTLBlockVectorDiscreteFunction
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::ISTLBlockVectorDiscreteFunction<Space>;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::ISTLLinearOperator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType =
    typename std::conditional<symmetric,
                              Dune::Fem::ISTLInverseOperator< DiscreteFunctionType<DSpace>, SolverParameter::cg>,
                              Dune::Fem::ISTLInverseOperator< DiscreteFunctionType<DSpace>, SolverParameter::bicgstab>
                             >::type;
  };
#else
  template< bool dummy >
  struct AvailableSolvers< Solver::Enum::istl, dummy >
  {
    static_warning(false, "You have chosen the istl solver backend which is currently not installed. Falling back to standard solver!");
  };
#endif // HAVE_ISTL

#if USE_UMFPACKSOLVER
  template <bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::umfpack,symmetric>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction<Space>;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType = Dune::Fem::UMFPACKOp< DiscreteFunctionType<DSpace>, LinearOperatorType<DSpace,RSpace>, symmetric >;
  };
#else
  //template< bool dummy >
  //struct AvailableSolvers< Solver::Enum::umfpack, dummy >
  //{
  //  static_warning(false, "You have chosen the UMFPACK solver backend which is currently not installed. Falling back to standard solver!");
  //};
#endif
#undef USE_UMFPACKSOLVER

#if HAVE_PETSC
  template <bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::petsc,symmetric>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::PetscDiscreteFunction<Space>;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::PetscLinearOperator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType = Dune::Fem::PetscInverseOperator< DiscreteFunctionType<DSpace> >;
    // to switch between solvers for symmetric and non symmetric operators
    // use the parameter petsc.kspsolver.method
  };
#else
  //template< bool dummy >
  //struct AvailableSolvers< Solver::Enum::petsc, dummy >
  //{
  //  static_warning(false, "You have chosen the PetSc solver backend which is currently not installed. Falling back to standard solver!");
  //};
#endif

  ///////////////////////////////////////////////////////////////////////////
// MatrixFreeSolverSelector
///////////////////////////////////////////////////////////////////////////

  //template <Solver::Enum solver, bool symmetric>
  //struct MatrixFreeSolverSelector
  //{
  //  static const bool solverConfigured = false; // this implementation is used for not installed packages
  //  // choose type of discrete function, Matrix implementation and solver implementation
  //  // this should work with any discrete function implementation
  //  typedef Dune::Fem::DynamicVector<double>                                                                        DofVectorType;
  //  template<class Space>
  //  using DiscreteFunctionType = Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction<DSpace,DofVectorType > >;
  //  template<class DSpace, class RSpace = DSpace>
  //  using LinearOperatorType = Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
  //  template<class DSpace, class RSpace = DSpace>
  //  using LinearInverseOperatorType = Dune::Fem::CGInverseOperator< DiscreteFunctionType<DSpace> >;
  //};

  template <Solver::Enum solver, bool symmetric >
  struct MatrixFreeSolverSelector
  {
    static const bool solverConfigured = false;
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction<Space>;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::Operator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType
      = typename std::conditional<symmetric,
                                  Dune::Fem::CgInverseOperator< DiscreteFunctionType<DSpace> >,
                                  Dune::Fem::KrylovInverseOperator< DiscreteFunctionType<DSpace> > >::type;
  };

  template <bool symmetric>
  struct MatrixFreeSolverSelector<Solver::Enum::fem,symmetric>
  {
    static const bool solverConfigured = true;
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction< Space >;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::Operator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType = Dune::Fem::KrylovInverseOperator< DiscreteFunctionType<DSpace> >;
  };

#if HAVE_DUNE_ISTL
  template <bool symmetric>
  struct MatrixFreeSolverSelector<Solver::Enum::istl,symmetric>
  {
    static const bool solverConfigured = true;
    template<class Space>
    using DiscreteFunctionType = Dune::Fem::ISTLBlockVectorDiscreteFunction<Space>;
    template<class DSpace, class RSpace = DSpace>
    using LinearOperatorType = Dune::Fem::Operator< DiscreteFunctionType<DSpace>, DiscreteFunctionType<RSpace> >;
    template<class DSpace, class RSpace = DSpace>
    using LinearInverseOperatorType
      = typename std::conditional<symmetric,
                                  Dune::Fem::ISTLInverseOperator< DiscreteFunctionType<DSpace>, SolverParameter::cg>,
                                  Dune::Fem::ISTLInverseOperator< DiscreteFunctionType<DSpace>, SolverParameter::bicgstab>
                                 >::type;
  };
#endif

///////////////////////////////////////////////////////////////////////////
// SolverSelector
///////////////////////////////////////////////////////////////////////////

  template <Solver::Enum solver, bool symmetric, bool matrixfree>
  struct SolverSelector;

  template <Solver::Enum solver, bool symmetric >
  struct SolverSelector< solver, symmetric, true >
  {
    typedef MatrixFreeSolverSelector< solver, symmetric > type;
  };

  template <Solver::Enum solver, bool symmetric >
  struct SolverSelector< solver, symmetric, false >
  {
    typedef MatrixSolverSelector< solver, symmetric > type;
  };


///////////////////////////////////////////////////////////////////////////
// DiscreteFunctionSpaceSelector
///////////////////////////////////////////////////////////////////////////
#ifdef USE_BASEFUNCTIONSET_CODEGEN
  using Storage = Dune::Fem::CodegenStorage;
#else
  using Storage = Dune::Fem::CachingStorage;
#endif

  template< class FunctionSpaceImp, class GridPartImp, int polOrder, DiscreteFunctionSpaces::Enum dfType, Galerkin::Enum opType >
  struct DiscreteFunctionSpaceSelector;

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::lagrange, Galerkin::Enum::cg >
  {
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, Storage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::legendre, Galerkin::Enum::dg >
  {
    typedef LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, Storage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::orthonormal, Galerkin::Enum::dg >
  {
    typedef DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, Storage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::lagrange, Galerkin::Enum::dg >
  {
    typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, Storage > type;
  };


  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::hierarchic_legendre, Galerkin::Enum::dg >
  {
    typedef HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, Storage > type;
  };

#if HAVE_DUNE_LOCALFUNCTIONS
  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::gausslobatto, Galerkin::Enum::dg >
  {
    typedef FixedOrderDGLagrangeSpace< FunctionSpaceImp, GridPartImp, polOrder, Dune::GaussLobattoPointSet, Storage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::gausslegendre, Galerkin::Enum::dg >
  {
    typedef FixedOrderDGLagrangeSpace< FunctionSpaceImp, GridPartImp, polOrder, Dune::GaussLegendrePointSet, Storage > type;
  };
#endif

  template< class ModelImp, class DiscreteFunctionSpaceImp,
            DiffusionFlux::Enum diffFluxId, Formulation::Enum form >
  struct DiffusionFluxSelector;

  template< class ModelImp, class DiscreteFunctionSpaceImp, DiffusionFlux::Enum diffFluxId >
  struct DiffusionFluxSelector< ModelImp, DiscreteFunctionSpaceImp, diffFluxId, Formulation::Enum::local >
  {
    typedef DGLocalDiffusionFlux< DiscreteFunctionSpaceImp, ModelImp, diffFluxId > type;
  };

  template< class ModelImp, class DiscreteFunctionSpaceImp, DiffusionFlux::Enum diffFluxId >
  struct DiffusionFluxSelector< ModelImp, DiscreteFunctionSpaceImp, diffFluxId, Formulation::Enum::primal >
  {
    typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceImp, ModelImp, diffFluxId > type;
  };


/////////////////////////////////////////////////////////////////////////
// Extra Parameters
/////////////////////////////////////////////////////////////////////////

  namespace details
  {
    template< class ParameterSpacesImp, class AC, class polOrds >
    struct ExtraParameterSelectorImpl;

    template< class... ParameterSpaces, class... AC, int... polOrds >
    struct ExtraParameterSelectorImpl< std::tuple<ParameterSpaces...>, std::tuple<AC...>, std::tuple< std::integral_constant< int, polOrds> ...> >
    {
      typedef std::tuple< typename AC::template DiscreteFunctions< ParameterSpaces, polOrds >... > type;
    };
  }


  template< class ParameterSpacesImp, class AC, int... polOrds >
  class ExtraParameterSelector;

  template< class... ParameterSpaces, class... AC, int... polOrds >
  class ExtraParameterSelector< std::tuple<ParameterSpaces...>, std::tuple<AC...>, polOrds... >
  {
    static_assert( sizeof...(ParameterSpaces) <= sizeof...(AC), "We expect for each parameter space an algorithm configurator!" );
    static_assert( sizeof...(ParameterSpaces) <= sizeof...(polOrds), "We expect for each parameter space a polynomial order!" );

    typedef std::make_integer_sequence< int, sizeof...(ParameterSpaces) > SequenceType;

    typedef typename tuple_reducer< std::tuple<AC...>, SequenceType >::type RedAC;
    typedef typename tuple_reducer< std::tuple<std::integral_constant<int,polOrds>...>, SequenceType >::type RedPolOrds;

  public:
    typedef typename details::ExtraParameterSelectorImpl< std::tuple<ParameterSpaces...>, RedAC, RedPolOrds >::type type;
  };


} // end namespace Fem
} // end namespace Dune
#endif
