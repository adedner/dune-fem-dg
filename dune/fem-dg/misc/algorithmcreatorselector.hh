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
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/pardginverseoperators.hh>
#include <dune/fem/solver/oemsolver.hh>

#if HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#endif

#if HAVE_UMFPACK || HAVE_SUITESPARSE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#endif

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscsolver.hh>
#endif

#if HAVE_EIGEN
#include <dune/fem/storage/eigenvector.hh>
#include <dune/fem/operator/linear/eigenoperator.hh>
#include <dune/fem/solver/eigen.hh>
#endif

//include operators
#include <dune/fem-dg/operator/dg/primaloperator.hh>
#include <dune/fem-dg/operator/dg/fluxoperator.hh>

#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>

#include <dune/fem-dg/operator/dg/operatortraits.hh>

#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/fluxes.hh>




namespace Dune
{
namespace Fem
{

  /**
   *  \brief Namespace containing an Enum class to describe the limiting
   *  of advection operators.
   */
  namespace AdvectionLimiter
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! use default value
      default_,
      //! limitation of advection term
      unlimited,
      //! no limitation of advection term
      limited,
    };
  }

  /**
   *  \brief Namespace containing an Enum class to describe the discrete function
   *  and the local finite element.
   */
  namespace DiscreteFunctionSpaces
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! use default value
      default_,
      //! Discrete function space with Lagrange Finite Elements
      lagrange,
      //! Discrete function space with Legendre Finite Elements
      legendre,
      //! Discrete function space with hierarchic Legendre Finite Elements
      hierarchic_legendre,
      //! Discrete function space with hierarchic orthonormal monomial basis functions
      orthonormal
      //! p-adaptive space from dune-fem, implementing dg and lagrange
      // padaptive
    };
  }


  /**
   *  \brief Namespace containing an Enum class to describe the Galerkin type
   *  of the discretization scheme.
   */
  namespace Galerkin
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! use default value
      default_,
      //! Continuous Galerkin
      cg,
      //! Discontinuous Galerkin
      dg
    };
  }

  /**
   *  \brief Namespace containing an Enum class to describe whether adaptiv
   *  calculations should be possible or not.
   */
  namespace Adaptivity
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! use default value
      default_,
      //! no Adaptivity
      no,
      //! Allow Adaptivity
      yes
    };
  }

  /**
   *  \brief Namespace containing an Enum class to describe the solver backends.
   */
  namespace Solver
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! use default value
      default_,
      //! use the matrix based version of the dune-fem solvers
      fem,
      //! use the matrix based version of the dune-fem solvers with blas
      femoem,
      //! use the dune-istl solvers
      istl,
      //! use the direct solver umfpack
      umfpack,
      //! use the petsc package
      petsc,
      //! use the eigen package
      eigen
    };
  }

  /**
   *  \brief Namespace containing an Enum class to describe whether we want to
   *  assemble a matrix or not.
   */
  namespace Matrix
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! use default value
      default_,
      //! use matrix free operator
      matrixfree,
      //! use matrix based operator
      assembled
    };
  }

  /**
   *  \brief Namespace containing an Enum class to describe the explicit/implicit
   *  operator splitting.
   */
  namespace OperatorSplit
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! use default value
      default_,
      full,
      expl,
      impl,
      rhs
    };
  }

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
  class OperatorSelector;

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

  template <Solver::Enum solver, bool symmetric, class DomainDFSpace, class RangeDFSpace = DomainDFSpace>
  struct MatrixSolverSelector
  {
    static const bool solverConfigured = false; // this implementation is used for not installed packages
    // choose type of discrete function, Matrix implementation and solver implementation
    // this should work with any discrete function implementation
    typedef Dune::DynamicVector<double>                                                                        DofVectorType;
    typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DomainDFSpace, DofVectorType > > DomainDiscreteFunctionType;
    typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< RangeDFSpace, DofVectorType > >  RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                              DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType >             LinearOperatorType;
    typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType >                                                    LinearInverseOperatorType;
  };

  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::fem,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    typedef Dune::Fem::AdaptiveDiscreteFunction< DomainDFSpace >                                        DomainDiscreteFunctionType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< RangeDFSpace >                                         RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                  DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename std::conditional<symmetric,
            Dune::Fem::CGInverseOperator< DiscreteFunctionType >,
            Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType > > :: type          LinearInverseOperatorType;
  };

  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::femoem,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    // this work with a discrete function implementation based on a double* dof storage
    typedef Dune::Fem::AdaptiveDiscreteFunction< DomainDFSpace >                                        DomainDiscreteFunctionType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< RangeDFSpace >                                         RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                  DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename std::conditional<symmetric,
            Dune::Fem::OEMCGOp< DiscreteFunctionType, LinearOperatorType >,
            Dune::Fem::OEMBICGSTABOp< DiscreteFunctionType, LinearOperatorType > > :: type              LinearInverseOperatorType;
  };

#if HAVE_DUNE_ISTL
  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::istl,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    // here we need the special ISTLBlockVectorDiscreteFunction
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DomainDFSpace >                            DomainDiscreteFunctionType;
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< RangeDFSpace >                             RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                             DiscreteFunctionType;
    typedef Dune::Fem::ISTLLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename std::conditional<symmetric,
            Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType >,
            Dune::Fem::ISTLBICGSTABOp< DiscreteFunctionType, LinearOperatorType > > :: type        LinearInverseOperatorType;
  };
#else
  template< bool dummy >
  struct AvailableSolvers< Solver::Enum::istl, dummy >
  {
    static_warning(false, "You have chosen the istl solver backend which is currently not installed. Falling back to standard solver!");
  };
#endif // HAVE_ISTL

#if HAVE_UMFPACK || HAVE_SUITESPARSE_UMFPACK
  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::umfpack,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    typedef Dune::Fem::AdaptiveDiscreteFunction< DomainDFSpace >                                        DomainDiscreteFunctionType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< RangeDFSpace >                                         RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                  DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef Dune::Fem::UMFPACKOp< DiscreteFunctionType, LinearOperatorType, symmetric >                 LinearInverseOperatorType;
  };
#else
  template< bool dummy >
  struct AvailableSolvers< Solver::Enum::umfpack, dummy >
  {
    static_warning(false, "You have chosen the UMFPACK solver backend which is currently not installed. Falling back to standard solver!");
  };
#endif

#if HAVE_PETSC
  template <class DomainDFSpace, class RangeDFSpace,bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::petsc,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    typedef Dune::Fem::PetscDiscreteFunction< DomainDFSpace >                                       DomainDiscreteFunctionType;
    typedef Dune::Fem::PetscDiscreteFunction< RangeDFSpace >                                        RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                              DiscreteFunctionType;
    typedef Dune::Fem::PetscLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, LinearOperatorType >             LinearInverseOperatorType;
    // to switch between solvers for symmetric and non symmetric operators
    // use the parameter petsc.kspsolver.method
  };
#else
  template< bool dummy >
  struct AvailableSolvers< Solver::Enum::petsc, dummy >
  {
    static_warning(false, "You have chosen the PetSc solver backend which is currently not installed. Falling back to standard solver!");
  };
#endif

#if HAVE_EIGEN
#warning "Eigen solver is not working at the moment.!"
  template <class DomainDFSpace, class RangeDFSpace,bool symmetric>
  struct MatrixSolverSelector<Solver::Enum::eigen,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    typedef Dune::Fem::EigenVector<double>                                                                          DofVectorType;
    typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DomainDFSpace, DofVectorType > > DomainDiscreteFunctionType;
    typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< RangeDFSpace, DofVectorType > >  RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                              DiscreteFunctionType;
    typedef Dune::Fem::EigenLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef Dune::Fem::EigenCGInverseOperator< DiscreteFunctionType >                               LinearInverseOperatorType;
  };
#else
  template< bool dummy >
  struct AvailableSolvers< Solver::Enum::eigen, dummy >
  {
    static_warning(false, "You have chosen the Eigen solver backend which is currently not installed. Falling back to standard solver!");
  };
#endif


///////////////////////////////////////////////////////////////////////////
// MatrixFreeSolverSelector
///////////////////////////////////////////////////////////////////////////

  //template <Solver::Enum solver, bool symmetric, class DomainDFSpace, class RangeDFSpace = DomainDFSpace>
  //struct MatrixFreeSolverSelector
  //{
  //  static const bool solverConfigured = false; // this implementation is used for not installed packages
  //  // choose type of discrete function, Matrix implementation and solver implementation
  //  // this should work with any discrete function implementation
  //  typedef Dune::Fem::DynamicVector<double>                                                                        DofVectorType;
  //  typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DomainDFSpace, DofVectorType > > DomainDiscreteFunctionType;
  //  typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< RangeDFSpace, DofVectorType > >  RangeDiscreteFunctionType;
  //  typedef DomainDiscreteFunctionType                                                                              DiscreteFunctionType;
  //  typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType >             LinearOperatorType;
  //  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType >                                                    LinearInverseOperatorType;
  //};

  template <Solver::Enum solver, bool symmetric, class DomainDFSpace, class RangeDFSpace = DomainDFSpace>
  struct MatrixFreeSolverSelector
  {
    static const bool solverConfigured = false;
    typedef Dune::Fem::AdaptiveDiscreteFunction< DomainDFSpace >                                        DomainDiscreteFunctionType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< RangeDFSpace >                                         RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                  DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename std::conditional<symmetric,
            Dune::Fem::CGInverseOperator< DiscreteFunctionType >,
            Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType > > :: type          LinearInverseOperatorType;
  };

  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct MatrixFreeSolverSelector<Solver::Enum::fem,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    typedef Dune::Fem::AdaptiveDiscreteFunction< DomainDFSpace >                                        DomainDiscreteFunctionType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< RangeDFSpace >                                         RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                  DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename std::conditional<symmetric,
            Dune::Fem::CGInverseOperator< DiscreteFunctionType >,
            Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType > > :: type          LinearInverseOperatorType;
  };

///////////////////////////////////////////////////////////////////////////
// SolverSelector
///////////////////////////////////////////////////////////////////////////

  template <Solver::Enum solver, bool symmetric, bool matrixfree, class DomainDFSpace, class RangeDFSpace = DomainDFSpace>
  struct SolverSelector;

  template <Solver::Enum solver, bool symmetric, class DomainDFSpace, class RangeDFSpace >
  struct SolverSelector< solver, symmetric, true, DomainDFSpace, RangeDFSpace >
  {
    typedef MatrixFreeSolverSelector< solver, symmetric, DomainDFSpace, RangeDFSpace > type;
  };

  template <Solver::Enum solver, bool symmetric, class DomainDFSpace, class RangeDFSpace >
  struct SolverSelector< solver, symmetric, false, DomainDFSpace, RangeDFSpace >
  {
    typedef MatrixSolverSelector< solver, symmetric, DomainDFSpace, RangeDFSpace > type;
  };


///////////////////////////////////////////////////////////////////////////
// DiscreteFunctionSpaceSelector
///////////////////////////////////////////////////////////////////////////

  template< class FunctionSpaceImp, class GridPartImp, int polOrder, DiscreteFunctionSpaces::Enum dfType, Galerkin::Enum opType >
  struct DiscreteFunctionSpaceSelector;

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::lagrange, Galerkin::Enum::cg >
  {
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, CachingStorage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::legendre, Galerkin::Enum::dg >
  {
    typedef LegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, CachingStorage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::orthonormal, Galerkin::Enum::dg >
  {
    typedef DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, CachingStorage > type;
  };

  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::lagrange, Galerkin::Enum::dg >
  {
    typedef LagrangeDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, CachingStorage > type;
  };


  template< class FunctionSpaceImp, class GridPartImp, int polOrder >
  struct DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrder, DiscreteFunctionSpaces::Enum::hierarchic_legendre, Galerkin::Enum::dg >
  {
    typedef HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, CachingStorage > type;
  };

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

} // end namespace Fem
} // end namespace Dune
#endif
