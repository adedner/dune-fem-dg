#ifndef FEMDG_PROBLEMCREATOR_SELECTOR_HH
#define FEMDG_PROBLEMCREATOR_SELECTOR_HH

// iostream includes
#include <iostream>

#include <dune/fem-dg/misc/static_warning.hh>

// include gridpart
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include discrete function
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/storage/vector.hh>
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
      //! limitation of advection term
      unlimited = 0,
      //! no limitation of advection term
      limited = 1,
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
      //! Discrete function space with Lagrange Finite Elements
      lagrange = 0,
      //! Discrete function space with Legendre Finite Elements
      legendre = 1,
      //! Discrete function space with hierarchic Legendre Finite Elements
      hierarchic_legendre = 2
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
      //! Continuous Galerkin
      cg = 0,
      //! Discontinuous Galerkin
      dg = 1
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
      //! no Adaptivity
      no = 0,
      //! Allow Adaptivity
      yes = 1
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
      //! use the matrix based version of the dune-fem solvers
      fem = 0,
      //! use the matrix based version of the dune-fem solvers with blas
      femoem = 1,
      //! use the dune-istl solvers
      istl = 2,
      //! use the direct solver umfpack
      umfpack = 3,
      //! use the petsc package
      petsc = 4,
      //! use the eigen package
      eigen = 5
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
    typedef Dune::Fem::LeafGridPart< Grid > type;
  };

  template <class Grid>
  struct GridPartSelector<Grid, Galerkin::Enum::dg, Adaptivity::Enum::no>
  {
    typedef Dune::Fem::LeafGridPart< Grid > type;
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

  template< class Op, class DiffusionOp, class AdvectionOp, bool advection, bool diffusion >
  struct ImplExplOperatorSelector
  {
    typedef Op                   FullOperatorType;
    typedef DiffusionOp          ImplicitOperatorType;
    typedef AdvectionOp          ExplicitOperatorType;
  };
  template< class Op, class DiffusionOp, class AdvectionOp, bool advection >
  struct ImplExplOperatorSelector< Op, DiffusionOp, AdvectionOp, advection, false >
  {
    typedef AdvectionOp          FullOperatorType;
    typedef FullOperatorType     ImplicitOperatorType;
    typedef FullOperatorType     ExplicitOperatorType;
  };
  template<class Op, class DiffusionOp, class AdvectionOp, bool diffusion >
  struct ImplExplOperatorSelector< Op, DiffusionOp, AdvectionOp, false, diffusion >
  {
    typedef DiffusionOp          FullOperatorType;
    typedef FullOperatorType     ImplicitOperatorType;
    typedef FullOperatorType     ExplicitOperatorType;
  };


///////////////////////////////////////////////////////////////////////////
// AdvectionDiffusionOperatorSelector
///////////////////////////////////////////////////////////////////////////

  template< class OperatorTraits, AdvectionLimiter::Enum op >
  class AdvectionDiffusionOperatorSelector;

  template< class OperatorTraits >
  class AdvectionDiffusionOperatorSelector< OperatorTraits, AdvectionLimiter::Enum::unlimited >
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
  class AdvectionDiffusionOperatorSelector< OperatorTraits, AdvectionLimiter::Enum::limited >
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

  template< class OperatorTraits, AdvectionLimiter::Enum op, OperatorSplit::Enum split, Matrix::Enum ass >
  class OperatorSelector;

  //matrixfree
  template< class OperatorTraits, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, op, OperatorSplit::Enum::expl, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, op >::ExplicitOperatorType type;
  };

  template< class OperatorTraits, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, op, OperatorSplit::Enum::impl, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, op >::ImplicitOperatorType type;
  };

  template< class OperatorTraits, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, op, OperatorSplit::Enum::full, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, op >::FullOperatorType type;
  };

  template< class OperatorTraits, AdvectionLimiter::Enum op >
  struct OperatorSelector< OperatorTraits, op, OperatorSplit::Enum::rhs, Matrix::Enum::matrixfree >
  {
    typedef typename AdvectionDiffusionOperatorSelector< OperatorTraits, op >::FullOperatorType type;
  };

  //assembled
  //TODO improve DGPrimalMatrixAssembly for correct splitting
  template< class AssemblerTraitsImp, AdvectionLimiter::Enum op >
  struct OperatorSelector< AssemblerTraitsImp, op, OperatorSplit::Enum::full, Matrix::Enum::assembled >
  {
    typedef DGPrimalMatrixAssembly< AssemblerTraitsImp > type;
  };

  template< class AssemblerTraitsImp, AdvectionLimiter::Enum op >
  struct OperatorSelector< AssemblerTraitsImp, op, OperatorSplit::Enum::rhs, Matrix::Enum::assembled >
  {
    typedef DGPrimalMatrixAssembly< AssemblerTraitsImp > type;
  };


///////////////////////////////////////////////////////////////////////////
// SolverSelector
///////////////////////////////////////////////////////////////////////////

  template< Solver::Enum solver, bool dummy = false >
  struct AvailableSolvers
  {};

  template <Solver::Enum solver, bool symmetric, class DomainDFSpace, class RangeDFSpace = DomainDFSpace>
  struct SolverSelector
  {
    static const bool solverConfigured = false; // this implementation is used for not installed packages
    // choose type of discrete function, Matrix implementation and solver implementation
    // this should work with any discrete function implementation
    typedef Dune::Fem::DynamicVector<double>                                                                        DofVectorType;
    typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DomainDFSpace, DofVectorType > > DomainDiscreteFunctionType;
    typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< RangeDFSpace, DofVectorType > >  RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                              DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType >             LinearOperatorType;
    typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType >                                                    LinearInverseOperatorType;
  };

  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct SolverSelector<Solver::Enum::fem,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    typedef Dune::Fem::AdaptiveDiscreteFunction< DomainDFSpace >                                        DomainDiscreteFunctionType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< RangeDFSpace >                                         RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                  DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename Dune::conditional<symmetric,
            Dune::Fem::CGInverseOperator< DiscreteFunctionType >,
            Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType > > :: type          LinearInverseOperatorType;
  };

  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct SolverSelector<Solver::Enum::femoem,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    // this work with a discrete function implementation based on a double* dof storage
    typedef Dune::Fem::AdaptiveDiscreteFunction< DomainDFSpace >                                        DomainDiscreteFunctionType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< RangeDFSpace >                                         RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                                  DiscreteFunctionType;
    typedef Dune::Fem::SparseRowLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename Dune::conditional<symmetric,
            Dune::Fem::OEMCGOp< DiscreteFunctionType, LinearOperatorType >,
            Dune::Fem::OEMBICGSTABOp< DiscreteFunctionType, LinearOperatorType > > :: type              LinearInverseOperatorType;
  };

#if HAVE_DUNE_ISTL
  template <class DomainDFSpace, class RangeDFSpace, bool symmetric>
  struct SolverSelector<Solver::Enum::istl,symmetric,DomainDFSpace,RangeDFSpace>
  {
    static const bool solverConfigured = true;
    // choose type of discrete function, Matrix implementation and solver implementation
    // here we need the special ISTLBlockVectorDiscreteFunction
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DomainDFSpace >                            DomainDiscreteFunctionType;
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< RangeDFSpace >                             RangeDiscreteFunctionType;
    typedef DomainDiscreteFunctionType                                                             DiscreteFunctionType;
    typedef Dune::Fem::ISTLLinearOperator< DomainDiscreteFunctionType, RangeDiscreteFunctionType > LinearOperatorType;
    typedef typename Dune::conditional<symmetric,
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
  struct SolverSelector<Solver::Enum::umfpack,symmetric,DomainDFSpace,RangeDFSpace>
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
  struct SolverSelector<Solver::Enum::petsc,symmetric,DomainDFSpace,RangeDFSpace>
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
  struct SolverSelector<Solver::Enum::eigen,symmetric,DomainDFSpace,RangeDFSpace>
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


}
}
#endif
