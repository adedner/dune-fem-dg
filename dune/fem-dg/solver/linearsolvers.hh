#ifndef DUNE_FEMDG_LINEARSOLVERS_HH
#define DUNE_FEMDG_LINEARSOLVERS_HH

// iostream includes
#include <iostream>

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

#if HAVE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#endif

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscsolver.hh>
#endif

/*********************************************************/

enum SolverType
{
  matrixFree,  // use the matrix free version of the dune-fem solvers
  fem,         // use the matrix based version of the dune-fem solvers
  femoem,      // use the matrix based version of the dune-fem solvers with blas
  istl,        // use the dune-istl solvers
  umfpack,     // use the direct solver umfpack
  petsc        // use the petsc package
};

enum GalerkinType
{
  cg, dg
};

template <class Grid, GalerkinType op>
struct GridPartChooser;
template <class Grid>
struct GridPartChooser<Grid,cg>
{
  typedef Dune::Fem::AdaptiveLeafGridPart< Grid, Dune::InteriorBorder_Partition > Type;
};
template <class Grid>
struct GridPartChooser<Grid,dg>
{
  // typedef Dune::Fem::LeafGridPart< Grid, Dune::All_Partition > Type;
  typedef Dune::Fem::LeafGridPart< Grid > Type;
};

template <class DFSpace, SolverType solver, bool symmetric>
struct Solvers
{
  static const bool solverConfigured = false; // this implementation is used for not installed packages
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  // this should work with any discrete function implementation
  typedef Dune::Fem::DynamicVector<double> DofVectorType;
  typedef Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DiscreteFunctionSpaceType, DofVectorType > > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;

};

template <class DFSpace, bool symmetric>
struct Solvers<DFSpace,fem,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef typename Dune::conditional<symmetric,
          Dune::Fem::CGInverseOperator< DiscreteFunctionType >,
          Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType > > :: type
          LinearInverseOperatorType;

};
template <class DFSpace,bool symmetric>
struct Solvers<DFSpace,femoem,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  // this work with a discrete function implementation based on a double* dof storage
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef typename Dune::conditional<symmetric,
          Dune::Fem::OEMCGOp< DiscreteFunctionType, LinearOperatorType >,
          Dune::Fem::OEMBICGSTABOp< DiscreteFunctionType, LinearOperatorType > > :: type
          LinearInverseOperatorType;
};

#if HAVE_DUNE_ISTL
template <class DFSpace,bool symmetric>
struct Solvers<DFSpace,istl,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  // here we need the special ISTLBlockVectorDiscreteFunction
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef typename Dune::conditional<symmetric,
          Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType >,
          Dune::Fem::ISTLBICGSTABOp< DiscreteFunctionType, LinearOperatorType > > :: type
          LinearInverseOperatorType;
};
#endif // HAVE_ISTL
#if HAVE_UMFPACK
template <class DFSpace, bool symmetric>
struct Solvers<DFSpace,umfpack,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::UMFPACKOp< DiscreteFunctionType, LinearOperatorType, symmetric > LinearInverseOperatorType;
};
#endif
#if HAVE_PETSC
template <class DFSpace,bool symmetric>
struct Solvers<DFSpace,petsc,symmetric>
{
  static const bool solverConfigured = true;
  typedef DFSpace DiscreteFunctionSpaceType;
  // choose type of discrete function, Matrix implementation and solver implementation
  typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::PetscInverseOperator< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
  // to switch between solvers for symmetric and non symmetric operators
  // use the parameter petsc.kspsolver.method
};
#endif

#endif // end #if DUNE_FEMDG_LINEARSOLVERS_HH
