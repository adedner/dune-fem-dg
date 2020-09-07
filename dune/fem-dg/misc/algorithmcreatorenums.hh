#ifndef FEMDG_ALGORITHMCREATOR_ENUMS_HH
#define FEMDG_ALGORITHMCREATOR_ENUMS_HH

// iostream includes
#include <iostream>
#include <type_traits>

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
      //! no limitation of advection term
      unlimited,
      //! limitation of advection term
      limited,
      //! scaling limitation of advection term
      scalinglimited
    };
  }

  /**
   *  \brief Namespace containing an Enum class to describe the limiting
   *  function used in the limited advection operators.
   */
  namespace AdvectionLimiterFunction
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      default_,
      //! MinMod limiter
      minmod,
      //! Superbee limiter
      superbee,
      //! van Leer limiter
      vanleer,
      //! no limiter
      none
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
      default_, // default is orthonormal which is the most general space
      //! Discrete function space with hierarchic orthonormal monomial basis functions
      orthonormal,
      //! Discrete function space with Legendre Finite Elements
      legendre,
      //! Discrete function space with hierarchic Legendre Finite Elements
      hierarchic_legendre,
      //! Discrete function space with Lagrange Finite Elements
      lagrange,
      //! p-adaptive space from dune-fem, implementing dg and lagrange
      padaptive,
      //! Lagrange space with GaussLobatto interpolation points
      gausslobatto,
      //! Lagrange space with GaussLegendre interpolation points
      gausslegendre
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

} // end namespace Fem
} // end namespace Dune
#endif
