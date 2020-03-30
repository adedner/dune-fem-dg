#ifndef DUNE_OPTIM_LP_HH
#define DUNE_OPTIM_LP_HH

#include <dune/common/densematrix.hh>

#include <dune/optim/activeindexmapper.hh>
#include <dune/optim/common/smallobject.hh>
#include <dune/optim/std/subarray.hh>

namespace Dune
{

  namespace Optim
  {

    namespace __LinearProgramming
    {

      // ContraintMatrix
      // ---------------

      template< class Constraints, class ActiveIndexMapper >
      struct ConstraintMatrix
      {
        ConstraintMatrix ( const Constraints &constraints, const ActiveIndexMapper &active )
          : constraints_( constraints ), active_( active )
        {}

        const Constraints &constraints () const noexcept { return constraints_; }
        const ActiveIndexMapper &active () const { return active_; }

        const std::size_t N () const { return active().size(); }
        const std::size_t M () const { return active().size(); }

      private:
        const Constraints &constraints_;
        const ActiveIndexMapper &active_;
      };

    } // namespace __LinearProgramming



    // LinearProgramming
    // -----------------

    template< class LinearSolver, bool verbose = false >
    struct LinearProgramming
    {
      typedef typename LinearSolver::Field Field;

      explicit LinearProgramming ( const Field &epsilon, const LinearSolver &linearSolver = LinearSolver() ) : epsilon_( epsilon ), linearSolver_( linearSolver ) {}

      /**
       * \brief solve inequality-constrained LP problem
       *
       * \note The constraints are such that \f$n * x \le c\f$.
       **/
      template< class DomainVector, class ConstraintArray, class ActiveIndexMapper >
      void operator() ( const DomainVector &descent, const ConstraintArray &constraints, DomainVector &x, ActiveIndexMapper &active ) const;

    private:
      Field epsilon_;
      LinearSolver linearSolver_;
    };



    // Implementation of LinearProgramming
    // -----------------------------------

    template< class LinearSolver, bool verbose >
    template< class DomainVector, class ConstraintArray, class ActiveIndexMapper >
    inline void
    LinearProgramming< LinearSolver, verbose >::operator() ( const DomainVector &descent, const ConstraintArray &constraints, DomainVector &x, ActiveIndexMapper &active ) const
    {
      typedef typename ConstraintArray::value_type Constraint;
      typedef typename ActiveIndexMapper::InactiveIterator InactiveIterator;

      __LinearProgramming::ConstraintMatrix< ConstraintArray, ActiveIndexMapper > constraintMatrix( constraints, active );

      assert( active.size() == x.size() );
      auto linearInverse = linearSolver_( constraintMatrix );

      DomainVector dx( x ), lambda( x );
      while( true )
      {
        // compute Langrange multipliers
        linearInverse.mtv( descent, lambda );

        // find a constraint with negative Lagrange multiplier
        int q = -1;
        Field lambda_min = -epsilon_;
        for( unsigned int i = 0; i < active.size(); ++i )
        {
          const Field &lambda_i = lambda[ i ];
          if( lambda_i >= lambda_min )
            continue;
          lambda_min = lambda_i;
          q = i;
        }
        if( q == -1 )
          return;
        if( verbose )
          std::cout << "Releasing constraint " << active[ q ] << std::endl;

        // compute search direction
        lambda = Field( 0 );
        lambda[ q ] = Field( -1 );
        linearInverse.mv( lambda, dx );

        // find the nearest constraint in the descent direction
        Field alpha = std::numeric_limits< Field >::infinity();
        int p = 0;
        const InactiveIterator end = active.endInactive();
        for( InactiveIterator it = active.beginInactive(); it != end; ++it )
        {
          const Constraint &constraint = constraints[ *it ];
          const Field ndx = constraint.normal() * dx;
          if( ndx < epsilon_ )
            continue;
          const Field beta = -constraint.evaluate( x ) / ndx;
          if( beta >= alpha )
            continue;

          alpha = beta;
          p = *it;
        }
        if( verbose )
          std::cout << "Adding Constraint: " << p << " (alpha = " << alpha << ")" << std::endl;

        // walk towards this constraint
        x.axpy( alpha, dx );

        // update active indices and linear solver
        active.update( q, p );
        linearInverse.updateRow( q, constraintMatrix );
      }
    }

  } // namespace Optim



  // DenseMatrixAssigner for Optim::__LinearProgramming::ConstraintMatrix
  // ---------------------------------------------------------------------

  template< class DenseMatrix, class Constraints, class ActiveIndexMapper >
  struct DenseMatrixAssigner< DenseMatrix, Optim::__LinearProgramming::ConstraintMatrix< Constraints, ActiveIndexMapper > >
  {
    static void apply ( DenseMatrix &denseMatrix, const Optim::__LinearProgramming::ConstraintMatrix< Constraints, ActiveIndexMapper > &constraintMatrix )
    {
      const Constraints &constraints = constraintMatrix.constraints();
      const ActiveIndexMapper &active = constraintMatrix.active();

      const std::size_t size = active.size();
      assert( (denseMatrix.N() == size) && (denseMatrix.M() == size) );
      for( std::size_t i = 0; i < size; ++i )
        denseMatrix[ i ] = constraints[ active[ i ] ].normal();
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_OPTIM_LP_HH
