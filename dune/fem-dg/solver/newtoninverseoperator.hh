// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_METSTROEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_METSTROEM_NEWTONINVERSEOPERATOR_HH

#include <cfloat>
#include <iostream>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>

namespace Dune
{

  namespace Fem
  {

    // NewtonParameter
    // ---------------

    struct NSNewtonParameter
#ifndef DOXYGEN
    : public LocalParameter< NSNewtonParameter, NSNewtonParameter >
#endif
    {
      NSNewtonParameter () {}

      virtual double toleranceParameter () const
      {
        return Parameter::getValue< double >( "fem.solver.newton.tolerance", 1e-6 );
      }

      virtual double linAbsTolParameter ( const double &tolerance )  const
      {
        return Parameter::getValue< double >( "fem.solver.newton.linabstol", tolerance / 8 );
      }

      virtual double linReductionParameter ( const double &tolerance ) const
      {
        return Parameter::getValue< double >( "fem.solver.newton.linreduction", tolerance / 8 );
      }

      virtual bool verbose () const
      {
        const bool v = Parameter::getValue< bool >( "fem.solver.verbose", false );
        return Parameter::getValue< bool >( "fem.solver.newton.verbose", v );
      }

      virtual bool linearSolverVerbose () const
      {
        const bool v = Parameter::getValue< bool >( "fem.solver.verbose", false );
        return Parameter::getValue< bool >( "fem.solver.newton.linear.verbose", v );
      }

      virtual int maxIterationsParameter () const
      {
        return Parameter::getValue< int >( "fem.solver.newton.maxiterations", std::numeric_limits< int >::max() );
      }

      virtual int maxLinearIterationsParameter () const
      {
        return Parameter::getValue< int >( "fem.solver.newton.maxlineariterations", std::numeric_limits< int >::max() );
      }

      virtual bool eisenstatWalker() const
      {
        return Parameter::getValue< bool >("fem.solver.newton.eisenstatwalker", false );
      }
    };



    // NewtonInverseOperator
    // ---------------------

    /** \class NewtonInverseOperator
     *  \brief inverse operator based on a newton scheme
     *
     *  \tparam  Op      operator to invert (must be a DifferentiableOperator)
     *  \tparam  LInvOp  linear inverse operator
     *
     *  \note Verbosity of the NewtonInverseOperator is controlled via the
     *        paramter <b>fem.solver.newton.verbose</b>; it defaults to
     *        <b>fem.solver.verbose</b>.
     */
    template< class Op, class LInvOp >
    class NSNewtonInverseOperator
    : public Operator< typename Op::JacobianOperatorType::RangeFunctionType, typename Op::JacobianOperatorType::DomainFunctionType >
    {
      typedef NSNewtonInverseOperator< Op, LInvOp > ThisType;
      typedef Operator< typename Op::JacobianOperatorType::RangeFunctionType, typename Op::JacobianOperatorType::DomainFunctionType > BaseType;

    public:
      //! type of operator's Jacobian
      typedef typename Op::JacobianOperatorType JacobianOperatorType;
      //typedef typename Op::SpaceOperatorType::PreconditionOperatorType PreconditionOperatorType;

      //! type of operator to invert
      typedef DifferentiableOperator< JacobianOperatorType > OperatorType;

      //! type of linear inverse operator
      typedef LInvOp LinearInverseOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef NSNewtonParameter ParametersType;

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */
      explicit NSNewtonInverseOperator ( const Op &op,
                                         const NSNewtonParameter &parameter = NSNewtonParameter() )
      : op_( op ),
        tolerance_( parameter.toleranceParameter() ),
        linAbsTol_( parameter.linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameter.linReductionParameter( tolerance_ ) ),
        verbose_( parameter.verbose() && MPIManager::rank () == 0 ),
        linVerbose_( parameter.linearSolverVerbose() ),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() ),
        //preconditioner_(op.spaceOperator().preconditioner()),
        useEisenstatWalker_( parameter.eisenstatWalker() )
      {
        actTolerance_ = tolerance_;
      }

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  epsilon  tolerance for norm of residual
       */
      NSNewtonInverseOperator ( const OperatorType &op, const DomainFieldType &epsilon,
                                const NSNewtonParameter &parameter = NSNewtonParameter() )
      : op_( op ),
        tolerance_( epsilon ),
        linAbsTol_( parameter.linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameter.linReductionParameter( tolerance_ ) ),
        verbose_( parameter.verbose() ),
        linVerbose_( parameter.linearSolverVerbose() ),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() ),
        //preconditioner_(op.spaceOperator().preconditioner())
        useEisenstatWalker_( parameter.eisenstatWalker() )
      {
        actTolerance_ = tolerance_;
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const;

      int iterations () const { return iterations_; }
      int linearIterations () const { return linearIterations_; }
      int spaceOperatorCalls () const { return spaceOperatorCalls_; }

      bool converged () const
      {
        // check for finite |residual| - this also works for -ffinite-math-only (gcc)
        const bool finite = (delta_ < std::numeric_limits< DomainFieldType >::max());
        return finite && (iterations_ < maxIterations_) && (linearIterations_ < maxLinearIterations_);
      }

    protected:
      double solve( const double actTol, const DomainFunctionType &u, RangeFunctionType &w ) const;

      const Op &op_;
      const double tolerance_, linAbsTol_, linReduction_;
      mutable double actTolerance_;
      const bool verbose_;
      const bool linVerbose_;
      const int maxIterations_;
      const int maxLinearIterations_;

      mutable DomainFieldType delta_;
      mutable int iterations_;
      mutable int linearIterations_;
      mutable int spaceOperatorCalls_;

      //const PreconditionOperatorType *preconditioner_;

      const bool useEisenstatWalker_;
    };



    // Implementation of NewtonInverseOperator
    // ---------------------------------------

    template< class JacobianOperator, class LInvOp >
    inline void NSNewtonInverseOperator< JacobianOperator, LInvOp >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      /*
      static int newtonCall = 0 ;
      ++newtonCall;

      if( false && newtonCall % 100 == 1 )
      {
        const RangeFunctionType initW( w );

        actTolerance_ = 0.5 ;
        double ratio = actTolerance_;
        double residuum = solve( actTolerance_, u, w, ratio );
        w.assign( initW );
        actTolerance_ *= 0.5;
        double newResiduum  = solve( actTolerance_, u, w, ratio );
        std::cout << "Check toler = " << actTolerance_ << " " << residuum << " " << newResiduum << std::endl;
        if( std::abs( newResiduum - residuum ) < 1e-14 ) return ;
        if( ( newResiduum / residuum ) < tolerance_ )
        {
          std::cout << "Check toler = " << actTolerance_ << " " << residuum << " " << newResiduum << std::endl;
          return ;
        }

        while( (actTolerance_ > 1e-10 ) )
        {
          std::cout << "Check toler = " << actTolerance_ << " " << residuum << " " << newResiduum << std::endl;
          if( std::abs( newResiduum - residuum ) < 1e-14 ) return ;
          if( ( newResiduum / residuum ) < tolerance_ )
          {
            std::cout << "Choose actTol = " << actTolerance_ << std::endl;
            return ;
          }
          residuum = newResiduum ;
          actTolerance_ *= 0.5;
          w.assign( initW );
          newResiduum = solve( ratio, u, w, ratio );
        }
        std::cout << "Choose actTol = " << actTolerance_ << std::endl;
      }
      else
      {
        //double ratio = 0;
      */
        solve( tolerance_, u, w );
      //}
    }


    template< class JacobianOperator, class LInvOp >
    inline double NSNewtonInverseOperator< JacobianOperator, LInvOp >::
    solve( const double actTol, const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      DomainFunctionType residual( "Newton-residual", u.space() );
      RangeFunctionType dw("Newton-dw", w.space() );
      JacobianOperatorType jOp( "jacobianOperator", dw.space(), u.space() );

      //preconditioner_->setUtilde(u);

      // compute initial residual
      op_( w, residual );
      residual -= u;
      double normRes = std::sqrt( residual.scalarProductDofs( residual ) );
      delta_ = normRes;

      //const double tolerance = tolerance_ * normRes;
      //delta_ = 1;
      //
      const double etaMax = 0.9;
      double deltaold = 0.5 * normRes ;
      const double gamma = 0.5 ;
      double linAbsTol = etaMax;
      //const double alpha2 = ( 1.0 + std::sqrt( 5.0 ) ) * 0.5;
      //const double alpha  = 1.1;
      const double alpha = 2.0;
      // const double alpha  = ( 1.0 + std::sqrt( 5.0 ) ) * 0.5;

      //std::cout << "Newton start delta: " << normRes << std::endl;

      for( iterations_ = 1, linearIterations_ = 0; converged(); ++iterations_ )
      {
        // evaluate operator's jacobian
        op_.jacobian( w, jOp );

        const int remLinearIts = maxLinearIterations_ - linearIterations_;
        dw.clear();

        // update linAbsTol
        if( iterations_ > 1 )
        {
          //std::cout << normRes << " n | o " << deltaold << std::endl;
          double newLinAbsTol = std::min( etaMax, (gamma * std::pow( normRes / deltaold , alpha )) );
          if( gamma * std::pow( deltaold, alpha ) > 0.1 )
            linAbsTol = std::min( etaMax, std::max( newLinAbsTol, gamma * ( deltaold * deltaold ) ) );
          else
            linAbsTol = newLinAbsTol;
        }

        double linReduction = ( useEisenstatWalker_ ) ? linAbsTol : linReduction_ ;
        //std::cout << "Using lin abs tol(" << iterations_ << ") : " << linAbsTol << std::endl;
        //linAbsTol = ( ! useEisenstatWalker_ ) ? linAbsTol_ : tolerance_ * linAbsTol;
        linAbsTol = ( ! useEisenstatWalker_ ) ? linAbsTol_ : actTol * linAbsTol;
        int linIterations;
        /*
        if (preconditioner_)
        {
          const LinearInverseOperatorType jInv( jOp, *preconditioner_, linReduction, linAbsTol, remLinearIts, linVerbose_ );
          jInv( residual, dw );
          linIterations = jInv.iterations();
        }
        else
        */
        {
          const LinearInverseOperatorType jInv( jOp, linReduction, linAbsTol, remLinearIts, linVerbose_ );
          jInv( residual, dw );
          linIterations = jInv.iterations();
        }
        linearIterations_ += linIterations;
        if( verbose_ )
          std::cerr << "Newton iteration " << iterations_
                    << ": linear iterations used = " << linIterations
                    << ": |residual| = " << delta_ << std::endl;

        w -= dw;

        op_( w, residual );
        residual -= u;

        deltaold = normRes ;
        normRes  = std::sqrt( residual.scalarProductDofs( residual ) );
        delta_   = normRes;

        //std::cout << "new/old " << normRes / deltaold << " actTol = " << actTol << std::endl;
        if (delta_ < actTol ) break;
      }
      if( verbose_ )
        std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;
      return delta_;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
