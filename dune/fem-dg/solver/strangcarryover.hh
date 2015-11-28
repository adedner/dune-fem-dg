#ifndef DUNE_FEM_SOLVER_STRANGCARRYOVER_HH
#define DUNE_FEM_SOLVER_STRANGCARRYOVER_HH

//- system includes
#include <cassert>
#include <limits>
#include <sstream>
#include <vector>

//- dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

//- dune-fem includes
#include <dune/fem/solver/odesolver.hh>

#include <dune/fem-dg/solver/newtoninverseoperator.hh>

namespace DuneODE
{

  /** \brief Implicit RungeKutta ODE solver. */
  template< class ExplicitOperator, class HelmholtzOperator, class LinearSolver >
  class StrangCarryoverSolver
  : public OdeSolverInterface< typename HelmholtzOperator::DomainFunctionType >
  {
    typedef StrangCarryoverSolver< ExplicitOperator, HelmholtzOperator, LinearSolver > ThisType;
    typedef OdeSolverInterface< typename HelmholtzOperator::DomainFunctionType > BaseType;
    typedef Dune::Fem::TimeProviderBase TimeProviderType;
    typedef ImplicitRungeKuttaTimeStepControl TimeStepControlType;
    typedef typename TimeStepControlType::ParametersType ParametersType;
    typedef typename HelmholtzOperator::JacobianOperatorType JacobianOperatorType;
    //typedef typename HelmholtzOperator::SpaceOperatorType::PreconditionOperatorType PreconditionOperatorType;

  public:
    typedef typename BaseType::MonitorType MonitorType;
    typedef typename BaseType::DestinationType DestinationType;

    typedef HelmholtzOperator HelmholtzOperatorType;
    typedef ExplicitOperator ExplicitOperatorType;
    typedef LinearSolver LinearSolverType;

    /** \brief constructor
     *
     *  \param[in]  helmholtzOp      Helmholtz operator \f$L\f$
     *  \param[in]  butcherTable     butcher table to use
     *  \param[in]  timeStepControl  time step controller
     */
    StrangCarryoverSolver ( ExplicitOperatorType &explicitOp,
                            HelmholtzOperatorType &helmholtzOp,
                            TimeProviderType &timeProvider,
                            int order = 1,
                            const NSNewtonParameter &parameter = NSNewtonParameter(),
                            const ParametersType &parameterTS = ParametersType() )
    : helmholtzOp_( helmholtzOp ),
      explicitOp_( explicitOp_ ),
      timeStepControl_( timeProvider, parameterTS ),
      rhs_( "RK rhs", helmholtzOp_.space() ),
      G_( "carry over", helmholtzOp_.space() ),
      q_( "tmp", helmholtzOp_.space() ),
      jOp_ ( "jacobianOperator", helmholtzOp.space(), helmholtzOp.space() ),
      //preconditioner_(helmholtzOp_.spaceOperator().preconditioner()),
      tolerance_( parameter.toleranceParameter() ),
      linAbsTol_( parameter.linAbsTolParameter( tolerance_ ) ),
      linReduction_( parameter.linReductionParameter( tolerance_ ) ),
      verbose_( parameter.verbose() && MPIManager::rank () == 0 ),
      linVerbose_( parameter.linearSolverVerbose() ),
      maxIterations_( parameter.maxIterationsParameter() ),
      maxLinearIterations_( parameter.maxLinearIterationsParameter() )
    {
    }

    /** \brief destructor */
    ~StrangCarryoverSolver ()
    {
    }

    //! apply operator once to get dt estimate
    void initialize ( const DestinationType &U0 )
    {
      const double time = timeStepControl_.time();

      helmholtzOp_.setTime( time );
      helmholtzOp_.initializeTimeStepSize( U0 );
      const double helmholtzEstimate = helmholtzOp_.timeStepEstimate();

#if 0
      explicitOp_.setTime( time );
      explicitOp_.initializeTimeStepSize( U0 );
      double sourceTermEstimate = explicitOp_.timeStepEstimate();
#else
      double sourceTermEstimate = 0;
#endif

      // negative time step is given by the empty source term
      if( sourceTermEstimate < 0.0 ) sourceTermEstimate = helmholtzEstimate ;

      timeStepControl_.initialTimeStepSize( helmholtzEstimate, sourceTermEstimate );

      helmholtzOp_.setLambda( 0.5 );
      helmholtzOp_.spaceOperator()(U0,rhs_);
      helmholtzOp_.jacobian( U0, jOp_ );
      /*
      if (preconditioner_)
      {
        const LinearSolver jInv( jOp_, *preconditioner_, linReduction_, linAbsTol_, maxLinearIterations_, linVerbose_ );
        jInv(rhs_,G_);
      }
      else
      */
      {
        const LinearSolver jInv( jOp_, linReduction_, linAbsTol_, maxLinearIterations_, linVerbose_ );
        jInv(rhs_,G_);
      }
    }

    using BaseType::solve;

    //! solve the system
    void solve ( DestinationType &U, MonitorType &monitor )
    {
      monitor.reset();
      //helmholtzOp_.spaceOperator().resetCalls();

      const double time = timeStepControl_.time();
      const double timeStepSize = timeStepControl_.timeStepSize();
      assert( timeStepSize > 0.0 );

      // explicitOp_.setTime( time );
      helmholtzOp_.setTime( time );

      U.axpy(timeStepSize/2.,G_);
#if 0
      // add semi implicit stuff if required
      // q2: q=Hf(u,1) = u + dt f(u)
      // q3: q=3/4u + 1/4 Hf(q,1) = 1/4q + 3/4u +  dt/4dt f(u)
      // q4: u=1/3u + 2/3 Hf(q,1) = 1/3u + 2/3q + 2dt/3dt f(u)
      explicitOp_(U,rhs_);
      q_.assign(U);
      q_.axpy(dt,rhs_);
      explicitOp_(q_,rhs_);
      q_ *= 1./4.;
      q_.axpy(3./4.,U);
      q_.axpy(dt/4.,rhs_);
      explicitOp_(q_,rhs_);
      U *= 1./3.;
      q_.axpy(1./4.,U);
      q_.axpy(dt*2./3.,rhs_);
#endif
      helmholtzOp_.setLambda( 0.5 );
      helmholtzOp_.spaceOperator()(U,rhs_);
      helmholtzOp_.jacobian( U, jOp_ );
      /*
      if (preconditioner_)
      {
        const LinearSolver jInv( jOp_, *preconditioner_, linReduction_, linAbsTol_, maxLinearIterations_, linVerbose_ );
        jInv(rhs_,G_);
      }
      else
      */
      {
        const LinearSolver jInv( jOp_, linReduction_, linAbsTol_, maxLinearIterations_, linVerbose_ );
        jInv(rhs_,G_);
      }
      int linIterations = 0;//helmholtzOp_.spaceOperator().calls();
      monitor.linearSolverIterations_ += linIterations;
      U.axpy(timeStepSize/2.,G_);

      // update time step size
      // timeStepControl_.timeStepEstimate( helmholtzOp_.timeStepEstimate(), explicitOp_.timeStepEstimate(), monitor );
      timeStepControl_.timeStepEstimate( helmholtzOp_.timeStepEstimate(), 0., monitor );
    }

    void description ( std::ostream &out ) const
    {
      out << "Strang Carryover solver.\\\\" << std::endl;
    }

  protected:
    HelmholtzOperatorType &helmholtzOp_;
    ExplicitOperator &explicitOp_;
    TimeStepControlType timeStepControl_;

    DestinationType rhs_, G_, q_;

    JacobianOperatorType jOp_;
    //const PreconditionOperatorType *preconditioner_;

    const double tolerance_, linAbsTol_, linReduction_;
    const bool verbose_;
    const bool linVerbose_;
    const int maxIterations_;
    const int maxLinearIterations_;
  };

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BASICIMPLICIT_HH
