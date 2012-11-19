#ifndef DUNE_FEM_DG_SMARTODESOLVER_HH
#define DUNE_FEM_DG_SMARTODESOLVER_HH

#include <limits>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/solver/odesolver.hh>

namespace Dune {


struct SmartOdeSolverParameters : public DuneODE :: ODEParameters 
{
  virtual double explicitFactor() const
  {
    return Fem::Parameter::getValue< double >( "fem.ode.explicitfactor" , 1.0 );
  }

  SmartOdeSolverParameters* clone() const 
  {
    return new SmartOdeSolverParameters( *this );
  }

  virtual int obtainOdeSolverType () const 
  {
    // we need this choice of explicit or implicit ode solver
    // defined here, so that it can be used later in two different
    // methods
    static const std::string odeSolver[]  = { "EX", "IM" ,"IMEX", "IMEX+"  };
    std::string key( "fem.ode.odesolver" );
    if( Fem :: Parameter :: exists( key ) )
      return Fem::Parameter::getEnum( "fem.ode.odesolver", odeSolver, 0 );
    else 
    {
      std::cerr << "WARNING: deprecated key, use `fem.ode.odesolver' instread!" << std::endl;
      return Fem::Parameter::getEnum( "femhowto.odesolver", odeSolver, 0 );
    }
  }

  virtual int obtainRungeKuttaSteps( const int defaultRKOrder ) const
  {
    std::string key("fem.ode.order");
    if ( Fem :: Parameter :: exists( key ) )
      return Fem::Parameter::getValue< int > ( key );
    else
      return defaultRKOrder ;
  }

};

template <class Operator, 
          class AdvectionOperator, 
          class DiffusionOperator> 
class SmartOdeSolver : 
  public DuneODE :: OdeSolverInterface< typename Operator :: DestinationType >
{
  typedef Operator           OperatorType;
  typedef AdvectionOperator  AdvectionOperatorType;
  typedef DiffusionOperator  DiffusionOperatorType;

public:
  typedef typename OperatorType :: DestinationType DestinationType ;
  typedef DestinationType  DiscreteFunctionType;

  // The ODE Solvers                                                  
  typedef DuneODE :: OdeSolverInterface< DestinationType >           OdeSolverInterfaceType;
  typedef DuneODE :: ImplicitOdeSolver< DiscreteFunctionType >       ImplicitOdeSolverType; 
  typedef DuneODE :: ExplicitOdeSolver< DiscreteFunctionType >       ExplicitOdeSolverType; 
  typedef DuneODE :: SemiImplicitOdeSolver< DiscreteFunctionType >   SemiImplicitOdeSolverType;

  typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

protected:
  OperatorType&          operator_;
  Fem::TimeProviderBase&      timeProvider_; 

  AdvectionOperatorType& advectionOperator_;
  DiffusionOperatorType& diffusionOperator_;
  const SmartOdeSolverParameters* param_;

  OdeSolverInterfaceType* explicitSolver_;
  OdeSolverInterfaceType* odeSolver_;

  const double explFactor_ ;
  const int verbose_ ;
  const int rkSteps_ ; 
  const int odeSolverType_;
  int imexCounter_ , exCounter_;
  int minIterationSteps_, maxIterationSteps_ ;
  bool imex_ ;
public:
  SmartOdeSolver( Fem::TimeProviderBase& tp,
                  OperatorType& op, 
                  AdvectionOperatorType& advOp,
                  DiffusionOperatorType& diffOp,
                  const SmartOdeSolverParameters &parameter = SmartOdeSolverParameters() )
   : operator_( op ), 
     timeProvider_( tp ),
     advectionOperator_( advOp ),
     diffusionOperator_( diffOp ),
     param_( parameter.clone() ),
     explicitSolver_( 0 ),
     odeSolver_( 0 ),
     explFactor_( param_->explicitFactor() ),
     verbose_( param_->verbose() ),
     rkSteps_( param_->obtainRungeKuttaSteps( operator_.space().order() + 1 ) ),
     odeSolverType_( param_->obtainOdeSolverType() ),
     imexCounter_( 0 ), exCounter_ ( 0 ),
     minIterationSteps_( std::numeric_limits< int > :: max() ),
     maxIterationSteps_( 0 ),
     imex_( odeSolverType_ > 1 )
  {
    // create implicit or explicit ode solver
    if( odeSolverType_ == 0 )
    {
      odeSolver_ = new ExplicitOdeSolverType( operator_, tp, rkSteps_);
    }
    else if (odeSolverType_ == 1)
    {
      odeSolver_ = new ImplicitOdeSolverType( operator_, tp, rkSteps_);
    }
    else if( odeSolverType_ == 2 )
    {
      odeSolver_ = new SemiImplicitOdeSolverType
        (advectionOperator_, diffusionOperator_, tp, rkSteps_);
    }
    else if( odeSolverType_ == 3 ) 
    {
      odeSolver_ = new SemiImplicitOdeSolverType
          (advectionOperator_, diffusionOperator_, tp, rkSteps_); 
      explicitSolver_ = new ExplicitOdeSolverType( operator_, tp, rkSteps_);
    }
    else 
    {
      DUNE_THROW(NotImplemented,"Wrong ODE solver selected");
    }
  }

  //! destructor 
  ~SmartOdeSolver() 
  {
    delete param_;           param_ = 0;
    delete odeSolver_;       odeSolver_ = 0;
    delete explicitSolver_ ; explicitSolver_ = 0;
  }

  //! initialize method 
  void initialize( const DestinationType& U )
  {
    if( explicitSolver_ )
    {
      explicitSolver_->initialize( U );
    }
    assert( odeSolver_ );
    odeSolver_->initialize( U );

  }

  //! solver the ODE 
  void solve( DestinationType& U , 
              MonitorType& monitor ) 
  {
    // take CPU time of solution process 
    Timer timer ;

    operator_.switchupwind();
    advectionOperator_.switchupwind();
    diffusionOperator_.switchupwind();

    // reset compute time counter 
    resetComputeTime();

    double maxAdvStep  = 0;
    double maxDiffStep = 0;

    if( explicitSolver_ && ! imex_ ) 
    {
      explicitSolver_->solve( U, monitor );

      maxAdvStep  = operator_.maxAdvectionTimeStep();
      maxDiffStep = operator_.maxDiffusionTimeStep();
      ++exCounter_ ;
    }
    else 
    {
      assert( odeSolver_ );
      odeSolver_->solve( U, monitor );

      maxAdvStep  = advectionOperator_.maxAdvectionTimeStep();
      maxDiffStep = diffusionOperator_.maxDiffusionTimeStep();
      ++imexCounter_ ;

      const int iterationSteps = monitor.newtonIterations_ * monitor.linearSolverIterations_ ;
      minIterationSteps_ = std::min( minIterationSteps_, iterationSteps );
      maxIterationSteps_ = std::max( maxIterationSteps_, iterationSteps );
    }

    if( explicitSolver_ ) 
    {
      const int averageIterationSteps = (minIterationSteps_ + maxIterationSteps_) / 2;
      double factor = explFactor_ ;
      if( averageIterationSteps > 0 ) 
        factor *= averageIterationSteps / (rkSteps_ + 1 ) ;

      // if true solve next time step with semi implicit solver
      imex_ = ( maxDiffStep > (factor * maxAdvStep) ) ;

      if( verbose_ == 2 ) 
      {
        std::cout << maxAdvStep << " a | d " << maxDiffStep << "  factor: " << factor
          << "  " << minIterationSteps_ << " min | max " << maxIterationSteps_ << "  use imex = " << imex_ << std::endl;
      }

      // make sure the correct time step is used for the explicit solver 
      if( ! imex_ ) 
        timeProvider_.provideTimeStepEstimate( operator_.timeStepEstimate() ) ;
    }

    // store needed time 
    monitor.odeSolveTime_     = timer.elapsed();
    monitor.operatorTime_     = operatorTime();
    monitor.numberOfElements_ = numberOfElements();
  }

  //! return CPU time needed for the operator evaluation
  double operatorTime() const 
  {
    if( imex_ ) 
    {
      return advectionOperator_.computeTime() +
             diffusionOperator_.computeTime() ;
    }
    else 
      return operator_.computeTime();
  }

  //! return number of elements meat during operator evaluation 
  size_t numberOfElements() const  
  {
    if( imex_ ) 
      return advectionOperator_.numberOfElements();
    else 
      return operator_.numberOfElements();
  }

  void description(std::ostream&) const {}

  // gather information from the space operator, the time integratior
  // and the problem to output before each table in tex file
  std::string description() const
  {
    std::string latexInfo;

    if ((odeSolverType_==0) || (odeSolverType_==1))
      latexInfo = operator_.description();
    else
      latexInfo = advectionOperator_.description()
                  + diffusionOperator_.description();

    std::stringstream odeInfo; 
    odeSolver_->description( odeInfo );

    latexInfo += odeInfo.str() + "\n";
    std::stringstream info; 
    info << "Regular  Solver used: " << imexCounter_ << std::endl;
    info << "Explicit Solver sued: " << exCounter_ << std::endl;
    latexInfo += info.str();

    return latexInfo;
  }

protected:
  void resetComputeTime() const 
  {
    // this will reset the internal time counters 
    operator_.computeTime() ;
    advectionOperator_.computeTime();
    diffusionOperator_.computeTime();
  }
}; // end SmartOdeSolver 

} // end namespace Dune 
#endif
