#ifndef DUNE_FEM_DG_SMARTODESOLVER_HH
#define DUNE_FEM_DG_SMARTODESOLVER_HH

#include <limits>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/solver/rungekutta/explicit.hh>
#include <dune/fem/solver/odesolver.hh>

namespace Dune {


struct SmartOdeSolverParameters : public DuneODE :: ODEParameters
{
  using DuneODE :: ODEParameters :: keyPrefix_;

  virtual double explicitFactor() const
  {
    return Fem::Parameter::getValue< double >( keyPrefix_ + "explicitfactor" , 1.0 );
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
    std::string key( keyPrefix_ + "odesolver" );
    if( Fem :: Parameter :: exists( key ) )
      return Fem::Parameter::getEnum( key, odeSolver, 0 );
    else
    {
      std::cerr << "WARNING: deprecated key, use `fem.ode.odesolver' instread!" << std::endl;
      return Fem::Parameter::getEnum( "femhowto.odesolver", odeSolver, 0 );
    }
  }

  virtual int obtainRungeKuttaSteps( const int defaultRKOrder ) const
  {
    std::string key( keyPrefix_ + "order");
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

  typedef DuneODE :: OdeSolverInterface< typename Operator :: DestinationType > BaseType;

  template <class AdvOp, class DiffOp>
  struct IsSame
  {
    static bool check(const AdvOp&, const DiffOp& )
    {
      return false;
    }
  };

  template <class AdvOp>
  struct IsSame< AdvOp, AdvOp>
  {
    static bool check(const AdvOp& a, const AdvOp& d)
    {
      return true;
    }
  };

public:
  typedef typename OperatorType :: DestinationType DestinationType ;
  typedef DestinationType  DiscreteFunctionType;

  template < class Op, class DF >
  struct OdeSolverSelection
  {
    typedef DuneODE :: ExplicitRungeKuttaSolver< DiscreteFunctionType >  ExplicitOdeSolverType;
    typedef ExplicitOdeSolverType ImplicitOdeSolverType ;
    typedef ExplicitOdeSolverType SemiImplicitOdeSolverType ;

    static ExplicitOdeSolverType*
    createExplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps )
    {
      return new ExplicitOdeSolverType( op, tp, rkSteps );
    }

    static ImplicitOdeSolverType*
    createImplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps )
    {
      DUNE_THROW(NotImplemented," ");
      return 0;
    }

    template <class ExplOp, class ImplOp>
    static SemiImplicitOdeSolverType*
    createSemiImplicitSolver( ExplOp& explOp, ImplOp& implOp, Fem::TimeProviderBase& tp, const int rkSteps )
    {
      DUNE_THROW(NotImplemented," ");
      return 0;
    }

  };

  template < class Op >
  struct OdeSolverSelection< Op, Fem :: AdaptiveDiscreteFunction< typename Op::SpaceType > >
  {
    typedef Fem :: AdaptiveDiscreteFunction< typename Operator::SpaceType >  DiscreteFunctionType ;
    typedef DuneODE :: ImplicitOdeSolver< DiscreteFunctionType >       ImplicitOdeSolverType;
    typedef DuneODE :: ExplicitOdeSolver< DiscreteFunctionType >       ExplicitOdeSolverType;
    typedef DuneODE :: SemiImplicitOdeSolver< DiscreteFunctionType >   SemiImplicitOdeSolverType;

    static ExplicitOdeSolverType*
    createExplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps )
    {
      return new ExplicitOdeSolverType( op, tp, rkSteps );
    }

    static ImplicitOdeSolverType*
    createImplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps )
    {
      return new ImplicitOdeSolverType( op, tp, rkSteps );
    }

    template <class ExplOp, class ImplOp>
    static SemiImplicitOdeSolverType*
    createSemiImplicitSolver( ExplOp& explOp, ImplOp& implOp, Fem::TimeProviderBase& tp, const int rkSteps )
    {
      return new SemiImplicitOdeSolverType( explOp, implOp, tp, rkSteps );
    }

  };

  // The ODE Solvers
  typedef DuneODE :: OdeSolverInterface< DestinationType >           OdeSolverInterfaceType;
  typedef OdeSolverSelection< OperatorType, DestinationType >  OdeSolversType ;
  typedef typename OdeSolversType :: ImplicitOdeSolverType  ImplicitOdeSolverType ;
  typedef typename OdeSolversType :: ExplicitOdeSolverType  ExplicitOdeSolverType ;
  typedef typename OdeSolversType :: SemiImplicitOdeSolverType  SemiImplicitOdeSolverType ;

  typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

  using BaseType :: solve;
protected:
  OperatorType&           operator_;
  Fem::TimeProviderBase&  timeProvider_;

  AdvectionOperatorType&  advectionOperator_;
  DiffusionOperatorType&  diffusionOperator_;
  const SmartOdeSolverParameters* param_;

  OdeSolverInterfaceType* explicitSolver_;
  OdeSolverInterfaceType* odeSolver_;

  const double explFactor_ ;
  const int verbose_ ;
  const int rkSteps_ ;
  const int odeSolverType_ ;
  int imexCounter_ , exCounter_;
  int minIterationSteps_, maxIterationSteps_ ;
  bool useImex_ ;
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
     useImex_( odeSolverType_ > 1 )
  {
    // create implicit or explicit ode solver
    if( odeSolverType_ == 0 )
    {
      odeSolver_ = OdeSolversType :: createExplicitSolver( operator_, tp, rkSteps_);
    }
    else if (odeSolverType_ == 1)
    {
      odeSolver_ = OdeSolversType :: createImplicitSolver( operator_, tp, rkSteps_);
    }
    else if( odeSolverType_ > 1 )
    {
      // make sure that advection and diffusion operator are different
      if( IsSame< AdvectionOperatorType, DiffusionOperatorType >::check( advectionOperator_, diffusionOperator_ ) )
      {
        DUNE_THROW(Dune::InvalidStateException,"Advection and Diffusion operator are the same, therefore IMEX cannot work!");
      }
      odeSolver_ = OdeSolversType :: createSemiImplicitSolver( advectionOperator_, diffusionOperator_, tp, rkSteps_);

      // IMEX+
      if( odeSolverType_ == 3 )
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

  void getAdvectionDiffsionTimeSteps( double& advStep, double& diffStep ) const
  {
    if( odeSolverType_ > 1 )
    {
      double steps[ 2 ] = { advectionOperator_.timeStepEstimate(),
                            diffusionOperator_.timeStepEstimate() };
      // get global min
      Dune :: Fem :: MPIManager :: comm().min( &steps[ 0 ], 2 );

      advStep  = steps[ 0 ];
      diffStep = steps[ 1 ];
    }
  }

  //! solver the ODE
  void solve( DestinationType& U ,
              MonitorType& monitor )
  {
    // take CPU time of solution process
    Timer timer ;

    // switch upwind direction
    operator_.switchupwind();
    advectionOperator_.switchupwind();
    diffusionOperator_.switchupwind();

    // reset compute time counter
    resetComputeTime();

    double maxAdvStep  = std::numeric_limits< double > :: max();
    double maxDiffStep = std::numeric_limits< double > :: max();

    if( explicitSolver_ && ! useImex_ )
    {
      explicitSolver_->solve( U, monitor );
      ++exCounter_ ;
    }
    else
    {
      assert( odeSolver_ );
      odeSolver_->solve( U, monitor );

      ++imexCounter_ ;

      const int iterationSteps = monitor.newtonIterations_ * monitor.linearSolverIterations_ ;
      minIterationSteps_ = std::min( minIterationSteps_, iterationSteps );
      maxIterationSteps_ = std::max( maxIterationSteps_, iterationSteps );

      if( verbose_ == 3 )
      {
        // get advection and diffusion time step
        getAdvectionDiffsionTimeSteps( maxAdvStep, maxDiffStep );

        double factor = explFactor_ ;
        //if( averageIterationSteps > 0 )
        //  factor *= averageIterationSteps / (rkSteps_ + 1 ) ;
        std::cout << maxAdvStep << " a | d " << maxDiffStep << "  factor: " << factor
          << "  " << minIterationSteps_ << " min | max " << maxIterationSteps_ << "  use imex = " << useImex_ << std::endl;
      }
    }

    if( explicitSolver_ )
    {
      // get advection and diffusion time step
      getAdvectionDiffsionTimeSteps( maxAdvStep, maxDiffStep );

      const int averageIterationSteps = (minIterationSteps_ + maxIterationSteps_) / 2;
      double factor = explFactor_ ;
      if( averageIterationSteps > 0 )
        factor *= double(rkSteps_ + 1) / double(averageIterationSteps) ;

      // if true solve next time step with semi implicit solver
      useImex_ = ( maxDiffStep < (factor * maxAdvStep) ) ;

      if( verbose_ == 3 )
      {
        std::cout << maxAdvStep << " a | d " << maxDiffStep << "  factor: " << factor
          << "  " << minIterationSteps_ << " min | max " << maxIterationSteps_
          << "  use imex = " << useImex_ << "  ex steps: " << exCounter_ << std::endl;
      }

      // make sure the correct time step is used for the explicit solver
      //if( ! useImex_ )
      //  timeProvider_.provideTimeStepEstimate( operator_.timeStepEstimate() ) ;
    }

    // store needed time
    monitor.odeSolveTime_     = timer.elapsed();
    monitor.operatorTime_     = operatorTime();
    monitor.numberOfElements_ = numberOfElements();
  }

  //! return CPU time needed for the operator evaluation
  double operatorTime() const
  {
    if( useImex_ )
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
    if( useImex_ )
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
    info << "Explicit Solver used: " << exCounter_ << std::endl;
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
