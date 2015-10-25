#ifndef DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH
#define DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include <dune/fem-dg/stepper/stepperbase.hh>
#include "../stokes/stokesalgorithm.hh"

template <class GridImp,
          class ProblemTraits,
          int polynomialOrder,
          class ExtraParameterTuple = std::tuple<> >
struct AdvectionDiffusionStepper
  : public StepperBase< GridImp, ProblemTraits, polynomialOrder >
{
  typedef StepperBase< GridImp, ProblemTraits, polynomialOrder, ExtraParameterTuple > BaseType ;

  typedef typename BaseType :: Traits     Traits ;
  typedef typename Traits :: StokesAlgorithmType  StokesAlgorithmType;

  // type of Grid
  typedef typename BaseType :: GridType                 GridType;

  // Choose a suitable GridView
  typedef typename BaseType :: GridPartType             GridPartType;

  // initial data type
  typedef typename BaseType :: InitialDataType          InitialDataType;

  // An analytical version of our model
  typedef typename BaseType :: ModelType                 ModelType;

  // The flux for the discretization of advection terms
  typedef typename BaseType :: FluxType                  FluxType;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef typename BaseType :: FullOperatorType          FullOperatorType;
  typedef typename BaseType :: ExplicitOperatorType      ExplicitOperatorType;
  typedef typename BaseType :: ImplicitOperatorType      ImplicitOperatorType;

  typedef typename Traits :: VelocityFunctionType        VelocityFunctionType;

  typedef typename Traits :: RhsOperatorType             RhsOperatorType;
  typedef typename Traits :: RhsStokesOperatorType       RhsStokesOperatorType;

  typedef typename BaseType :: LinearInverseOperatorType LinearInverseOperatorType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename BaseType :: DiscreteFunctionType      DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename BaseType :: DiscreteSpaceType         DiscreteSpaceType;

  // The ODE Solvers
  typedef typename BaseType :: OdeSolverType     OdeSolverType;

  typedef typename BaseType :: TimeProviderType       TimeProviderType;
  typedef typename BaseType :: AdaptationManagerType  AdaptationManagerType;
  typedef typename BaseType :: AdaptationHandlerType  AdaptationHandlerType;

  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId =
    BaseType::Traits::DiffusionFluxId ;

  typedef typename BaseType::OperatorTraits OperatorTraits;
  typedef typename BaseType::SolverMonitorType  SolverMonitorType;

  // advection = true , diffusion = true
  typedef Dune :: DGAdaptationIndicatorOperator< OperatorTraits, true, true >  DGIndicatorType;

  // gradient estimator
  typedef Estimator< DiscreteFunctionType, InitialDataType > GradientIndicatorType ;

  // type of 64bit unsigned integer
  typedef typename BaseType :: UInt64Type  UInt64Type;

  typedef typename OperatorTraits :: ExtraParameterTupleType  ExtraParameterTupleType;


  using BaseType :: grid_;
  using BaseType :: gridPart_;
  using BaseType :: space;
  using BaseType :: problem;
  using BaseType :: adaptationHandler_ ;
  using BaseType :: adaptParam_;
  using BaseType :: adaptive ;
  using BaseType :: doEstimateMarkAdapt ;
  using BaseType :: name ;
  using BaseType :: solution;
  using BaseType :: overallTimer_;
  using BaseType :: odeSolver_;
  using BaseType :: odeSolverMonitor_;

  AdvectionDiffusionStepper( GridType& grid,
                             const std::string name = "",
                             ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
    BaseType( grid, name ),
    subTimeProvider_( 0.0, 1.0, grid ),
    stokes_( grid, name ),
    velocity_( "velocity", stokes_.space() ),
    tmp_( "tmp", stokes_.space() ),
    rhs_( "rhs", stokes_.space() ),
    tuple_( &velocity_, &rhs_ ),
    tupleV_( &velocity_ ),
    dgOperator_( gridPart_, problem(), tuple_, name ),
    rhsOperator_( gridPart_, problem(), tupleV_, name ),
    rhsStokesOperator_( gridPart_, problem(), std::tuple<>(), name )
  {
  }

  //! return overal number of grid elements
  virtual UInt64Type gridSize() const
  {
    // is adaptation handler exists use the information to avoid global comm
    if( adaptationHandler_ )
    {
      UInt64Type globalElements = adaptationHandler_->globalNumberOfElements() ;
      if( Dune::Fem::Parameter::verbose () )
      {
        std::cout << "grid size (sum,min,max) = ( "
          << globalElements << " , "
          << adaptationHandler_->minNumberOfElements() << " , "
          << adaptationHandler_->maxNumberOfElements() << ")" << std::endl;
      }
      return globalElements;
    }

    // one of them is not zero,
    size_t dgSize      = dgOperator_.numberOfElements();
    UInt64Type grSize  = dgSize ;
    double minMax[ 2 ] = { double(grSize), 1.0/double(grSize) } ;
    grid_.comm().max( &minMax[ 0 ], 2 );
    if( Dune::Fem::Parameter :: verbose () )
    {
      std::cout << "grid size (min,max) = ( " << size_t(1.0/minMax[ 1 ]) << " , " << size_t(minMax[ 0 ]) << ")" << std::endl;
    }
    return grid_.comm().sum( grSize );
  }

  virtual OdeSolverType* createOdeSolver(TimeProviderType& tp)
  {
    // create ODE solver
    typedef RungeKuttaSolver< FullOperatorType, FullOperatorType, FullOperatorType,
                              LinearInverseOperatorType > OdeSolverImpl;
    return new OdeSolverImpl( subTimeProvider_,
                              dgOperator_,
                              dgOperator_,
                              dgOperator_,
                              name() );
  }

  // before first step, do data initialization
  void initializeStep( TimeProviderType& tp, const int loop )
  {
    const double theta = problem().theta();
    BaseType :: initializeStep( tp, loop );
    // pass dt estimate to time provider
    const double dtEstimate = subTimeProvider_.timeStepEstimate();
    tp.provideTimeStepEstimate( dtEstimate / (1.0 - 2.0 * theta) );
  }

  void step(TimeProviderType& tp,
            SolverMonitorType& monitor )
  {
    const double theta = problem().theta();
    const double time  = tp.time();
    const double dt    = tp.deltaT();

    DiscreteFunctionType& U = solution();
    std::cout << "Sizes: ";
    std::cout << velocity_.size() << "  " << stokes_.solution().size() << "  " <<
      U.size() << std::endl;

/*
    typedef typename InitialDataType :: TimeDependentFunctionType
      TimeDependentFunctionType;

    // communication is needed when blocking communication is used
    // but has to be avoided otherwise (because of implicit solver)
    const bool doCommunicate = ! NonBlockingCommParameter :: nonBlockingCommunication ();

    // create L2 projection
    Fem :: L2Projection< TimeDependentFunctionType,
        DiscreteFunctionType > l2pro( 2 * U.space().order(), doCommunicate );

    // L2 project initial data
    l2pro( problem().fixedTimeFunction( tp.time() ), U );
    tp.provideTimeStepEstimate( 1e-1 );
    return ;
*/

    assert( velocity_.size() == stokes_.solution().size() );
    assert( U.size() == velocity_.size() );

    stokes_.assemble() ;

    // reset overall timer
    overallTimer_.reset();

    // u^* = u^n
    velocity_.assign( U );

    // set time for problem
    problem().setTime( time );

    typedef typename RhsOperatorType::DestinationType RhsDestinationType;
    RhsDestinationType tmp1( "tmp1", space() );
    RhsDestinationType tmp2( "tmp2", space() );
    tmp1.assign( U );

    // set operator time
    rhsOperator_.setTime( time );
    // compute right hand side
    rhsOperator_( tmp1, tmp2 );
    rhs_.assign( tmp2 );

    // stokes solve (step 1)
    stokes_.solve( rhs_ );

    // update solution
    U.assign( stokes_.solution() );

    return ;

    // u^* = (2 theta - 1)/theta * u^n + (1 - theta)/theta u^n+theta
    velocity_ *= ( 2.0*theta - 1.0 ) / theta ;
    velocity_.axpy( (1.0-theta)/theta, stokes_.solution() );

    // set time for problem
    problem().setTime( time + theta * dt );

    tmp1.assign( stokes_.solution() );
    // compute laplace U
    rhsStokesOperator_( tmp1, tmp2 );
    rhs_.assign( tmp2 );

    // compute grad p
    stokes_.pressureGradient( tmp_ );

    // rhs = \alpha\mu\laplace U - \grad P
    rhs_ -= tmp_;

    // set time for ode solve
    subTimeProvider_.init();
    subTimeProvider_.provideTimeStepEstimate( time + dt * theta );
    subTimeProvider_.next();
    subTimeProvider_.next( (1.0 - 2.0 *theta) * dt );

    // set time for problem
    problem().setTime( time + (1.0 - 2.0*theta) * dt );

    std::cout << time << "  t | dt  " << dt << "  " << std::endl;
    std::cout << subTimeProvider_.time() << "  t | dt  " << subTimeProvider_.deltaT() << std::endl;
    // subTimeProvider_.setTime( time + dt * theta );
    // compute time step for advection-diffusion step
    // subTimeProvider_.provideTimeStepEstimate( (1.0 - 2.0 *theta) * dt / subTimeProvider_.factor() );
    // TODO: set time step correctly

    // solve advection-diffusion step (step 2)
    assert(odeSolver_);
    odeSolver_->solve( U, odeSolverMonitor_ );

    // set time for problem
    problem().setTime( time + (1.0-theta) * dt );

    // set operator time
    rhsOperator_.setTime( time + (1.0-theta) * dt );
    tmp1.assign( U );
    // compute right hand side
    rhsOperator_( tmp1, tmp2 );
    rhs_.assign( tmp2 );

    // stokes solve (step 3)
    stokes_.solve( rhs_ );
    U.assign( stokes_.solution() );

    // limit solution if necessary
    // limitSolution ();

    // copy information to solver monitor
    monitor.newton_iterations     = odeSolverMonitor_.newtonIterations_;
    monitor.ils_iterations        = odeSolverMonitor_.linearSolverIterations_;
    monitor.max_newton_iterations = odeSolverMonitor_.maxNewtonIterations_ ;
    monitor.max_ils_iterations    = odeSolverMonitor_.maxLinearSolverIterations_;
    monitor.operator_calls        = odeSolverMonitor_.spaceOperatorCalls_;

    // set time step size to monitor
    monitor.setTimeStepInfo( tp );

    const double dtEstimate = subTimeProvider_.timeStepEstimate();
    tp.provideTimeStepEstimate( dtEstimate / (1.0 - 2.0 * theta) );
  }


  //! estimate and mark solution
  virtual void initialEstimateMarkAdapt( )
  {
  }

  //! estimate and mark solution
  virtual void estimateMarkAdapt( )
  {
  }

  const ModelType& model() const { return dgOperator_.model(); }

protected:
  //typename OperatorTraits::SpaceType vSpace_;
  //typename OperatorTraits::VeloType  velo_;
  TimeProviderType        subTimeProvider_;
  StokesAlgorithmType     stokes_;

  VelocityFunctionType    velocity_;
  VelocityFunctionType    tmp_;
  VelocityFunctionType    rhs_;

  std::tuple< VelocityFunctionType*, VelocityFunctionType * > tuple_;
  std::tuple< VelocityFunctionType* > tupleV_;

  FullOperatorType        dgOperator_;
  RhsOperatorType         rhsOperator_;
  RhsStokesOperatorType   rhsStokesOperator_;
};
#endif // FEMHOWTO_STEPPER_HH
