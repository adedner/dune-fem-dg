#ifndef DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH
#define DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include "stepperbase.hh"
#include "../stokes/stokes.hh"

template <class GridImp,
          class ProblemTraits,
          int polynomialOrder,
          class ExtraParameterTuple = std::tuple<> >
struct AdvectionDiffusionStepper
  : public StepperBase< GridImp, ProblemTraits, polynomialOrder >
{
  typedef StepperBase< GridImp, ProblemTraits, polynomialOrder, ExtraParameterTuple > BaseType ;

  typedef typename ProblemTraits :: StokesProblemTraits StokesProblemTraits;
  typedef StokesAlgorithm< GridImp, StokesProblemTraits, polynomialOrder >  StokesAlgorithm;

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
  typedef typename BaseType :: FullOperatorType               FullOperatorType;
  typedef typename BaseType :: ExplicitOperatorType           ExplicitOperatorType;
  typedef typename BaseType :: ImplicitOperatorType           ImplicitOperatorType;

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

  AdvectionDiffusionStepper( GridType& grid,
                             const std::string name = "",
                             ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
    BaseType( grid, name ),
    stokes_( grid, name ),
    velocity_( stokes_.solution() ),
    rhs_( stokes_.solution() ),
    tuple_( &velo_, &rhs_ ),
    tupleV_( &velo_ ),
    subTimeProvider_( 0.0, grid ),
    dgOperator_( gridPart_, problem(), tuple_, name ),
    rhsOperator_( gridPart_, problem(), tupleV_, name )
    // rhsOperator2_( gridPart_, problem(), tuple_, name )
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
    size_t advSize     = dgAdvectionOperator_.numberOfElements();
    size_t diffSize    = dgDiffusionOperator_.numberOfElements();
    size_t dgIndSize   = gradientIndicator_.numberOfElements();
    size_t dgSize      = dgOperator_.numberOfElements();
    UInt64Type grSize  = std::max( std::max(advSize, dgSize ), std::max( diffSize, dgIndSize ) );
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
    // create adaptation handler in case of apost indicator
    if( adaptive() )
    {
      if( ! adaptationHandler_ && adaptParam_.aposterioriIndicator() )
      {
        adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
        dgIndicator_.setAdaptation( *adaptationHandler_ );
      }
    }

    // create ODE solver
    typedef RungeKuttaSolver< FullOperatorType, FullOperatorType, FullOperatorType,
                              LinearInverseOperatorType > OdeSolverImpl;
    return new OdeSolverImpl( subTimeProvider_,
                              dgOperator_,
                              dgOperator_,
                              dgOperator_,
                              name() );
  }

  void step(TimeProviderType& tp,
            SolverMonitorType& monitor )
  {
    const double theta = 1.0 - 1.0/M_SQRT2 ;
    const double time  = tp.time();
    const double dt    = tp.deltaT();

    DiscreteFunctionType& U = solution();

    // reset overall timer
    overallTimer_.reset();

    // u^* = u^n
    velocity_.assign( stokes_.solution() );

    // set operator time
    rhsOperator_.setTime( time );
    // compute right hand side
    rhsOperator_( U, rhs_ );

    // stokes solve (step 1)
    stokes_.solve( &rhs_ );
    // update solution
    U.assign( stokes_.solution() );

    // u^* = (2 theta - 1)/theta * u^n + (1 - theta)/theta u^n+theta
    velocity_ *= ( 2.0*theta - 1.0 ) / theta ;
    velocity_.axpy( (1.0-theta)/theta, U );

    rhsOperator2_( U, rhs_ );

    // set time for ode solve
    subTimeProvider_.setTime( time + dt * theta );
    // compute time step for advection-diffusion step
    subTimeProvider_.provideTimeStepEstimate( (1.0 - 2.0 *theta) * dt / subTimeProvider_.factor() );
    // TODO: set time step correctly

    // solve advection-diffusion step (step 2)
    assert(odeSolver_);
    odeSolver_->solve( U, odeSolverMonitor_ );

    // set operator time
    rhsOperator_.setTime( time + (1.0-theta) * dt );
    // compute right hand side
    rhsOperator_( U, rhs_ );

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

  }


  //! estimate and mark solution
  virtual void initialEstimateMarkAdapt( )
  {
    doEstimateMarkAdapt( dgIndicator_, gradientIndicator_, true );
  }

  //! estimate and mark solution
  virtual void estimateMarkAdapt( )
  {
    doEstimateMarkAdapt( dgIndicator_, gradientIndicator_, false );
  }

  const ModelType& model() const { return dgOperator_.model(); }

protected:
  //typename OperatorTraits::SpaceType vSpace_;
  //typename OperatorTraits::VeloType  velo_;
  ExtraParameterTupleType tuple_;
  StokesAlgorithm         stokes_;
  TimeProviderType        subTimeProvider_;

  FullOperatorType        dgOperator_;
  ExplicitOperatorType    dgAdvectionOperator_;
  ImplicitOperatorType    dgDiffusionOperator_;
  DGIndicatorType         dgIndicator_;
  GradientIndicatorType   gradientIndicator_;
};
#endif // FEMHOWTO_STEPPER_HH
