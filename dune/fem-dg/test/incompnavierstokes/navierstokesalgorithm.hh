#ifndef DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH
#define DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include <dune/fem-dg/algorithm/evolution.hh>

#include <dune/fem-dg/operator/dg/primaloperator.hh>

namespace Dune
{
namespace Fem
{

template <class GridImp,
          class ProblemTraits,
            int polynomialOrder >
  struct NavierStokesStepper
    : public EvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder >
{
  typedef EvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder > BaseType;

  typedef typename BaseType::Traits::StokesAlgorithmType    StokesAlgorithmType;

  // type of Grid
  typedef typename BaseType::GridType                       GridType;

  // Choose a suitable GridView
  typedef typename BaseType::GridPartType                   GridPartType;

  // initial data type
  typedef typename BaseType::ProblemType                    ProblemType;

  // An analytical version of our model
  typedef typename BaseType::ModelType                      ModelType;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef typename BaseType::OperatorType::FullType         FullOperatorType;
  typedef typename BaseType::OperatorType::ExplicitType     ExplicitOperatorType;
  typedef typename BaseType::OperatorType::ImplicitType     ImplicitOperatorType;

  typedef typename Traits::VelocityFunctionType             VelocityFunctionType;

  typedef typename Traits::RhsOperatorType                  RhsOperatorType;
  typedef typename Traits::RhsLaplaceOperatorType           RhsLaplaceOperatorType;

  typedef typename BaseType::BasicLinearSolverType          BasicLinearSolverType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename BaseType::DiscreteFunctionType           DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename BaseType::DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;

  // The ODE Solvers
  typedef typename BaseType::OdeSolverType                  OdeSolverType;

  typedef typename BaseType::TimeProviderType               TimeProviderType;

  //typedef typename BaseType::SolverMonitorType            SolverMonitorType;

  // type of 64bit unsigned integer
  typedef typename BaseType::UInt64Type                     UInt64Type;

  typedef typename BaseType::ExtraParameterTupleType        ExtraParameterTupleType;

  using BaseType::grid_;
  using BaseType::gridPart_;
  using BaseType::space;
  using BaseType::problem;
  using BaseType::name;
  using BaseType::adaptHandler_
  using BaseType::solution;
  using BaseType::overallTimer_;
  using BaseType::odeSolver_;
  //using BaseType::odeSolverMonitor_;

  NavierStokesStepper( GridType& grid,
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
    rhsLaplaceOperator_( gridPart_, problem(), std::tuple<>(), name )
  {
    adaptHandler_.setIndicator( problem(), tuple );
  }

  //! return overal number of grid elements
  virtual UInt64Type gridSize() const
  {
    int globalElements = adaptHandler_.globalNumberOfElements();
    if( globalElements > 0 )
      return globalElements;

    // one of them is not zero,
    size_t dgSize      = dgOperator_.numberOfElements();
    UInt64Type grSize  = dgSize;
    double minMax[ 2 ] = { double(grSize), 1.0/double(grSize) } ;
    grid_.comm().max( &minMax[ 0 ], 2 );
    if( Dune::Fem::Parameter::verbose () )
    {
      std::cout << "grid size (min,max) = ( " << size_t(1.0/minMax[ 1 ]) << " , " << size_t(minMax[ 0 ]) << ")" << std::endl;
    }
    return grid_.comm().sum( grSize );
  }

  virtual OdeSolverType* createOdeSolver(TimeProviderType& tp)
  {
    // create ODE solver
    typedef RungeKuttaSolver< FullOperatorType, ExplicitOperatorType, ImplicitOperatorType,
                              BasicLinearSolverType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgOperator_,
                              dgAdvectionOperator_,
                              dgDiffusionOperator_,
                              name() );
  }

  // before first step, do data initialization
  void initializeStep( TimeProviderType& tp, const int loop )
  {
    const double theta = problem().theta();
    BaseType::initializeStep( tp, loop );
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

    //----------- STEP 1 -------------------------

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

    ////compute right hand side////
    // set operator time
    rhsOperator_.setTime( time );
    tmp1.assign( U );
    rhsOperator_( tmp1, tmp2 );
    rhs_.assign( tmp2 );

    // stokes solve (step 1)
    stokes_.solve( rhs_ );

    // update solution
    U.assign( stokes_.solution() );

    return ;
    // ---------- STEP 2 -------------------------

    // u^* = (2 theta - 1)/theta * u^n + (1 - theta)/theta u^n+theta
    velocity_ *= ( 2.0*theta - 1.0 ) / theta ;
    velocity_.axpy( (1.0-theta)/theta, stokes_.solution() );

    // set time for problem
    problem().setTime( time + theta * dt );

    tmp1.assign( stokes_.solution() );

    //// compute rhs////
    // compute laplace U
    rhsLaplaceOperator_( tmp1, tmp2 );
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

    // ---------- STEP 3 -------------------------

    // set time for problem
    problem().setTime( time + (1.0-theta) * dt );

    ////compute rhs////
    // set operator time
    rhsOperator_.setTime( time + (1.0-theta) * dt );
    tmp1.assign( U );
    // compute right hand side
    rhsOperator_( tmp1, tmp2 );
    rhs_.assign( tmp2 );

    // stokes solve (step 3)
    stokes_.solve( rhs_ );
    U.assign( stokes_.solution() );

    const double dtEstimate = subTimeProvider_.timeStepEstimate();
    tp.provideTimeStepEstimate( dtEstimate / (1.0 - 2.0 * theta) );
  }


  const ModelType& model() const { return dgOperator_.model(); }

protected:
  TimeProviderType        subTimeProvider_;
  StokesAlgorithmType     stokes_;

  VelocityFunctionType    velocity_;
  VelocityFunctionType    tmp_;
  VelocityFunctionType    rhs_;

  std::tuple< VelocityFunctionType*, VelocityFunctionType * > tuple_;
  std::tuple< VelocityFunctionType* > tupleV_;

  FullOperatorType        dgOperator_;
  RhsOperatorType         rhsOperator_;
  RhsLaplaceOperatorType   rhsLaplaceOperator_;
};

}
}
#endif // FEMHOWTO_STEPPER_HH
