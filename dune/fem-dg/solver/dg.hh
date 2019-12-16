#ifndef DUNE_FEM_DG_SOLVER_OPERATOR_HH
#define DUNE_FEM_DG_SOLVER_OPERATOR_HH

// system includes
#include <string>

// dune-fem includes
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#include <dune/fem-dg/algorithm/evolution.hh>

// dune-fem-dg includes
#include <dune/fem-dg/algorithm/evolution.hh>
#include <dune/fem-dg/operator/dg/dgpyoperator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>
#include <dune/fem-dg/models/modelwrapper.hh>
#include <dune/fem-dg/misc/algorithmcreatorselector.hh>

#ifdef EULER_WRAPPER_TEST
#include <dune/fem-dg/models/additional.hh>
#endif

#if HAVE_DUNE_FEMPY
#include <dune/fempy/quadrature/fempyquadratures.hh>
#endif

namespace Dune
{
namespace Fem
{

  // DG solver operator
  //---------------------

  template < class DestinationImp,
             class AdvectionModel,
             class DiffusionModel,
             class Additional>
#ifdef EULER_WRAPPER_TEST
#error
  class DGSolver : public DuneODE :: OdeSolverInterface< DestinationImp >
#else
  class DGSolver : public Fem::SpaceOperatorInterface< DestinationImp >
#endif
  {
  public:
    typedef DestinationImp   DestinationType;
    typedef typename DestinationType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    static const int polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder ;

    typedef typename DiscreteFunctionSpaceType :: GridPartType    GridPartType;
    typedef typename GridPartType::GridType                       GridType;

    typedef typename GridType :: CollectiveCommunication          CollectiveCommunicationType;
    typedef TimeProvider< CollectiveCommunicationType >           TimeProviderType;

    typedef DGOperator< DestinationType, AdvectionModel, DiffusionModel, Additional > DGOperatorType;

    typedef typename DGOperatorType :: FullOperatorType      FullOperatorType;
    typedef typename DGOperatorType :: ExplicitOperatorType  ExplicitOperatorType;
    typedef typename DGOperatorType :: ImplicitOperatorType  ImplicitOperatorType;

    typedef DuneODE::OdeSolverInterface< DestinationType >      OdeSolverInterfaceType;

    static constexpr bool symmetric  =  false ;
    static constexpr bool matrixfree =  true  ;
    static constexpr bool threading  = Additional::threading;

    static const Solver::Enum solverId  = Additional::solverId;

    // solver selection, available fem, istl, petsc, ...
    typedef typename MatrixFreeSolverSelector< solverId, symmetric > :: template LinearInverseOperatorType< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType >  LinearSolverType ;

    // type of runge kutta solver
    typedef RungeKuttaSolver< FullOperatorType, ExplicitOperatorType, ImplicitOperatorType,
                              LinearSolverType > RKSolverType;

    typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

    static std::string name() { return std::string(""); }

    DGSolver( const DiscreteFunctionSpaceType& space,
              const AdvectionModel &advectionModel,
              const DiffusionModel &diffusionModel,
              const TimeSteppingParameters& param = TimeSteppingParameters() )
      : dgOperator_( space, advectionModel, diffusionModel ),
        extra_(),
        tpPtr_( new TimeProviderType(space.gridPart().comm()) ),
        tp_( *tpPtr_ ),
        rkSolver_( tp_, dgOperator_.fullOperator(), dgOperator_.explicitOperator(), dgOperator_.implicitOperator(), name() ),
        initialized_( false )
    {
      //Dune::Fem::Parameter::append("fem.parallel.numberofthreads", std::to_string( Additional::nThreads ) );

      const double maxTimeStep = param.maxTimeStep();
      fixedTimeStep_ = param.fixedTimeStep();

      // start first time step with prescribed fixed time step
      // if it is not 0 otherwise use the internal estimate
      tp_.provideTimeStepEstimate(maxTimeStep);

      // adjust fixed time step with timeprovider.factor()
      fixedTimeStep_ /= tp_.factor() ;
      if ( fixedTimeStep_ > 1e-20 )
        tp_.init( fixedTimeStep_ );
      else
        tp_.init();

      std::cout << "cfl = " << double(tp_.factor()) << " " << tp_.time() << std::endl;
    }

    DGSolver( const DiscreteFunctionSpaceType& space,
              const AdvectionModel &advectionModel,
              const DiffusionModel &diffusionModel,
              const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : dgOperator_( space, advectionModel, diffusionModel ),
        extra_(),
        tpPtr_( new TimeProviderType(space.gridPart().comm(),parameter) ),
        tp_( *tpPtr_ ),
        rkSolver_( tp_, dgOperator_.fullOperator(), dgOperator_.explicitOperator(), dgOperator_.implicitOperator(), name() ),
        initialized_( false )
    {
      //Dune::Fem::Parameter::append("fem.parallel.numberofthreads", std::to_string( Additional::nThreads ) );

      const TimeSteppingParameters param("femdg.stepper.",parameter);
      const double maxTimeStep = param.maxTimeStep();
      fixedTimeStep_ = param.fixedTimeStep();

      // start first time step with prescribed fixed time step
      // if it is not 0 otherwise use the internal estimate
      tp_.provideTimeStepEstimate(maxTimeStep);

      // adjust fixed time step with timeprovider.factor()
      fixedTimeStep_ /= tp_.factor() ;
      if ( fixedTimeStep_ > 1e-20 )
        tp_.init( fixedTimeStep_ );
      else
        tp_.init();

      std::cout << "cfl = " << double(tp_.factor()) << " " << tp_.time() << std::endl;
    }

    DGSolver( TimeProviderType& tp,
              const DiscreteFunctionSpaceType& space,
              const AdvectionModel &advectionModel,
              const DiffusionModel &diffusionModel,
              const TimeSteppingParameters& param = TimeSteppingParameters())
      : dgOperator_( space, advectionModel, diffusionModel ),
        extra_(),
        tpPtr_(),
        tp_( tp ),
        rkSolver_( tp_, dgOperator_.fullOperator(), dgOperator_.explicitOperator(), dgOperator_.implicitOperator(), name() ),
        initialized_( false )
    {
      //Dune::Fem::Parameter::append("fem.parallel.numberofthreads", std::to_string( Additional::nThreads ) );
      std::cout << "cfl = " << double(tp_.factor()) << " " << tp_.time() << std::endl;
    }

    virtual void initialize( const DestinationType& dest )
    {
      checkInitialize( dest );
    }

    virtual void description( std::ostream& out) const { dgOperator_.description( out ); }

    const DiscreteFunctionSpaceType& space () const { return dgOperator_.space(); }
    const DiscreteFunctionSpaceType& domainSpace () const { return space(); }
    const DiscreteFunctionSpaceType& rangeSpace () const { return space(); }

    //! evaluate the operator
    void operator()( const DestinationType& arg, DestinationType& dest ) const
    {
      dest.assign( arg );
      solve( dest );
    }

    void limit( DestinationType &u) const { dgOperator_.limit(u); }

#ifdef EULER_WRAPPER_TEST
    void solve( DestinationType& dest, MonitorType& )
    {
      solve( dest );
    }
#endif

    void checkInitialize( const DestinationType& dest ) const
    {
      if( !initialized_ )
      {
        rkSolver_.initialize( dest );
        // initialize TimeProvider with the estimate obtained form operator
        if ( fixedTimeStep_ > 1e-20 )
          tp_.init( fixedTimeStep_ );
        else
          tp_.init();
        initialized_ = true;
      }
    }

    //! deprecated method
    void solve( DestinationType& dest ) const
    {
      step( dest );
    }

    void step( DestinationType& dest ) const
    {
      // check if initialization needs to be done
      checkInitialize( dest );

      // make sure the current time step is valid
      assert( tp_.timeStepValid() );

      // solve ODE
      rkSolver_.solve( dest, monitor_ );

      if( tpPtr_ )
      {
        // next time step is prescribed by fixedTimeStep
        if ( fixedTimeStep_ > 1e-20 )
          tp_.next( fixedTimeStep_ );
        else
          tp_.next();
      }

      //std::cout << "t = " << tp_.time() << "  dt = " << tp_.deltaT() << std::endl;

#ifndef EULER_WRAPPER_TEST
      // return limited solution, to be discussed
      // this is only enabled when AdvectionLimiter is not unlimited
      limit( dest );
#endif
    }

    void setTimeStepSize( const double dt )
    {
      fixedTimeStep_  = dt ;
      fixedTimeStep_ /= tp_.factor() ;
      tp_.provideTimeStepEstimate( dt );
    }

    double deltaT() const { return tp_.deltaT(); }

  protected:
    DGOperatorType dgOperator_;

    std::tuple<>                          extra_;
    std::unique_ptr< TimeProviderType >   tpPtr_;
    TimeProviderType&                     tp_;
    mutable MonitorType                   monitor_;
    mutable RKSolverType                  rkSolver_;
    mutable double                        fixedTimeStep_ ;
    mutable bool                          initialized_;
  };

} // end namespace Fem
} // end namespace Dune
#endif
