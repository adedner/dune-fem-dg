#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH


#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

// helm holz operator wrapper
//#include <dune/fem/operator/dg/helmholtz.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/operator/dg/passtraits.hh>
#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>
#include <dune/fem/misc/femtimer.hh>

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem-dg/pass/threadpass.hh>
#include <dune/common/timer.hh>
#include "monitor.hh"


#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>

namespace Dune
{

namespace Fem
{


  class StepperParameters
  : public Dune::Fem::LocalParameter< StepperParameters, StepperParameters >
  {
    protected:
    const std::string keyPrefix_;

    public:
    StepperParameters( const std::string keyPrefix = "femdg.stepper." )
      : keyPrefix_( keyPrefix )
    {}

    virtual double fixedTimeStep() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "fixedtimestep" , 0.0 );
    }

    virtual double fixedTimeStepEocLoopFactor() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "fixedtimestepeocloopfactor" , 1.0 );
    }

    virtual double startTime() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "starttime" , 0.0 );
    }

    virtual double endTime() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "endtime"/*, 1.0 */);
    }

    virtual int printCount() const
    {
      return Dune::Fem::Parameter::getValue< int >( keyPrefix_ + "printcount" , -1 );
    }

    virtual double maxTimeStep() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "maxtimestep", std::numeric_limits<double>::max());
    }

    virtual int maximalTimeSteps () const
    {
      return Dune::Fem::Parameter::getValue< int >(  keyPrefix_ + "maximaltimesteps", std::numeric_limits<int>::max());
    }

    virtual bool stopAtEndTime() const
    {
      return Dune::Fem::Parameter::getValue< bool >( keyPrefix_ + "stopatendtime", bool(false) );
    }

  };



  // EvolutionAlgorithmTraits
  // -------------------------

  template< class Grid,
            class ProblemTraits,
            int polOrder,
            class HandlerTraits >
  struct EvolutionAlgorithmTraits
  {
    static const int polynomialOrder = polOrder;

    typedef ProblemTraits                                         ProblemTraitsType;

    // type of Grid
    typedef Grid                                          GridType;
    typedef typename ProblemTraits :: HostGridPartType    HostGridPartType;
    typedef typename ProblemTraits :: GridPartType        GridPartType;

    typedef typename ProblemTraits::template AnalyticalTraits< HostGridPartType >  AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< HostGridPartType, polynomialOrder >  DiscreteTraits;

    // obtain the problem dependent types, analytical context
    typedef typename AnalyticalTraits::ModelType                   ModelType;
    typedef typename AnalyticalTraits::ProblemType                 ProblemType;
    typedef typename AnalyticalTraits::InitialDataType             InitialDataType;

    // type of discrete function space and discrete function
    typedef typename DiscreteTraits::InitialProjectorType          InitialProjectorType;

    // type of dg operator
    typedef typename DiscreteTraits::OperatorType                  OperatorType;

    typedef typename ProblemTraits::FunctionSpaceType              FunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;

    typedef typename DiscreteTraits::ExtraParameterTuple           ExtraParameterTupleType;
    typedef typename DiscreteTraits::IOTupleType                   IOTupleType;

    // wrap operator
    typedef GridTimeProvider< GridType >                           TimeProviderType;

    typedef typename DiscreteTraits::OdeSolverType                 OdeSolverType;

    typedef typename DiscreteTraits::BasicLinearSolverType         BasicLinearSolverType;

    typedef SolverMonitor<1>                                       SolverMonitorType;

    typedef typename HandlerTraits::DiagnosticsHandlerType         DiagnosticsHandlerType;
    typedef typename HandlerTraits::SolverMonitorHandlerType       SolverMonitorHandlerType;
    typedef typename HandlerTraits::CheckPointHandlerType          CheckPointHandlerType;
    typedef typename HandlerTraits::DataWriterHandlerType          DataWriterHandlerType;
    typedef typename HandlerTraits::AdditionalOutputHandlerType    AdditionalOutputHandlerType;
    typedef typename HandlerTraits::SolutionLimiterHandlerType     SolutionLimiterHandlerType;
    typedef typename HandlerTraits::AdaptHandlerType               AdaptHandlerType;
  };


  template< class Grid, class ProblemTraits, int polOrder >
  struct NoHandlerTraits
  {
    typedef NoDiagnosticsHandler                                                         DiagnosticsHandlerType;
    typedef NoSolverMonitorHandler                                                       SolverMonitorHandlerType;
    typedef NoCheckPointHandler< Grid >                                                  CheckPointHandlerType;
    typedef NoDataWriterHandler                                                          DataWriterHandlerType;
    typedef NoAdditionalOutputHandler                                                    AdditionalOutputHandlerType;
    typedef NoSolutionLimiterHandler                                                     SolutionLimiterHandler;
    typedef NoAdaptHandler                                                               AdaptHandlerType;
  };


  template< class EvolutionAlgorithmTraits >
  class EvolutionAlgorithmBase;

  template< class Grid, class ProblemTraits, int polOrder, class HandlerTraits = typename ProblemTraits::template DiscreteTraits< Grid, polOrder >::HandlerTraits >
  class EvolutionAlgorithm
    : public EvolutionAlgorithmBase< EvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder, HandlerTraits > >
  {
    typedef EvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder, HandlerTraits >     Traits;
    typedef EvolutionAlgorithmBase< Traits >                                             BaseType;
  public:
    EvolutionAlgorithm( Grid &grid, const std::string name = "" )
    : BaseType( grid, name  )
    {}
  };

  template< class Grid, class ProblemTraits, int polOrder >
  struct EvolutionAlgorithmNoHandler
  : public EvolutionAlgorithm< Grid, ProblemTraits, polOrder, NoHandlerTraits< Grid, ProblemTraits, polOrder > >
  {
    typedef EvolutionAlgorithm< Grid, ProblemTraits, polOrder, NoHandlerTraits< Grid, ProblemTraits, polOrder > >
                                                                                           BaseType;

    EvolutionAlgorithmNoHandler( Grid& grid, const std::string name = "" )
      : BaseType( grid, name )
    {}

  };

  // EvolutionAlgorithm
  // ------------------

  template< class EvolutionAlgorithmTraits >
  class EvolutionAlgorithmBase
    : public AlgorithmBase< EvolutionAlgorithmTraits >
  {
    typedef EvolutionAlgorithmTraits                                                    Traits;
    typedef AlgorithmBase< Traits >                                                     BaseType;

  public:
    typedef typename BaseType::GridType                          GridType;
    typedef typename BaseType::IOTupleType                       IOTupleType;
    typedef typename BaseType::SolverMonitorType                 SolverMonitorType;
    typedef typename Traits::TimeProviderType                    TimeProviderType;

    typedef typename Traits::HostGridPartType                    HostGridPartType;

    // An analytical version of our model
    typedef typename Traits::ModelType                           ModelType;
    typedef typename Traits::ProblemType                         ProblemType;

    typedef typename Traits::GridPartType                        GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType           DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType                DiscreteFunctionType;

    typedef typename Traits::BasicLinearSolverType               BasicLinearSolverType;

    // The DG space operator
    typedef typename Traits::OperatorType                        OperatorType;

    // The ODE Solvers
    typedef typename Traits::OdeSolverType                       OdeSolverType;

    // type of initial interpolation
    typedef typename Traits::InitialProjectorType                InitialProjectorType;

    // analytical Tratis
    typedef typename Traits::AnalyticalTraits                    AnalyticalTraits;

    // discrete Traits
    typedef typename Traits::DiscreteTraits                      DiscreteTraits;

    typedef typename Traits::DiscreteTraits::GridExactSolutionType
                                                                 GridExactSolutionType;

    typedef typename Traits::DiagnosticsHandlerType                     DiagnosticsHandlerType;
    typedef typename Traits::SolverMonitorHandlerType                   SolverMonitorHandlerType;
    typedef typename Traits::CheckPointHandlerType                      CheckPointHandlerType;
    typedef typename Traits::DataWriterHandlerType                      DataWriterHandlerType;
    typedef typename Traits::AdditionalOutputHandlerType                AdditionalOutputHandlerType;
    typedef typename Traits::SolutionLimiterHandlerType                 SolutionLimiterHandlerType;
    typedef typename Traits::AdaptHandlerType                           AdaptHandlerType;

    typedef typename Traits::ExtraParameterTupleType  ExtraParameterTupleType;

    typedef uint64_t                                              UInt64Type ;
    typedef StepperParameters                                     StepperParametersType;

    using BaseType::grid_;
    using BaseType::eocParam_;
    using BaseType::grid;

    EvolutionAlgorithmBase ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  ),
      param_( StepperParametersType( Dune::ParameterKey::generate( "", "femdg.stepper." ) ) ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      solution_( "U_"+name, space() ),
      problem_( Traits::ProblemTraitsType::problem() ),
      model_( *problem_ ),
      checkPointHandler_( grid_, "" ),
      dataWriterHandler_( grid_, "" ),
      diagnosticsHandler_(),
      solverMonitorHandler_( "" ),
      additionalOutputHandler_( space_ ),
      solutionLimiterHandler_( name ),
      adaptHandler_( gridPart_, solution(), solutionLimiterHandler_, name ),
      overallTimer_(),
      odeSolver_(),
      timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
      fixedTimeStep_( param_.fixedTimeStep() )
    {
      if( problem().hasExactSolution() )
      {
        exact_.reset( new DiscreteFunctionType("exact solution", space() ) );
      }
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type grSize = grid().size( 0 );
      return grid().comm().sum( grSize );
    }

    virtual bool checkDofsValid( TimeProviderType& tp, const int loop  ) const
    {
      return solution_.dofsValid();
    }

    //! default time loop implementation, overload for changes in derived classes !!!
    void solve ( const int loop )
    {
      // get start and end time from parameter file
      const double startTime = param_.startTime();
      const double endTime   = param_.endTime();

      // Initialize TimeProvider
      TimeProviderType tp( startTime, this->grid() );

      // call solve implementation taking start and end time
      solve( loop, tp, endTime );
    }

    //! default time loop implementation, overload for changes in derived classes !!!
    virtual void solve ( const int loop, TimeProviderType& tp, const double endTime )
    {
      // get grid reference
      GridType& grid = this->grid();

      // print info on each printCount step
      const int printCount = param_.printCount();

      double maxTimeStep = param_.maxTimeStep();

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      // in codegen modus make endTime large and only compute one timestep
      const int maximalTimeSteps = 1;
#else
      // if this variable is set then only maximalTimeSteps timesteps will be computed
      const int maximalTimeSteps = param_.maximalTimeSteps();
#endif

      // set initial data (and create ode solver)
      // -> calls initializeStep()
      preInitializeStep( tp, loop );

      // start first time step with prescribed fixed time step
      // if it is not 0 otherwise use the internal estimate
      tp.provideTimeStepEstimate(maxTimeStep);

      // adjust fixed time step with timeprovider.factor()
      const double fixedTimeStep = fixedTimeStep_/tp.factor() ;
      if ( fixedTimeStep > 1e-20 )
        tp.init( fixedTimeStep );
      else
        tp.init();

      // true if last time step should match end time
      const bool stopAtEndTime =  param_.stopAtEndTime();

      //******************************
      //*  Time Loop                 *
      //******************************
      for( ; tp.time() < endTime; )
      {
        // write data for current time
        dataWriterHandler_.step( tp );

        // possibly write check point (default is disabled)
        checkPointHandler_.step( tp );

        // reset time step estimate
        tp.provideTimeStepEstimate( maxTimeStep );

        // current time step size
        const double deltaT = tp.deltaT();

        //************************************************
        //* Compute an ODE timestep                      *
        //************************************************
        Dune::FemTimer::start( timeStepTimer_ );
        overallTimer_.reset();

        // estimate, mark, adapt
        adaptHandler_.step( tp );

        // perform the solve for one time step, i.e. solve ODE
        step( tp );

        // stop FemTimer for this time step
        Dune::FemTimer::stop(timeStepTimer_,Dune::FemTimer::max);

        // Check that no NAN have been generated
        if( !checkDofsValid( tp, loop ) )
        {
          dataWriterHandler_.finalize( tp );
          std::abort();
        }

        if( (printCount > 0) && (((tp.timeStep()) % printCount) == 0))
        {
          UInt64Type grSize = gridSize();
          if( grid.comm().rank() == 0 )
          {
            std::cout << "step: " << tp.timeStep()-1 << "  time = " << tp.time()+tp.deltaT() << ", dt = " << deltaT
                      <<",  grid size: " << grSize << ", elapsed time: ";
            Dune::FemTimer::print(std::cout,timeStepTimer_);
            solverMonitorHandler_.stepPrint();
          }
        }

        // next advance should not exceed endtime
        if( stopAtEndTime )
          tp.provideTimeStepEstimate( (endTime - tp.time()) );

        // next time step is prescribed by fixedTimeStep
        if ( fixedTimeStep > 1e-20 )
          tp.next( fixedTimeStep );
        else
          tp.next();

        // for debugging and codegen only
        if( tp.timeStep() >= maximalTimeSteps )
        {
          if( Fem::Parameter::verbose() )
            std::cerr << "ABORT: time step count reached max limit of " << maximalTimeSteps << std::endl;
          break ;
        }

        if (tp.timeStep()<2)
        {
          // write parameters used (before simulation starts)
          Fem::Parameter::write("parameter.log");
        }
      } /****** END of time loop *****/

      // finalize eoc step
      finalizeStep( tp );

      // prepare the fixed time step for the next eoc loop
      fixedTimeStep_ /= param_.fixedTimeStepEocLoopFactor();
    }

    //! finalize problem, i.e. calculated EOC ...
    virtual void finalize ( const int eocloop )
    {}

    virtual SolverMonitorType& monitor()
    {
      return solverMonitorHandler_.monitor();
    }

    // function creating the ode solvers
    virtual OdeSolverType* createOdeSolver( TimeProviderType& ) = 0;

    // return reference to the discrete function space
    const DiscreteFunctionSpaceType &space () const
    {
      return space_;
    }

    // return reference to discrete function holding solution
    DiscreteFunctionType &solution ()
    {
      return solution_;
    }

    virtual IOTupleType dataTuple ()
    {
      return std::make_tuple( &solution(), exact_.operator->() );
    }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const
    {
      return problem().dataPrefix();
    }

    // before data initialization
    void preInitializeStep ( TimeProviderType &tp, int loop )
    {
      // restoreData if checkpointing is enabled (default is disabled)
      bool newStart = ( eocParam_.steps() == 1) ? checkPointHandler_.restoreData( tp ) : false;

      initializeStep( tp, loop );

      if( adaptHandler_.adaptive() && newStart )
      {
        // adapt the grid to the initial data
        for( int startCount = 0; startCount < adaptHandler_.finestLevel(); ++ startCount )
        {
          // call initial adaptation
          adaptHandler_.init();

          // setup problem again
          initializeStep( tp, loop );

          // some info in verbose mode
          if( Fem::Parameter::verbose() )
          {
            std::cout << "Start adaptation: step " << startCount << ",  dt = " << tp.deltaT() << ",  grid size: " << gridSize()
                      << std::endl;
          }
        }
      }
      dataWriterHandler_.init( tp, dataTuple(), eocParam_.dataOutputParameters( loop, problem_->dataPrefix() ) );

      // register data functions to check pointer
      checkPointHandler_.registerData( solution() );
    }

    // before first step, do data initialization
    virtual void initializeStep ( TimeProviderType &tp, int loop )
    {
      DiscreteFunctionType& U = solution();

      odeSolver_.reset( this->createOdeSolver( tp ) );
      assert( odeSolver_ );

      // communication is needed when blocking communication is used
      // but has to be avoided otherwise (because of implicit solver)
      const bool doCommunicate = ! NonBlockingCommParameter :: nonBlockingCommunication ();

      // create projection
      InitialProjectorType projection( 2 * U.space().order(), doCommunicate );

      // project initial data
      projection( problem().fixedTimeFunction( tp.time() ), U );

      // ode.initialize applies the DG Operator once to get an initial
      // estimate on the time step. This does not change the initial data u.
      odeSolver_->initialize( U );
    }


    virtual void step ( TimeProviderType &tp )
    {
      DiscreteFunctionType& U = solution();

      // solve ODE
      assert(odeSolver_);

      typename OdeSolverType::MonitorType odeSolverMonitor;
      odeSolver_->solve( U, odeSolverMonitor );

      // limit solution if necessary
      solutionLimiterHandler_.step( U );

      // copy information to solver monitor and update time step information
      solverMonitorHandler_.step( odeSolverMonitor, tp );

      diagnosticsHandler_.step( tp, solution(), odeSolverMonitor, overallTimer_, adaptHandler_ );
    }

    void finalizeStep ( TimeProviderType &tp )
    {
      DiscreteFunctionType& u = solution();

      // flush diagnostics data
      diagnosticsHandler_.finalize();

      // setup exact solution
      if( exact_ )
      {
        // TODO: Add this projection to addtional output writer, later
        Dune::Fem::DGL2ProjectionImpl::project( problem().exactSolution( tp.time() ), *exact_ );
      }

      // write last time step
      dataWriterHandler_.finalize( tp );

      // adjust average time step size
      solverMonitorHandler_.finalize( gridWidth(), gridSize() );

      AnalyticalTraits::addEOCErrors( tp, u, model(), problem() );

      // delete ode solver
      odeSolver_.reset();

      adaptHandler_.finalize();
    }

    std::string description () const { return problem().description(); }

    ProblemType &problem ()
    {
      assert( problem_ );
      return *problem_;
    }

    const ProblemType &problem () const
    {
      assert( problem_ );
      return *problem_;
    }

    const ModelType &model () const { return model_; }
    ModelType &model () { return model_; }

    GridPartType &gridPart () { return gridPart_; }


  protected:
    OdeSolverType &odeSolver ()
    {
      assert( odeSolver_ );
      return *odeSolver_;
    }

    StepperParametersType          param_;
    GridPartType                   gridPart_;
    DiscreteFunctionSpaceType      space_;

    // the solution
    DiscreteFunctionType           solution_;

    // InitialDataType evaluates to $u_0$
    std::unique_ptr< ProblemType > problem_;
    ModelType                      model_;
    std::unique_ptr< DiscreteFunctionType > exact_;

    CheckPointHandlerType          checkPointHandler_;
    DataWriterHandlerType          dataWriterHandler_;
    DiagnosticsHandlerType         diagnosticsHandler_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    AdditionalOutputHandlerType    additionalOutputHandler_;
    SolutionLimiterHandlerType     solutionLimiterHandler_;
    AdaptHandlerType               adaptHandler_;

    Dune::Timer                    overallTimer_;
    std::unique_ptr< OdeSolverType > odeSolver_;
    unsigned int                   timeStepTimer_;
    double                         fixedTimeStep_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
