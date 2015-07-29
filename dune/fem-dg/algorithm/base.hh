#ifndef DUNE_FEMDG_ALGORITHM_BASE_HH
#define DUNE_FEMDG_ALGORITHM_BASE_HH

// system includes
#include <sstream>
#include <sys/resource.h>

// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// local includes
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/misc/femeoc.hh>

namespace Dune
{

namespace Fem
{
  struct EocDataOutputParameters :
         public Fem::LocalParameter<Fem::DataWriterParameters,EocDataOutputParameters>
  {
    protected:
    std::string loop_;
    public:
    EocDataOutputParameters(int loop, const std::string& name) {
      std::stringstream ss;
      ss << name << loop;
      loop_ = ss.str();
    }
    EocDataOutputParameters(const EocDataOutputParameters& other)
    : loop_(other.loop_) {}

    std::string path() const {
      return loop_;
    }
  };

  struct EocParameters :
         public Fem::LocalParameter<EocParameters,EocParameters>
  {
    protected:
    std::string keyPrefix_;

    public:
    EocParameters( const std::string keyPrefix = "fem.eoc.")
      : keyPrefix_( keyPrefix )
    {}

    const EocDataOutputParameters dataOutputParameters( int loop, const std::string& name ) const
    {
      return EocDataOutputParameters( loop, name );
    }

    virtual int steps() const
    {
      return Fem::Parameter::getValue<int>( keyPrefix_ + "steps", 1);
    }

    virtual std::string outputPath() const
    {
      return Fem::Parameter::getValue<std::string>( keyPrefix_ + "outputpath", "./");
    }

    virtual std::string fileName() const
    {
      return Fem::Parameter::getValue<std::string>( keyPrefix_ + "filename", "eoc" );
    }

    virtual int quadOrder() const
    {
      return Fem::Parameter::getValue< int >( keyPrefix_ + "quadorder", -1 );
    }

  };


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


  //! get memory in MB
  inline double getMemoryUsage()
  {
    struct rusage info;
    getrusage( RUSAGE_SELF, &info );
    return (info.ru_maxrss / 1024.0);
  }


  // AlgorithmBase
  // -------------

  template< class AlgorithmTraits >
  class AlgorithmBase
  {
    typedef AlgorithmTraits Traits;

  public:
    // type of Grid
    typedef typename Traits::GridType            GridType;

    // type of IOTuple
    typedef typename Traits::IOTupleType         IOTupleType;

    // type of statistics monitor
    typedef typename Traits::SolverMonitorType   SolverMonitorType;

    //! constructor
    AlgorithmBase ( GridType &grid, const std::string algorithmName = "" )
      : grid_( grid ),
        algorithmName_( algorithmName )
    {}

    virtual ~AlgorithmBase () {}

    //! return reference to hierarchical grid
    GridType &grid () const { return grid_; }

    // return size of grid
    virtual std::size_t gridSize () const
    {
      std::size_t grSize = grid_.size( 0 );
      return grid_.comm().sum( grSize );
    }

    //! return default data tuple for data output
    virtual IOTupleType dataTuple() = 0;

    // solve the problem for eoc loop 'loop' and return statistics
    virtual SolverMonitorType solve ( const int loop ) = 0;

    //! finalize problem, i.e. calculated EOC ...
    virtual void finalize ( const int eocloop ) = 0;

    //! name of the scheme
    virtual std::string description () const = 0;

    virtual const std::string name() { return algorithmName_; }

  protected:
    GridType &grid_;
    const std::string algorithmName_;
  };


  // EvolutionAlgorithmBase
  // ----------------------

  template< class AlgorithmTraits >
  class EvolutionAlgorithmBase
    : public AlgorithmBase< AlgorithmTraits >
  {
    typedef AlgorithmBase< AlgorithmTraits > BaseType;
    typedef EvolutionAlgorithmBase< AlgorithmTraits > ThisType;

  public:
    typedef typename BaseType::GridType                  GridType;
    typedef typename BaseType::IOTupleType               IOTupleType;
    typedef typename BaseType::SolverMonitorType         SolverMonitorType;
    typedef GridTimeProvider< GridType >                 TimeProviderType;

    typedef uint64_t                                     UInt64Type ;
    typedef StepperParameters                            StepperParametersType;
    typedef EocParameters                                EocParametersType;
    typedef AdaptationParameters                         AdaptationParametersType;

    using BaseType::grid;

    //! constructor
    EvolutionAlgorithmBase ( GridType &grid, const std::string algorithmName = "" )
      : BaseType( grid, algorithmName ),
        param_( StepperParametersType( Dune::ParameterKey::generate( "", "femdg.stepper." ) ) ),
        eocParam_( EocParametersType( Dune::ParameterKey::generate( "", "fem.eoc." ) ) ),
        adaptParam_( AdaptationParametersType( Dune::ParameterKey::generate( "", "fem.adaptation." ) ) ),
        timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
        fixedTimeStep_( param_.fixedTimeStep() )
    {}

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth() const { return 0.0; }

    // return size of grid
    virtual UInt64Type gridSize() const
    {
      UInt64Type grSize = grid().size( 0 );
      return grid().comm().sum( grSize );
    }

    // --------- CheckPointing -----------
    virtual bool checkPointing() const { return false; }
    virtual bool doCheckPointRestoreData( const GridType& grid, TimeProviderType& tp ) const { return false; }
    virtual void doCheckPointWriteData( TimeProviderType& tp ) const { }
    void checkPointWriteData( TimeProviderType& tp ) const { if( checkPointing() ) doCheckPointWriteData( tp ); }
    bool checkPointRestoreData( const GridType& grid, TimeProviderType& tp ) const
    {
      if( checkPointing() )
        return doCheckPointRestoreData( grid, tp );
      return false;
    }
    // -----------------------------------

    // ---------- Adaptation ------------
    virtual bool adaptive () const { return false ; }
    //! estimate and mark cell for refinement/coarsening (default does nothing)
    virtual void estimateMarkAdapt( ) { }
    //! initial estimate and mark cell for refinement/coarsening
    //! the default redirects to estimateMarkAdapt
    virtual void initialEstimateMarkAdapt( ) { estimateMarkAdapt( ); }

    template <class IndicatorOperator, class GradientIndicator>
    void doEstimateMarkAdapt( const IndicatorOperator& dgIndicator,
                              GradientIndicator& gradientIndicator,
                              const bool initialAdaptation = false ){}
    // -----------------------------------

    // ---------- DataWriter ---------------
    //! write data, this should be overloaded in the derived implementation (default does nothing)
    virtual void writeData ( TimeProviderType &tp, const bool reallyWrite = false ) {}
    // -----------------------------------

    virtual bool checkDofsValid( TimeProviderType& tp, const int loop  ) const { return true; }

    //! initialize method for time loop, i.e. L2-project initial data
    virtual void initializeStep( TimeProviderType &tp, const int loop, SolverMonitorType &monitor ) = 0;

    //! solve one time step
    virtual void step( TimeProviderType &tp, SolverMonitorType &monitor ) = 0;

    //! call limiter if necessary
    virtual void limitSolution () {}

    //! finalize this EOC loop and possibly calculate EOC ...
    virtual void finalizeStep( TimeProviderType &tp ) = 0;

    //! default time loop implementation, overload for changes in derived classes !!!
    virtual SolverMonitorType solve ( const int loop )
    {
      // get start and end time from parameter file
      const double startTime = param_.startTime();
      const double endTime   = param_.endTime();

      // Initialize TimeProvider
      TimeProviderType tp( startTime, this->grid() );

      // call solve implementation taking start and end time
      return solve( loop, tp, endTime );
    }

    //! default time loop implementation, overload for changes in derived classes !!!
    virtual SolverMonitorType solve ( const int loop, TimeProviderType& tp, const double endTime  )
    {
      // get grid reference
      GridType& grid = this->grid();

      // verbosity
      const bool verbose = Fem::Parameter::verbose ();
      // print info on each printCount step
      const int printCount = param_.printCount();

      double maxTimeStep = param_.maxTimeStep();

  #ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      // in codegen modus make endTime large and only compute one timestep
      // const double endTime = startTime + 1e8;
      const int maximalTimeSteps = 1;
  #else
      // if this variable is set then only maximalTimeSteps timesteps will be computed
      const int maximalTimeSteps = param_.maximalTimeSteps();
  #endif

      // create monitor object (initialize all varialbes)
      SolverMonitorType monitor;

      // restoreData if checkpointing is enabled (default is disabled)
      const bool newStart = checkPointRestoreData( grid, tp );

      // set initial data (and create ode solver)
      initializeStep( tp, loop, monitor );

      // start first time step with prescribed fixed time step
      // if it is not 0 otherwise use the internal estimate
      tp.provideTimeStepEstimate(maxTimeStep);

      // adjust fixed time step with timeprovider.factor()
      const double fixedTimeStep = fixedTimeStep_/tp.factor() ;

      if ( fixedTimeStep > 1e-20 )
        tp.init( fixedTimeStep );
      else
        tp.init();

      // for simulation new start do start adaptation
      if( adaptive() && newStart )
      {
        // adapt the grid to the initial data
        for( int startCount = 0; startCount < adaptParam_.finestLevel(); ++ startCount )
        {
          // call initial adaptation
          initialEstimateMarkAdapt( );

          // setup problem again
          initializeStep( tp, loop, monitor );

          // get grid size (outside of verbose if)
          UInt64Type grSize = gridSize();
          // some info in verbose mode
          if( verbose )
          {
            std::cout << "Start adaptation: step " << startCount << ",  dt = " << tp.deltaT() << ",  grid size: " << grSize
                      << std::endl;
          }
        }
      }

      // true if last time step should match end time
      const bool stopAtEndTime =  param_.stopAtEndTime();

      //******************************
      //*  Time Loop                 *
      //******************************
      for( ; tp.time() < endTime; )
      {
        // write data for current time
        writeData( tp );

        // possibly write check point (default is disabled)
        checkPointWriteData( tp );

        // reset time step estimate
        tp.provideTimeStepEstimate( maxTimeStep );

        // current time step size
        const double deltaT = tp.deltaT();

        // current time step number
        const int timeStep  = tp.timeStep();

        //************************************************
        //* Compute an ODE timestep                      *
        //************************************************
        Dune::FemTimer::start( timeStepTimer_ );

        if( adaptive() )
        {
          // grid adaptation (including marking of elements)
          if( timeStep % adaptParam_.adaptCount() == 0 )
            estimateMarkAdapt();
        }

        // perform the solve for one time step, i.e. solve ODE
        step( tp, monitor );

        // stop FemTimer for this time step
        Dune::FemTimer::stop(timeStepTimer_,Dune::FemTimer::max);

        // Check that no NAN have been generated
        if( !checkDofsValid( tp, loop ) )
        {
          writeData( tp, true );
          std::abort();
        }

        if( (printCount > 0) && (((timeStep+1) % printCount) == 0))
        {
          UInt64Type grSize = gridSize();
          if( grid.comm().rank() == 0 )
          {
            std::cout << "step: " << timeStep << "  time = " << tp.time()+tp.deltaT() << ", dt = " << deltaT
                      <<",  grid size: " << grSize << ", elapsed time: ";
            Dune::FemTimer::print(std::cout,timeStepTimer_);
            std::cout << ",  Newton: " << *monitor.newton_iterations
                      << "  ILS: " << *monitor.ils_iterations
                      << "  OC: " << *monitor.operator_calls << std::endl;

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

      // write last time step
      writeData( tp, true );

      // finalize eoc step
      finalizeStep( tp );

      // prepare the fixed time step for the next eoc loop
      fixedTimeStep_ /= param_.fixedTimeStepEocLoopFactor();

      // adjust average time step size
      monitor.finalize( gridWidth(), gridSize() );

      return monitor;
    }

    //! finalize problem, i.e. calculated EOC ...
    virtual void finalize ( const int eocloop )
    {}

  protected:
    StepperParametersType param_;
    EocParametersType eocParam_;
    AdaptationParametersType adaptParam_;
    unsigned int timeStepTimer_;
    double fixedTimeStep_;
  };

  /////////////////////////////////////////////////////////////////////////////
  //
  //  compute algorithm
  //
  /////////////////////////////////////////////////////////////////////////////
  template <class Algorithm>
  void compute(Algorithm& algorithm)
  {
    const std::string name = "";

    typedef typename Algorithm::GridType GridType;
    GridType& grid = algorithm.grid();

    typename Algorithm::IOTupleType dataTup = algorithm.dataTuple() ;

    typedef Fem::DataOutput<GridType, typename Algorithm::IOTupleType> DataOutputType;
    DataOutputType dataOutput( grid, dataTup );

    typedef typename Algorithm::SolverMonitorType SolverMonitorType;

    // initialize FemEoc if eocSteps > 1
    EocParameters eocParam( Dune::ParameterKey::generate( "", "fem.eoc." ) );
    FemEoc::clear();
    FemEoc::initialize( eocParam.outputPath(), eocParam.fileName(), algorithm.description() );

    const unsigned int femTimerId = FemTimer::addTo("timestep");
    for(int eocloop=0; eocloop < eocParam.steps(); ++eocloop )
    {
      // start fem timer (enable with -DFEMTIMER)
      FemTimer :: start(femTimerId);

      // additional timer
      Dune::Timer timer ;

      // call algorithm and return solver statistics and some info
      SolverMonitorType monitor = algorithm.solve( eocloop );

      // get run time
      const double runTime = timer.elapsed();

      // also get times for FemTimer if enabled
      FemTimer::stop(femTimerId);
      FemTimer::printFile("./timer.out");
      FemTimer::reset(femTimerId);

      // finalize the algorithm
      algorithm.finalize( eocloop );

      std::stringstream eocInfo ;
      // generate EOC information
      Fem::FemEoc::write( *monitor.gridWidth, *monitor.elements, runTime, *monitor.timeSteps,
                          monitor.doubleValues(), monitor.intValues(), eocInfo );

      // in verbose mode write EOC info to std::cout
      if( Fem::Parameter::verbose() )
      {
        std::cout << std::endl << "EOC info: " << std::endl << eocInfo.str() << std::endl;
      }

      // write eoc step
      dataOutput.writeData( eocloop );

      // Refine the grid for the next EOC Step. If the scheme uses adaptation,
      // the refinement level needs to be set in the algorithms' initialize method.
      if( eocParam.steps() > 1 )
        Fem::GlobalRefine::apply(grid,Dune::DGFGridInfo<GridType>::refineStepsForHalf());

    } /***** END of EOC Loop *****/

    // FemTimer cleanup
    FemTimer::removeAll();
  }




} // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_BASE_HH
