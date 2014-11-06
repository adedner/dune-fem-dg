#ifndef FEMHOWTO_BASEEVOL_HH
#define FEMHOWTO_BASEEVOL_HH

// system includes
#include <sstream>

// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/gridwidth.hh>

// local includes
#include <dune/fem-dg/stepper/base.hh>
#include <dune/fem-dg/misc/cons2prim.hh>

/////////////////////////////////////////////////////////////////////////////
//
//  EOC output parameter class 
//
/////////////////////////////////////////////////////////////////////////////

struct EocDataOutputParameters :   /*@LST1S@*/
       public Dune::Fem::LocalParameter<Dune::Fem::DataWriterParameters,EocDataOutputParameters> 
{
  std::string loop_;
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

//! solver statistics and further info 
//! such as newton iterations and iterations of the linear solver 
class SolverMonitor 
{
public:
  double gridWidth;
  double avgTimeStep;
  double minTimeStep;
  double maxTimeStep;
  size_t timeSteps;
  int newton_iterations;
  int ils_iterations;
  int total_newton_iterations;
  int total_ils_iterations;
  int max_newton_iterations;
  int max_ils_iterations;

  unsigned long elements ;

  SolverMonitor() 
  {
    gridWidth = 0;
    avgTimeStep = 0;
    minTimeStep = std::numeric_limits<double>::max();
    maxTimeStep = 0;
    timeSteps = 0;
    newton_iterations = 0;
    ils_iterations = 0;
    total_newton_iterations = 0;
    total_ils_iterations = 0;
    max_newton_iterations = 0;
    max_ils_iterations = 0;
    elements = 0;
  }

  void setTimeStepInfo( const Dune::Fem::TimeProviderBase& tp ) 
  { 
    const double ldt = tp.deltaT();
    // calculate time step info 
    minTimeStep  = std::min( ldt, minTimeStep );
    maxTimeStep  = std::max( ldt, maxTimeStep );
    avgTimeStep += ldt;

    timeSteps = tp.timeStep() + 1;

    // accumulate iteration information 
    total_newton_iterations += newton_iterations;
    total_ils_iterations += ils_iterations;
  }

  void finalize( const double h = 0, const unsigned long el = 0 )
  {
    avgTimeStep /= double( timeSteps );
    gridWidth = h;
    elements = el;
  }

  void dump( std::ostream& out ) const
  {
    out << "SolverMonitorInfo (timesteps):" << std::endl;
    out << "min  dt: " << minTimeStep << std::endl;
    out << "max  dt: " << maxTimeStep << std::endl;
    out << "aver dt: " << avgTimeStep << std::endl;

    if( max_newton_iterations > 0 || max_ils_iterations > 0 )
    {
      out << "SolverMonitorInfo (iterations):" << std::endl;
      out << "newton max: " << max_newton_iterations << std::endl;
      out << "newton tot: " << total_newton_iterations << std::endl;
      out << "ils    max: " << max_ils_iterations << std::endl;
      out << "ils    tot: " << total_ils_iterations << std::endl;
    }
    out << std::endl;
  }
};

// ostream operator for SolverMonitor
inline std::ostream& operator << ( std::ostream& out, const SolverMonitor& monitor ) 
{
  monitor.dump( out );
  return out;
}

/////////////////////////////////////////////////////////////////////////////
//
//  Base class for evolutionary problems 
//
/////////////////////////////////////////////////////////////////////////////
template <class TraitsImp> 
class AlgorithmBase
{
  typedef TraitsImp Traits ;
public:
  // type of Grid
  typedef typename Traits :: GridType                  GridType;

  // choose a suitable GridView
  typedef typename Traits :: GridPartType              GridPartType;

  // the indicator function type (for limiting only)
  typedef typename Traits :: IndicatorType             IndicatorType;

  // type of IOTuple 
  typedef typename Traits :: IOTupleType  IOTupleType ;

  // type of time provider organizing time for time loops 
  typedef Dune::Fem::GridTimeProvider< GridType >                 TimeProviderType;

  // type of 64bit unsigned integer 
  typedef uint64_t UInt64Type ;

  // type of statistics monitor 
  typedef SolverMonitor  SolverMonitorType ;

  // type of parameter class 
  typedef Dune::Fem::Parameter  ParameterType ;

  //! constructor 
  AlgorithmBase( GridType& grid ) 
   : grid_( grid ),
     // Initialize Timer for CPU time measurements
     timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
     fixedTimeStep_( ParameterType::getValue<double>("fixedTimeStep",0) ),
     fixedTimeStepEocLoopFactor_( ParameterType::getValue<double>("fixedTimeStepEocLoopFactor",1.) )
  {
  }

  // nothing to do here
  virtual ~AlgorithmBase() {} 

  //! return default data tuple for data output
  virtual IOTupleType dataTuple() = 0 ;

  //! return reference to hierarchical grid 
  GridType& grid() { return grid_ ; }

  // return grid width of grid (overload in derived classes)
  virtual double gridWidth() const { return 0.0; }

  // return size of grid 
  virtual UInt64Type gridSize() const 
  { 
    UInt64Type grSize = grid_.size( 0 );
    return grid_.comm().sum( grSize );
  }

  // return reference to additional variables 
  virtual IndicatorType* indicator()  { return 0; }

  //! initialize method for time loop, i.e. L2-project initial data 
  virtual void initializeStep( TimeProviderType& tp, const int loop ) = 0;

  //! solve one time step 
  virtual void step(TimeProviderType& tp,
                    SolverMonitorType& monitor ) = 0;

  //! call limiter if necessary 
  virtual void limitSolution () {} 

  //! finalize this EOC loop and possibly calculate EOC ...
  virtual void finalizeStep(TimeProviderType& tp) = 0;

  //! restore all data from check point (overload to do something)
  //! return true if the simulation was newly started 
  virtual bool restoreFromCheckPoint( TimeProviderType& tp ) { abort(); return true;  } 

  //! write a check point (overload to do something)
  virtual void writeCheckPoint(TimeProviderType& tp) const {}

  //! estimate and mark cell for refinement/coarsening (default does nothing)
  virtual void estimateMarkAdapt( ) { }

  //! initial estimate and mark cell for refinement/coarsening
  //! the default redirects to estimateMarkAdapt
  virtual void initialEstimateMarkAdapt( )
  {
    estimateMarkAdapt( );  
  }

  virtual void checkDofsValid( TimeProviderType& tp, const int loop  ) const {}

  virtual bool adaptive () const { return false ; }

  virtual void checkAdditionalVariables( TimeProviderType& tp ) {} 

  //! write data, this should be overloaded in the derived implementation (default does nothing)
  virtual void writeData( TimeProviderType& tp, 
                          const bool reallyWrite = false ) {}

  //! default time loop implementation, overload for changes in derived classes !!!
  virtual SolverMonitorType solve ( const int loop )
  {
    // get start and end time from parameter file
    const double startTime = ParameterType::getValue<double>("femhowto.startTime", 0.0);
    const double endTime   = ParameterType::getValue<double>("femhowto.endTime");

    // Initialize TimeProvider
    TimeProviderType tp( startTime, this->grid() );

    // call solve implementation taking start and end time
    return solve( loop, tp, endTime );
  }

  //! time loop implementation taking TimeProvider and endTime
  virtual SolverMonitorType solve ( const int loop,
                                    TimeProviderType& tp,
                                    const double endTime )
  {
    // get grid reference 
    GridType& grid = this->grid();

    // verbosity 
    const bool verbose = ParameterType :: verbose ();
    // print info on each printCount step 
    const int printCount = ParameterType::getValue<int>("femhowto.printCount", -1);

    double maxTimeStep =
      ParameterType::getValue("femhowto.maxTimeStep", std::numeric_limits<double>::max());

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
    // in codegen modus make endTime large and only compute one timestep
    // const double endTime = startTime + 1e8;
    const int maximalTimeSteps = 1;
#else 
    // if this variable is set then only maximalTimeSteps timesteps will be computed 
    const int maximalTimeSteps =
      ParameterType::getValue("femhowto.maximaltimesteps", std::numeric_limits<int>::max());
#endif

    // create monitor object (initialize all varialbes) 
    SolverMonitorType monitor; 

    // if adaptCount is 0 then no dynamics grid adaptation
    int adaptCount = 0;
    int maxAdaptationLevel = 0;
    if( adaptive() )
    {
      adaptCount = ParameterType::getValue<int>("fem.adaptation.adaptcount");
      maxAdaptationLevel = ParameterType::getValue<int>("fem.adaptation.finestLevel");
    }

    // only do checkpointing when number of EOC steps is 1 
    const bool doCheckPointing = ( InitFemEoc :: eocSteps() == 1 ); 

    // restoreData if checkpointing is enabled (default is disabled)
    const bool newStart = ( doCheckPointing ) ? restoreFromCheckPoint( tp ) : false ;

    // set initial data (and create ode solver)
    initializeStep( tp, loop );

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
    if( newStart && adaptCount > 0 ) 
    {
      // adapt the grid to the initial data
      for( int startCount = 0; startCount < maxAdaptationLevel; ++ startCount )
      {
        // call initial adaptation
        initialEstimateMarkAdapt( );

        // setup problem again 
        initializeStep( tp, loop );

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
    const bool stopAtEndTime =  ParameterType::getValue("femdg.stopatendtime", bool(false) );

    //******************************    
    //*  Time Loop                 *
    //******************************
    for( ; tp.time() < endTime; ) 
    {
      // write data for current time 
      writeData( tp );

      // only do checkpointing for loop 0, i.e. not in EOC calculations 
      if( doCheckPointing ) 
      {
        // possibly write check point (default is disabled)
        writeCheckPoint( tp );
      }

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

      // grid adaptation (including marking of elements)
      if( (adaptCount > 0) && (timeStep % adaptCount) == 0 )
        estimateMarkAdapt();

      // perform the solve for one time step, i.e. solve ODE 
      step( tp, monitor ); 

      // stop FemTimer for this time step
      Dune::FemTimer::stop(timeStepTimer_,Dune::FemTimer::max);

      // Check that no NAN have been generated
      checkDofsValid( tp, loop );

      if( (printCount > 0) && (timeStep % printCount == 0)) 
      {
        UInt64Type grSize = gridSize();
        if( grid.comm().rank() == 0 )
        {
          std::cout << "step: " << timeStep << "  time = " << tp.time() << ", dt = " << deltaT
                    <<",  grid size: " << grSize << ",  Newton: " << monitor.newton_iterations << "  ILS: " << monitor.ils_iterations << std::endl;
        }
      }

      // next advance should not exceed endtime
      if( stopAtEndTime ) 
        tp.provideTimeStepEstimate( (endTime - tp.time()) );

      // next time step is prescribed by fixedTimeStep
      // it fixedTimeStep is not 0
      if ( fixedTimeStep > 1e-20 )
        tp.next( fixedTimeStep );
      else
        tp.next();

      // for debugging and codegen only 
      if( tp.timeStep() >= maximalTimeSteps ) 
      {
        if( ParameterType :: verbose() )
          std::cerr << "ABORT: time step count reached max limit of " << maximalTimeSteps << std::endl;
        break ;
      }

      if (tp.timeStep()<2)
      {
        // write parameters used (before simulation starts)
        Dune::Fem::Parameter::write("parameter.log");
      }
    } /****** END of time loop *****/

    // write last time step  
    writeData( tp, true );

    // finalize eoc step 
    finalizeStep( tp ); 
    
    // prepare the fixed time step for the next eoc loop
    fixedTimeStep_ /= fixedTimeStepEocLoopFactor_; 

    // adjust average time step size 
    monitor.finalize( gridWidth(),  // h 
                      gridSize() ); // elements

    return monitor;
  }

  //! finalize problem, i.e. calculated EOC ...
  virtual void finalize( const int eocloop ) 
  {
  }

protected:      
  // reference to hierarchical grid 
  GridType& grid_ ;

  unsigned int timeStepTimer_; 

  // use fixed time step if fixedTimeStep>0
  double fixedTimeStep_;        
  double fixedTimeStepEocLoopFactor_;        
};

#endif
