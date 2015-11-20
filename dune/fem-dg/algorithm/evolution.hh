#ifndef DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH


#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem/misc/femtimer.hh>

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem-dg/pass/threadpass.hh>
#include <dune/common/timer.hh>

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

  // internal forward declarations
  // -----------------------------

  template< class Traits >
  class EvolutionAlgorithmBase;


  // StepperParametersType
  // ---------------------

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
  template< int polOrder, class ... ProblemTraits >
  struct EvolutionAlgorithmTraits
  {
    // type of Grid
    typedef typename std::tuple_element<0, std::tuple< ProblemTraits... > >::type::GridType
                                                                        GridType;

    // wrap operator
    typedef GridTimeProvider< GridType >                                TimeProviderType;

    //typedef ...
    typedef std::tuple< typename std::add_pointer< typename ProblemTraits::template Stepper<polOrder>::Type >::type... >
                                                                        StepperTupleType;

    //typedef typename Std::make_index_sequence_impl< std::tuple_size< StepperTupleType >::value >::type
    //                                                                  IndexSequenceType;

    typedef Dune::Fem::AdaptHandler< StepperTupleType >                 AdaptHandlerType;
    typedef Dune::Fem::CheckPointHandler< StepperTupleType >            CheckPointHandlerType;
    typedef Dune::Fem::SolverMonitorHandler< StepperTupleType >         SolverMonitorHandlerType;
    typedef Dune::Fem::DataWriterHandler< StepperTupleType >            DataWriterHandlerType;
    typedef Dune::Fem::DiagnosticsHandler< StepperTupleType >           DiagnosticsHandlerType;
    typedef Dune::Fem::SolutionLimiterHandler< StepperTupleType >       SolutionLimiterHandlerType;

    typedef typename DataWriterHandlerType::IOTupleType                 IOTupleType;

  };


  // EvolutionAlgorithm
  // ------------------

  template< int polOrder, class ... ProblemTraits >
  class EvolutionAlgorithm
  : public EvolutionAlgorithmBase< EvolutionAlgorithmTraits< polOrder, ProblemTraits ... > >
  {
    typedef EvolutionAlgorithmTraits< polOrder, ProblemTraits... > Traits;
    typedef EvolutionAlgorithmBase< Traits > BaseType;
  public:
    typedef typename BaseType::GridType GridType;

    EvolutionAlgorithm ( GridType &grid, const std::string name = "" )
      : BaseType( grid, name )
    {}
  };



  // EvolutionAlgorithmBase
  // ----------------------

  template< class Traits >
  class EvolutionAlgorithmBase
    : public AlgorithmInterface< Traits >
  {
    typedef AlgorithmInterface< Traits >                         BaseType;
  public:
    typedef typename BaseType::GridType                          GridType;
    typedef typename BaseType::IOTupleType                       IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType          SolverMonitorHandlerType;

    typedef typename Traits::StepperTupleType                    StepperTupleType;
    typedef typename Traits::TimeProviderType                    TimeProviderType;

    typedef typename Traits::DiagnosticsHandlerType              DiagnosticsHandlerType;
    typedef typename Traits::CheckPointHandlerType               CheckPointHandlerType;
    typedef typename Traits::DataWriterHandlerType               DataWriterHandlerType;
    typedef typename Traits::SolutionLimiterHandlerType          SolutionLimiterHandlerType;
    typedef typename Traits::AdaptHandlerType                    AdaptHandlerType;

    typedef uint64_t                                             UInt64Type ;

    typedef StepperParameters                                    StepperParametersType;

    using BaseType::eocParams;
    using BaseType::grid;

    static const int numSteppers = std::tuple_size< StepperTupleType >::value;

  private:
    struct Initialize {
    private:
      template<class T, class AdaptHandler, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::DiagnosticsHandlerType >::value >::type
      getDiagnostics( T, AdaptHandler&, Args&& ... ){}
      template<class T, class AdaptHandler, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::DiagnosticsHandlerType >::value >::type
      getDiagnostics( T e, AdaptHandler& handler, Args &&... a )
      {
        if( e->diagnostics() )
        {
          e->diagnostics()->registerData( "AdaptationTime", &handler.adaptationTime() );
          e->diagnostics()->registerData( "LoadBalanceTime", &handler.loadBalanceTime() );
        }
      }
    public:
      template< class T, class AdaptHandler, class ... Args > static void apply ( T& e, AdaptHandler& handler, Args && ... a )
      {
        e->initialize( std::forward<Args>(a)... );
        getDiagnostics( e, handler, std::forward<Args>(a)... );
      }
    };
    struct PreSolve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->preSolve( std::forward<Args>(a)... ); }
    };
    struct Solve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->solve( std::forward<Args>(a)... ); }
    };
    struct PostSolve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->postSolve( std::forward<Args>(a)... ); }
    };
    struct Finalize {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->finalize( std::forward<Args>(a)... ); }
    };
    struct GridWidth {
      template< class T, class ... Args > static void apply ( T& e, double& res, Args && ... a )
      { res = std::max( res, e->gridWidth(std::forward<Args>(a)... ) ); }
    };
    struct GridSize {
      template<class T, class... Args > static void apply ( T& e, UInt64Type& res, Args && ... a )
      { res = std::max( res, e->gridSize(std::forward<Args>(a)... ) ); }
    };
    struct CheckSolutionValid {
      template<class T, class... Args > static void apply( T e, bool& res, Args&& ... a )
      { res &= e->checkSolutionValid( std::forward<Args>(a)... ); }
    };
    template< class Caller >
    class LoopCallee
    {
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... Args >
        static void apply ( Tuple &tuple, Args&& ... a )
        {
          Caller::apply( std::get<i>( tuple ), std::forward<Args>(a)... );
        }
      };
    };

    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0,  numSteppers-1 >;

    template< int i >
    using ElementType=typename std::remove_pointer< typename std::tuple_element< i, StepperTupleType >::type >::type;

    template< std::size_t ...i >
    static StepperTupleType createStepper ( Std::index_sequence< i... >, GridType &grid )
    {
      static auto tuple = std::make_tuple( new ElementType<i>( grid, ElementType<i>::name() ) ... );
      return tuple;
    }

    // create Tuple of contained sub algorithms
    static StepperTupleType createStepper( GridType &grid )
    {
      return createStepper( typename Std::make_index_sequence_impl< std::tuple_size< StepperTupleType >::value >::type(), grid );
    }

  public:

    EvolutionAlgorithmBase ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  ),
      tuple_( createStepper( grid ) ),
      param_( StepperParametersType( ParameterKey::generate( "", "femdg.stepper." ) ) ),
      checkPointHandler_( tuple_ ),
      dataWriterHandler_( tuple_ ),
      diagnosticsHandler_( tuple_ ),
      solverMonitorHandler_( tuple_ ),
      solutionLimiterHandler_( tuple_ ),
      adaptHandler_( tuple_ ),
      fixedTimeStep_( param_.fixedTimeStep() )
    {}

    // return grid width of grid (overload in derived classes)
    double gridWidth () const
    {
      double res=0.0;
      ForLoopType< GridWidth >::apply( tuple_, res );
      return res;
    }

    // return size of grid
    UInt64Type gridSize () const
    {
      UInt64Type res=0;
      ForLoopType< GridSize >::apply( tuple_, res );
      return res;
    }

    bool checkSolutionValid( const int loop, TimeProviderType& tp ) const
    {
      bool res = true;
      ForLoopType< CheckSolutionValid >::apply( tuple_, res, loop, &tp );
      return res;
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
    void solve ( const int loop, TimeProviderType& tp, const double endTime )
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

      // restoreData if checkpointing is enabled (default is disabled)
      bool newStart = ( eocParams().steps() == 1) ? checkPointHandler_.restoreData( tp ) : false;

      initialize( loop, tp );

      if( adaptHandler_.adaptive() && newStart )
      {
        // adapt the grid to the initial data
        for( int startCount = 0; startCount < adaptHandler_.finestLevel(); ++ startCount )
        {
          // call initial adaptation
          adaptHandler_.init();

           // setup problem again
           initialize( loop, tp );

           // some info in verbose mode
           if( Fem::Parameter::verbose() )
           {
             std::cout << "Start adaptation: step " << startCount << ",  dt = " << tp.deltaT() << ",  grid size: " << gridSize()
                       << std::endl;
          }
        }
      }

      dataWriterHandler_.init( tp, eocParams().dataOutputParameters( loop, dataPrefix() ) );

      // register data functions to check pointer
      checkPointHandler_.registerData();

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

        preStep( loop, tp );

        // estimate, mark, adapt
        adaptHandler_.step( tp );

        // perform the solve for one time step, i.e. solve ODE
        step( loop, tp );

        // limit solution if necessary
        solutionLimiterHandler_.step();

        postStep( loop, tp );

        // update time step information
        solverMonitorHandler_.step( tp );

        //write diagnostics
        diagnosticsHandler_.step( tp );


        // stop FemTimer for this time step
        Dune::FemTimer::stop(timeStepTimer_,Dune::FemTimer::max);

        // Check that no NAN have been generated
        if( !checkSolutionValid( loop, tp ) )
        {
          dataWriterHandler_.finalize( tp );
          std::abort();
        }

        if( (printCount > 0) && (((tp.timeStep()) % printCount) == 0))
        {
          if( grid.comm().rank() == 0 )
          {
            std::cout << "step: " << tp.timeStep() << "  time = " << tp.time()+tp.deltaT() << ", dt = " << deltaT
                      <<",  grid size: " << gridSize() << ", elapsed time: ";
            Dune::FemTimer::print(std::cout,timeStepTimer_);
            solverMonitorHandler_.print( "Newton", "ILS", "OC" );
            std::cout << std::endl;
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
      finalize( loop, tp );

      // prepare the fixed time step for the next eoc loop
      fixedTimeStep_ /= param_.fixedTimeStepEocLoopFactor();
    }

    virtual SolverMonitorHandlerType& monitor()
    {
      return solverMonitorHandler_;
    }

    virtual IOTupleType dataTuple ()
    {
      return dataWriterHandler_.dataTuple();
    }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const
    {
      //only dataPrefix from first tuple element
      return std::get<0>( tuple_ )->problem().dataPrefix();
    }

    // before first step, do data initialization
    virtual void initialize ( int loop, TimeProviderType &tp )
    {
      ForLoopType< Initialize >::apply( tuple_, adaptHandler_, loop, &tp );
    }

    virtual void preStep ( int loop, TimeProviderType &tp )
    {
      ForLoopType< PreSolve >::apply( tuple_, loop, &tp );
    }

    //Needs to be overridden to enable fancy steps
    virtual void step ( int loop, TimeProviderType &tp )
    {
      ForLoopType< Solve >::apply( tuple_, loop, &tp );
    }

    virtual void postStep ( int loop, TimeProviderType &tp )
    {
      ForLoopType< PostSolve >::apply( tuple_, loop, &tp );
    }


    void finalize ( int loop, TimeProviderType &tp )
    {
      // flush diagnostics data
      diagnosticsHandler_.finalize();

      // write last time step
      dataWriterHandler_.finalize( tp );

      // adjust average time step size
      solverMonitorHandler_.finalize( gridWidth(), gridSize() );

      adaptHandler_.finalize();

      ForLoopType< Finalize >::apply( tuple_, loop, &tp );
    }

    StepperTupleType &stepperTuple () { return tuple_; }
    const StepperTupleType &stepperTuple () const { return tuple_; }

  protected:
    StepperTupleType               tuple_;
    StepperParametersType          param_;

    CheckPointHandlerType          checkPointHandler_;
    DataWriterHandlerType          dataWriterHandler_;
    DiagnosticsHandlerType         diagnosticsHandler_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    SolutionLimiterHandlerType     solutionLimiterHandler_;
    AdaptHandlerType               adaptHandler_;

    Dune::Timer                    overallTimer_;
    unsigned int                   timeStepTimer_;
    double                         fixedTimeStep_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
