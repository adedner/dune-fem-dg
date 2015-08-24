#ifndef DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH


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
#include <dune/fem-dg/algorithm/monitor.hh>


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
            int polOrder,
            class HandlerTraits,
            class ProblemTraits,
            class... Tail  >
  struct CombinedEvolutionAlgorithmTraits
  {
    static const int polynomialOrder = polOrder;

    typedef ProblemTraits                                         ProblemTraitsType;

    // type of Grid
    typedef Grid                                          GridType;
    typedef typename ProblemTraits :: HostGridPartType    HostGridPartType;
    typedef typename ProblemTraits :: GridPartType        GridPartType;

    typedef typename ProblemTraits::AnalyticalTraits               AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polynomialOrder >  DiscreteTraits;

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

    typedef typename HandlerTraits::DiagnosticsHandlerType         DiagnosticsHandlerType;
    typedef typename HandlerTraits::SolverMonitorHandlerType       SolverMonitorHandlerType;
    typedef typename HandlerTraits::CheckPointHandlerType          CheckPointHandlerType;
    typedef typename HandlerTraits::DataWriterHandlerType          DataWriterHandlerType;
    typedef typename HandlerTraits::AdditionalOutputHandlerType    AdditionalOutputHandlerType;
    typedef typename HandlerTraits::SolutionLimiterHandlerType     SolutionLimiterHandlerType;
    typedef typename HandlerTraits::AdaptHandlerType               AdaptHandlerType;
  };


  // EvolutionAlgorithm
  // ------------------

  template< int polOrder, class HandlerTraits, class ProblemHead, class... ProblemTail >
  class CombinedEvolutionAlgorithm
    : public AlgorithmBase< CombinedEvolutionAlgorithmTraits< typename ProblemHead::GridType,
                                                              polOrder,
                                                              HandlerTraits,
                                                              ProblemHead,
                                                              ProblemTail... > >
  {
    typedef CombinedEvolutionAlgorithmTraits< typename ProblemHead::GridType,
                                                              polOrder,
                                                              HandlerTraits,
                                                              ProblemHead,
                                                              ProblemTail... >                 Traits;
    typedef AlgorithmBase< Traits >                                                     BaseType;
  public:
    typedef typename BaseType::GridType                          GridType;
    typedef typename BaseType::IOTupleType                       IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType          SolverMonitorHandlerType;

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
    typedef typename Traits::CheckPointHandlerType                      CheckPointHandlerType;
    typedef typename Traits::DataWriterHandlerType                      DataWriterHandlerType;
    typedef typename Traits::AdditionalOutputHandlerType                AdditionalOutputHandlerType;
    typedef typename Traits::SolutionLimiterHandlerType                 SolutionLimiterHandlerType;
    typedef typename Traits::AdaptHandlerType                           AdaptHandlerType;

    typedef typename Traits::ExtraParameterTupleType  ExtraParameterTupleType;

    typedef uint64_t                                              UInt64Type ;
    typedef StepperParameters                                     StepperParametersType;

    typedef std::tuple< typename std::add_pointer< typename ProblemHead::template Stepper<polOrder>::Type >::type,
                        typename std::add_pointer< typename ProblemTail::template Stepper<polOrder>::Type >::type... >     StepperTupleType;

    using BaseType::grid_;
    using BaseType::eocParam_;
    using BaseType::grid;

    template< int i >
    struct Constructor
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        typedef typename std::remove_pointer< typename std::tuple_element< i, Tuple >::type >::type Element;
        std::get< i >( tuple ) = new Element( std::forward< Args >( args ) ... );
      }
    };
    template< int i >
    struct InitializeStep
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->initializeStep( args... );
      }
    };
    template< int i >
    struct Step
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->step( args... );
      }
    };
    template< int i >
    struct FinalizeStep
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->finalizeStep( args... );
      }
    };
    template< int i >
    struct GridWidth
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, double& res, Args && ... args )
      {
        res = std::max( res, std::get< i >( tuple )->gridWidth( args... ) );
      }
    };
    template< int i >
    struct GridSize
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, UInt64Type& res, Args && ... args )
      {
        res = std::max( res, std::get< i >( tuple )->gridSize( args... ) );
      }
    };
    template< int i >
    struct CheckDofsValid
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, bool& res, Args && ... args )
      {
        res &= std::get< i >( tuple )->checkDofsValid( args... );
      }
    };

    CombinedEvolutionAlgorithm ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  ),
      tuple_( createStepper( grid, name ) ),
      param_( StepperParametersType( Dune::ParameterKey::generate( "", "femdg.stepper." ) ) ),
      checkPointHandler_( tuple_ ),
      solverMonitorHandler_( tuple_ ),
      dataWriterHandler_( tuple_ ),
      diagnosticsHandler_( tuple_ ),
      additionalOutputHandler_( tuple_ ),
      solutionLimiterHandler_( tuple_ ),
      adaptHandler_( tuple_ ),
      fixedTimeStep_( param_.fixedTimeStep() )
    {}

    // create Tuple of contained subspaces
    static StepperTupleType createStepper( GridType &grid, const std::string name = "" )
    {
      StepperTupleType tuple;
      ForLoop< Constructor, 0, sizeof ... ( ProblemTail ) >::apply( tuple, grid, name );
      return tuple;
    }

    // return grid width of grid (overload in derived classes)
    double gridWidth () const
    {
      double res=0.0;
      ForLoop< GridWidth, 0, sizeof ... ( ProblemTail ) >::apply( tuple_, res );
      return res;
    }

    // return size of grid
    UInt64Type gridSize () const
    {
      UInt64Type res=0;
      ForLoop< GridSize, 0, sizeof ... ( ProblemTail ) >::apply( tuple_, res );
      return res;
    }

    bool checkDofsValid( TimeProviderType& tp, const int loop  ) const
    {
      bool res = true;
      ForLoop< CheckDofsValid, 0, sizeof ... ( ProblemTail ) >::apply( tuple_, res, tp, loop );
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

      dataWriterHandler_.init( tp, eocParam_.dataOutputParameters( loop, dataPrefix() ) );

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

        // estimate, mark, adapt
        adaptHandler_.step( tp );

        // perform the solve for one time step, i.e. solve ODE
        step( tp );


        // limit solution if necessary
        solutionLimiterHandler_.step();

        // update time step information
        solverMonitorHandler_.step( tp );

        //write diagnostics
        diagnosticsHandler_.step( tp, 123456789 /*odeSolverMonitor.numberOfElements_*/ );


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
          if( grid.comm().rank() == 0 )
          {
            std::cout << "step: " << tp.timeStep()-1 << "  time = " << tp.time()+tp.deltaT() << ", dt = " << deltaT
                      <<",  grid size: " << gridSize() << ", elapsed time: ";
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

    std::string description () const { return "global algorithm"; }

    //! finalize problem, i.e. calculated EOC ...
    virtual void finalize ( const int eocloop )
    {}

    virtual SolverMonitorHandlerType& monitor()
    {
      return solverMonitorHandler_;
    }

    virtual IOTupleType* dataTuple ()
    {
      return &dataWriterHandler_.dataTuple();
    }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const
    {
      //only dataPrefix from first tuple element
      return std::get<0>( tuple_ )->dataPrefix();
    }

    // before first step, do data initialization
    virtual void initializeStep ( TimeProviderType &tp, int loop )
    {
      ForLoop< InitializeStep, 0, sizeof ... ( ProblemTail ) >::apply( tuple_, tp, loop );
    }

    //Needs to be overridden to enable fancy steps
    virtual void step ( TimeProviderType &tp )
    {
      ForLoop< Step, 0, sizeof ... ( ProblemTail ) >::apply( tuple_, tp );
    }

    void finalizeStep ( TimeProviderType &tp )
    {
      // flush diagnostics data
      diagnosticsHandler_.finalize();

      // write last time step
      dataWriterHandler_.finalize( tp );

      // adjust average time step size
      solverMonitorHandler_.finalize( gridWidth(), gridSize() );

      adaptHandler_.finalize();

      ForLoop< FinalizeStep, 0, sizeof ... ( ProblemTail ) >::apply( tuple_, tp );
    }

  protected:
    StepperTupleType               tuple_;
    StepperParametersType          param_;

    CheckPointHandlerType          checkPointHandler_;
    DataWriterHandlerType          dataWriterHandler_;
    DiagnosticsHandlerType         diagnosticsHandler_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    AdditionalOutputHandlerType    additionalOutputHandler_;
    SolutionLimiterHandlerType     solutionLimiterHandler_;
    AdaptHandlerType               adaptHandler_;

    Dune::Timer                    overallTimer_;
    unsigned int                   timeStepTimer_;
    double                         fixedTimeStep_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
