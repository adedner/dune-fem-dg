#ifndef DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH


#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/coupling.hh>
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
#include <dune/fem-dg/algorithm/handler/postprocessing.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>

namespace Dune
{
namespace Fem
{

  // internal forward declarations
  // -----------------------------

  template< class Traits, template<class, class > class EvolutionCreatorType >
  class EvolutionAlgorithmBase;


  /**
   *  \brief Parameter class for time stepping in instationary equations.
   *
   *  \ingroup ParameterClass
   */
  class TimeSteppingParameters
  : public Dune::Fem::LocalParameter< TimeSteppingParameters, TimeSteppingParameters >
  {
    protected:
    const std::string keyPrefix_;

    public:
    /**
     * \brief Constructor
     *
     * \param keyPrefix the key prefix for the parameter file.
     */
    TimeSteppingParameters( const std::string keyPrefix = "femdg.stepper." )
      : keyPrefix_( keyPrefix )
    {}

    /**
     * \brief returns a fixed time step \f$ \Delta t=\text{const} \f$
     *
     * \note To choose the time step size adaptively (i.e. non fixed time steps),
     * set this value to \f$ \leq0 \f$.
     */
    virtual double fixedTimeStep() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "fixedtimestep" , 0.0 );
    }

    /**
     * \brief return an additional scaling of fixedTimeStep() which is applied
     * for each eoc loop.
     *
     * Example: A value of \f$ 2 \f$ would half the fixed time step size for each eoc loop.
     */
    virtual double fixedTimeStepEocLoopFactor() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "fixedtimestepeocloopfactor" , 1.0 );
    }

    /**
     * \brief returns the end time \f$ t_\text{start} \f$ of the time stepping.
     */
    virtual double startTime() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "starttime" , 0.0 );
    }

    /**
     * \brief returns the end time \f$ t_\text{end} \f$ of the time stepping.
     */
    virtual double endTime() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "endtime"/*, 1.0 */);
    }

    /**
     * \brief returns how often to print additional information about time stepping.
     *
     * \f$ 1 \f$ means every time step, \f$ 2 \f$ every second and so on.
     *
     * \note a value \f$ <1 \f$ disables printing.
     */
    virtual int printCount() const
    {
      return Dune::Fem::Parameter::getValue< int >( keyPrefix_ + "printcount" , -1 );
    }

    /**
     * \brief returns the maximal time step \f$ \Delta_\text{max}} t \f$ for the
     * time stepping.
     */
    virtual double maxTimeStep() const
    {
      return Dune::Fem::Parameter::getValue< double >( keyPrefix_ + "maxtimestep", std::numeric_limits<double>::max());
    }

    /**
     * \brief returns the maximal number of time steps.
     */
    virtual int maximalTimeSteps () const
    {
      return Dune::Fem::Parameter::getValue< int >(  keyPrefix_ + "maximaltimesteps", std::numeric_limits<int>::max());
    }

    /**
     * \brief returns true whether the last time step should reach \f$ t_\text{end} \f$ exactly or not.
     *
     * Usually, it is the case, that the last time step
     *
     * To avoid this inaccuracy (needed for eoc measurements etc.)
     * you can set this option to true.
     *
     * \note The last time steps may be time consuming because
     * the last time steps may become very small.
     */
    virtual bool stopAtEndTime() const
    {
      return Dune::Fem::Parameter::getValue< bool >( keyPrefix_ + "stopatendtime", bool(false) );
    }

  };


  /**
   *  \brief Traits class
   */
  template< int polOrder, class ... ProblemTraits >
  struct EvolutionAlgorithmTraits
  {
    // type of Grid
    typedef typename std::tuple_element<0, std::tuple< ProblemTraits... > >::type::GridType
                                                                        GridType;

    // wrap operator
    typedef GridTimeProvider< GridType >                                TimeProviderType;

    //typedef ...
    typedef std::tuple< typename std::add_pointer< typename ProblemTraits::template Algorithm<polOrder>::Type >::type... >
                                                                        SubAlgorithmTupleType;

    //typedef typename Std::make_index_sequence_impl< std::tuple_size< SubAlgorithmTupleType >::value >::type
    //                                                                  IndexSequenceType;

    typedef Dune::Fem::AdaptHandler< SubAlgorithmTupleType >            AdaptHandlerType;
    typedef Dune::Fem::CheckPointHandler< SubAlgorithmTupleType >       CheckPointHandlerType;
    typedef Dune::Fem::SolverMonitorHandler< SubAlgorithmTupleType >    SolverMonitorHandlerType;
    typedef Dune::Fem::DataWriterHandler< SubAlgorithmTupleType >       DataWriterHandlerType;
    typedef Dune::Fem::DiagnosticsHandler< SubAlgorithmTupleType >      DiagnosticsHandlerType;
    typedef Dune::Fem::PostProcessingHandler< SubAlgorithmTupleType >   PostProcessingHandlerType;

    typedef typename DataWriterHandlerType::IOTupleType                 IOTupleType;

  };


  /**
   *  \brief A global algorithm class
   *
   * \ingroup Algorithms
   */
  template< int polOrder, template<class, class > class EvolutionCreatorType, class ... ProblemTraits >
  class EvolutionAlgorithm
  : public EvolutionAlgorithmBase< EvolutionAlgorithmTraits< polOrder, ProblemTraits ... >, EvolutionCreatorType >
  {
    typedef EvolutionAlgorithmTraits< polOrder, ProblemTraits... > Traits;
    typedef EvolutionAlgorithmBase< Traits, EvolutionCreatorType > BaseType;
  public:
    typedef typename BaseType::GridType GridType;

    EvolutionAlgorithm ( GridType &grid, const std::string name = "" )
      : BaseType( grid, name )
    {}
  };

  /**
   *  \brief A global algorithm class
   *
   * \ingroup Algorithms
   */
  template< class Traits, template<class, class > class EvolutionCreatorType >
  class EvolutionAlgorithmBase
    : public AlgorithmInterface< Traits >
  {
    typedef AlgorithmInterface< Traits >                         BaseType;
  public:
    typedef typename BaseType::GridType                          GridType;
    typedef typename BaseType::IOTupleType                       IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType          SolverMonitorHandlerType;

    typedef typename Traits::SubAlgorithmTupleType               SubAlgorithmTupleType;
    typedef typename Traits::TimeProviderType                    TimeProviderType;

    typedef typename Traits::DiagnosticsHandlerType              DiagnosticsHandlerType;
    typedef typename Traits::CheckPointHandlerType               CheckPointHandlerType;
    typedef typename Traits::DataWriterHandlerType               DataWriterHandlerType;
    typedef typename Traits::PostProcessingHandlerType           PostProcessingHandlerType;
    typedef typename Traits::AdaptHandlerType                    AdaptHandlerType;

    typedef uint64_t                                             UInt64Type ;

    typedef TimeSteppingParameters                               TimeSteppingParametersType;

    using BaseType::eocParams;
    using BaseType::grid;

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
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0, std::tuple_size< SubAlgorithmTupleType >::value-1 >;

  public:

    /**
     * \brief Constructor
     *
     * \param grid
     * \param name the name of the algorithm
     */
    EvolutionAlgorithmBase ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  ),
      tuple_( EvolutionCreatorType< SubAlgorithmTupleType, GridType >::apply( grid ) ),
      checkPointHandler_( tuple_ ),
      dataWriterHandler_( tuple_ ),
      diagnosticsHandler_( tuple_ ),
      solverMonitorHandler_( tuple_ ),
      postProcessingHandler_( tuple_ ),
      adaptHandler_( tuple_ ),
      param_( TimeSteppingParametersType( ParameterKey::generate( "", "femdg.stepper." ) ) ),
      overallTimer_(),
      timeStepTimer_( Dune::FemTimer::addTo("sum time for timestep") ),
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

    /**
     *  \brief checks whether the computed solution is physically valid, i.e. contains NaNs.
     */
    bool checkSolutionValid( const int loop, TimeProviderType& tp ) const
    {
      bool res = true;
      ForLoopType< CheckSolutionValid >::apply( tuple_, res, loop, &tp );
      return res;
    }

    /**
     *  \brief Solves the whole problem
     *
     *  \param loop the number of the eoc loop
     */
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

    /**
     *  \brief Solves the whole problem
     *
     *
     *
     *  \param loop the number of the eoc loop
     *  \param tp time provider
     *  \param endTime end Time of the simulation
     */
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
      bool newStart = ( eocParams().steps() == 1) ? checkPointHandler_.initialize_pre( this, loop, tp ) : false;

      initialize( loop, tp );

      if( newStart )
        adaptHandler_.initialize_post( this, loop, tp );
      dataWriterHandler_.initialize_post( this, loop, tp );
      checkPointHandler_.initialize_post( this, loop, tp );

      //HandlerCaller handler_( this, loop, tp );
      //handler_.initialize_post( adaptHandler_, dataWriterHandler_, checkPointHandler_ )

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
        dataWriterHandler_.preSolve_pre( this, loop, tp );
        checkPointHandler_.preSolve_pre( this, loop, tp );

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
        adaptHandler_.solve_pre( this, loop, tp );

        // perform the solve for one time step, i.e. solve ODE
        step( loop, tp );

        // do post processing, i.e. limit solution if necessary
        postProcessingHandler_.solve_post( this, loop, tp );

        postStep( loop, tp );

        solverMonitorHandler_.postSolve_post( this, loop, tp );
        diagnosticsHandler_.postSolve_post( this, loop, tp );

        // stop FemTimer for this time step
        Dune::FemTimer::stop(timeStepTimer_,Dune::FemTimer::sum);

        // Check that no NAN have been generated
        if( !checkSolutionValid( loop, tp ) )
        {
          dataWriterHandler_.finalize_pre( this, loop, tp );
          std::cerr << "Solution is not valid. Aborting." << std::endl;
          std::abort();
        }

        const int timeStep = tp.timeStep() + 1 ;
        if( (printCount > 0) && (timeStep % printCount) == 0)
        {
          if( grid.comm().rank() == 0 )
          {
            std::cout << "step: " << timeStep << "  time = " << tp.time()+tp.deltaT() << ", dt = " << deltaT
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
        if( timeStep >= maximalTimeSteps )
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
      diagnosticsHandler_.finalize_pre( this, loop, tp );

      // write last time step
      dataWriterHandler_.finalize_pre( this, loop, tp );

      // adjust average time step size
      solverMonitorHandler_.finalize_pre( this, loop, tp );

      adaptHandler_.finalize_pre( this, loop, tp );

      ForLoopType< Finalize >::apply( tuple_, loop, &tp );
    }

    SubAlgorithmTupleType &subAlgorithmTuple () { return tuple_; }
    const SubAlgorithmTupleType &subAlgorithmTuple () const { return tuple_; }

  protected:
    SubAlgorithmTupleType          tuple_;

    CheckPointHandlerType          checkPointHandler_;
    DataWriterHandlerType          dataWriterHandler_;
    DiagnosticsHandlerType         diagnosticsHandler_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    PostProcessingHandlerType      postProcessingHandler_;
    AdaptHandlerType               adaptHandler_;

    TimeSteppingParametersType     param_;
    Dune::Timer                    overallTimer_;
    unsigned int                   timeStepTimer_;
    double                         fixedTimeStep_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
