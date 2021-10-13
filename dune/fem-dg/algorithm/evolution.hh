#ifndef DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_COMBINED_EVOLUTION_HH


#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/createsubalgorithms.hh>
#include <dune/fem/misc/femtimer.hh>

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem-dg/pass/threadpass.hh>
#include <dune/common/timer.hh>

#include <dune/fem-dg/algorithm/caller/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/checkpoint.hh>
#include <dune/fem-dg/algorithm/caller/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/eocwriter.hh>
#include <dune/fem-dg/algorithm/caller/postprocessing.hh>
#include <dune/fem-dg/algorithm/caller/adapt.hh>
#include <dune/fem-dg/misc/timesteppingparameter.hh>

#include "containers.hh"

namespace Dune
{
namespace Fem
{

  // internal forward declarations
  // -----------------------------

  template< class Traits >
  class EvolutionAlgorithmBase;

  /**
   *  \brief Traits class
   */
  template< int polOrder, class ... ProblemTraits >
  struct EvolutionAlgorithmTraits
  {
    // type of Grid
    typedef typename std::tuple_element_t<0, std::tuple< ProblemTraits... > >::GridType
                                                                        GridType;

    typedef std::tuple< typename std::add_lvalue_reference<typename ProblemTraits::GridType>::type ... >
                                                                        GridTypes;


    // wrap operator
    typedef GridTimeProvider< GridType >                                TimeProviderType;

    typedef CreateSubAlgorithms< typename ProblemTraits::template Algorithm<polOrder>...  >
                                                                        CreateSubAlgorithmsType;

    typedef typename CreateSubAlgorithmsType::SubAlgorithmTupleType     SubAlgorithmTupleType;

    //typedef typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >
    //                                                                  IndexSequenceType;

    typedef Dune::Fem::AdaptCaller< SubAlgorithmTupleType >             AdaptCallerType;
    typedef Dune::Fem::CheckPointCaller< SubAlgorithmTupleType >        CheckPointCallerType;
    typedef Dune::Fem::SolverMonitorCaller< SubAlgorithmTupleType >     SolverMonitorCallerType;
    typedef Dune::Fem::DataWriterCaller< SubAlgorithmTupleType >        DataWriterCallerType;
    typedef Dune::Fem::EocWriterCaller< SubAlgorithmTupleType >         EocWriterCallerType;
    typedef Dune::Fem::DiagnosticsCaller< SubAlgorithmTupleType >       DiagnosticsCallerType;
    typedef Dune::Fem::PostProcessingCaller< SubAlgorithmTupleType >    PostProcessingCallerType;
  };


  /**
   *  \brief A global algorithm class
   *
   * \ingroup Algorithms
   */
  template< int polOrder, class ... ProblemTraits >
  class EvolutionAlgorithm
  : public EvolutionAlgorithmBase< EvolutionAlgorithmTraits< polOrder, ProblemTraits ... > >
  {
    typedef EvolutionAlgorithmTraits< polOrder, ProblemTraits... > Traits;
    typedef EvolutionAlgorithmBase< Traits > BaseType;
  public:
    typedef typename BaseType::GridType GridType;

    template< class GlobalContainerImp >
    EvolutionAlgorithm ( const std::string name, const std::shared_ptr<GlobalContainerImp>& cont )
    : BaseType( name, cont )
    {}

  };

  /**
   *  \brief A global algorithm class
   *
   * \ingroup Algorithms
   */
  template< class Traits >
  class EvolutionAlgorithmBase
    : public EOCAlgorithm< Traits >
  {
    typedef EOCAlgorithm< Traits >                               BaseType;
  public:
    typedef typename BaseType::GridType                          GridType;
    typedef typename BaseType::SolverMonitorCallerType           SolverMonitorCallerType;

    typedef typename Traits::SubAlgorithmTupleType               SubAlgorithmTupleType;
    typedef typename Traits::TimeProviderType                    TimeProviderType;

    typedef typename Traits::DiagnosticsCallerType               DiagnosticsCallerType;
    typedef typename Traits::CheckPointCallerType                CheckPointCallerType;
    typedef typename Traits::DataWriterCallerType                DataWriterCallerType;
    typedef typename Traits::EocWriterCallerType                 EocWriterCallerType;
    typedef typename Traits::PostProcessingCallerType            PostProcessingCallerType;
    typedef typename Traits::AdaptCallerType                     AdaptCallerType;

    typedef typename Traits::CreateSubAlgorithmsType             CreateSubAlgorithmsType;

    typedef typename BaseType::UInt64Type                        UInt64Type;

    typedef TimeSteppingParameters                               TimeSteppingParametersType;

    using BaseType::eocParams;
    using BaseType::grid;
    using BaseType::tuple_;
    using BaseType::dataWriterCaller_;
    using BaseType::solverMonitorCaller_;

  private:
    struct Initialize {
    private:
      template<class T, class AdaptCaller, class... Args >
      static typename std::enable_if< std::is_void< typename T::element_type::DiagnosticsType >::value >::type
      getDiagnostics( T&, AdaptCaller&, Args&& ... ){}
      template<class T, class AdaptCaller, class... Args >
      static typename std::enable_if< !std::is_void< typename T::element_type::DiagnosticsType >::value >::type
      getDiagnostics( T& e, AdaptCaller& caller, Args &&... a )
      {
        if( e->diagnostics() )
        {
          e->diagnostics()->registerData( "AdaptationTime", &caller.adaptationTime() );
          e->diagnostics()->registerData( "LoadBalanceTime", &caller.loadBalanceTime() );
        }
      }
    public:
      template< class T, class AdaptCaller, class ... Args > static void apply ( T& e, AdaptCaller& caller, Args && ... a )
      {
        e->initialize( std::forward<Args>(a)... );
        getDiagnostics( e, caller, std::forward<Args>(a)... );
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
      template<class T, class... Args > static void apply( T& e, bool& res, Args&& ... a )
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

    template< class GlobalContainerImp >
    EvolutionAlgorithmBase ( const std::string name, const std::shared_ptr<GlobalContainerImp>& cont )
    : BaseType( name, cont ),
      checkPointCaller_( tuple_ ),
      diagnosticsCaller_( tuple_ ),
      postProcessingCaller_( tuple_ ),
      adaptCaller_( tuple_ ),
      param_( TimeSteppingParametersType( ParameterKey::generate( "", "femdg.stepper." ) ) ),
      overallTimer_(),
      timeStepTimer_( Dune::FemTimer::addTo("sum time for timestep") ),
      fixedTimeStep_( param_.fixedTimeStep() )
    {}

    /**
     * \brief Returns grid width of grid (overload in derived classes)
     */
    double gridWidth () const
    {
      double res=0.0;
      ForLoopType< GridWidth >::apply( tuple_, res );
      return res;
    }

    /**
     * \brief Returns size of grid.
     */
    UInt64Type gridSize () const
    {
      UInt64Type res=0;
      ForLoopType< GridSize >::apply( tuple_, res );
      return res;
    }

    /**
     * \brief checks whether the computed solution is physically valid, i.e. contains NaNs.
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
    void eocSolve ( const int loop )
    {
      // get start and end time from parameter file
      const double startTime = param_.startTime();
      const double endTime   = param_.endTime();

      // Initialize TimeProvider
      TimeProviderType tp( startTime, this->grid() );

      // call solve implementation taking start and end time
      eocSolve( loop, tp, endTime );
    }

    /**
     *  \brief Solves the whole problem
     *
     *  \param loop the number of the eoc loop
     *  \param tp time provider
     *  \param endTime end Time of the simulation
     */
    void eocSolve ( const int loop, TimeProviderType& tp, const double endTime )
    {
      double maxTimeStep = param_.maxTimeStep();

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      // in codegen modus make endTime large and only compute one timestep
      const int maximalTimeSteps = 1;
#else
      // if this variable is set then only maximalTimeSteps timesteps will be computed
      const int maximalTimeSteps = param_.maximalTimeSteps();
#endif

      // CALLER
      checkPointCaller_.initializeStart( this, loop, tp );

      initialize( loop, tp );

      // CALLER
      if( !checkPointCaller_.checkPointRestored() )
        adaptCaller_.initializeEnd( this, loop, tp );
      dataWriterCaller_.initializeEnd( this, loop, tp );
      checkPointCaller_.initializeEnd( this, loop, tp );

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

      dataWriterCaller_.postSolveEnd( this, loop, tp );
      //******************************
      //*  Time Loop                 *
      //******************************
      for( ; tp.time() < endTime; )
      {
        // CALLER
        dataWriterCaller_.preSolveStart( this, loop, tp );
        checkPointCaller_.preSolveStart( this, loop, tp );

        // reset time step estimate
        tp.provideTimeStepEstimate( maxTimeStep );

        //************************************************
        //* Compute an ODE timestep                      *
        //************************************************
        Dune::FemTimer::start( timeStepTimer_ );
        overallTimer_.reset();

        preSolve( loop, tp );

        // CALLER
        adaptCaller_.solveStart( this, loop, tp );

        // perform the solve for one time step, i.e. solve ODE
        solve( loop, tp );

        //go to next time step, if time step was invalidated by solver
        if( tp.timeStepValid() )
        {
          // CALLER
          postProcessingCaller_.solveEnd( this, loop, tp );

          postSolve( loop, tp );

          // CALLER
          solverMonitorCaller_.postSolveEnd( this, loop, tp );
          diagnosticsCaller_.postSolveEnd( this, loop, tp );
          dataWriterCaller_.postSolveEnd( this, loop, tp );
        }

        // stop FemTimer for this time step
        Dune::FemTimer::stop(timeStepTimer_,Dune::FemTimer::sum);

        // next advance should not exceed endtime
        if( stopAtEndTime )
          tp.provideTimeStepEstimate( (endTime - tp.time()) );

        // next time step is prescribed by fixedTimeStep
        if ( fixedTimeStep > 1e-20 )
          tp.next( fixedTimeStep );
        else
          tp.next();

        const int timeStep = tp.timeStep();
        if( tp.timeStepValid() )
          printTimeStepInformation( timeStep, tp );

        // for debugging and codegen only
        if( timeStep >= maximalTimeSteps )
        {
          if( Fem::Parameter::verbose() )
            std::cerr << "ABORT: time step count reached max limit of " << maximalTimeSteps << std::endl;
          break;
        }

        if (tp.timeStep()<2)
        {
          // write parameters used (before simulation starts)
          Fem::Parameter::write("parameter.log");
        }

      } /****** END of time loop *****/

      // CALLER
      diagnosticsCaller_.finalizeStart( this, loop, tp );
      dataWriterCaller_.finalizeStart( this, loop, tp );
      solverMonitorCaller_.finalizeStart( this, loop, tp );
      adaptCaller_.finalizeStart( this, loop, tp );
      checkPointCaller_.finalizeStart( this, loop, tp );

      finalize( loop, tp );

      // prepare the fixed time step for the next eoc loop
      fixedTimeStep_ /= param_.fixedTimeStepEocLoopFactor();
    }

    /**
     * \brief Returns type of post processing caller
     */
    virtual PostProcessingCallerType& postProcessing()
    {
      return postProcessingCaller_;
    }

    /**
     * \brief Initializes all sub-algorithms.
     *
     * This method calls consecutively `initialize()` for all sub-algorithms.
     * \code
     * alg1_.initialize( loop, tp );
     * alg2_.initialize( loop, tp );
     * alg3_.initialize( loop, tp );
     * // and so on.
     * \endcode
     *
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    virtual void initialize ( int loop, TimeProviderType &tp )
    {
      ForLoopType< Initialize >::apply( tuple_, adaptCaller_, loop, &tp );
    }

    /**
     * \brief Prepare solve of all sub-algorithms.
     *
     * This method calls consecutively `preSolve()` for all sub-algorithms.
     * \code
     * alg1_.preSolve( loop, tp );
     * alg2_.preSolve( loop, tp );
     * alg3_.preSolve( loop, tp );
     * // and so on.
     * \endcode
     *
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    virtual void preSolve ( int loop, TimeProviderType &tp )
    {
      ForLoopType< PreSolve >::apply( tuple_, loop, &tp );
      //forEach( tuple_, []( auto& e, auto&&... a ){ e->preSolve( std::forward<decltype(a)>(a)... ) } );
    }

    /**
     * \brief Solves all sub-algorithms.
     *
     * This method calls consecutively `solve()` for all sub-algorithms.
     * \code
     * alg1_.solve( loop, tp );
     * alg2_.solve( loop, tp );
     * alg3_.solve( loop, tp );
     * // and so on.
     * \endcode
     *
     * \note For coupled algorithms, this method has to be overloaded.
     *
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    virtual void solve ( int loop, TimeProviderType &tp )
    {
      ForLoopType< Solve >::apply( tuple_, loop, &tp );
    }

    /**
     * \brief Finish solves of all sub-algorithms.
     *
     * This method calls consecutively `postSolve()` for all sub-algorithms.
     * \code
     * alg1_.postSolve( loop, tp );
     * alg2_.postsolve( loop, tp );
     * alg3_.postsolve( loop, tp );
     * // and so on.
     * \endcode
     *
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    virtual void postSolve ( int loop, TimeProviderType &tp )
    {
      ForLoopType< PostSolve >::apply( tuple_, loop, &tp );
    }

    /**
     * \brief Finish solves of all sub-algorithms.
     *
     * This method calls consecutively `finalize()` for all sub-algorithms.
     * \code
     * alg1_.finalize( loop, tp );
     * alg2_.finalize( loop, tp );
     * alg3_.finalize( loop, tp );
     * // and so on.
     * \endcode
     *
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    void finalize ( int loop, TimeProviderType &tp )
    {
      ForLoopType< Finalize >::apply( tuple_, loop, &tp );
    }

  private:

    void printTimeStepInformation( int timeStep, TimeProviderType& tp )
    {
      const int printCount = param_.printCount();
      if( (printCount > 0) && (timeStep % printCount) == 0)
      {
        // obtain grid size (requires all reduce)
        const UInt64Type grdsize = gridSize();
        if( grid().comm().rank() == 0 )
        {
          double dt = tp.deltaT();
          std::cout << "step: " << timeStep << "  time = " << tp.time()+dt << ", dt = " << dt
                    <<",  grid size: " << grdsize << ", elapsed time: ";
          Dune::FemTimer::print(std::cout,timeStepTimer_);
          std::cout << "Newton:  " << solverMonitorCaller_.getData( "Newton" ) << ", ";
          std::cout << "ILS:  " << solverMonitorCaller_.getData( "ILS" ) << ", ";
          std::cout << "OC:  " << solverMonitorCaller_.getData( "OC" ) << ", ";
          std::cout << std::endl;
        }
      }
    }

  protected:
    CheckPointCallerType           checkPointCaller_;
    DiagnosticsCallerType          diagnosticsCaller_;
    PostProcessingCallerType       postProcessingCaller_;
    AdaptCallerType                adaptCaller_;

    TimeSteppingParametersType     param_;
    Dune::Timer                    overallTimer_;
    unsigned int                   timeStepTimer_;
    double                         fixedTimeStep_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
