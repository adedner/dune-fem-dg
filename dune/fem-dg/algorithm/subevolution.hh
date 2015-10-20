#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH

#include <iostream>
#include <string>
#include <type_traits>

#include <dune/fem-dg/pass/threadpass.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>

#include <dune/fem/misc/femtimer.hh>
#include <dune/common/timer.hh>
#include <dune/fem/space/common/interpolate.hh>

#include <dune/fem-dg/algorithm/interface.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>

namespace Dune
{

namespace Fem
{



  // SubEvolutionAlgorithmTraits
  // -------------------------

  template< class Grid,
            class ProblemTraits,
            int polOrder >
  struct SubEvolutionAlgorithmTraits
  {
    typedef Grid GridType;
    typedef ProblemTraits ProblemTraitsType;
    enum { polynomialOrder = polOrder };


    typedef typename ProblemTraits::AnalyticalTraits                    AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polOrder > DiscreteTraits;

    typedef typename DiscreteTraits::InitialProjectorType               InitialProjectorType;

    typedef typename DiscreteTraits::Operator                           OperatorType;

    typedef typename DiscreteTraits::ExtraParameterTuple                ExtraParameterTupleType;

    typedef typename DiscreteTraits::Solver                             SolverType;

  };



  template< class SubEvolutionAlgorithmTraits >
  class SubEvolutionAlgorithmBase;

  template< class Grid, class ProblemTraits, int polOrder >
  class SubEvolutionAlgorithm
    : public SubEvolutionAlgorithmBase< SubEvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder > >
  {
    typedef SubEvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder >                    Traits;
    typedef SubEvolutionAlgorithmBase< Traits >                                             BaseType;
  public:
    SubEvolutionAlgorithm( Grid &grid, const std::string name = "" )
    : BaseType( grid, name  )
    {}
  };

  // SubEvolutionAlgorithm
  // ------------------

  template< class SubTraits >
  class SubEvolutionAlgorithmBase
    : public SubAlgorithmInterface< typename SubTraits::GridType,
                                    typename SubTraits::ProblemTraitsType,
                                    SubTraits::polynomialOrder >
  {
    typedef SubTraits                          Traits;

    typedef SubAlgorithmInterface< typename Traits::GridType,
                                   typename Traits::ProblemTraitsType,
                                   Traits::polynomialOrder > BaseType;

  public:
    typedef typename BaseType::GridType                              GridType;
    typedef typename BaseType::GridPartType                          GridPartType;
    typedef typename BaseType::HostGridPartType                      HostGridPartType;

    typedef typename BaseType::TimeProviderType                      TimeProviderType;

    typedef typename BaseType::ModelType                             ModelType;
    typedef typename BaseType::ProblemType                           ProblemType;

    typedef typename BaseType::DiscreteFunctionType                  DiscreteFunctionType;

    typedef typename BaseType::UInt64Type                            UInt64Type;

    typedef typename BaseType::CheckPointDiscreteFunctionType        CheckPointDiscreteFunctionType;
    typedef typename BaseType::LimitDiscreteFunctionType             LimitDiscreteFunctionType;
    typedef typename BaseType::AdaptationDiscreteFunctionType        AdaptationDiscreteFunctionType;

    typedef typename BaseType::IOTupleType                           IOTupleType;
    typedef typename BaseType::AdaptIndicatorType                    AdaptIndicatorType;
    typedef typename BaseType::DiagnosticsHandlerType                DiagnosticsHandlerType;
    typedef typename BaseType::SolverMonitorHandlerType              SolverMonitorHandlerType;
    typedef typename BaseType::AdditionalOutputHandlerType           AdditionalOutputHandlerType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename Traits::InitialProjectorType                    InitialProjectorType;
    typedef typename Traits::ExtraParameterTupleType                 ExtraParameterTupleType;

    // The DG space operator
    typedef typename Traits::OperatorType                            OperatorType;

    // type of steady state solver
    typedef typename Traits::SolverType                              SolverType;

    // type of analytical traits
    typedef typename Traits::AnalyticalTraits                        AnalyticalTraits;

    // type of discrete traits
    typedef typename Traits::DiscreteTraits                          DiscreteTraits;


    using BaseType::grid_;
    using BaseType::algorithmName_;
    using BaseType::problem;
    using BaseType::model;
    using BaseType::gridSize;

    SubEvolutionAlgorithmBase ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      solution_( "U_"+name, space_ ),
      exactSolution_( "U_exact"+name, space_ ),
      diagnosticsHandler_( name ),
      solverMonitorHandler_(  name ),
      additionalOutputHandler_( nullptr ),
      odeSolverMonitor_(),
      overallTimer_(),
      solver_(),
      overallTime_( 0 ),
      ioTuple_( std::make_tuple( &solution(), &exactSolution_ ) )
    {}

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    // function creating the ode solvers
    virtual typename SolverType::type* createSolver ( TimeProviderType& ) = 0;

    // return reference to discrete function holding solution
    DiscreteFunctionType& solution () { return solution_; }

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType* monitor() { return solverMonitorHandler_.value(); }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType* diagnostics() { return diagnosticsHandler_.value(); }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType* additionalOutput() { return additionalOutputHandler_.value(); }

    //LIMITING
    virtual void limit(){}
    virtual LimitDiscreteFunctionType* limitSolution () { return &solution_; }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() { return nullptr; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return &solution_; }

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () { return &solution_; }

    //DATAWRITING
    virtual IOTupleType& dataTuple () { return ioTuple_; }

    virtual bool checkSolutionValid ( const int loop, TimeProviderType* tp ) const { return solution_.dofsValid(); }

    virtual void initialize ( const int loop, TimeProviderType* tp )
    {
      // project initial data
      const bool doCommunicate = ! NonBlockingCommParameter :: nonBlockingCommunication ();
      InitialProjectorType projection( 2 * solution().space().order(), doCommunicate );
      projection( problem().fixedTimeFunction( tp->time() ), solution() );
      //TODO check whether this version would also work (goal: avoid InitialProjectorType typedef :D
      //auto ftp = problem().fixedTimeFunction( tp.time() );
      //GridFunctionAdapter< typename ProblemType::InstationaryFunctionType, GridPartType >
      //  adapter( "-exact", ftp, gridPart() );
      //interpolate( adapter, solution() );
      //if( NonBlockingCommParameter::nonBlockingCommunication() )
      //  solution().communicate();

      // setup ode solver
      solver_.reset( this->createSolver( *tp ) );

      // initialize ode solver
      solver().initialize( solution() );

      //initialize solverMonitor
      if( solverMonitorHandler_ )
      {
        solverMonitorHandler_.registerData( "GridWidth", &solverMonitorHandler_.monitor().gridWidth, nullptr, true );
        solverMonitorHandler_.registerData( "Elements", &solverMonitorHandler_.monitor().elements, nullptr, true );
        solverMonitorHandler_.registerData( "TimeSteps", &solverMonitorHandler_.monitor().timeSteps, nullptr, true );
        solverMonitorHandler_.registerData( "AvgTimeStep", &solverMonitorHandler_.monitor().avgTimeStep );
        solverMonitorHandler_.registerData( "MinTimeStep", &solverMonitorHandler_.monitor().minTimeStep );
        solverMonitorHandler_.registerData( "MaxTimeStep", &solverMonitorHandler_.monitor().maxTimeStep );
        solverMonitorHandler_.registerData( "Newton", &solverMonitorHandler_.monitor().newton_iterations,       &odeSolverMonitor_.newtonIterations_);
        solverMonitorHandler_.registerData( "ILS", &solverMonitorHandler_.monitor().ils_iterations,             &odeSolverMonitor_.linearSolverIterations_);
        //solverMonitorHandler_.registerData( "OC", &solverMonitorHandler_.monitor().operator_calls,              &odeSolverMonitor_.spaceOperatorCalls_);
        solverMonitorHandler_.registerData( "MaxNewton",&solverMonitorHandler_.monitor().max_newton_iterations, &odeSolverMonitor_.maxNewtonIterations_ );
        solverMonitorHandler_.registerData( "MaxILS",&solverMonitorHandler_.monitor().max_ils_iterations,    &odeSolverMonitor_.maxLinearSolverIterations_ );
      }

      //initialize diagnosticsHandler
      if( diagnosticsHandler_ )
      {
        diagnosticsHandler_.registerData( "OperatorTime", &odeSolverMonitor_.operatorTime_ );
        diagnosticsHandler_.registerData( "OdeSolveTime", &odeSolverMonitor_.odeSolveTime_ );
        diagnosticsHandler_.registerData( "OverallTimer", &overallTime_ );
        diagnosticsHandler_.registerData( "NumberOfElements", &odeSolverMonitor_.numberOfElements_ );
      }
    }

    virtual void preSolve( const int loop, TimeProviderType* tp )
    {}

    virtual void solve ( const int loop, TimeProviderType* tp )
    {
      overallTimer_.reset();
      odeSolverMonitor_.reset();

      // solve ODE
      solver().solve( solution(), odeSolverMonitor_ );

      overallTime_ = overallTimer_.stop();
    }

    virtual void postSolve( const int loop, TimeProviderType* tp )
    {}

    virtual void finalize ( const int loop, TimeProviderType* tp )
    {
      // add eoc errors
      AnalyticalTraits::addEOCErrors( *tp, solution(), model(), problem() );

      // delete ode solver
      solver_.reset();
    }

  protected:
    typename SolverType::type &solver ()
    {
      assert( solver_ );
      return *solver_;
    }

    GridPartType                   gridPart_;
    DiscreteFunctionSpaceType      space_;

    // the solution
    DiscreteFunctionType           solution_;
    DiscreteFunctionType           exactSolution_;

    DiagnosticsHandlerOptional< DiagnosticsHandlerType >           diagnosticsHandler_;
    SolverMonitorHandlerOptional< SolverMonitorHandlerType >       solverMonitorHandler_;
    AdditionalOutputHandlerOptional< AdditionalOutputHandlerType > additionalOutputHandler_;
    typename SolverType::type::MonitorType odeSolverMonitor_;

    Dune::Timer                    overallTimer_;
    std::unique_ptr< typename SolverType::type > solver_;
    double overallTime_;
    IOTupleType ioTuple_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
