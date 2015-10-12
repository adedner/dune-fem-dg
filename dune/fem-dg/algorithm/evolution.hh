#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH

#include <iostream>
#include <string>
#include <type_traits>

#include <dune/fem-dg/pass/threadpass.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem-dg/misc/typedefcheck.hh>

#include <dune/fem/misc/femtimer.hh>
#include <dune/common/timer.hh>
#include <dune/fem/space/common/interpolate.hh>

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
  private:
    CHECK_TYPEDEF_EXISTS( AdaptIndicatorType )
    CHECK_TYPEDEF_EXISTS( AdditionalOutputHandlerType )
    CHECK_TYPEDEF_EXISTS( SolverMonitorHandlerType )
    CHECK_TYPEDEF_EXISTS( DiagnosticsHandlerType )

  public:
    static const int polynomialOrder = polOrder;

    typedef ProblemTraits                                          ProblemTraitsType;

    // type of Grid
    typedef Grid                                                   GridType;
    typedef typename ProblemTraits::HostGridPartType               HostGridPartType;
    typedef typename ProblemTraits::GridPartType                   GridPartType;

    typedef typename ProblemTraits::AnalyticalTraits               AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polynomialOrder >  DiscreteTraits;

    // obtain the problem dependent types, analytical context
    typedef typename AnalyticalTraits::ModelType                   ModelType;
    typedef typename AnalyticalTraits::ProblemType                 ProblemType;
    typedef typename AnalyticalTraits::InitialDataType             InitialDataType;

    // type of discrete function space and discrete function
    typedef typename DiscreteTraits::InitialProjectorType          InitialProjectorType;

    // type of dg operator
    typedef typename DiscreteTraits::Operator                      OperatorType;

    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;

    typedef typename DiscreteTraits::ExtraParameterTuple           ExtraParameterTupleType;
    typedef typename DiscreteTraits::IOTupleType                   IOTupleType;

    // wrap operator
    typedef GridTimeProvider< GridType >                           TimeProviderType;

    typedef typename DiscreteTraits::Solver                        SolverType;

    typedef typename AdaptIndicatorTypes< DiscreteTraits >::type            AdaptIndicatorType;
    typedef typename AdditionalOutputHandlerTypes< DiscreteTraits >::type   AdditionalOutputHandlerType;
    typedef typename SolverMonitorHandlerTypes< DiscreteTraits >::type      SolverMonitorHandlerType;
    typedef typename DiagnosticsHandlerTypes< DiscreteTraits >::type        DiagnosticsHandlerType;
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

  template< class SubEvolutionAlgorithmTraits >
  class SubEvolutionAlgorithmBase
  {
    typedef SubEvolutionAlgorithmTraits                          Traits;

    //static_assert( !std::is_void< typename Traits::AdditionalOutputHandlerType >::value,
    //               "SubEvolutionAlgorithmBase wants to create an object of type AdditionalOuputHandlerType."
    //               "Please define 'AdditionalOutputHandlerType' in your SubProblemCreator creating this class." );
  public:
    typedef typename Traits::GridType                            GridType;
    typedef typename Traits::IOTupleType                         IOTupleType;
    typedef typename Traits::TimeProviderType                    TimeProviderType;

    typedef typename Traits::HostGridPartType                    HostGridPartType;

    // An analytical version of our model
    typedef typename Traits::ModelType                           ModelType;
    typedef typename Traits::ProblemType                         ProblemType;

    typedef typename Traits::GridPartType                        GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType           DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType                DiscreteFunctionType;

    typedef typename Traits::SolverType                          SolverType;

    // The DG space operator
    typedef typename Traits::OperatorType                        OperatorType;

    // type of initial interpolation
    typedef typename Traits::InitialProjectorType                InitialProjectorType;

    // analytical Tratis
    typedef typename Traits::AnalyticalTraits                    AnalyticalTraits;

    // discrete Traits
    typedef typename Traits::DiscreteTraits                      DiscreteTraits;

    typedef typename Traits::ExtraParameterTupleType  ExtraParameterTupleType;

    typedef uint64_t                                              UInt64Type ;

    typedef DiscreteFunctionType                                  CheckPointDiscreteFunctionType;
    typedef DiscreteFunctionType                                  LimitDiscreteFunctionType;
    typedef DiscreteFunctionType                                  AdaptationDiscreteFunctionType;

    typedef typename Traits::AdaptIndicatorType                   AdaptIndicatorType;
    typedef typename Traits::DiagnosticsHandlerType               DiagnosticsHandlerType;
    typedef typename Traits::SolverMonitorHandlerType             SolverMonitorHandlerType;
    typedef typename Traits::AdditionalOutputHandlerType          AdditionalOutputHandlerType;


    SubEvolutionAlgorithmBase ( GridType &grid, const std::string name = "" )
    : grid_( grid ),
      algorithmName_( name ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      solution_( "U_"+name, space() ),
      exactSolution_( "U_exact"+name, space() ),
      problem_( Traits::ProblemTraitsType::problem() ),
      model_( *problem_ ),
      diagnosticsHandler_( name ),
      solverMonitorHandler_(  name ),
      additionalOutputHandler_( nullptr ),
      odeSolverMonitor_(),
      overallTimer_(),
      solver_(),
      overallTime_( 0 )
    {}
    virtual const std::string name () { return algorithmName_; }

    GridType& grid () const { return grid_; }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    // return size of grid
    virtual UInt64Type gridSize () const { UInt64Type grSize = grid().size(0); return grid().comm().sum( grSize); }

    virtual bool checkSolutionValid ( const int loop, TimeProviderType& tp ) const { return solution_.dofsValid(); }

    // function creating the ode solvers
    virtual typename SolverType::type* createSolver ( TimeProviderType& ) = 0;

    // return reference to the discrete function space
    const DiscreteFunctionSpaceType& space () const { return space_; }

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
    virtual IOTupleType dataTuple () { return std::make_tuple( &solution(), &exactSolution_ ); }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const {  return problem().dataPrefix(); }

    virtual void initialize ( int loop, TimeProviderType &tp )
    {
      // project initial data
      const bool doCommunicate = ! NonBlockingCommParameter :: nonBlockingCommunication ();
      InitialProjectorType projection( 2 * solution().space().order(), doCommunicate );
      projection( problem().fixedTimeFunction( tp.time() ), solution() );
      //TODO check whether this version would also work (goal: avoid InitialProjectorType typedef :D
      //auto ftp = problem().fixedTimeFunction( tp.time() );
      //GridFunctionAdapter< typename ProblemType::InstationaryFunctionType, GridPartType >
      //  adapter( "-exact", ftp, gridPart() );
      //interpolate( adapter, solution() );
      //if( NonBlockingCommParameter::nonBlockingCommunication() )
      //  solution().communicate();

      // setup ode solver
      solver_.reset( this->createSolver( tp ) );

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

    virtual void preSolve( int loop, TimeProviderType& tp )
    {}

    virtual void solve ( int loop, TimeProviderType &tp )
    {
      overallTimer_.reset();
      odeSolverMonitor_.reset();

      // solve ODE
      solver().solve( solution(), odeSolverMonitor_ );

      overallTime_ = overallTimer_.stop();
    }

    virtual void postSolve( int loop, TimeProviderType& tp )
    {}

    void finalize ( int loop, TimeProviderType &tp )
    {
      // add eoc errors
      AnalyticalTraits::addEOCErrors( tp, solution(), model(), problem() );

      // delete ode solver
      solver_.reset();
    }

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
    typename SolverType::type &solver ()
    {
      assert( solver_ );
      return *solver_;
    }

    GridType&                      grid_;
    std::string                    algorithmName_;
    GridPartType                   gridPart_;
    DiscreteFunctionSpaceType      space_;

    // the solution
    DiscreteFunctionType           solution_;
    DiscreteFunctionType           exactSolution_;

    // InitialDataType evaluates to $u_0$
    std::unique_ptr< ProblemType > problem_;
    ModelType                      model_;

    DiagnosticsHandlerOptional< DiagnosticsHandlerType >           diagnosticsHandler_;
    SolverMonitorHandlerOptional< SolverMonitorHandlerType >       solverMonitorHandler_;
    AdditionalOutputHandlerOptional< AdditionalOutputHandlerType > additionalOutputHandler_;
    typename SolverType::type::MonitorType odeSolverMonitor_;

    Dune::Timer                    overallTimer_;
    std::unique_ptr< typename SolverType::type > solver_;
    double overallTime_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
