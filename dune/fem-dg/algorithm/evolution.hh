#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH

#include <iostream>
#include <string>

#include <dune/fem-dg/pass/threadpass.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/misc/femtimer.hh>
#include <dune/common/timer.hh>
#include <dune/fem/space/common/interpolate.hh>

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
    typedef typename DiscreteTraits::OperatorType                  OperatorType;

    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;

    typedef typename DiscreteTraits::ExtraParameterTuple           ExtraParameterTupleType;
    typedef typename DiscreteTraits::IOTupleType                   IOTupleType;

    // wrap operator
    typedef GridTimeProvider< GridType >                           TimeProviderType;

    typedef typename DiscreteTraits::OdeSolverType                 OdeSolverType;

    typedef typename DiscreteTraits::BasicLinearSolverType         BasicLinearSolverType;

    typedef typename DiscreteTraits::AdaptIndicatorType            AdaptIndicatorType;
    typedef typename DiscreteTraits::SolverMonitorHandlerType      SolverMonitorHandlerType;
    typedef typename DiscreteTraits::DiagnosticsHandlerType        DiagnosticsHandlerType;
    typedef typename DiscreteTraits::AdditionalOutputHandlerType   AdditionalOutputHandlerType;
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
      solverMonitorHandler_( name ),
      additionalOutputHandler_( exactSolution_, name ),
      odeSolverMonitor_(),
      overallTimer_(),
      odeSolver_(),
      overallTime_( 0 )
    {}
    virtual const std::string name () { return algorithmName_; }

    GridType& grid () const { return grid_; }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    // return size of grid
    virtual UInt64Type gridSize () const { UInt64Type grSize = grid().size(0); return grid().comm().sum( grSize); }

    virtual bool checkDofsValid ( TimeProviderType& tp, const int loop  ) const { return solution_.dofsValid(); }

    // function creating the ode solvers
    virtual OdeSolverType* createOdeSolver ( TimeProviderType& ) = 0;

    // return reference to the discrete function space
    const DiscreteFunctionSpaceType& space () const { return space_; }

    // return reference to discrete function holding solution
    DiscreteFunctionType& solution () { return solution_; }

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType& monitor() { return solverMonitorHandler_; }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType& diagnostics() { return diagnosticsHandler_; }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType& additionalOutput() { return additionalOutputHandler_; }

    //LIMITING
    virtual void limit(){}
    virtual LimitDiscreteFunctionType* limitSolution () { return &solution_; }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() { return (AdaptIndicatorType*)0; }
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
      odeSolver_.reset( this->createOdeSolver( tp ) );

      // initialize ode solver
      odeSolver().initialize( solution() );

      //initialize solverMonitor
      solverMonitorHandler_.registerData( "GridWidth", solverMonitorHandler_.monitor().gridWidth, nullptr, true );
      solverMonitorHandler_.registerData( "Elements", solverMonitorHandler_.monitor().elements, nullptr, true );
      solverMonitorHandler_.registerData( "TimeSteps", solverMonitorHandler_.monitor().timeSteps, nullptr, true );
      solverMonitorHandler_.registerData( "AvgTimeStep", solverMonitorHandler_.monitor().avgTimeStep );
      solverMonitorHandler_.registerData( "MinTimeStep", solverMonitorHandler_.monitor().minTimeStep );
      solverMonitorHandler_.registerData( "MaxTimeStep", solverMonitorHandler_.monitor().maxTimeStep );
      solverMonitorHandler_.registerData( "Newton", solverMonitorHandler_.monitor().newton_iterations,       &odeSolverMonitor_.newtonIterations_);
      solverMonitorHandler_.registerData( "ILS", solverMonitorHandler_.monitor().ils_iterations,             &odeSolverMonitor_.linearSolverIterations_);
      //solverMonitorHandler_.registerData( "OC", solverMonitorHandler_.monitor().operator_calls,              &odeSolverMonitor_.spaceOperatorCalls_);
      solverMonitorHandler_.registerData( "MaxNewton",solverMonitorHandler_.monitor().max_newton_iterations, &odeSolverMonitor_.maxNewtonIterations_ );
      solverMonitorHandler_.registerData( "MaxILS",solverMonitorHandler_.monitor().max_ils_iterations,    &odeSolverMonitor_.maxLinearSolverIterations_ );

      //initialize diagnosticsHandler
      diagnosticsHandler_.registerData( "OperatorTime", &odeSolverMonitor_.operatorTime_ );
      diagnosticsHandler_.registerData( "OdeSolveTime", &odeSolverMonitor_.odeSolveTime_ );
      diagnosticsHandler_.registerData( "OverallTimer", &overallTime_ );
      diagnosticsHandler_.registerData( "NumberOfElements", &odeSolverMonitor_.numberOfElements_ );
    }

    virtual void preSolve( int loop, TimeProviderType& tp )
    {}

    virtual void solve ( int loop, TimeProviderType &tp )
    {
      overallTimer_.reset();
      odeSolverMonitor_.reset();

      // solve ODE
      odeSolver().solve( solution(), odeSolverMonitor_ );

      overallTime_ = overallTimer_.stop();
    }

    virtual void postSolve( int loop, TimeProviderType& tp )
    {}

    void finalize ( int loop, TimeProviderType &tp )
    {
      // add eoc errors
      AnalyticalTraits::addEOCErrors( tp, solution(), model(), problem() );

      // delete ode solver
      odeSolver_.reset();
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

    DiagnosticsHandlerType         diagnosticsHandler_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    AdditionalOutputHandlerType    additionalOutputHandler_;
    typename OdeSolverType::MonitorType odeSolverMonitor_;

    Dune::Timer                    overallTimer_;
    std::unique_ptr< OdeSolverType > odeSolver_;
    double overallTime_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
