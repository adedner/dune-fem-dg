#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH

#include <iostream>
#include <string>
#include <type_traits>

#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>

#include <dune/fem/misc/femtimer.hh>
#include <dune/common/timer.hh>
#include <dune/fem/space/common/interpolate.hh>

#include <dune/fem-dg/algorithm/sub/interface.hh>
#include <dune/fem-dg/algorithm/handler/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/sub/adapt.hh>

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

    typedef typename DiscreteTraits::Operator                           OperatorType;

    typedef typename DiscreteTraits::Solver                             SolverType;

  };

  template <class DiscreteFunctionImp >
  struct SubEvolutionContainer
  {
    typedef DiscreteFunctionImp                                          DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridType                 GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType             GridPartType;

    using DiscreteFunction = DiscreteFunctionType;
    using DiscreteFunctionSpace = DiscreteFunctionSpaceType;

  public:

    SubEvolutionContainer( GridType& grid, const std::string name = "" )
    : grid_( grid ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      solution_( new DiscreteFunctionType( name + "u", space() ) ),
      exactSolution_( new DiscreteFunctionType( name + "u-exact", space() ) )
    {}

    //grid
    const GridType& grid() const
    {
      return grid_;
    }
    GridType& grid()
    {
      return grid_;
    }

    //grid part
    const GridPartType& gridPart() const
    {
      return gridPart_;
    }
    GridPartType& gridPart()
    {
      return gridPart_;
    }

    //spaces
    const DiscreteFunctionSpaceType& space() const
    {
      return space_;
    }
    DiscreteFunctionSpaceType& space()
    {
      return space_;
    }

    //solution
    std::shared_ptr< DiscreteFunction > solution() const
    {
      return solution_;
    }
    void setSolution( std::shared_ptr< DiscreteFunction > solution )
    {
      solution_ = solution;
    }

    //exact solution
    std::shared_ptr< DiscreteFunction > exactSolution() const
    {
      return exactSolution_;
    }
    void setExactSolution( std::shared_ptr< DiscreteFunction > exactSolution )
    {
      exactSolution_ = exactSolution;
    }

  private:
    GridType&                 grid_;
    GridPartType              gridPart_;
    DiscreteFunctionSpace     space_;

    std::shared_ptr< DiscreteFunction > solution_;
    std::shared_ptr< DiscreteFunction > exactSolution_;
    std::shared_ptr< DiscreteFunction > rhs_;
  };



  template< class SubEvolutionAlgorithmTraits >
  class SubEvolutionAlgorithmBase;

  template< class Grid, class ProblemTraits, int polOrder >
  class SubEvolutionAlgorithm
    : public SubEvolutionAlgorithmBase< SubEvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder > >
  {
    typedef SubEvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder > Traits;
    typedef SubEvolutionAlgorithmBase< Traits >                          BaseType;
  public:
    typedef typename BaseType::ContainerType                             ContainerType;
    SubEvolutionAlgorithm( Grid &grid, ContainerType& container )
    : BaseType( grid, container )
    {}
  };

  /**
   *  \brief Algorithm for solving an instationary PDE.
   *
   *  \ingroup SubAlgorithms
   */
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
    typedef typename BaseType::DiagnosticsType                DiagnosticsType;
    typedef typename BaseType::SolverMonitorType              SolverMonitorType;
    typedef typename BaseType::AdditionalOutputType           AdditionalOutputType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // The DG space operator
    typedef typename Traits::OperatorType                            OperatorType;

    // type of steady state solver
    typedef typename Traits::SolverType                              SolverType;

    // type of analytical traits
    typedef typename Traits::AnalyticalTraits                        AnalyticalTraits;

    // type of discrete traits
    typedef typename Traits::DiscreteTraits                          DiscreteTraits;

    using BaseType::grid;
    using BaseType::name;
    using BaseType::problem;
    using BaseType::model;
    using BaseType::gridSize;

    typedef SubEvolutionContainer< DiscreteFunctionType >               ContainerType;

    SubEvolutionAlgorithmBase ( GridType &grid, const ContainerType& container )
    : BaseType( grid ),
      container_( container ),
      overallTime_( 0 ),
      overallTimer_(),
      gridPart_( grid ),
      space_( gridPart_ ),
      solution_( container_.solution() ),
      exactSolution_( container_.exactSolution() ),
      solver_( nullptr ),
      ioTuple_( new IOTupleType( std::make_tuple( solution_.get(), exactSolution_.get() ) ) ),
      diagnostics_( name() ),
      solverMonitor_( name() ),
      additionalOutput_( nullptr ),
      odeSolverMonitor_()
    {}

    typename SolverType::type* solver()
    {
      return solver_.get();
    }

    DiscreteFunctionType& solution ()
    {
      assert( solution_ );
      return *solution_;
    }

    DiscreteFunctionType& exactSolution ()
    {
      assert( exactSolution_ );
      return *exactSolution_;
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    //SOLVERMONITOR
    virtual SolverMonitorType* monitor() { return solverMonitor_.value(); }

    //DIAGNOSTICS
    virtual DiagnosticsType* diagnostics() { return diagnostics_.value(); }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputType* additionalOutput() { return additionalOutput_.value(); }

    //LIMITING
    virtual void limit(){}
    virtual LimitDiscreteFunctionType* limitSolution () { return solution_.get(); }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() { return nullptr; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return solution_.get(); }

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () { return solution_.get(); }

    //DATAWRITING
    virtual IOTupleType& dataTuple () { assert( ioTuple_ ); return *ioTuple_; }

  private:
    virtual std::shared_ptr< typename SolverType::type > doCreateSolver( TimeProviderType& tp )
    {
      return nullptr;
    }

    virtual bool doCheckSolutionValid ( const int loop, TimeProviderType& tp ) const override
    {
      assert( solution_ );
      return solution_->dofsValid();
    }

    virtual void doInitialize ( const int loop, TimeProviderType& tp ) override
    {
      // project initial data
      //TODO check whether this version work
      auto ftp = problem().fixedTimeFunction( tp.time() );
      GridFunctionAdapter< typename ProblemType::InstationaryFunctionType, GridPartType >
        adapter( "-exact", ftp, gridPart_, space_.order() + 2 );
      interpolate( adapter, solution() );
      if( NonBlockingCommParameter::nonBlockingCommunication() )
        solution().communicate();

      // setup ode solver
      solver_ = this->doCreateSolver( tp );

      // initialize ode solver
      solver()->initialize( solution() );

      //initialize solverMonitor
      if( solverMonitor_ )
      {
        solverMonitor_.registerData( "GridWidth", &solverMonitor_.monitor().gridWidth, nullptr, true );
        solverMonitor_.registerData( "Elements", &solverMonitor_.monitor().elements, nullptr, true );
        solverMonitor_.registerData( "TimeSteps", &solverMonitor_.monitor().timeSteps, nullptr, true );
        solverMonitor_.registerData( "AvgTimeStep", &solverMonitor_.monitor().avgTimeStep );
        solverMonitor_.registerData( "MinTimeStep", &solverMonitor_.monitor().minTimeStep );
        solverMonitor_.registerData( "MaxTimeStep", &solverMonitor_.monitor().maxTimeStep );
        solverMonitor_.registerData( "Newton", &solverMonitor_.monitor().newton_iterations,       &odeSolverMonitor_.newtonIterations_);
        solverMonitor_.registerData( "ILS", &solverMonitor_.monitor().ils_iterations,             &odeSolverMonitor_.linearSolverIterations_);
        //solverMonitor_.registerData( "OC", &solverMonitor_.monitor().operator_calls,              &odeSolverMonitor_.spaceOperatorCalls_);
        solverMonitor_.registerData( "MaxNewton",&solverMonitor_.monitor().max_newton_iterations, &odeSolverMonitor_.maxNewtonIterations_ );
        solverMonitor_.registerData( "MaxILS",&solverMonitor_.monitor().max_ils_iterations,    &odeSolverMonitor_.maxLinearSolverIterations_ );
      }

      //initialize diagnostics
      if( diagnostics_ )
      {
        diagnostics_.registerData( "OperatorTime", &odeSolverMonitor_.operatorTime_ );
        diagnostics_.registerData( "OdeSolveTime", &odeSolverMonitor_.odeSolveTime_ );
        diagnostics_.registerData( "OverallTimer", &overallTime_ );
        diagnostics_.registerData( "NumberOfElements", &odeSolverMonitor_.numberOfElements_ );
      }
    }

    virtual void doSolve ( const int loop, TimeProviderType& tp ) override
    {
      overallTimer_.reset();
      odeSolverMonitor_.reset();

      // solve ODE
      solver()->solve( solution(), odeSolverMonitor_ );

      overallTime_ = overallTimer_.stop();
    }

    virtual void doFinalize ( const int loop, TimeProviderType& tp ) override
    {
      // add eoc errors
      AnalyticalTraits::addEOCErrors( tp, solution(), model(), problem() );

      // delete ode solver
      solver_ = nullptr;
    }

  protected:
    const ContainerType&                    container_;
    double                                  overallTime_;
    Dune::Timer                             overallTimer_;
    GridPartType                            gridPart_;
    DiscreteFunctionSpaceType               space_;

    // the solution
    std::shared_ptr< DiscreteFunctionType > solution_;
    std::shared_ptr< DiscreteFunctionType > exactSolution_;

    std::shared_ptr< typename SolverType::type > solver_;
    std::unique_ptr< IOTupleType >               ioTuple_;

    DiagnosticsOptional< DiagnosticsType >           diagnostics_;
    SolverMonitorOptional< SolverMonitorType >       solverMonitor_;
    AdditionalOutputOptional< AdditionalOutputType > additionalOutput_;
    typename SolverType::type::MonitorType           odeSolverMonitor_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
