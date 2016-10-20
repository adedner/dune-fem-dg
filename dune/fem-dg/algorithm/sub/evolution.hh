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
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>

#include "container.hh"

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

    //short cut for tuple containing shared_ptr
    template< class... Args >
    using shared_tuple = std::tuple< std::shared_ptr<Args>... >;
  public:
    template< int i >
    using DiscreteFunction = DiscreteFunctionType;
    template< int i >
    using DiscreteFunctionSpace = typename DiscreteFunction<i>::DiscreteFunctionSpaceType;
    template< int i >
    using GridPart = typename DiscreteFunctionSpace<i>::GridPartType;
    template< int i >
    using Grid = typename GridPart<i>::GridType;

    // owning container
    SubEvolutionContainer( Grid<0>& grid, const std::string name = "" )
    : stringId_( "123" ),
      grid_( &grid ),
      gridPart_( std::make_tuple( std::make_shared< GridPart<0> >( *grid_ ) ) ),
      space_( std::make_tuple( std::make_shared<DiscreteFunctionSpace<0> >( *std::get<0>(gridPart_) ) ) ),
      solution_( std::make_tuple( std::make_shared<DiscreteFunction<0> >( name + "u" + stringId_, *space<0>() ) ) ),
      exactSolution_( std::make_tuple( std::make_shared<DiscreteFunction<0> >( name + "u-exact" + stringId_, *space<0>() ) ) )
    {}

    // non owning container, for coupling
    SubEvolutionContainer( const std::string name = "" )
    : stringId_( "123" ),
      grid_( nullptr ),
      gridPart_( std::make_tuple( nullptr ) ),
      space_( std::make_tuple( nullptr ) ),
      solution_( std::make_tuple( nullptr ) ),
      exactSolution_( std::make_tuple( nullptr ) )
    {}

    //spaces
    template< int i >
    std::shared_ptr< DiscreteFunctionSpace<i> > space() const
    {
      assert( std::get<i>(space_) );
      return std::get<i>(space_);
    }

    //solution
    template< int i >
    std::shared_ptr< DiscreteFunction<i> > solution() const
    {
      return std::get<i>(solution_);
    }
    template< int i >
    void setSolution( std::shared_ptr< DiscreteFunction<i> > solution )
    {
      std::get<i>(solution_) = solution;
    }

    //exact solution
    template< int i >
    std::shared_ptr< DiscreteFunction<i> > exactSolution() const
    {
      return std::get<i>( exactSolution_);
    }
    template< int i >
    void setExactSolution( std::shared_ptr< DiscreteFunction<i> > exactSolution )
    {
     std::get<i>( exactSolution_ )= exactSolution;
    }

  private:
    const std::string                         stringId_;
    Grid<0>*                                 grid_;
    shared_tuple< GridPart<0> >              gridPart_;
    shared_tuple< DiscreteFunctionSpace<0> > space_;

    shared_tuple< DiscreteFunction<0> >      solution_;
    shared_tuple< DiscreteFunction<0> >      exactSolution_;
    //shared_tuple< DiscreteFunction<0> >      rhs_;
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
    typedef typename BaseType::DiagnosticsType                       DiagnosticsType;
    typedef typename BaseType::SolverMonitorType                     SolverMonitorType;
    typedef typename BaseType::AdditionalOutputType                  AdditionalOutputType;

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
      solution_( container_.template solution<0>() ),
      exactSolution_( container_.template exactSolution<0>() ),
      solver_( nullptr ),
      ioTuple_( new IOTupleType( std::make_tuple( solution_.get(), exactSolution_.get() ) ) ),
      diagnostics_( name() ),
      solverMonitor_( name() ),
      additionalOutput_( exactSolution() ),
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
    virtual double gridWidth () const { return GridWidth::calcGridWidth( solution_->space().gridPart() ); }

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
        adapter( "-exact", ftp, solution_->space().gridPart(), solution_->space().order() + 2 );
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
