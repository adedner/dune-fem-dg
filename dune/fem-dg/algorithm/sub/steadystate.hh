#ifndef DUNE_FEMDG_ALGORITHM_STEADYSTATEALGORITHM_HH
#define DUNE_FEMDG_ALGORITHM_STEADYSTATEALGORITHM_HH

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/operator/linear/spoperator.hh>

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>
#include <dune/fem-dg/misc/uniquefunctionname.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/sub/evolution.hh>
#include <dune/fem-dg/algorithm/sub/interface.hh>

#include <dune/fem-dg/algorithm/handler/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/sub/additionaloutput.hh>


namespace Dune
{
namespace Fem
{

  template <class DiscreteFunctionImp >
  struct SubSteadyStateContainer
  {

    typedef DiscreteFunctionImp                                           DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridType                  GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType              GridPartType;

    using DiscreteFunction = DiscreteFunctionType;
    using DiscreteFunctionSpace = DiscreteFunctionSpaceType;

  public:

    SubSteadyStateContainer( GridType& grid, const std::string name = "" )
    : stringId_( FunctionIDGenerator::instance().nextId() ),
      grid_( grid ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      solution_( new DiscreteFunctionType( name + "u" + stringId_, space() ) ),
      exactSolution_( new DiscreteFunctionType( name + "u-exact" + stringId_, space() ) ),
      rhs_( new DiscreteFunctionType( name + "rhs-u" + stringId_, space() ) )
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

    //rhs
    std::shared_ptr< DiscreteFunction > rhs() const
    {
      return rhs_;
    }
    void setRhs( std::shared_ptr< DiscreteFunction > rhs )
    {
      rhs_ = rhs;
    }

  private:
    const std::string                   stringId_;
    GridType&                           grid_;
    GridPartType                        gridPart_;
    DiscreteFunctionSpaceType           space_;

    std::shared_ptr< DiscreteFunction > solution_;
    std::shared_ptr< DiscreteFunction > exactSolution_;
    std::shared_ptr< DiscreteFunction > rhs_;
  };



  template< class Obj >
  class RhsOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    RhsOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };


  class EmptyOperator
  {
  public:
    template< class... Args >
    EmptyOperator( Args&&... ){}

    template< class... Args >
    void operator()( Args&&... args )
    {}
  };

  template<>
  class RhsOptional< void >
    : public OptionalNullPtr< EmptyOperator >
  {
    typedef OptionalNullPtr< EmptyOperator >    BaseType;
  public:
    template< class... Args >
    RhsOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };


  template< class Grid,
            class ProblemTraits,
            int polOrd >
  struct SubSteadyStateTraits
  {
  private:

    CHECK_TYPEDEF_EXISTS( RhsType )

  public:
    typedef typename ProblemTraits::AnalyticalTraits                   AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polOrd >  DiscreteTraits;

    typedef typename DiscreteTraits::Solver                            SolverType;

    typedef typename DiscreteTraits::Operator                          OperatorType;

    typedef typename RhsTypes< OperatorType >::type                    RhsType;
  };

  /**
   *  \brief Algorithm for solving a stationary PDE.
   *
   *  \ingroup Algorithms
   */
  template< class Grid, class ProblemTraits, int polOrder >
  class SubSteadyStateAlgorithm
    : public SubAlgorithmInterface< Grid, ProblemTraits, polOrder >
  {
    typedef SubAlgorithmInterface< Grid, ProblemTraits, polOrder >     BaseType;
    typedef SubSteadyStateTraits< Grid, ProblemTraits, polOrder >    Traits;

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

    // The DG space operator
    typedef typename Traits::OperatorType                            OperatorType;

    // type of steady state solver
    typedef typename Traits::SolverType                              SolverType;

    // type of analytical traits
    typedef typename Traits::AnalyticalTraits                        AnalyticalTraits;

    // type of discrete traits
    typedef typename Traits::DiscreteTraits                          DiscreteTraits;

    typedef typename Traits::RhsType                                 RhsType;

  public:
    using BaseType::grid;
    using BaseType::name;
    using BaseType::problem;
    using BaseType::model;
    using BaseType::gridSize;

  public:

    typedef SubSteadyStateContainer< DiscreteFunctionType >          ContainerType;

    SubSteadyStateAlgorithm ( GridType &grid, ContainerType& container  )
      : BaseType( grid ),
        container_( container ),
        gridPart_( container_.gridPart() ),
        space_( container_.space() ),
        solverIterations_( 0 ),
        solution_( container_.solution() ),
        exactSolution_( container_.exactSolution() ),
        rhs_( container_.rhs() ),
        rhsOperator_( doCreateRhsOperator() ),
        solver_( nullptr ),
        ioTuple_( new IOTupleType( std::make_tuple( solution_.get(), exactSolution_.get() ) ) ),
        solverMonitorHandler_( name() ),
        diagnosticsHandler_( name() ),
        additionalOutputHandler_( name() )
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
    const DiscreteFunctionType& solution () const
    {
      assert( solution_ );
      return *solution_;
    }

    DiscreteFunctionType& rhs ()
    {
      assert( rhs_ );
      return *rhs_;
    }

    DiscreteFunctionType& exactSolution ()
    {
      assert( exactSolution_ );
      return *exactSolution_;
    }

    virtual void setTime( const double time ){}

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    //ADAPTATION
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return solution_.get(); }

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType* monitor() { return solverMonitorHandler_.value(); }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType* additionalOutput() { return additionalOutputHandler_.value(); }

    //DATAWRITING
    IOTupleType& dataTuple () { return *ioTuple_; }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType* diagnostics()
    {
      return diagnosticsHandler_.value();
    }

  protected:

    virtual std::shared_ptr< RhsOptional< RhsType > > doCreateRhsOperator()
    {
      return std::make_shared< RhsOptional< RhsType > >(gridPart_, problem(), std::tuple<>(), name() );
    }

    virtual std::shared_ptr< typename SolverType::type > doCreateSolver()
    {
      return nullptr;
    }

    virtual bool doCheckSolutionValid ( const int loop ) const override
    {
      assert( solution_ );
      return solution_->dofsValid();
    }

    virtual void doInitialize ( const int loop )
    {
      //initialize solverMonitor
      solverMonitorHandler_.registerData( "GridWidth", &solverMonitorHandler_.monitor().gridWidth, nullptr, true );
      solverMonitorHandler_.registerData( "Elements", &solverMonitorHandler_.monitor().elements, nullptr, true );
      solverMonitorHandler_.registerData( "TimeSteps", &solverMonitorHandler_.monitor().timeSteps, nullptr, true );
      solverMonitorHandler_.registerData( "ILS", &solverMonitorHandler_.monitor().ils_iterations, &solverIterations_ );
      solverMonitorHandler_.registerData( "MaxILS", &solverMonitorHandler_.monitor().max_ils_iterations );
    }

    virtual void doPreSolve ( const int loop )
    {
      if( *rhsOperator_ ) //rhs by external rhs operator
        (*rhsOperator_)( solution(), rhs() );
      solution().clear();
      solver_ = this->doCreateSolver();
    }

    virtual void doSolve ( const int loop )
    {
      Dune::Timer timer;
      double time = 0;
      timer.reset();
      (*solver_)( rhs(), solution() );
      solverIterations_ = solver_->iterations();
      time = timer.stop();
      std::cout << "Solve time: " << time << std::endl;
    }

    virtual void doPostSolve( const int loop )
    {
      monitor()->finalize( gridWidth(), gridSize() );
    }

    virtual void doFinalize ( const int loop )
    {
      // add eoc errors
      //AnalyticalTraits::addEOCErrors( solution(), model(), problem() );

      solver_ = nullptr;
    }
  protected:

    ContainerType&             container_;

    GridPartType&                                gridPart_; // reference to grid part, i.e. the leaf grid
    const DiscreteFunctionSpaceType&             space_;    // the discrete function space
    int                                          solverIterations_;

    std::shared_ptr< DiscreteFunctionType >      solution_;
    std::shared_ptr< DiscreteFunctionType >      exactSolution_;

    std::shared_ptr< DiscreteFunctionType >      rhs_;
    std::shared_ptr< RhsOptional< RhsType > >    rhsOperator_;

    std::unique_ptr< IOTupleType >               ioTuple_;

    std::shared_ptr< typename SolverType::type > solver_;

    SolverMonitorHandlerOptional< SolverMonitorHandlerType >       solverMonitorHandler_;
    DiagnosticsHandlerOptional< DiagnosticsHandlerType >           diagnosticsHandler_;
    AdditionalOutputHandlerOptional< AdditionalOutputHandlerType > additionalOutputHandler_;
  };


}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
