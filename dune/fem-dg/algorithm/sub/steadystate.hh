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

#include <dune/fem-dg/misc/integral_constant.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>
#include <dune/fem-dg/misc/uniquefunctionname.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/sub/evolution.hh>
#include <dune/fem-dg/algorithm/sub/interface.hh>

#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>

#include "containers.hh"

namespace Dune
{
namespace Fem
{

  template <class DiscreteFunctionImp >
  struct SubSteadyStateContainerItem
  {
    using CItem = ContainerItem< DiscreteFunctionImp >;
  public:
    using DiscreteFunction = DiscreteFunctionImp;

    // owning container
    template< class SameObject >
    SubSteadyStateContainerItem( const std::shared_ptr<SameObject>& obj, const std::string name = "" )
    : stringId_( FunctionIDGenerator::instance().nextId() ),
      solution_(      std::make_shared< CItem >( name + "u" + stringId_, obj ) ),
      exactSolution_( std::make_shared< CItem >( name + "u-exact" + stringId_, *solution_ ) ),
      rhs_(           std::make_shared< CItem >( name + "u-rhs" + stringId_, *solution_ ) )
    {}

    //solution
    shared_ptr< DiscreteFunction > solution() const
    {
      return solution_->shared();
    }

    //exact solution
    shared_ptr< DiscreteFunction > exactSolution() const
    {
      return exactSolution_->shared();
    }

    //rhs
    shared_ptr< DiscreteFunction > rhs() const
    {
      return rhs_->shared();
    }
  private:
    const std::string          stringId_;

    std::shared_ptr< CItem > solution_;
    std::shared_ptr< CItem > exactSolution_;
    std::shared_ptr< CItem > rhs_;
  };


  template <class... DiscreteFunctions >
  struct SubSteadyStateContainer
  : public OneArgContainer< SubSteadyStateContainerItem, std::tuple< DiscreteFunctions... > >
  {
    typedef OneArgContainer< SubSteadyStateContainerItem, std::tuple< DiscreteFunctions...> > BaseType;
  public:
    using BaseType::operator();

    // constructor: do not touch/delegate everything
    template< class ... Args>
    SubSteadyStateContainer( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}
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
    typedef typename ProblemTraits::template DiscreteTraits< polOrd >  DiscreteTraits;

    typedef typename DiscreteTraits::Solver                            SolverType;

    typedef typename DiscreteTraits::Operator                          OperatorType;

    typedef typename RhsTypes< OperatorType >::type                    RhsType;
  };

  /**
   *  \brief Algorithm for solving a stationary PDE.
   *
   *  \ingroup SubAlgorithms
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
    typedef typename BaseType::DiagnosticsType                       DiagnosticsType;
    typedef typename BaseType::SolverMonitorType                     SolverMonitorType;
    typedef typename BaseType::AdditionalOutputType                  AdditionalOutputType;


    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // The DG space operator
    typedef typename Traits::OperatorType                            OperatorType;

    // type of steady state solver
    typedef typename Traits::SolverType                              SolverType;

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

    //type for a standalone container
    typedef SubSteadyStateContainer< DiscreteFunctionType >          ContainerType;

    template< class ContainerImp, class ExtraArgsImp >
    SubSteadyStateAlgorithm( const std::shared_ptr< ContainerImp >& cont,
                             const std::shared_ptr< ExtraArgsImp >& extra )
      : BaseType( const_cast< GridType& >( (*cont)(_0)->solution()->gridPart().grid() ) ),
        solverIterations_( 0 ),
        //newContainer_(  cont ? nullptr : std::make_shared< ContainerType >( grid, "name" ) )
        //solution_( cont ? : (*cont)(_0)->solution() : (*newContainer_)(_0)->solution() ),
        solution_( (*cont)(_0)->solution() ),
        exactSolution_( (*cont)(_0)->exactSolution() ),
        rhs_( (*cont)(_0)->rhs() ),
        rhsOperator_( std::make_shared< RhsOptional< RhsType > >( solution().space().gridPart(), model(), extra, name() ) ),
        ioTuple_( std::make_unique<IOTupleType>( std::make_tuple( solution_.get(), exactSolution_.get() ) ) ),
        solver_( nullptr ),
        solverMonitor_( name() ),
        diagnostics_( name() ),
        additionalOutput_( name() )
    {}

    SubSteadyStateAlgorithm ( GridType &grid  )
      : SubSteadyStateAlgorithm( grid, std::make_shared< ContainerType >( grid, "name" ) )
    {}

    typename SolverType::type* solver() const
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
    virtual double gridWidth () const { return GridWidth::calcGridWidth( solution().space().gridPart() ); }

    //ADAPTATION
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return solution_.get(); }

    //SOLVERMONITOR
    virtual SolverMonitorType* monitor() { return solverMonitor_.value(); }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputType* additionalOutput() { return additionalOutput_.value(); }

    //DATAWRITING
    IOTupleType& dataTuple () { return *ioTuple_; }

    //DIAGNOSTICS
    virtual DiagnosticsType* diagnostics()
    {
      return diagnostics_.value();
    }

  protected:
    virtual std::shared_ptr< typename SolverType::type > doCreateSolver()
    {
      return nullptr;
    }

    virtual bool doCheckSolutionValid ( const int loop ) const override
    {
      return solution().dofsValid();
    }

    virtual void doInitialize ( const int loop )
    {
      //initialize solverMonitor
      solverMonitor_.registerData( "GridWidth", &solverMonitor_.monitor().gridWidth, nullptr, true );
      solverMonitor_.registerData( "Elements", &solverMonitor_.monitor().elements, nullptr, true );
      solverMonitor_.registerData( "TimeSteps", &solverMonitor_.monitor().timeSteps, nullptr, true );
      solverMonitor_.registerData( "ILS", &solverMonitor_.monitor().ils_iterations, &solverIterations_ );
      solverMonitor_.registerData( "MaxILS", &solverMonitor_.monitor().max_ils_iterations );
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
      //ProblemTraits::addEOCErrors( solution(), model(), problem() );

      solver_ = nullptr;
    }
  protected:

    int                                          solverIterations_;

    std::shared_ptr< DiscreteFunctionType >      solution_;
    std::shared_ptr< DiscreteFunctionType >      exactSolution_;

    std::shared_ptr< DiscreteFunctionType >      rhs_;
    std::shared_ptr< RhsOptional< RhsType > >    rhsOperator_;

    std::unique_ptr< IOTupleType >               ioTuple_;

    std::shared_ptr< typename SolverType::type > solver_;

    SolverMonitorOptional< SolverMonitorType >       solverMonitor_;
    DiagnosticsOptional< DiagnosticsType >           diagnostics_;
    AdditionalOutputOptional< AdditionalOutputType > additionalOutput_;
  };


}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
