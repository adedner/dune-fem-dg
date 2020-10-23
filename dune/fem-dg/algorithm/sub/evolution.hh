#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH

#include <iostream>
#include <string>
#include <type_traits>
#include <memory>

#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>

#include <dune/fem/misc/femtimer.hh>
#include <dune/common/timer.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem-dg/misc/uniquefunctionname.hh>

#include <dune/fem-dg/algorithm/sub/interface.hh>
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>

#include "containers.hh"

#include <dune/fem-dg/misc/interpolatethreaded.hh>

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


    typedef typename ProblemTraits::template DiscreteTraits< polOrder > DiscreteTraits;

    typedef typename DiscreteTraits::Operator                           OperatorType;

    typedef typename DiscreteTraits::Solver                             SolverType;

  };

  template <class DiscreteFunctionImp >
  struct SubEvolutionContainerItem
  {
    using CItem = ContainerItem< DiscreteFunctionImp >;
  public:
    using DiscreteFunction = DiscreteFunctionImp;

    // owning container
    template< class SameObject >
    SubEvolutionContainerItem( const std::shared_ptr<SameObject>& obj, const std::string name = "" )
    : stringId_( FunctionIDGenerator::instance().nextId() ),
      solution_(      std::make_shared< CItem >( name + "u" + stringId_, obj ) ),
      exactSolution_( std::make_shared< CItem >( name + "u-exact" + stringId_, *solution_ ) )
    {}

    //solution
    std::shared_ptr< DiscreteFunction > solution() const
    {
      return solution_->shared();
    }

    //exact solution
    std::shared_ptr< DiscreteFunction > exactSolution() const
    {
      return exactSolution_->shared();
    }
  private:
    const std::string          stringId_;

    std::shared_ptr< CItem > solution_;
    std::shared_ptr< CItem > exactSolution_;
  };





  template <class... DiscreteFunctions >
  struct SubEvolutionContainer
  : public OneArgContainer< SubEvolutionContainerItem, std::tuple< DiscreteFunctions... > >
    {
    typedef OneArgContainer< SubEvolutionContainerItem, std::tuple< DiscreteFunctions...> > BaseType;
  public:
    using BaseType::operator();

    // constructor: do not touch/delegate everything
    template< class ... Args>
    SubEvolutionContainer( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}
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

    template< class ContainerImp, class ExtraArgsImp >
    SubEvolutionAlgorithm( const std::shared_ptr< ContainerImp >& cont,
                           const std::shared_ptr< ExtraArgsImp >& extra )
    : BaseType( cont, extra )
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
  protected:
    using BaseType::model_;
  public:
    typedef typename BaseType::GridType                              GridType;
    typedef typename BaseType::GridPartType                          GridPartType;
    typedef typename BaseType::HostGridPartType                      HostGridPartType;

    typedef typename BaseType::TimeProviderType                      TimeProviderType;

    typedef typename BaseType::DiscreteFunctionType                  DiscreteFunctionType;

    typedef typename BaseType::UInt64Type                            UInt64Type;

    typedef typename BaseType::CheckPointDiscreteFunctionType        CheckPointDiscreteFunctionType;
    typedef typename BaseType::LimitDiscreteFunctionType             LimitDiscreteFunctionType;
    typedef typename BaseType::AdaptationDiscreteFunctionType        AdaptationDiscreteFunctionType;

    typedef typename BaseType::AdaptIndicatorType                    AdaptIndicatorType;
    typedef typename BaseType::DiagnosticsType                       DiagnosticsType;
    typedef typename BaseType::SolverMonitorType                     SolverMonitorType;
    typedef typename BaseType::DataWriterType                        DataWriterType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // The DG space operator
    typedef typename Traits::OperatorType                            OperatorType;

    // type of steady state solver
    typedef typename Traits::SolverType                              SolverType;

    // type of discrete traits
    typedef typename Traits::DiscreteTraits                          DiscreteTraits;

    using BaseType::grid;
    using BaseType::name;
    using BaseType::gridSize;

    typedef SubEvolutionContainer< DiscreteFunctionType >               ContainerType;

    template< class ContainerImp, class ExtraArgsImp >
    SubEvolutionAlgorithmBase( const std::shared_ptr< ContainerImp >& cont,
                               const std::shared_ptr< ExtraArgsImp >& extra )
    : BaseType( const_cast< GridType& >( (*cont)(_0)->solution()->gridPart().grid() ) ),
      overallTime_( 0 ),
      solution_( (*cont)(_0)->solution() ),
      exactSolution_( (*cont)(_0)->exactSolution() ),
      solver_( nullptr ),
      diagnostics_( name() ),
      solverMonitor_( name() ),
      dataWriter_( name() ),
      odeSolverMonitor_()
    {}

    typename SolverType::type* solver()
    {
      return solver_.get();
    }

    DiscreteFunctionType& solution () override
    {
      assert( solution_ );
      return *solution_;
    }
    const DiscreteFunctionType& solution () const
    {
      assert( solution_ );
      return *solution_;
    }

    //DiscreteFunctionType& exactSolution ()
    //{
    //  assert( exactSolution_ );
    //  return *exactSolution_;
    //}
    //const DiscreteFunctionType& exactSolution () const
    //{
    //  assert( exactSolution_ );
    //  return *exactSolution_;
    //}

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const override { return GridWidth::calcGridWidth( solution_->space().gridPart() ); }

    //SOLVERMONITOR
    virtual SolverMonitorType* monitor() override { return solverMonitor_.value(); }

    //DIAGNOSTICS
    virtual DiagnosticsType* diagnostics() override { return diagnostics_.value(); }

    //DATAWRITER
    virtual DataWriterType* dataWriter() override { return dataWriter_.value(); }

    //LIMITING
    virtual void limit() override {}
    virtual LimitDiscreteFunctionType* limitSolution () override { return solution_.get(); }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() override { return nullptr; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () override { return solution_.get(); }

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () override { return solution_.get(); }

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
      auto ftp = model_.problem().fixedTimeFunction( tp.time() );
      interpolateThreaded( gridFunctionAdapter( ftp, solution_->space().gridPart(), solution_->space().order() + 2 ), solution() );
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

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      {
        const int polOrder = solution().space().order();
        const int maxOrder = 2 * polOrder + 1 ;
        const int dim = GridType :: dimension;
        std::vector< int > quadOrder( maxOrder );
        for( int i=2; i<maxOrder+2; ++i )
          quadOrder.push_back( i );
        std::string filename( "autogeneratedcode_"+std::to_string( dim ) + "d_k" + std::to_string( polOrder ) + ".hh" );
        Dune::Fem::generateCode( solution().space(), quadOrder, quadOrder, "./", filename );

      }
#endif
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
      // set current time to model for evaluation of exact solution
      model_.setTime( tp.time() );

      // add eoc errors
      model_.eocErrors( solution() );

      // delete ode solver
      solver_ = nullptr;
    }

  protected:
    double                                  overallTime_;
    Dune::Timer                             overallTimer_;

    // the solution
    std::shared_ptr< DiscreteFunctionType > solution_;
    std::shared_ptr< DiscreteFunctionType > exactSolution_;

    std::shared_ptr< typename SolverType::type > solver_;

    DiagnosticsOptional< DiagnosticsType >           diagnostics_;
    SolverMonitorOptional< SolverMonitorType >       solverMonitor_;
    DataWriterOptional< DataWriterType >         dataWriter_;
    typename SolverType::type::MonitorType           odeSolverMonitor_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
