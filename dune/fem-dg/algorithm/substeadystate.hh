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

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/subevolution.hh>
#include <dune/fem-dg/algorithm/interface.hh>

#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>


namespace Dune
{

namespace Fem
{
  class IDGenerator
  {
    public:
      static IDGenerator& instance ()
      {
        static IDGenerator generator;
        return generator;
      }

      std::string nextId()
      {
        id_++;
        if( id_ == 0 )
          return "";
        std::stringstream s;
        s << "[" << id_ << "]";
        return s.str();
      }
      std::string id()
      {
        if( id_ == 0 )
          return "";
        std::stringstream s;
        s << "[" << id_ << "]";
        return s.str();
      }
    private:
      IDGenerator () : id_(-1) {}

      int id_;
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


  template<class Grid,
            class ProblemTraits,
            int polOrd >
  struct SubSteadyStateTraits
  {
  private:

    CHECK_TYPEDEF_EXISTS( RhsType )

  public:
    typedef typename ProblemTraits::AnalyticalTraits                   AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polOrd >  DiscreteTraits;

    typedef typename DiscreteTraits::ExtraParameterTuple               ExtraParameterTuple;

    typedef typename DiscreteTraits::Solver                            SolverType;

    typedef typename DiscreteTraits::Operator                          OperatorType;

    typedef typename RhsTypes< OperatorType >::type                    RhsType;
  };


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

    using BaseType::grid_;
    using BaseType::algorithmName_;
    using BaseType::problem;
    using BaseType::model;
    using BaseType::gridSize;

  public:
    using BaseType::checkSolutionValid;
    using BaseType::initialize;
    using BaseType::preSolve;
    using BaseType::solve;
    using BaseType::postSolve;
    using BaseType::finalize;

    SubSteadyStateAlgorithm ( GridType &grid, const std::string name = ""  )
      : BaseType( grid, name ),
        solverMonitorHandler_( name ),
        diagnosticsHandler_( name ),
        additionalOutputHandler_( name ),
        gridPart_( grid ),
        space_( gridPart_ ),
        solution_( "solution-" + name + IDGenerator::instance().nextId(), space_ ),
        exactSolution_( "exact-solution-" + name + IDGenerator::instance().id(), space_ ),
        rhsOperator_( gridPart_, problem(), std::tuple<>(), name ),
        ioTuple_( std::make_tuple( &solution(), &exactSolution_ ) ),
        solverIterations_()
    {
      solution().clear();
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    virtual DiscreteFunctionType& solution () { return solution_; }

    virtual DiscreteFunctionType& rhs() = 0;

    virtual typename SolverType::type* createSolver( DiscreteFunctionType* rhs ) = 0;

    virtual bool checkSolutionValid ( const int loop ) const { return solution_.dofsValid(); }

    virtual void initialize ( const int loop )
    {
      //initialize solverMonitor
      solverMonitorHandler_.registerData( "GridWidth", &solverMonitorHandler_.monitor().gridWidth, nullptr, true );
      solverMonitorHandler_.registerData( "Elements", &solverMonitorHandler_.monitor().elements, nullptr, true );
      solverMonitorHandler_.registerData( "TimeSteps", &solverMonitorHandler_.monitor().timeSteps, nullptr, true );
      solverMonitorHandler_.registerData( "ILS", &solverMonitorHandler_.monitor().ils_iterations, &solverIterations_ );
      solverMonitorHandler_.registerData( "MaxILS", &solverMonitorHandler_.monitor().max_ils_iterations );
    }

    virtual void preSolve( const int loop )
    {
      if( rhsOperator_ )
        rhsOperator_( solution(), rhs() );
      solution().clear();
      solver_.reset( this->createSolver( &rhs() ) );
    }

    virtual void solve ( const int loop )
    {
      (*solver_)( rhs(), solution() );
      solverIterations_ = solver_->iterations();
    }

    virtual void postSolve( const int loop )
    {
      monitor()->finalize( gridWidth(), gridSize() );
    }

    virtual void finalize ( const int loop )
    {
      // add eoc errors
      //AnalyticalTraits::addEOCErrors( solution(), model(), problem() );

      solver_.reset( nullptr );
    }

    //ADAPTATION
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return &solution_; }

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType* monitor() { return solverMonitorHandler_.value(); }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType* additionalOutput() { return additionalOutputHandler_.value(); }

    //DATAWRITING
    IOTupleType& dataTuple () { return ioTuple_; }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType* diagnostics() { return diagnosticsHandler_.value(); }

  protected:

    SolverMonitorHandlerOptional< SolverMonitorHandlerType >       solverMonitorHandler_;
    DiagnosticsHandlerOptional< DiagnosticsHandlerType >           diagnosticsHandler_;
    AdditionalOutputHandlerOptional< AdditionalOutputHandlerType > additionalOutputHandler_;

    GridPartType               gridPart_; // reference to grid part, i.e. the leaf grid
    DiscreteFunctionSpaceType  space_;    // the discrete function space
    DiscreteFunctionType       solution_;
    DiscreteFunctionType       exactSolution_;

    RhsOptional<RhsType >      rhsOperator_;

    IOTupleType                ioTuple_;

    std::unique_ptr< typename SolverType::type > solver_;
    int solverIterations_;
  };


}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
