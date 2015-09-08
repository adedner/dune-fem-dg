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

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/evolution.hh>

#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>


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

  template<class Grid,
            class ProblemTraits,
            int polOrd >
  struct SubSteadyStateTraits
  {
     private:
    // if not AdaptIndicatorType is defined, use void as "marker"
    template< class T, bool exists >
    struct AdaptIndicatorHelper
    { typedef void type; };

    template< class T >
    struct AdaptIndicatorHelper< T, true >
    { typedef typename T::AdaptIndicatorType type; };

    template< class T >
    struct AdaptIndicators
    {
      template <typename TT>
      static auto apply(TT const&) -> decltype( typename TT::AdaptIndicatorType(), std::true_type()) {
        return std::true_type();
      }
      static std::false_type apply(...) { return std::false_type(); }

      typedef typename AdaptIndicatorHelper< T, decltype( apply( std::declval<T>() ) )::value >::type type;
    };

    // if AdditionalOutputHandlerType is not defined, use void as "marker"
    template< class T, bool exists >
    struct AdditionalOutputHandlerHelper
    { typedef void type; };

    template< class T >
    struct AdditionalOutputHandlerHelper< T, true >
    { typedef typename T::AdditionalOutputHandlerType type; };

    template< class T >
    struct AdditionalOutputHandlers
    {
      template <typename TT>
      static auto apply(TT const&) -> decltype( typename TT::AdditionalOutputHandlerType(), std::true_type()) {
        return std::true_type();
      }
      static std::false_type apply(...) { return std::false_type(); }

      typedef typename AdditionalOutputHandlerHelper< T, decltype( apply( std::declval<T>() ) )::value >::type type;
    };


  public:
    enum { polynomialOrder = polOrd };

    // type of Grid
    typedef Grid                                                   GridType;
    typedef typename ProblemTraits::HostGridPartType               HostGridPartType;
    typedef typename ProblemTraits::GridPartType                   GridPartType;

    typedef typename ProblemTraits::AnalyticalTraits               AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polynomialOrder >  DiscreteTraits;

    // obtain the problem dependent types, analytical context
    typedef typename AnalyticalTraits::ModelType                   ModelType;
    typedef typename AnalyticalTraits::ProblemType                 ProblemType;

    // type of discrete function space and discrete function
    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;

    typedef typename DiscreteTraits::ExtraParameterTuple           ExtraParameterTuple;

    typedef typename DiscreteTraits::Solver                        SolverType;

    typedef typename DiscreteTraits::Operator                      OperatorType;

    // type of IOTuple
    typedef typename DiscreteTraits::IOTupleType                   IOTupleType;

    typedef typename AdaptIndicators< DiscreteTraits >::type       AdaptIndicatorType;
    typedef typename AdditionalOutputHandlers< DiscreteTraits >::type AdditionalOutputHandlerType;
    typedef typename DiscreteTraits::SolverMonitorHandlerType      SolverMonitorHandlerType;
    typedef typename DiscreteTraits::DiagnosticsHandlerType        DiagnosticsHandlerType;
  };


  template< class Grid, class ProblemTraits, int polOrder >
  class SubSteadyStateAlgorithm
  {
    typedef SubSteadyStateTraits< Grid, ProblemTraits, polOrder >    Traits;

  public:
    typedef typename Traits::GridType                           GridType;
    typedef typename Traits::IOTupleType                        IOTupleType;

    typedef typename Traits::HostGridPartType                     HostGridPartType;

    // initial data type
    typedef typename Traits::ProblemType                          ProblemType;

    // An analytical version of our model
    typedef typename Traits::ModelType                            ModelType;

    typedef typename Traits::GridPartType                         GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType            DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType                 DiscreteFunctionType;

    // The DG space operator
    typedef typename Traits::OperatorType                         OperatorType;

    // type of steady state solver
    typedef typename Traits::SolverType                           SolverType;

    // type of analytical traits
    typedef typename Traits::AnalyticalTraits                     AnalyticalTraits;

    // type of discrete traits
    typedef typename Traits::DiscreteTraits                       DiscreteTraits;

    typedef uint64_t                                              UInt64Type ;

    typedef DiscreteFunctionType                                  CheckPointDiscreteFunctionType;
    typedef DiscreteFunctionType                                  LimitDiscreteFunctionType;
    typedef DiscreteFunctionType                                  AdaptationDiscreteFunctionType;

    typedef typename Traits::AdaptIndicatorType                   AdaptIndicatorType;
    typedef typename Traits::AdditionalOutputHandlerType          AdditionalOutputHandlerType;
    typedef typename Traits::SolverMonitorHandlerType             SolverMonitorHandlerType;
    typedef typename Traits::DiagnosticsHandlerType               DiagnosticsHandlerType;

    SubSteadyStateAlgorithm ( GridType &grid, const std::string name = ""  )
      : grid_( grid ),
        algorithmName_( name ),
        problem_( ProblemTraits::problem() ),
        model_( problem() ),
        solverMonitorHandler_( "" ),
        diagnosticsHandler_( "" ),
        gridPart_( grid ),
        space_( gridPart_ ),
        solution_( "solution-" + name + IDGenerator::instance().nextId(), space_ ),
        exactSolution_( "exact-solution-" + name + IDGenerator::instance().id(), space_ ),
        solverIterations_()
    {
      solution().clear();
    }

    virtual const std::string name () { return algorithmName_; }

    GridType& grid () const { return grid_; }

    //! return reference to discrete space
    DiscreteFunctionSpaceType& space () { return space_; }
    const DiscreteFunctionSpaceType& space () const { return space_; }

    virtual DiscreteFunctionType& solution () { return solution_; }

    virtual DiscreteFunctionType& rhs() = 0;

    std::string description () const { return problem().description(); }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const
    {
      return problem().dataPrefix();
    }

    virtual AdaptIndicatorType* adaptIndicator() { return nullptr; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return &solution_; }

    virtual SolverMonitorHandlerType& monitor()
    {
      return solverMonitorHandler_;
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type grSize = grid().size( 0 );
      return grid().comm().sum( grSize );
    }

    virtual typename SolverType::type* createSolver( DiscreteFunctionType* rhs ) = 0;

    virtual bool checkSolutionValid ( const int loop  ) const { return solution_.dofsValid(); }

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
      solverMonitorHandler_.finalize( gridWidth(), gridSize() );
    }

    void finalize ( const int loop )
    {
      // add eoc errors
      AnalyticalTraits::addEOCErrors( solution(), model(), problem() );

      solver_.reset( nullptr );
    }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType* additionalOutput() { return nullptr; }

    //DATAWRITING
    IOTupleType dataTuple () { return std::make_tuple( &solution(), nullptr ); }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType& diagnostics() { return diagnosticsHandler_; }

    const ProblemType &problem () const
    {
      assert( problem_ );
      return *problem_;
    }

    ProblemType &problem ()
    {
      assert( problem_ );
      return *problem_;
    }

    ModelType &model () { return model_; }
    const ModelType &model () const { return model_; }

  protected:
    GridType&   grid_;
    std::string algorithmName_;
    ProblemType *problem_;
    ModelType model_;

    SolverMonitorHandlerType       solverMonitorHandler_;
    DiagnosticsHandlerType         diagnosticsHandler_;

    GridPartType gridPart_;      // reference to grid part, i.e. the leaf grid
    DiscreteFunctionSpaceType space_;    // the discrete function space
    DiscreteFunctionType solution_;
    DiscreteFunctionType exactSolution_;

    std::unique_ptr< typename SolverType::type > solver_;
    int solverIterations_;
  };

  // SubEvolutionAlgorithmTraits
  // -------------------------

  template< class Grid,
            class ProblemTraits,
            int polOrder >
  struct SubEvolutionAlgorithmTraits2
  {
  private:
    // if AdaptIndicatorType is not defined, use "void" as "marker"
    template< class T, bool exists >
    struct AdaptIndicatorHelper
    { typedef void type; };

    template< class T >
    struct AdaptIndicatorHelper< T, true >
    { typedef typename T::AdaptIndicatorType type; };

    template< class T >
    struct AdaptIndicators
    {
      template <typename TT>
      static auto apply(TT const&) -> decltype( std::declval<typename TT::AdaptIndicatorType>(), std::true_type()) {
        return std::true_type();
      }
      static std::false_type apply(...) { return std::false_type(); }

      typedef typename AdaptIndicatorHelper< T, decltype( apply( std::declval<T>() ) )::value >::type type;
    };

    // if AdditionalOutputHandlerType is not defined, use "void" as "marker"
    template< class T, bool exists >
    struct AdditionalOutputHandlerHelper
    { typedef void type; };

    template< class T >
    struct AdditionalOutputHandlerHelper< T, true >
    { typedef typename T::AdditionalOutputHandlerType type; };

    template< class T >
    struct AdditionalOutputHandlers
    {
      template <typename TT>
      static auto apply(TT const&) -> decltype( std::declval<typename TT::AdditionalOutputHandlerType>(), std::true_type()) {
        return std::true_type();
      }
      static std::false_type apply(...) { return std::false_type(); }

      typedef typename AdditionalOutputHandlerHelper< T, decltype( apply( std::declval<T>() ) )::value >::type type;
    };

  public:
    static const int polynomialOrder = polOrder;

    typedef ProblemTraits                                          ProblemTraitsType;

    // type of Grid
    typedef Grid                                                   GridType;
    typedef typename ProblemTraits::HostGridPartType               HostGridPartType;
    typedef typename ProblemTraits::GridPartType                   GridPartType;

    typedef typename ProblemTraits::AnalyticalTraits               AnalyticalTraits;
    typedef typename ProblemTraits::DiscreteTraits                 DiscreteTraits;

    // obtain the problem dependent types, analytical context
    typedef typename AnalyticalTraits::ModelType                   ModelType;
    typedef typename AnalyticalTraits::ProblemType                 ProblemType;
    typedef typename AnalyticalTraits::InitialDataType             InitialDataType;

    // type of dg operator
    typedef typename DiscreteTraits::Operator                      OperatorType;

    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;

    typedef typename DiscreteTraits::ExtraParameterTuple           ExtraParameterTupleType;
    typedef typename DiscreteTraits::IOTupleType                   IOTupleType;

    // wrap operator
    typedef GridTimeProvider< GridType >                           TimeProviderType;

    typedef typename DiscreteTraits::Solver                        SolverType;

    typedef typename AdaptIndicators< DiscreteTraits >::type       AdaptIndicatorType;
    typedef typename AdditionalOutputHandlers< DiscreteTraits >::type AdditionalOutputHandlerType;


    typedef typename DiscreteTraits::SolverMonitorHandlerType      SolverMonitorHandlerType;
    typedef typename DiscreteTraits::DiagnosticsHandlerType        DiagnosticsHandlerType;
  };



  template< class Grid, class SteadyStateProblemTraits, int polOrder >
  class SubSteadyState2EvolutionAlgorithm
    : public SteadyStateProblemTraits
  {
    typedef SubEvolutionAlgorithmTraits2< Grid, SteadyStateProblemTraits, polOrder > Traits;
    typedef SteadyStateProblemTraits                                                 BaseType;

  public:
    typedef typename Traits::GridType                           GridType;
    typedef typename Traits::IOTupleType                        IOTupleType;

    typedef typename Traits::HostGridPartType                     HostGridPartType;

    // initial data type
    typedef typename Traits::ProblemType                          ProblemType;

    // An analytical version of our model
    typedef typename Traits::ModelType                            ModelType;

    typedef typename Traits::GridPartType                         GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType            DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType                 DiscreteFunctionType;

    // The DG space operator
    typedef typename Traits::OperatorType                         OperatorType;

    // type of steady state solver
    typedef typename Traits::SolverType                           SolverType;

    // type of analytical traits
    typedef typename Traits::AnalyticalTraits                     AnalyticalTraits;

    // type of discrete traits
    typedef typename Traits::DiscreteTraits                       DiscreteTraits;

    typedef uint64_t                                              UInt64Type ;

    typedef DiscreteFunctionType                                  CheckPointDiscreteFunctionType;
    typedef DiscreteFunctionType                                  LimitDiscreteFunctionType;
    typedef DiscreteFunctionType                                  AdaptationDiscreteFunctionType;

    typedef typename Traits::AdaptIndicatorType                   AdaptIndicatorType;
    typedef typename Traits::AdditionalOutputHandlerType          AdditionalOutputHandlerType;
    typedef typename Traits::SolverMonitorHandlerType             SolverMonitorHandlerType;
    typedef typename Traits::DiagnosticsHandlerType               DiagnosticsHandlerType;

    typedef typename Traits::TimeProviderType                     TimeProviderType;

    SubSteadyState2EvolutionAlgorithm ( GridType &grid, const std::string name = ""  )
      : BaseType( grid, name )
    {}

    virtual const std::string name () { return BaseType::name(); }

    GridType& grid () const { return BaseType::grid(); }

    //! return reference to discrete space
    DiscreteFunctionSpaceType& space () { return BaseType::space(); }
    const DiscreteFunctionSpaceType& space () const { return BaseType::space(); }

    virtual DiscreteFunctionType& solution () { return BaseType::solution(); }

    virtual DiscreteFunctionType& rhs() { return BaseType::rhs(); };

    std::string description () const { return BaseType::description(); }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const
    {
      return BaseType::dataPrefix();
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return BaseType::gridWidth(); }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      return BaseType::gridSize();
    }

    virtual bool checkSolutionValid ( const int loop, TimeProviderType& tp ) const { return BaseType::checkSolutionValid( loop ); }

    virtual void initialize ( const int loop, TimeProviderType& tp )
    {
      BaseType::initialize( loop );
    }

    virtual void preSolve( const int loop, TimeProviderType& tp  )
    {
      BaseType::preSolve( loop );
    }

    virtual void solve ( const int loop, TimeProviderType& tp  )
    {
      BaseType::solve( loop );
    }

    virtual void postSolve( const int loop, TimeProviderType& tp  )
    {
      BaseType::postSolve( loop );
    }

    void finalize ( const int loop, TimeProviderType& tp  )
    {
      BaseType::finalize( loop );
    }

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType& monitor() { return BaseType::monitor(); }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType& diagnostics() { return BaseType::diagnostics(); }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType* additionalOutput() { return nullptr; }

    //LIMITING
    virtual void limit(){}
    virtual LimitDiscreteFunctionType* limitSolution () { return nullptr; }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() { return nullptr; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return nullptr; }

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () { return &BaseType::solution(); }

    //DATAWRITING
    virtual IOTupleType dataTuple () { return IOTupleType(); }

    const ProblemType &problem () const
    {
      return BaseType::problem();
    }

    ProblemType &problem ()
    {
      return BaseType::problem();
    }

    ModelType &model () { return BaseType::model(); }
    const ModelType &model () const { return BaseType::model(); }

  };



}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
