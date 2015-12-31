#ifndef DUNE_FEMDG_ALGORITHM_COMBINEDSTEADYSTATEALGORITHM_HH
#define DUNE_FEMDG_ALGORITHM_COMBINEDSTEADYSTATEALGORITHM_HH

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem-dg/algorithm/base.hh>

#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>

#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>

namespace Dune
{
namespace Fem
{

  template< int polOrder, class ... ProblemTraits >
  struct SteadyStateTraits
  {
    // type of Grid
    typedef typename std::tuple_element<0, std::tuple< ProblemTraits... > >::type::GridType
                                                                           GridType;

    // wrap operator
    typedef GridTimeProvider< GridType >                                   TimeProviderType;

    //typedef ...
    typedef std::tuple< typename std::add_pointer< typename ProblemTraits::template Stepper<polOrder>::Type >::type... >
                                                                           StepperTupleType;

    //typedef typename Std::make_index_sequence_impl< std::tuple_size< StepperTupleType >::value >::type
    //                                                                     IndexSequenceType;

    typedef Dune::Fem::SolverMonitorHandler< StepperTupleType >            SolverMonitorHandlerType;
    typedef Dune::Fem::DataWriterHandler< StepperTupleType >               DataWriterHandlerType;

    typedef typename DataWriterHandlerType::IOTupleType                    IOTupleType;
  };

  //template< int polOrder, class... ProblemTraits >
  //class SteadyStateAlgorithm;




  template< int polOrder/*, class ProblemTraitsHead*/, class... ProblemTraits >
  class SteadyStateAlgorithm/*< polOrder, ProblemTraitsHead, ProblemTraits... >*/
    : public AlgorithmInterface< SteadyStateTraits< polOrder/*, ProblemTraitsHead*/, ProblemTraits... > >
  {

    typedef SteadyStateTraits< polOrder/*, ProblemTraitsHead*/, ProblemTraits... > Traits;
    typedef AlgorithmInterface< Traits >                                BaseType;
  public:
    typedef typename BaseType::GridType                                 GridType;
    typedef typename BaseType::IOTupleType                              IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType                 SolverMonitorHandlerType;
    typedef typename Traits::DataWriterHandlerType                      DataWriterHandlerType;

    typedef uint64_t                                                    UInt64Type ;

    typedef std::tuple< /*typename std::add_pointer< typename ProblemTraitsHead::template Stepper<polOrder>::Type >::type,*/
                        typename std::add_pointer< typename ProblemTraits::template Stepper<polOrder>::Type >::type... > StepperTupleType;

  private:
    struct Initialize {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->initialize( std::forward<Args>(a)... ); }
    };
    struct PreSolve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->preSolve( std::forward<Args>(a)... ); }
    };
    struct Solve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->solve( std::forward<Args>(a)... ); }
    };
    struct PostSolve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->postSolve( std::forward<Args>(a)... ); }
    };
    struct Finalize {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->finalize( std::forward<Args>(a)... ); }
    };
    struct GridWidth {
      template< class T, class ... Args > static void apply ( T& e, double& res, Args && ... a )
      { res = std::max( res, e->gridWidth(std::forward<Args>(a)... ) ); }
    };
    struct GridSize {
      template<class T, class... Args > static void apply ( T& e, UInt64Type& res, Args && ... a )
      { res = std::max( res, e->gridSize(std::forward<Args>(a)... ) ); }
    };
    struct CheckSolutionValid {
      template<class T, class... Args > static void apply( T e, bool& res, Args&& ... a )
      { res &= e->checkSolutionValid( std::forward<Args>(a)... ); }
    };
    struct Init {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->init(); }
    };
    template< class Caller >
    class LoopCallee
    {
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... Args >
        static void apply ( Tuple &tuple, Args&& ... a )
        {
          Caller::apply( std::get<i>( tuple ), std::forward<Args>(a)... );
        }
      };
    };

    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0,  sizeof ... ( ProblemTraits )-1 >;

    template< std::size_t ...i >
    static StepperTupleType createStepper ( Std::index_sequence< i... >, GridType &grid )
    {
      static auto tuple = std::make_tuple( new typename std::remove_pointer< typename std::tuple_element< i, StepperTupleType >::type >::type( grid ) ... );
      ForLoopType< Init >::apply( tuple );
      return tuple;
    }

    // create Tuple of contained sub algorithms
    static StepperTupleType createStepper( GridType &grid )
    {
      return createStepper( typename Std::make_index_sequence_impl< std::tuple_size< StepperTupleType >::value >::type(), grid );
    }

  public:
    using BaseType::grid;

    SteadyStateAlgorithm ( GridType &grid )
    : BaseType( grid ),
      tuple_( createStepper( grid ) ),
      solverMonitorHandler_( tuple_ ),
      dataWriterHandler_( tuple_ )
    {}


    virtual IOTupleType dataTuple ()
    {
      return dataWriterHandler_.dataTuple();
    }

    virtual SolverMonitorHandlerType& monitor()
    {
      return solverMonitorHandler_;
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const
    {
      double res=0.0;
      ForLoopType< GridWidth >::apply( tuple_, res );
      return res;
    }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type res=0;
      ForLoopType< GridSize >::apply( tuple_, res );
      return res;
    }

    virtual void initialize ( const int loop )
    {
      ForLoopType< Initialize >::apply( tuple_, loop );
    }

    virtual void preSolve( const int loop )
    {
      ForLoopType< PreSolve >::apply( tuple_, loop );
    }

    virtual void solve ( const int loop )
    {
      initialize( loop );
      preSolve( loop );
      ForLoopType< Solve >::apply( tuple_, loop );
      postSolve( loop );
      finalize( loop );
    }

    virtual void postSolve( const int loop )
    {
      ForLoopType< PostSolve >::apply( tuple_, loop );
    }

    void finalize ( const int loop )
    {
      ForLoopType< Finalize >::apply( tuple_, loop );
    }

  protected:

    StepperTupleType               tuple_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    DataWriterHandlerType          dataWriterHandler_;
  };


#if 0
  //specialization
  template< int polOrder, class ProblemTraits >
  class SteadyStateAlgorithm< polOrder, ProblemTraits >
    : public AlgorithmInterface< SteadyStateTraits< polOrder, ProblemTraits > >
  {

    typedef SteadyStateTraits< polOrder, ProblemTraits >                Traits;
    typedef AlgorithmInterface< Traits >                                BaseType;
  public:
    typedef typename BaseType::GridType                                 GridType;
    typedef typename BaseType::IOTupleType                              IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType                 SolverMonitorHandlerType;
    typedef typename Traits::DataWriterHandlerType                      DataWriterHandlerType;

    typedef uint64_t                                                    UInt64Type ;

    typedef typename ProblemTraits::template Stepper<polOrder>::Type    AlgType;

    using BaseType::grid;

    SteadyStateAlgorithm ( GridType &grid )
    : BaseType( grid  ),
      alg_( grid ),
      solverMonitorHandler_( std::make_tuple( &alg_ ) ),
      dataWriterHandler_( std::make_tuple( &alg_ ) )
    {}

    virtual IOTupleType dataTuple ()
    {
      return dataWriterHandler_.dataTuple();
    }

    virtual SolverMonitorHandlerType& monitor()
    {
      return solverMonitorHandler_;
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const
    {
      return alg_.gridWidth();
    }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type res = alg_.gridSize();
      return res;
    }

    virtual void initialize ( const int loop )
    {
      alg_.initialize( loop );
    }

    virtual void preSolve( const int loop )
    {
      alg_.preSolve( loop );
    }

    virtual void solve ( const int loop )
    {
      initialize( loop );
      preSolve( loop );
      alg_.solve( loop );
      postSolve( loop );
      finalize( loop );
    }

    virtual void postSolve( const int loop )
    {
      alg_.solve( loop );
    }

    void finalize ( const int loop )
    {
      alg_.finalize( loop );
    }

  protected:

    AlgType                        alg_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    DataWriterHandlerType          dataWriterHandler_;
  };
#endif

}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
