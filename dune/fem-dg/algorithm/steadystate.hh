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
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
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
    typedef typename std::tuple_element<0, std::tuple< ProblemTraits... > >::type::GridType  GridType;

    // wrap operator
    typedef GridTimeProvider< GridType >                                   TimeProviderType;

    //typedef ...
    typedef std::tuple< typename std::add_pointer< typename ProblemTraits::template Stepper<polOrder>::Type >::type... > StepperTupleType;
    typedef typename Std::make_index_sequence_impl< std::tuple_size< StepperTupleType >::value >::type                   IndexSequenceType;

    typedef Dune::Fem::SolverMonitorHandler< StepperTupleType >   SolverMonitorHandlerType;
    typedef Dune::Fem::DataWriterHandler< StepperTupleType >      DataWriterHandlerType;

    typedef typename DataWriterHandlerType::IOTupleType                                                                      IOTupleType;
  };



  template< int polOrder, class... ProblemTraits >
  class SteadyStateAlgorithm
    : public AlgorithmInterface< SteadyStateTraits< polOrder, ProblemTraits... > >
  {

    typedef SteadyStateTraits< polOrder, ProblemTraits... >             Traits;
    typedef AlgorithmInterface< Traits >                                BaseType;
  public:
    typedef typename BaseType::GridType                                 GridType;
    typedef typename BaseType::IOTupleType                              IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType                 SolverMonitorHandlerType;
    typedef typename Traits::DataWriterHandlerType                      DataWriterHandlerType;

    typedef uint64_t                                                    UInt64Type ;

    typedef std::tuple< typename std::add_pointer< typename ProblemTraits::template Stepper<polOrder>::Type >::type... > StepperTupleType;

    template< int i >
    struct Constructor
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        typedef typename std::remove_pointer< typename std::tuple_element< i, Tuple >::type >::type Element;
        std::get< i >( tuple ) = new Element( std::forward< Args >( args ) ... );
      }
    };
    template< int i >
    struct Initialize
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->initialize( args... );
      }
    };
    template< int i >
    struct PreSolve
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->preSolve( args... );
      }
    };
    template< int i >
    struct Solve
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->solve( args... );
      }
    };
    template< int i >
    struct PostSolve
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->postSolve( args... );
      }
    };

    template< int i >
    struct Finalize
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->finalize( args... );
      }
    };
    template< int i >
    struct GridWidth
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, double& res, Args && ... args )
      {
        res = std::max( res, std::get< i >( tuple )->gridWidth( args... ) );
      }
    };
    template< int i >
    struct GridSize
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, UInt64Type& res, Args && ... args )
      {
        res = std::max( res, std::get< i >( tuple )->gridSize( args... ) );
      }
    };
    template< int i >
    struct CheckSolutionValid
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, bool& res, Args && ... args )
      {
        res &= std::get< i >( tuple )->checkSolutionValid( args... );
      }
    };

    using BaseType::grid;

    SteadyStateAlgorithm ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  ),
      tuple_( createStepper( grid, name ) ),
      solverMonitorHandler_( tuple_ ),
      dataWriterHandler_( tuple_ )
    {}

    // create Tuple of contained subspaces
    static StepperTupleType createStepper( GridType &grid, const std::string name = "" )
    {
      StepperTupleType tuple;
      ForLoop< Constructor, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple, grid, name );
      return tuple;
    }

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
      ForLoop< GridWidth, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple_, res );
      return res;
    }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type res=0;
      ForLoop< GridSize, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple_, res );
      return res;
    }

    virtual void initialize ( const int loop )
    {
      ForLoop< Initialize, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple_, loop );
    }

    virtual void preSolve( const int loop )
    {
      ForLoop< PreSolve, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple_, loop );
    }

    virtual void solve ( const int loop )
    {
      initialize( loop );
      preSolve( loop );
      ForLoop< Solve, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple_, loop );
      postSolve( loop );
      finalize( loop );
    }

    virtual void postSolve( const int loop )
    {
      ForLoop< PostSolve, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple_, loop );
    }

    void finalize ( const int loop )
    {
      ForLoop< Finalize, 0, sizeof ... ( ProblemTraits )-1 >::apply( tuple_, loop );
    }

  protected:

    StepperTupleType               tuple_;
    SolverMonitorHandlerType       solverMonitorHandler_;
    DataWriterHandlerType          dataWriterHandler_;
  };

}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
