#ifndef FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH
#define FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH

#include <string>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class SolverMonitorHandler;


  /**
   * \brief Handler class to collect all data from solver monitor from sub algorithms.
   *
   * \ingroup Handler
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class SolverMonitorHandler< AlgTupleImp, Std::index_sequence< Ints... > >
    : public HandlerInterface
  {
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef Std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    enum CombinationType { max, min, sum, avg };

    template< class Caller >
    class LoopCallee
    {
      template<class C, class T, class... A >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::SolverMonitorHandlerType >::value >::type
      getMonitor( T, A&& ... ){}
      template<class C, class T, class... A >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::SolverMonitorHandlerType >::value >::type
      getMonitor( T elem, A &&... a )
      {
        if( elem->monitor() )
          C::apply(elem->monitor(), std::forward<A>(a)... );
      }
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... A >
        static void apply ( Tuple &tuple, A&& ... a )
        { getMonitor< Caller >( std::get<i>( tuple ), std::forward<A>(a)... ); }
      };
    };

    struct Step {
      template<class T, class... A > static void apply( T e, A&& ... a )
      { e->step( std::forward<A>(a)... ); }
    };

    struct Finalize {
      template<class T, class... A > static void apply( T e, A&& ... a )
      { e->finalize( std::forward<A>(a)... ); }
    };

    struct GetData {
      template<class T, class... A > static void apply( T e, double& res, const std::string name, CombinationType comb, A&& ... a )
      {
        switch( comb )
        {
          case CombinationType::max:
            res = std::max( res, e->getData( name ) );
          case CombinationType::min:
            res = std::min( res, e->getData( name ) );
          case CombinationType::sum:
            res += e->getData( name );
          case CombinationType::avg:
            res += e->getData( name );
        }
      }
    };

    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0, numAlgs - 1 >;

  public:

    SolverMonitorHandler( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) )
    {}

    template< class SubAlgImp, class TimeProviderImp >
    void postSolve_post( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Step >::apply( tuple_, tp );
    }

    template< class SubAlgImp, class TimeProviderImp >
    void finalize_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Finalize >::apply( tuple_, alg->gridWidth(), alg->gridSize() );
    }

    template< class... StringType >
    void print( std::string str, StringType& ...tail )
    {
      std::cout << str << ":  " << getData( str ) << ", ";
      print( tail... );
    }
    void print() {}

    const double getData( const std::string name, CombinationType comb = CombinationType::max ) const
    {
      double res = (comb == CombinationType::max) ? std::numeric_limits<double>::max() : 0.0;
      ForLoopType< GetData >::apply( tuple_, res, name, comb );
      return res;
    }

  private:
    TupleType  tuple_;
  };

  template< class AlgTupleImp >
  class SolverMonitorHandler< AlgTupleImp, Std::index_sequence<> >
    : public HandlerInterface
  {
  public:
    template <class ... Args>
    SolverMonitorHandler(Args&& ... )
    {}

    template <class ... Args>
    void print(Args&& ... ) const {}

    template <class ... Args>
    const double getData( Args&& ... ) const
    { return 0.0; }

  };
}
}

#endif
