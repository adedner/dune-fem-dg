#ifndef FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH
#define FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH

#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/common/forloop.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class PostProcessingHandler;


  template< class AlgTupleImp, std::size_t... Ints >
  class PostProcessingHandler< AlgTupleImp, Std::index_sequence< Ints... > >
    : public HandlerInterface
  {
  public:
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef Std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef typename std::remove_pointer< typename std::tuple_element< 0, TupleType >::type >::type::GridType
                                                                                    GridType;

    template< template<int> class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;



    template< int i >
    struct LimitSolution
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args&& ... a )
      {
        std::get< i >( tuple )->limit( std::forward<Args>(a)... );
      }
    };

    PostProcessingHandler( const AlgTupleType& tuple )
      : tuple_( tuple )
    {}

    template< class SubAlgImp, class TimeProviderImp >
    void solve_post( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< LimitSolution >::apply( tuple_ );
    }

  private:
    TupleType tuple_;
  };


  template< class TupleImp >
  class PostProcessingHandler< TupleImp, Std::index_sequence<> >
    : public HandlerInterface
  {
  public:
    template< class ... Args >
    PostProcessingHandler ( Args && ... ) {}

  };

}
}

#endif
