#ifndef FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH
#define FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH

#include <string>
#include <dune/fem-dg/misc/tupleutility.hh>

namespace Dune
{
namespace Fem
{

  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class SolutionLimiterHandler;


  template< class AlgTupleImp, std::size_t... Ints >
  class SolutionLimiterHandler< AlgTupleImp, Std::index_sequence< Ints... > >
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

    SolutionLimiterHandler( const AlgTupleType& tuple )
      : tuple_( tuple )
    {}

    void step()
    {
      ForLoopType< LimitSolution >::apply( tuple_ );
    }

  private:
    TupleType tuple_;
  };


  template< class TupleImp >
  class SolutionLimiterHandler< TupleImp, Std::index_sequence<> >
  {
  public:
    template< class ... Args >
    SolutionLimiterHandler ( Args && ... ) {}

    void step () {}
  };

}
}

#endif
