#ifndef FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH
#define FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH

#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/common/forloop.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Caller class managing the post processing of data.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class PostProcessingCaller;


  template< class AlgTupleImp, std::size_t... Ints >
  class PostProcessingCaller< AlgTupleImp, Std::index_sequence< Ints... > >
    : public CallerInterface
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

  private:
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

  public:

    /**
     * \brief Constructor
     *
     * \param[in] tuple Tuple of pointer to sub-algorithm.
     */
    PostProcessingCaller( const AlgTupleType& tuple )
      : tuple_( tuple )
    {}

    /**
     * \brief Apply post processing to the solution.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void solveEnd( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< LimitSolution >::apply( tuple_ );
    }

  private:
    TupleType tuple_;
  };

 /**
   * \brief Caller class doing no post processing.
   *
   * \ingroup Callers
   */
  template< class TupleImp >
  class PostProcessingCaller< TupleImp, Std::index_sequence<> >
    : public CallerInterface
  {
  public:
    template< class ... Args >
    PostProcessingCaller ( Args && ... ) {}

  };

}
}

#endif