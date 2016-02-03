#ifndef FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH
#define FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH

#include <dune/fem-dg/misc/diagnostics.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{
  /**
   * \brief Handler class managing the data writing of diagnostic data.
   *
   * \ingroup Handlers
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class DiagnosticsHandler;

 /**
   * \brief Specialization of a handler class managing the writing of diagnostic data.
   *
   * \ingroup Handlers
   *
   * This class manages the writing of diagnostic data for a tuple of sub-algorithms.
   * For each sub-algorithm diagnostic writing can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef DataWriterHandler< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                            Std::index_sequence< 0, 2 > >
   *                                           MyHandler;
   * \endcode
   * This would enable data writing for `Alg1` and `Alg3`;
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam Std::index_sequence< Ints... > Index sequence for enabling the data writing feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class DiagnosticsHandler< AlgTupleImp, Std::index_sequence< Ints... > >
    : public HandlerInterface
  {
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef Std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef typename std::remove_pointer< typename std::tuple_element< 0, TupleType >::type >::type::GridType
                                                                                    GridType;

    template< class Caller >
    class LoopCallee
    {
      template<class C, class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::DiagnosticsType >::value >::type
      getDiagnostics( T, Args&& ... ){}
      template<class C, class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::DiagnosticsType >::value >::type
      getDiagnostics( T elem, Args &&... a )
      {
        if( elem->diagnostics() )
          C::apply( elem->diagnostics(), std::forward<Args>(a)... );
      }
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... Args >
        static void apply ( Tuple &tuple, Args&& ... a )
        {
          getDiagnostics< Caller >( std::get<i>( tuple ), std::forward<Args>(a)... );
        }
      };
    };

    struct Write {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->step( std::forward<Args>(a)... ); }
    };
    struct Finalize {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->finalize( std::forward<Args>(a)... ); }
    };

    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0, numAlgs - 1 >;

  public:

    /**
     * \brief Constructor.
     *
     * \param[in] tuple Tuple of all sub-algorithms.
     */
    DiagnosticsHandler( const AlgTupleType& tuple )
      : tuple_( tuple )
    {}

    /**
     * \brief Finalizes the writing of diagnostic data.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void finalize_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Finalize >::apply( tuple_ );
    }

    /**
     * \brief Writes diagnostic data.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void postSolve_post( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Write >::apply( tuple_, tp );
    }

  private:
    TupleType tuple_;
  };

  /**
   * \brief Specialization of a handler class doing no diagnostic writing.
   *
   * \ingroup Handlers
   */
  template< class TupleImp >
  class DiagnosticsHandler< TupleImp, Std::index_sequence<> >
    : public HandlerInterface
  {
  public:
    template< class ... Args >
    DiagnosticsHandler ( Args && ... ) {}
  };

}
}

#endif
