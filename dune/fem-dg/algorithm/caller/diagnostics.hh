#ifndef FEMDG_ALGORITHM_DIAGNOSTICSCALLER_HH
#define FEMDG_ALGORITHM_DIAGNOSTICSCALLER_HH

#include <dune/fem-dg/misc/diagnostics.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{
  /**
   * \brief Caller class managing the data writing of diagnostic data.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename std::make_index_sequence< std::tuple_size< AlgTupleImp >::value > >
  class DiagnosticsCaller;

 /**
   * \brief Specialization of a caller class managing the writing of diagnostic data.
   *
   * \ingroup Callers
   *
   * This class manages the writing of diagnostic data for a tuple of sub-algorithms.
   * For each sub-algorithm diagnostic writing can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef DataWriterCaller< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                            std::index_sequence< 0, 2 > >
   *                                           MyCaller;
   * \endcode
   * This would enable data writing for `Alg1` and `Alg3`;
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam std::index_sequence< Ints... > Index sequence for enabling the data writing feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class DiagnosticsCaller< AlgTupleImp, std::index_sequence< Ints... > >
    : public CallerInterface
  {
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef typename std::tuple_element_t< 0, TupleType >::element_type::GridType    GridType;

    template< class Caller >
    class LoopCallee
    {
      template<class C, class T, class... Args >
      static typename std::enable_if< std::is_void< typename T::element_type::DiagnosticsType >::value >::type
      getDiagnostics( T&, Args&& ... ){}
      template<class C, class T, class... Args >
      static typename std::enable_if< !std::is_void< typename T::element_type::DiagnosticsType >::value >::type
      getDiagnostics( T& elem, Args &&... a )
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
    explicit DiagnosticsCaller( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) )
    {}

    /**
     * \brief Finalizes the writing of diagnostic data.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void finalizeStart( AlgImp* alg, int loop, TimeProviderImp& tp )
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
    template< class AlgImp, class TimeProviderImp >
    void postSolveEnd( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Write >::apply( tuple_, tp );
    }

  private:
    TupleType tuple_;
  };

  /**
   * \brief Specialization of a caller class doing no diagnostic writing.
   *
   * \ingroup Callers
   */
  template< class TupleImp >
  class DiagnosticsCaller< TupleImp, std::index_sequence<> >
    : public CallerInterface
  {
  public:
    template< class ... Args >
    DiagnosticsCaller ( Args && ... ) {}
  };

}
}

#endif
