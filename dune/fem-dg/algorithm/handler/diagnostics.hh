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

  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class DiagnosticsHandler;


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
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::DiagnosticsHandlerType >::value >::type
      getDiagnostics( T, Args&& ... ){}
      template<class C, class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::DiagnosticsHandlerType >::value >::type
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

    DiagnosticsHandler( const AlgTupleType& tuple )
      : tuple_( tuple )
    {}

    template< class SubAlgImp, class TimeProviderImp >
    void postSolve_post( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Write >::apply( tuple_, tp );
    }

    template< class SubAlgImp, class TimeProviderImp >
    void finalize_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Finalize >::apply( tuple_ );
    }

  private:
    TupleType tuple_;
  };

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
