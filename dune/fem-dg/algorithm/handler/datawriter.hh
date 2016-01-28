#ifndef FEMDG_DATAWRITERHANDLER_HH
#define FEMDG_DATAWRITERHANDLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class DataWriterHandler;

  template< class AlgTupleImp, std::size_t... Ints >
  class DataWriterHandler< AlgTupleImp, Std::index_sequence< Ints... > >
    : public HandlerInterface
  {
    template< class TupleType > struct IOTupleExtractor;
    template< class ... Args > struct IOTupleExtractor< std::tuple< Args... > >
    { typedef typename tuple_concat< typename std::remove_pointer< Args >::type::IOTupleType::type... >::type type; };

    typedef AlgTupleImp                                                            AlgTupleType;

    typedef Std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef typename std::remove_pointer< typename std::tuple_element< 0, TupleType >::type >::type::GridType
                                                                                    GridType;

  public:
    typedef typename IOTupleExtractor< TupleType >::type                            IOTupleType;

    typedef DataWriter< GridType, IOTupleType >                                     DataWriterType;

    template <class Grid, class DataTuple>
    struct DataOutput
    {
      typedef Dune::Fem::DataWriter< Grid, DataTuple > Type;
    };


  private:
    template< int i >
    struct AdditionalOutput
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdditionalOutputType >::value >::type
      additionalOutput( T, Args&& ... ){}
      template<class T, class TimeProviderImp, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdditionalOutputType >::value >::type
      additionalOutput( T elem, TimeProviderImp& tp, Args && ... args )
      {
        if( elem->additionalOutput() )
          elem->additionalOutput()->step( tp, *elem, args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        additionalOutput( std::get<i>( tuple ), std::forward<Args>(args)... );
      }
    };

    template< template< int > class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;

  public:

    DataWriterHandler( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) ),
        dataWriter_(),
        dataTuple_( dataTuple( tuple_, IndexSequenceType() ) )
    {
    }

    template< class SubAlgImp, class TimeProviderImp >
    void initialize_post( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      dataWriter_.reset( new DataWriterType( std::get<0>(tuple_)->solution().space().grid(), dataTuple_, tp,
                                             alg->eocParams().dataOutputParameters( loop, alg->dataPrefix() ) ) );
    }

    template< class SubAlgImp, class TimeProviderImp >
    void preSolve_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        //update all additional Output
        ForLoopType< AdditionalOutput >::apply( tuple_, tp );

        //writeData
        dataWriter_->write( tp );
      }
    }

    template< class SubAlgImp, class TimeProviderImp >
    void finalize_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        //update all additional Output
        ForLoopType< AdditionalOutput >::apply( tuple_, tp );
        //writeData
        dataWriter_->write( tp );
      }
    }

    IOTupleType dataTuple()
    {
      return dataTuple_;
    }

  private:
    template< std::size_t ... i >
    IOTupleType dataTuple ( const TupleType &tuple, Std::index_sequence< i ... > )
    {
      return std::tuple_cat( (*std::get< i >( tuple )->dataTuple() )... );
    }

    TupleType                         tuple_;
    std::unique_ptr< DataWriterType > dataWriter_;
    IOTupleType                       dataTuple_;
  };


  template< class AlgTupleImp >
  class DataWriterHandler< AlgTupleImp, Std::index_sequence<> >
    : public HandlerInterface
  {
  public:
    template <class Grid, class DataTuple>
    struct DataOutput
    {
      typedef Dune::Fem::DataOutput< Grid, DataTuple > Type;
    };

    template< class ... Args >
    DataWriterHandler( Args&& ... ) {}

  };

}
}
#endif
