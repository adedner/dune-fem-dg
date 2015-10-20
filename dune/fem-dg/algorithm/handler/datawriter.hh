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

//#include <dune/common/classname.hh>
namespace Dune
{
namespace Fem
{

  template< class... StepperArg >
  class DataWriterHandler;

  template< class StepperHead, class ... StepperArg >
  class DataWriterHandler< StepperHead, StepperArg ... >
  {
  public:
    typedef std::tuple< typename std::add_pointer< StepperHead >::type, typename std::add_pointer< StepperArg >::type... > StepperTupleType;
    typedef typename std::remove_pointer< typename std::tuple_element<0,StepperTupleType>::type >::type  FirstStepperType;
    typedef typename StepperHead::GridType                                                               GridType;
    typedef typename tuple_concat< typename StepperHead::IOTupleType::type, typename StepperArg::IOTupleType::type... >::type IOTupleType;

    typedef DataWriter< GridType, IOTupleType >                                                          DataWriterType;

    static_assert( Std::are_all_same< GridType, typename StepperArg::GridType... >::value,
                   "DataWriterHandler: GridType has to be equal for all steppers" );

  private:
    template< int i >
    struct AdditionalOutput
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdditionalOutputHandlerType >::value >::type
      additionalOutput( T, Args&& ... ){}
      template<class T, class TimeProviderImp, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdditionalOutputHandlerType >::value >::type
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

//    struct DataTupleElement
//    {
//      template< class Tuple, class ... Args >
//      static void apply ( Tuple &tuple, Args && ... args )
//      {
//        auto data = std::get<i>( tuple );
//        //std::cout << i << ":" << className( data ) << std::end;
//        std::cout << data << ", ";
//      }
//    };
//
//    template< int i >
//    struct DataTuple
//    {
//      template< class Tuple, class ... Args >
//      static void apply ( Tuple &tuple, Args && ... args )
//      {
//        auto data = std::get<i>( tuple )->dataTuple();
//        //std::cout << "IOTuple" << i << ":" << className( data ) << std::endl << std::endl;
//        std::cout << "|";
//        ForLoop< DataTupleElement, 0, std::tuple_size<decltype( data )>::value -1 >::apply( data );
//      }
//    };
//
//    template< int i >
//    struct DataTupleConverted
//    {
//      template< class Tuple, class ... Args >
//      static void apply ( Tuple &tuple, Args && ... args )
//      {
//        auto data = std::get<i>( tuple );
//        std::cout << data << ", ";
//      }
//    };


  public:


    DataWriterHandler( const StepperTupleType& tuple )
      : tuple_( tuple ),
        dataWriter_(),
        dataTuple_( dataTuple( tuple, Std::index_sequence_for< StepperHead, StepperArg ... >() ) )
    {
      //std::cout << "============================" <<std::endl << "DataIOTuple";
      //ForLoop< DataTupleConverted, 0, tuple_size< decltype( dataTuple_ ) >::value - 1 >::apply( dataTuple_ );
      //std::cout << std::endl;
    }

    template< class TimeProviderImp, class ParameterType >
    void init( TimeProviderImp& tp, ParameterType& param  )
    {
      dataWriter_.reset( new DataWriterType( std::get<0>(tuple_)->solution().space().grid(), dataTuple_, tp, param ) );
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        //update all additional Output
        ForLoop< AdditionalOutput, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );

        //std::cout << "============================" <<std::endl;
        //ForLoop< DataTuple, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
        //std::cout << std::endl;
        //writeData
        dataWriter_->write( tp );
      }
    }

    template< class TimeProviderImp >
    void finalize( TimeProviderImp& tp )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        //update all additional Output
        ForLoop< AdditionalOutput, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
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
    IOTupleType dataTuple ( const StepperTupleType &tuple, Std::index_sequence< i ... > )
    {
      return std::tuple_cat( (*std::get< i >( tuple )->dataTuple() )... );
    }


    const StepperTupleType&           tuple_;
    std::unique_ptr< DataWriterType > dataWriter_;
    IOTupleType                       dataTuple_;
  };

  template<>
  class DataWriterHandler<>
  {
    public:

    template< class ... Args >
    DataWriterHandler( Args&& ... ) {}

    template< class ... Args >
    void init( Args&& ...  ) {}

    template< class ... Args >
    void step( Args&& ...  ) {}

    template< class ... Args >
    void finalize( Args&& ...  ) {}
  };

}
}
#endif
