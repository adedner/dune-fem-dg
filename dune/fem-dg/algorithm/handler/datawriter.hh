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


namespace Dune
{
namespace Fem
{

  template< class... StepperArg >
  class CombinedDefaultDataWriterHandler
  {
  public:
    typedef std::tuple< typename std::add_pointer< StepperArg >::type... >                               StepperTupleType;
    typedef typename std::remove_pointer< typename std::tuple_element<0,StepperTupleType>::type >::type  FirstStepperType;
    typedef typename FirstStepperType::GridType                                                          GridType;
   // typedef typename tuple_concat< typename StepperArg::IOTupleType... >::type                           IOTupleType;
    typedef typename FirstStepperType::IOTupleType                                                       IOTupleType;
    typedef DataWriter< GridType, IOTupleType >                                                          DataWriterType;


    static_assert( Std::are_all_same< typename StepperArg::GridType... >::value,
                   "CombinedDefaultDataWriterHandler: GridType has to be equal for all steppers" );
    static_assert( sizeof ... ( StepperArg ) > 0,
                   "CombinedDefaultDataWriterHandler: called with zero steppers" );


    CombinedDefaultDataWriterHandler( const StepperTupleType& tuple )
      : tuple_( tuple ),
        dataTuple_(),
        dataWriter_()
    {
      setDataTuple( tuple, Std::index_sequence_for< StepperArg ... >() );
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
        dataWriter_->write( tp );
      }
    }

    template< class TimeProviderImp >
    void finalize( TimeProviderImp& tp )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        dataWriter_->write( tp );
      }
    }

    IOTupleType& dataTuple()
    {
      return dataTuple_;
    }

    private:
    template< std::size_t ... i >
    void setDataTuple ( const StepperTupleType &tuple, Std::index_sequence< i ... > )
    {
      //TODO check tuple elements pointing to null pointer
      dataTuple_ = std::tuple_cat( *(std::get< i >( tuple )->dataTuple())... );
    }

    //template< class Element >
    //Element& returnDataTuple( Element&& elem )
    //{
    //  if( elem )
    //    return elem;
    //  return std::make_tuple( decltype( elem ) );
    //}

    const StepperTupleType&           tuple_;
    IOTupleType                       dataTuple_;
    std::unique_ptr< DataWriterType > dataWriter_;

  };


  class NoDataWriterHandler
  {
    public:

    template< class ... Args >
    NoDataWriterHandler( Args&& ... )
    {}

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
