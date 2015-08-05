#ifndef FEMDG_DATAWRITERHANDLER_HH
#define FEMDG_DATAWRITERHANDLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>

namespace Dune
{
namespace Fem
{


  template< class GridImp, class IOTupleImp >
  class DefaultDataWriterHandler
  {

  public:
    typedef DataWriter< GridImp, IOTupleImp >       DataWriterType;



    DefaultDataWriterHandler( const GridImp& grid, const std::string keyPrefix = "" )
      : dataTuple_(),
        dataWriter_(),
        grid_( grid ),
        keyPrefix_( keyPrefix )
    {}

    template< class TimeProviderImp, class DataTupleImp, class ParameterType >
    void init( TimeProviderImp& tp, DataTupleImp dataTup, ParameterType& param  )
    {
      dataTuple_ = dataTup;
      dataWriter_.reset( new DataWriterType( grid_, dataTuple_, tp, param ) );
    }


    //template< class InDiscreteFunction, class OutDiscreteFunction, class ... Args >
    //void preWriteData( const InDiscreteFunction& in, OutDiscreteFunction* out, Args& ... args )
    //{
    //  const bool reallyWrite = writeAnyway ? true : true /*dataWriter_->willWrite( tp )*/;
    //  if( reallyWrite && out )
    //  {
    //    DiscreteFunctionCalculatorImp::setup( in, *out, args... );
    //  }
    //}

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


    private:
    IOTupleImp             dataTuple_;
    std::unique_ptr< DataWriterType > dataWriter_;
    const GridImp&        grid_;
    const std::string keyPrefix_;

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
    void writeStep( Args&& ...  ) {}

    template< class ... Args >
    void finalize( Args&& ...  ) {}
  };


}
}
#endif
