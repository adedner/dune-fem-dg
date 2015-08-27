#ifndef FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH
#define FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH

#include <dune/fem-dg/misc/diagnostics.hh>

namespace Dune
{
namespace Fem
{

  template< class... StepperArg >
  class DiagnosticsHandler;


  template< class StepperHead, class ... StepperArg >
  class DiagnosticsHandler< StepperHead, StepperArg ... >
  {
    typedef std::tuple< typename std::add_pointer< StepperHead >::type,
      typename std::add_pointer< StepperArg >::type... >                                StepperTupleType;
    typedef typename StepperHead::GridType                                              GridType;

    template< int i >
    struct Write
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->diagnostics().step( args... );
      }
    };
    template< int i >
    struct Finalize
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->diagnostics().finalize( args... );
      }
    };

  public:

    DiagnosticsHandler( const StepperTupleType& tuple )
      : tuple_( tuple )
    {}

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      ForLoop< Write, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
    }

    void finalize() const
    {
      ForLoop< Finalize, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
    }

  private:
    StepperTupleType tuple_;
  };

  template<>
  class DiagnosticsHandler<>
  {
  public:
    template< class ... Args >
    DiagnosticsHandler ( Args && ... ) {}

    template <class ... Args>
    void step( Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };



  template< class... DiagnosticsImp >
  class SubDiagnosticsHandler;

  template< class DiagnosticsImp >
  class SubDiagnosticsHandler< DiagnosticsImp >
  {
  public:
    typedef DiagnosticsImp DiagnosticsType;
    typedef std::map< std::string, long unsigned int* > DataIntType;
    typedef std::map< std::string, double* > DataDoubleType;


    SubDiagnosticsHandler()
      : diagnostics_( true ),
        dataInt_(),
        dataDouble_()
    {}

    void registerData( const std::string name, double* diagnosticsData )
    {
      assert( diagnosticsData );
      dataDouble_.insert( std::make_pair(name, diagnosticsData ) );
    }

    void registerData( const std::string name, long unsigned int* diagnosticsData )
    {
      assert( diagnosticsData );
      dataInt_.insert( std::make_pair(name, diagnosticsData ) );
    }

    const double getData( const std::string name )
    {
      if( dataInt_.find(name) != dataInt_.end() )
      {
        assert( dataInt_[ name ] );
        return (double)*dataInt_[ name ];
      }
      assert( dataDouble_[ name ] );
      return *dataDouble_[ name ];
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      const double ldt = tp.deltaT();
      diagnostics_.write( tp.time() + ldt, ldt, getData( "Elements" ), std::vector<double>() );
    }

    void finalize() const
    {
      diagnostics_.flush();
    }

  private:
    DiagnosticsType diagnostics_;
    DataIntType        dataInt_;
    DataDoubleType     dataDouble_;
  };


  template<>
  class SubDiagnosticsHandler<>
  {
    struct NoDiagnosticsType
    {
      template <class ... Args>
      void step( Args&& ... ) const {};

      template <class ... Args>
      void finalize( Args&& ... ) const {};
    };

  public:
    typedef NoDiagnosticsType     DiagnosticsType;

    template <class ... Args>
    double getData( Args&& ... ) const { return 0.0; };

    template <class ... Args>
    void registerData( Args&& ... ) const {};

    template <class ... Args>
    void step( Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };

}
}

#endif
