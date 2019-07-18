#ifndef FEMDG_SUBALGORITHM_DIAGNOSTICSHANDLER_HH
#define FEMDG_SUBALGORITHM_DIAGNOSTICSHANDLER_HH

#include <dune/fem-dg/misc/diagnostics.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>

namespace Dune
{
namespace Fem
{

  template< class... DiagnosticsImp >
  class SubDiagnostics;


  template<>
  class SubDiagnostics<>
  {
    struct NoDiagnosticsType
    {
      template <class ... Args>
      void step( Args&& ... ) const {}

      template <class ... Args>
      void finalize( Args&& ... ) const {}
    };

  public:
    template <class ... Args>
    SubDiagnostics( Args&& ... )
    {}

    typedef NoDiagnosticsType     DiagnosticsType;

    template <class ... Args>
    double getData( Args&& ... ) const { return 0.0; }

    template <class ... Args>
    void registerData( Args&& ... ) const {}

    template <class ... Args>
    void step( Args&& ... ) const {}

    template <class ... Args>
    void finalize( Args&& ... ) const {}
  };


  template< class DiagnosticsImp >
  class SubDiagnostics< DiagnosticsImp >
  {
  public:
    typedef DiagnosticsImp DiagnosticsType;
    typedef std::map< std::string, long unsigned int* > DataIntType;
    typedef std::map< std::string, double* > DataDoubleType;


    SubDiagnostics( const std::string keyPrefix = "" )
      : keyPrefix_( keyPrefix ),
        diagnostics_( true, keyPrefix ),
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

    double getData( const std::string name )
    {
      if( dataInt_.find(name) != dataInt_.end() )
      {
        assert( dataInt_[ name ] );
        return (double)*dataInt_[ name ];
      }
      if( dataDouble_.find(name) != dataDouble_.end() )
      {
        assert( dataDouble_[ name ] );
        return *dataDouble_[ name ];
      }
      return 0.0;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      const double ldt = tp.timeStepValid() ? tp.deltaT() : 0.0;
      diagnostics_.write( tp.time() + ldt, ldt, getData( "Elements" ), std::vector<double>() );
    }

    void finalize() const
    {
      diagnostics_.flush();
    }

  private:
    std::string keyPrefix_;
    DiagnosticsType diagnostics_;
    DataIntType        dataInt_;
    DataDoubleType     dataDouble_;
  };


  template< class Obj >
  class DiagnosticsOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    DiagnosticsOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

  template<>
  class DiagnosticsOptional< void >
    : public OptionalNullPtr< SubDiagnostics<> >
  {
    typedef OptionalNullPtr< SubDiagnostics<> >    BaseType;
  public:
    template< class... Args >
    DiagnosticsOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };
}
}

#endif
