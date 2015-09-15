#ifndef FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH
#define FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH

#include <dune/fem-dg/misc/diagnostics.hh>
#include <dune/fem-dg/misc/optional.hh>

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
          C::applyImpl(elem->diagnostics(), std::forward<Args>(a)... );
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
      template<class T, class... Args > static void applyImpl( T e, Args&& ... a )
      { e->step( std::forward<Args>(a)... ); }
    };
    struct Finalize {
      template<class T, class... Args > static void applyImpl( T e, Args&& ... a )
      { e->finalize( std::forward<Args>(a)... ); }
    };

  public:

    DiagnosticsHandler( const StepperTupleType& tuple )
      : tuple_( tuple )
    {}

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      ForLoop< LoopCallee<Write>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
    }

    void finalize() const
    {
      ForLoop< LoopCallee<Finalize>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
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
    void finalize( Args&& ... ) const {};
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


    SubDiagnosticsHandler( const std::string keyPrefix = "" )
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

    const double getData( const std::string name )
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
      const double ldt = tp.deltaT();
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
    template <class ... Args>
    SubDiagnosticsHandler( Args&& ... )
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

  template< class Obj >
  class DiagnosticsHandlerOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    DiagnosticsHandlerOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

  template<>
  class DiagnosticsHandlerOptional< void >
    : public OptionalNullPtr< SubDiagnosticsHandler<> >
  {
    typedef OptionalNullPtr< SubDiagnosticsHandler<> >    BaseType;
  public:
    template< class... Args >
    DiagnosticsHandlerOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };
}
}

#endif
