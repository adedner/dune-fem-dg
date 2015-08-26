#ifndef FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH
#define FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH

#include <dune/fem-dg/misc/diagnostics.hh>

namespace Dune
{
namespace Fem
{

  template< class... StepperArg >
  class CombinedDefaultDiagnosticsHandler;


  template< class StepperHead, class ... StepperArg >
  class CombinedDefaultDiagnosticsHandler< StepperHead, StepperArg ... >
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

    CombinedDefaultDiagnosticsHandler( const StepperTupleType& tuple )
      : tuple_( tuple )
    {}

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp, int numElements )
    {
      ForLoop< Write, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp, numElements );
    }


    void finalize() const
    {
      ForLoop< Finalize, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
    }


  private:
    StepperTupleType tuple_;
  };


  template<>
  class CombinedDefaultDiagnosticsHandler<>
  {
  public:
    template< class ... Args >
    CombinedDefaultDiagnosticsHandler ( Args && ... ) {}

    template <class ... Args>
    void step( Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };


  class DefaultDiagnosticsHandler
  {
  public:

    DefaultDiagnosticsHandler()
      : diagnostics_( true ),
        timings_()
    {}

    template< class... DoubleArg >
    void addTimings( double head, DoubleArg... tail )
    {
      timings_.push_back( head );
      addTimings( tail... );
    }

    void addTimings()
    {}

    template< class... DoubleVectorArg >
    void addTimings( std::vector< double >& head, DoubleVectorArg&... tail )
    {
      timings_.insert( timings_.end(), head.begin(), head.end() );
      addTimings( tail... );
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp, int numElements )
    {
      const double ldt = tp.deltaT();
      diagnostics_.write( tp.time() + ldt, ldt, numElements, timings_ );
      timings_.clear();
    }


    void finalize() const
    {
      diagnostics_.flush();
    }


  private:
    Diagnostics diagnostics_;
    std::vector< double > timings_;
  };



  class NoDiagnosticsHandler
  {
  public:
    template <class ... Args>
    void step( Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };

}
}

#endif
