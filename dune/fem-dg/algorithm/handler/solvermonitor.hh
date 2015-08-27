#ifndef FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH
#define FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH

#include <string>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/algorithm/monitor.hh>

namespace Dune
{
namespace Fem
{

  template< class... StepperArg >
  class CombinedDefaultSolverMonitorHandler;

  template< >
  class CombinedDefaultSolverMonitorHandler<>
  {
  public:
    template <class ... Args>
    CombinedDefaultSolverMonitorHandler(Args&& ... )
    {}

    template <class ... Args>
    void print(Args&& ... ) const {};

    template <class ... Args>
    void step(Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };


  template< class StepperHead, class... StepperArg >
  class CombinedDefaultSolverMonitorHandler< StepperHead, StepperArg... >
  {
  public:
    typedef std::tuple< typename std::add_pointer< StepperHead >::type,
                        typename std::add_pointer< StepperArg >::type... >                     StepperTupleType;

    enum CombinationType { max, min, sum, avg };

    template< int i >
    struct Step
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
         std::get< i >( tuple )->monitor().step( args... );
      }
    };
    template< int i >
    struct Finalize
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
         std::get< i >( tuple )->monitor().finalize( args... );
      }
    };

    CombinedDefaultSolverMonitorHandler( const StepperTupleType& tuple )
      : tuple_( tuple )
    {}

    template< class... StringType >
    void print( std::string str, StringType& ...tail )
    {
      std::cout << str << ":  " << getData( str ) << ", ";
      print( tail... );
      }

    void print() {};

    const double getData( const std::string name, CombinationType comb = CombinationType::max ) const
    {
      return getData( tuple_, name, comb, Std::index_sequence_for< StepperHead, StepperArg ... >() );
      }

    template< std::size_t ... i >
    static double getData( const StepperTupleType& tuple, const std::string name, CombinationType comb, Std::index_sequence< i ... > )
    {
      switch( comb )
      {
        case CombinationType::max:
            return Std::max( std::get< i >( tuple )->monitor().getData( name )... );
        case CombinationType::min:
            //return Std::min( std::get< i >( tuple )->monitor().getData( name )... );
        case CombinationType::sum:
            return Std::sum( std::get< i >( tuple )->monitor().getData( name )... );
        case CombinationType::avg:
            return Std::sum( std::get< i >( tuple )->monitor().getData( name )... ) / std::tuple_size< StepperTupleType >::value;
        default: return 0;
      }
     }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      ForLoop< Step, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
    }

    void finalize( const double gridWidth, const double gridSize )
    {
      ForLoop< Finalize, 0, sizeof ... ( StepperArg ) >::apply( tuple_, gridWidth, gridSize );
    }

  private:
    StepperTupleType  tuple_;
  };


  template< class SolverMonitorImp >
  class DefaultSolverMonitorHandler
  {
  public:
    typedef SolverMonitorImp               SolverMonitorType;

    typedef std::map< std::string, std::tuple< double*, double*, bool > > DataDoubleType;
    typedef std::map< std::string, std::tuple< long unsigned int*, long unsigned int*, bool > >       DataIntType;

    DefaultSolverMonitorHandler( const std::string keyPrefix = "" )
      : solverMonitor_()
    {}

    // call inside stepper something like this:
    // registerData( tuple->monitor().newton_iterations_ )
    // externalMonitorData: if set, the externalMonitorData is copied to monitorData in each step() call
    void registerData( const std::string name, int* monitorData, int* externalMonitorData = nullptr, bool internal = false )
    {
      assert( monitorData );
      //dangerous cast
      dataInt_.insert( std::make_pair(name, std::make_tuple( reinterpret_cast<long unsigned int*>(monitorData), reinterpret_cast<long unsigned int*>(externalMonitorData), internal ) ) );
    }

    void registerData( const std::string name, long unsigned int* monitorData, long unsigned int* externalMonitorData = nullptr, bool internal = false )
    {
      assert( monitorData );
      dataInt_.insert( std::make_pair(name, std::make_tuple( monitorData, externalMonitorData, internal ) ) );
    }

    void registerData( const std::string name, double* monitorData, double* externalMonitorData = nullptr, bool internal = false )
    {
      assert( monitorData );
      dataDouble_.insert( std::make_pair(name, std::make_tuple( monitorData, externalMonitorData, internal ) ) );
    }

    const double getData( const std::string name )
    {
      if( dataInt_.find(name) != dataInt_.end() )
      {
        assert( *std::get<0>(dataInt_[ name ]) );
        return (double)*std::get<0>(dataInt_[ name ]);
      }
      assert( std::get<0>(dataDouble_[ name ]) );
      return *std::get<0>(dataDouble_[ name ]);
    }


    template< class... StringType >
    void print( std::string str, StringType&... tail )
  {
      std::cout << str << ":  " << getData( str ) << ", ";
      print( tail... );
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      for (auto it = std::begin(dataDouble_); it!=std::end(dataDouble_); ++it)
        if( std::get<1>(it->second) )
          *std::get<0>(it->second) = *std::get<1>(it->second);
      for (auto it = std::begin(dataInt_); it!=std::end(dataInt_); ++it)
        if( std::get<1>(it->second) )
          *std::get<0>(it->second) = *std::get<1>(it->second);
      solverMonitor_.setTimeStepInfo( tp );
    }

    void finalize( const double gridWidth, const double gridSize )
    {
      solverMonitor_.finalize( gridWidth, gridSize );
    }

    SolverMonitorType& monitor()
    {
      return solverMonitor_;
    }

  private:
    SolverMonitorType solverMonitor_;
    DataDoubleType            dataDouble_;
    DataIntType               dataInt_;
  };

  class NoSolverMonitorHandler
  {
  public:
    typedef SolverMonitor<1>     SolverMonitorType;

    NoSolverMonitorHandler( const std::string keyPrefix = "" )
      : monitor_()
    {}

    template <class ... Args>
    double registerData(Args&& ... ) const {};

    template <class ... Args>
    double getData(Args&& ... ) const {};

    template <class ... Args>
    void print(Args&& ... ) const {};

    template <class ... Args>
    void step(Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};

    template <class ... Args>
    SolverMonitorType& monitor(Args&& ... ) { return monitor_; };

  private:
    SolverMonitorType monitor_;

  };

}
}

#endif
