#ifndef FEMDG_SUBALGORITHM_SOLVERMONITORHANDLER_HH
#define FEMDG_SUBALGORITHM_SOLVERMONITORHANDLER_HH

#include <string>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>

namespace Dune
{
namespace Fem
{

  template< class... SolverMonitorImp >
  class SubSolverMonitorHandler;


  template< class SolverMonitorImp >
  class SubSolverMonitorHandler< SolverMonitorImp >
  {
  public:
    typedef SolverMonitorImp               SolverMonitorType;

    typedef std::map< std::string, std::tuple< double*, double*, bool > > DataDoubleType;
    typedef std::map< std::string, std::tuple< long unsigned int*, long unsigned int*, bool > >       DataIntType;

    SubSolverMonitorHandler( const std::string keyPrefix = "" )
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
        assert( std::get<0>(dataInt_[ name ]) );
        return (double)*std::get<0>(dataInt_[ name ]);
      }
      if( dataDouble_.find(name) != dataDouble_.end() )
      {
        assert( std::get<0>(dataDouble_[ name ]) );
        return *std::get<0>(dataDouble_[ name ]);
      }
      return 0.0;
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


  template<>
  class SubSolverMonitorHandler<>
  {
    //simple no solver monitor handler
    //capturing most famous issues
    struct NoSolverMonitorType
    {
      double gridWidth;
      unsigned long int elements;
      int timeSteps;
      int newton_iterations;
      int ils_iterations;
      int operator_calls;
      int max_newton_iterations;
      int max_ils_iterations;
      double avgTimeStep;
      double minTimeStep;
      double maxTimeStep;
    };
  public:

    typedef NoSolverMonitorType    SolverMonitorType;

    SubSolverMonitorHandler( const std::string keyPrefix = "" )
      : monitor_()
    {}

    template <class ... Args>
    void registerData(Args&& ... ) const {}

    template <class ... Args>
    double getData(Args&& ... ) const { return 0.0; }

    template <class ... Args>
    void print(Args&& ... ) const {}

    template <class ... Args>
    void step(Args&& ... ) const {}

    template <class ... Args>
    void finalize(Args&& ... ) const {}

    template <class ... Args>
    SolverMonitorType& monitor(Args&& ... ) { return monitor_; }

  private:
    SolverMonitorType monitor_;

  };


  template< class Obj >
  class SolverMonitorHandlerOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    SolverMonitorHandlerOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

  template<>
  class SolverMonitorHandlerOptional< void >
    : public OptionalNullPtr< SubSolverMonitorHandler<> >
  {
    typedef OptionalNullPtr< SubSolverMonitorHandler<> >    BaseType;
  public:
    template< class... Args >
    SolverMonitorHandlerOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

}
}

#endif
