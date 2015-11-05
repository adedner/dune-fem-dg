#ifndef FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH
#define FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH

#include <string>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>

namespace Dune
{
namespace Fem
{

  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class SolverMonitorHandler;

  template< class AlgTupleImp >
  class SolverMonitorHandler< AlgTupleImp, Std::index_sequence<> >
  {
  public:
    template <class ... Args>
    SolverMonitorHandler(Args&& ... )
    {}

    template <class ... Args>
    void print(Args&& ... ) const {};

    template <class ... Args>
    void step(Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};

    template <class ... Args>
    const double getData( Args&& ... ) const
    { return 0.0; }

  };

  template< class AlgTupleImp, std::size_t... Ints >
  class SolverMonitorHandler< AlgTupleImp, Std::index_sequence< Ints... > >
  {
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef Std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    enum CombinationType { max, min, sum, avg };

    template< class Caller >
    class LoopCallee
    {
      template<class C, class T, class... A >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::SolverMonitorHandlerType >::value >::type
      getMonitor( T, A&& ... ){}
      template<class C, class T, class... A >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::SolverMonitorHandlerType >::value >::type
      getMonitor( T elem, A &&... a )
      {
        if( elem->monitor() )
          C::applyImpl(elem->monitor(), std::forward<A>(a)... );
      }
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... A >
        static void apply ( Tuple &tuple, A&& ... a )
        { getMonitor< Caller >( std::get<i>( tuple ), std::forward<A>(a)... ); }
      };
    };

    struct Step {
      template<class T, class... A > static void applyImpl( T e, A&& ... a )
      { e->step( std::forward<A>(a)... ); }
    };

    struct Finalize {
      template<class T, class... A > static void applyImpl( T e, A&& ... a )
      { e->finalize( std::forward<A>(a)... ); }
    };

    struct GetData {
      template<class T, class... A > static void applyImpl( T e, double& res, const std::string name, CombinationType comb, A&& ... a )
      {
        switch( comb )
        {
          case CombinationType::max:
            res = std::max( res, e->getData( name ) );
          case CombinationType::min:
            res = std::min( res, e->getData( name ) );
          case CombinationType::sum:
            res += e->getData( name );
          case CombinationType::avg:
            res += e->getData( name );
        }
      }
    };

    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0, numAlgs - 1 >;

  public:

    SolverMonitorHandler( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) )
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
      double res = (comb == CombinationType::max) ? std::numeric_limits<double>::max() : 0.0;
      ForLoopType< GetData >::apply( tuple_, res, name, comb );
      return res;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      ForLoopType< Step >::apply( tuple_, tp );
    }

    void finalize( const double gridWidth, const double gridSize )
    {
      ForLoopType< Finalize >::apply( tuple_, gridWidth, gridSize );
    }

  private:
    TupleType  tuple_;
  };






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
