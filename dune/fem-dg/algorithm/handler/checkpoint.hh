#ifndef FEMDG_CHECKPOINTHANDLER_HH
#define FEMDG_CHECKPOINTHANDLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  template< class GridImp >
  class GridCheckPointHandler
  {
  public:
    typedef GridImp                                            GridType;
    typedef Dune::Fem::CheckPointer< GridType, std::tuple<> >  CheckPointerType;
    typedef Dune::Fem::CheckPointerParameters                  CheckPointerParametersType;

    static bool checkPointExists( const std::string keyPrefix = "" )
    {
      return (checkPointFile(keyPrefix).size() > 0 );
    }

    static Dune::GridPtr< GridType > restoreGrid( const std::string keyPrefix = "" )
    {
      if( 0 == Dune::Fem::MPIManager :: rank () )
      {
        std::cout << std::endl;
        std::cout << "********************************************************" << std::endl;
        std::cout << "**   Restart from checkpoint `" << checkPointFile() << "' " << std::endl;
        std::cout << "********************************************************" << std::endl;
        std::cout << std::endl;
      }
      Dune::GridPtr< GridType > gridptr ;
      gridptr = CheckPointerType::restoreGrid( checkPointFile( keyPrefix ), -1,
                                               CheckPointerParametersType( ParameterKey::generate( keyPrefix, "fem.io." ) ) );
      return gridptr;
    }

    static inline std::string checkPointFile ( const std::string keyPrefix = "" )
    {
      static std::string checkFileName ;
      const std::string key( ParameterKey::generate( keyPrefix, "fem.io.checkpointrestartfile" ) );
      if( Fem::Parameter::exists( key ) )
        checkFileName = Fem::Parameter::getValue<std::string> ( key );
      else
        checkFileName = "";
      return checkFileName;
    }


  };


  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class CheckPointHandler;


  template< class AlgTupleImp, std::size_t... Ints >
  class CheckPointHandler< AlgTupleImp, Std::index_sequence< Ints... > >
    : public HandlerInterface
  {

    template< int i >
    struct RegisterData
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        if( std::get< i >( tuple )->checkPointSolution() )
          Dune::Fem::persistenceManager << *(std::get< i >( tuple )->checkPointSolution());
      }
    };

  public:
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef Std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef GridCheckPointHandler< typename std::remove_pointer< typename std::tuple_element< 0, TupleType >::type >::type::GridType >
                                                                                BaseType;

    typedef typename BaseType::GridType                                         GridType;
    typedef typename BaseType::CheckPointerType                                 CheckPointerType;
    typedef typename BaseType::CheckPointerParametersType                       CheckPointerParametersType;

    template< template<int> class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;

  public:

    CheckPointHandler( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) ),
        checkPointer_(),
        keyPrefix_( std::get<0>( tuple_ )->name() ),
        checkParam_( ParameterKey::generate( keyPrefix_, "fem.io." ) )
    {}

    template< class SubAlgImp, class TimeProviderImp >
    bool initialize_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( BaseType::checkPointExists(keyPrefix_) )
      {
        // restore data
        checkPointer( tp ).restoreData( std::get<0>(tuple_)->solution().space().grid(), BaseType::checkPointFile( keyPrefix_ ) );
        return false;
      }
      return true;
    }

    template< class SubAlgImp, class TimeProviderImp >
    void initialize_post( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( BaseType::checkPointExists(keyPrefix_) )
        ForLoopType< RegisterData >::apply( tuple_ );
    }

    template< class SubAlgImp, class TimeProviderImp >
    void preSolve_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      checkPointer( tp ).write( tp );
    }

    template< class TimeProviderImp >
    CheckPointerType& checkPointer( const TimeProviderImp& tp  ) const
    {
      // create check point if not existent
      if( ! checkPointer_ )
        checkPointer_.reset( new CheckPointerType( std::get<0>(tuple_)->solution().space().grid(), tp, checkParam_) );
      return *checkPointer_;
    }

  protected:
    TupleType               tuple_;
    mutable std::unique_ptr< CheckPointerType > checkPointer_;
    const std::string keyPrefix_;
    CheckPointerParametersType checkParam_;
  };


  template< class AlgTupleImp >
  class CheckPointHandler< AlgTupleImp, Std::index_sequence<> >
    : public HandlerInterface
  {
    public:

    template< class ... Args >
    CheckPointHandler( Args&& ... ) {}

    template< class ... Args>
    static bool checkPointExists( Args&& ... ) {return false;}

    template< class GridImp, class ... Args>
    static Dune::GridPtr< GridImp > restoreGrid( Args&& ... )
    {
      Dune::GridPtr< GridImp > gridptr;
      return gridptr;
    }

  };

}
}
#endif
