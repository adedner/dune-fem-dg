#ifndef FEMDG_CHECKPOINTHANDLER_HH
#define FEMDG_CHECKPOINTHANDLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>

namespace Dune
{
namespace Fem
{
  template< class... StepperArg >
  class CombinedDefaultCheckPointHandler
  {
    typedef Dune::Fem::CheckPointerParameters                                                            CheckPointerParametersType;

    typedef std::tuple< typename std::add_pointer< StepperArg >::type... >                               StepperTupleType;
    typedef typename std::remove_pointer< typename std::tuple_element<0,StepperTupleType>::type >::type  FirstStepperType;
    typedef typename FirstStepperType::GridType                                                          GridType;
    typedef Dune::Fem::CheckPointer< GridType, std::tuple<> >                                             CheckPointerType;

    static_assert( Std::are_all_same< typename StepperArg::GridType... >::value,
                   "CombinedDefaultCheckPointHandler: GridType has to be equal for all steppers" );
    static_assert( sizeof ... ( StepperArg ) > 0,
                   "CombinedDefaultCheckPointHandler: called with zero steppers" );

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

    //explicit CombinedDefaultCheckPointHandler( StepperArg &&... arg )
    //  : tuple_( std::forward< StepperArg >( arg )... )
    //  ....
    //{}

    CombinedDefaultCheckPointHandler( const StepperTupleType& tuple )
      : tuple_( tuple ),
        checkPointer_(),
        keyPrefix_( std::get<0>( tuple_ )->name() ),
        checkParam_( Dune::ParameterKey::generate( keyPrefix_, "fem.io." ) )
    {
    }


    void registerData()
    {
      if( checkPointExists(keyPrefix_) )
        ForLoop< RegisterData, 0, sizeof ... ( StepperArg )-1 >::apply( tuple_ );
    }

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
                                               CheckPointerParametersType( Dune::ParameterKey::generate( keyPrefix, "fem.io." ) ) );
      return gridptr;
    }

    template< class TimeProviderImp >
    bool restoreData( TimeProviderImp& tp ) const
    {
      if( checkPointExists(keyPrefix_) )
      {
        // restore data
        checkPointer( tp ).restoreData( std::get<0>(tuple_)->space().grid(), checkPointFile( keyPrefix_ ) );
        return false;
      }
      return true;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp ) const
    {
      checkPointer( tp ).write( tp );
    }

  private:

    static inline std::string checkPointFile ( const std::string keyPrefix = "" )
    {
      static std::string checkFileName ;
      const std::string key( Dune::ParameterKey::generate( keyPrefix, "fem.io.checkpointrestartfile" ) );
      if( Fem::Parameter::exists( key ) )
        checkFileName = Fem::Parameter::getValue<std::string> ( key );
      else
        checkFileName = "" ;
      return checkFileName ;
    }

    template< class TimeProviderImp >
    CheckPointerType& checkPointer( const TimeProviderImp& tp  ) const
    {
      // create check point if not existent
      if( ! checkPointer_ )
        checkPointer_.reset( new CheckPointerType( std::get<0>(tuple_)->space().grid(), tp, checkParam_) );
      return *checkPointer_;
    }

    const StepperTupleType&           tuple_;
    mutable std::unique_ptr< CheckPointerType > checkPointer_;
    const std::string keyPrefix_;
    CheckPointerParametersType checkParam_;

  };



  template< class GridImp, class DataImp = std::tuple<> >
  class DefaultCheckPointHandler
  {
    typedef Dune::Fem::CheckPointer< GridImp, DataImp >   CheckPointerType;
    typedef Dune::Fem::CheckPointerParameters             CheckPointerParametersType;

    public:
    DefaultCheckPointHandler( const GridImp& grid, const std::string keyPrefix = "" )
      : checkPointer_(),
        keyPrefix_( keyPrefix ),
        grid_( grid ),
        checkParam_( Dune::ParameterKey::generate( keyPrefix, "fem.io." ) )
    {
    }

    static bool checkPointExists( const std::string keyPrefix = "" )
    {
      return (checkPointFile(keyPrefix).size() > 0 );
    }

    static Dune::GridPtr< GridImp > restoreGrid( const std::string keyPrefix = "" )
    {
      if( 0 == Dune::Fem::MPIManager :: rank () )
      {
        std::cout << std::endl;
        std::cout << "********************************************************" << std::endl;
        std::cout << "**   Restart from checkpoint `" << checkPointFile() << "' " << std::endl;
        std::cout << "********************************************************" << std::endl;
        std::cout << std::endl;
      }
      Dune::GridPtr< GridImp > gridptr ;
      gridptr = CheckPointerType::restoreGrid( checkPointFile( keyPrefix ), -1,
                                               CheckPointerParametersType( Dune::ParameterKey::generate( keyPrefix, "fem.io." ) ) );
      return gridptr;
    }

    template< class FirstArg, class ... Args>
    void registerData( FirstArg& head, Args& ... tail ) const
    {
      if( checkPointExists(keyPrefix_) )
      {
        Dune::Fem::persistenceManager << head;
        registerData( tail... );
      }
    }

    void registerData() const
    {}

    template< class TimeProviderImp >
    bool restoreData( TimeProviderImp& tp ) const
    {
      if( checkPointExists(keyPrefix_) )
      {
        // restore data
        checkPointer( tp ).restoreData( grid_, checkPointFile( keyPrefix_ ) );
        return false;
      }
      return true;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp ) const
    {
      checkPointer( tp ).write( tp );
    }

  private:

    static inline std::string checkPointFile ( const std::string keyPrefix = "" )
    {
      static std::string checkFileName ;
      const std::string key( Dune::ParameterKey::generate( keyPrefix, "fem.io.checkpointrestartfile" ) );
      if( Fem::Parameter::exists( key ) )
        checkFileName = Fem::Parameter::getValue<std::string> ( key );
      else
        checkFileName = "" ;
      return checkFileName ;
    }

    template< class TimeProviderImp >
    CheckPointerType& checkPointer( const TimeProviderImp& tp  ) const
    {
      // create check point if not existent
      if( ! checkPointer_ )
        checkPointer_.reset( new CheckPointerType( grid_, tp, checkParam_) );
      return *checkPointer_;
    }

    mutable std::unique_ptr< CheckPointerType > checkPointer_;
    const std::string keyPrefix_;
    const GridImp& grid_;
    CheckPointerParametersType checkParam_;

  };

  template< class GridImp>
  class NoCheckPointHandler
  {
    public:

    template< class ... Args >
    NoCheckPointHandler( Args&& ... )
    {}

    template< class ... Args>
    static bool checkPointExists( Args&& ... ) {return false;}

    template< class ... Args>
    static Dune::GridPtr< GridImp > restoreGrid( Args&& ... )
    {
      Dune::GridPtr< GridImp > gridptr;
      return gridptr;
    }

    template< class ... Args>
    void registerData( Args&& ... ) const {}

    template< class ... Args>
    bool restoreData( Args&& ... ) const { return false;}

    template< class ... Args>
    void step( Args&& ... ) const {}

  };


}
}
#endif
