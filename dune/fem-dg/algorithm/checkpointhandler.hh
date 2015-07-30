#ifndef FEMDG_CHECKPOINTHANDLER_HH
#define FEMDG_CHECKPGRIDHANDLER_HH

#include <memory>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>

namespace Dune
{
namespace Fem
{


  template< class GridImp, class DataImp = tuple<> >
  class DefaultCheckPointHandler
  {

    public:
    typedef Dune::Fem::CheckPointer< GridImp, DataImp >   CheckPointerType;

    DefaultCheckPointHandler( const std::string keyPrefix = "" )
      : checkPointer_(),
        keyPrefix_( keyPrefix )
    {}

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
                                               Dune::Fem::CheckPointerParameters( Dune::ParameterKey::generate( keyPrefix, "fem.io." ) ) );
      return gridptr;
    }


    template< class TimeProviderImp, class FirstArg, class ... Args>
    bool restoreData( const GridImp& grid, TimeProviderImp& tp, const FirstArg& head, const Args& ... tail ) const
    {
      addPersistentObject( head );
      return restoreData( grid, tp, tail... );
    };

    template< class TimeProviderImp >
    bool restoreData( const GridImp& grid, const TimeProviderImp& tp ) const
    {
      std::string checkPointRestartFile = checkPointFile( keyPrefix_ );

      // if check file is non-zero a restart is performed
      if( checkPointRestartFile.size() > 0 )
      {
        // restore data
        checkPointer( grid, tp ).restoreData( grid, checkPointRestartFile );
        return false;
      }
      return true;
    }

    template< class TimeProviderImp >
    bool writeData( const GridImp& grid, const TimeProviderImp& tp ) const
    {
      checkPointer( grid, tp ).write( tp );
      return true;
    }

  private:
    static inline std::string checkPointFile ( const std::string keyPrefix = "" )
    {
      static bool initialized = false ;
      static std::string checkFileName ;
      if( ! initialized )
      {
        const std::string key( Dune::ParameterKey::generate( keyPrefix, "fem.io.checkpointrestartfile" ) );
        if( Fem::Parameter::exists( key ) )
        {
          checkFileName = Fem::Parameter::getValue<std::string> ( key );
        }
        else
          checkFileName = "" ;
        initialized = true ;
      }

      return checkFileName ;
    }

    template< class TimeProviderImp >
    CheckPointerType& checkPointer( const GridImp& grid, const TimeProviderImp& tp  ) const
    {
      // create check point if not existent
      if( ! checkPointer_ )
        checkPointer_.reset( new CheckPointerType( grid, tp, Dune::Fem::CheckPointerParameters( Dune::ParameterKey::generate( keyPrefix_, "fem.io." ) )  ) );
      return *checkPointer_;
    }

    template< class ObjImp >
    void addPersistentObject( const ObjImp& obj ) const
    {
      Dune::Fem::persistenceManager << obj;
    }

    mutable std::unique_ptr< CheckPointerType > checkPointer_;
    const std::string keyPrefix_;

  };


  template< class GridImp, class DataImp = tuple<> >
  class NoCheckPointHandler
  {
    public:
    typedef Dune::Fem::CheckPointer< GridImp, DataImp >   CheckPointerType;

    NoCheckPointHandler()
    {}

    template< class ... Args>
    static bool checkPointExists( const Args& ... args )
    {
      return false;
    }

    template< class ... Args>
    static Dune::GridPtr< GridImp > restoreGrid( const Args& ... args )
    {
      Dune::GridPtr< GridImp > gridptr;
      return gridptr;
    }

    template< class ... Args>
    bool restoreData( const Args& ... args ) const
    {
      return false;
    }

    template< class ... Args>
    bool writeData( const Args& ... args ) const
    {
      return false;
    }

  };


}
}
#endif
