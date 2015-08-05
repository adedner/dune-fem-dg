#ifndef FEMDG_CHECKPOINTHANDLER_HH
#define FEMDG_CHECKPOINTHANDLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>

namespace Dune
{
namespace Fem
{


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
