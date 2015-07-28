#ifndef FEMDG_CHECKPOINTINITIALIZER_HH
#define FEMDG_CHECKPGRIDINITIALIZER_HH

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>

namespace Dune
{
namespace Fem
{


  template< class GridImp, class DataImp = tuple<> >
  class DefaultCheckPointInitializer
  {

    public:
    typedef Dune::Fem::CheckPointer< GridImp, DataImp >   CheckPointerType;

    DefaultCheckPointInitializer()
      : checkPointer_(0)
    {}


    ~DefaultCheckPointInitializer()
    {
      delete checkPointer_;
      checkPointer_ = 0;
    }



    static bool checkPointExists()
    {
      return (checkPointFile().size() > 0 );
    }

    static inline std::string checkPointFile ()
    {
      static bool initialized = false ;
      static std::string checkFileName ;
      if( ! initialized )
      {
        const std::string key( "fem.io.checkpointrestartfile" );
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

    static Dune::GridPtr< GridImp > restoreGrid()
    {
      Dune::GridPtr< GridImp > gridptr ;
      gridptr = CheckPointerType::restoreGrid( checkPointFile(), -1, Dune::Fem::CheckPointerParameters( Dune::ParameterKey::generate( "", "fem.io." ) ) );
      return gridptr;
    }


    template< class TimeProviderImp >
    bool restoreData( const GridImp& grid, const TimeProviderImp& tp ) const
    {
      std::string checkPointRestartFile = checkPointFile();

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
    }


    template< class TimeProviderImp >
    CheckPointerType& checkPointer( const GridImp& grid, TimeProviderImp& tp  ) const
    {
      // create check point if not existent
      if( ! checkPointer_ )
        checkPointer_ = new CheckPointerType( grid, tp, Dune::Fem::CheckPointerParameters( Dune::ParameterKey::generate( "", "fem.io." ) )  );
      return *checkPointer_;
    }

    template< class ObjImp >
    void addPersistentObject( const ObjImp& obj ) const
    {
      Dune::Fem::persistenceManager << obj;
    }

    private:
    mutable CheckPointerType* checkPointer_;

  };

}
}
#endif
