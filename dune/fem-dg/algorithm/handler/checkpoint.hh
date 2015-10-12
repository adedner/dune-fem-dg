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
                                               CheckPointerParametersType( Dune::ParameterKey::generate( keyPrefix, "fem.io." ) ) );
      return gridptr;
    }

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


  };


  template< class... StepperArg >
  class CheckPointHandler;

  template< class StepperHead, class ... StepperArg >
  class CheckPointHandler< StepperHead, StepperArg ... >
  {
    typedef GridCheckPointHandler< typename StepperHead::GridType >                                 BaseType;
    typedef typename BaseType::GridType                                                             GridType;
    typedef typename BaseType::CheckPointerType                                                     CheckPointerType;
    typedef typename BaseType::CheckPointerParametersType                                           CheckPointerParametersType;
    typedef std::tuple< typename std::add_pointer< StepperHead >::type, typename std::add_pointer< StepperArg >::type... > StepperTupleType;

    static_assert( Std::are_all_same< typename StepperHead::GridType, typename StepperArg::GridType... >::value,
                   "CheckPointHandler: GridType has to be equal for all steppers" );

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

    CheckPointHandler( const StepperTupleType& tuple )
      : tuple_( tuple ),
        checkPointer_(),
        keyPrefix_( std::get<0>( tuple_ )->name() ),
        checkParam_( Dune::ParameterKey::generate( keyPrefix_, "fem.io." ) )
    {}

    void registerData()
    {
      if( BaseType::checkPointExists(keyPrefix_) )
        ForLoop< RegisterData, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
    }

    template< class TimeProviderImp >
    bool restoreData( TimeProviderImp& tp ) const
    {
      if( BaseType::checkPointExists(keyPrefix_) )
      {
        // restore data
        checkPointer( tp ).restoreData( std::get<0>(tuple_)->solution().space().grid(), BaseType::checkPointFile( keyPrefix_ ) );
        return false;
      }
      return true;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp ) const
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

    const StepperTupleType&           tuple_;
    mutable std::unique_ptr< CheckPointerType > checkPointer_;
    const std::string keyPrefix_;
    CheckPointerParametersType checkParam_;
  };


  template<>
  class CheckPointHandler<>
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