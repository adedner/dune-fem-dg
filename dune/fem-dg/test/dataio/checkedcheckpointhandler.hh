#ifndef DUNE_FEMDG_DATAIOCHECKPOINTING_STEPPER_HH
#define DUNE_FEMDG_DATAIOCHECKPOINTING_STEPPER_HH


// local includes
#include <dune/fem-dg/algorithm/sub/evolution.hh>
#include <dune/fem/function/common/instationary.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>

namespace Dune
{
namespace Fem
{

  template< class... StepperArg >
  class CheckedCheckPointHandler;

  template< class StepperHead, class ... StepperArg >
  class CheckedCheckPointHandler< StepperHead, StepperArg ... >
    : public CheckPointHandler< StepperHead, StepperArg ... >
  {
    typedef CheckPointHandler< StepperHead, StepperArg ... >  BaseType;
    typedef typename BaseType::CheckPointerType               CheckPointerType;
    typedef typename BaseType::StepperTupleType               StepperTupleType;

    template< int i >
    struct RestoreData
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        if( std::get< i >( tuple )->checkPointSolution() )
          std::get< i >( tuple )->checkCheckPointSolutionValid( std::forward<Args>( args )... );
      }
    };

    using BaseType::tuple_;
    using BaseType::keyPrefix_;

  public:

    CheckedCheckPointHandler( const StepperTupleType& tuple )
      : BaseType( tuple )
    {}

    void registerData()
    {
      BaseType::registerData();
    }

    template< class TimeProviderImp >
    bool restoreData( TimeProviderImp& tp ) const
    {
      if( !BaseType::restoreData( tp ) )
      {
        ForLoop< RestoreData, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
        return false;
      }
      return true;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp ) const
    {
      BaseType::step( tp );
    }

    template< class TimeProviderImp >
    CheckPointerType& checkPointer( const TimeProviderImp& tp  ) const
    {
      return BaseType::checkPointer( tp );
    }
  };

}
}

#endif // #ifndef DUNE_FEMDG_CHECKPOINTING_STEPPER_HH
