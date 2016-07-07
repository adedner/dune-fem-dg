#ifndef DUNE_FEMDG_DATAIOCHECKPOINTING_STEPPER_HH
#define DUNE_FEMDG_DATAIOCHECKPOINTING_STEPPER_HH


// local includes
#include <dune/fem-dg/algorithm/sub/evolution.hh>
#include <dune/fem/function/common/instationary.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem-dg/algorithm/caller/checkpoint.hh>

namespace Dune
{
namespace Fem
{

  template< class AlgTupleImp,
            class IndexSequenceImp=typename std::make_index_sequence< std::tuple_size< AlgTupleImp >::value > >
  class CheckedCheckPointCaller;


  template< class AlgTupleImp, std::size_t... Ints >
  class CheckedCheckPointCaller< AlgTupleImp, std::index_sequence< Ints... > >
    : public CheckPointCaller< AlgTupleImp, std::index_sequence< Ints... > >
  {
    typedef CheckPointCaller< AlgTupleImp, std::index_sequence< Ints... > >  BaseType;
    typedef typename BaseType::CheckPointerType                              CheckPointerType;
    typedef typename BaseType::AlgTupleType                                  AlgTupleType;

    typedef std::index_sequence< Ints... >                                   IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();


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

    template< template<int> class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;

  public:

    CheckedCheckPointCaller( const AlgTupleType& tuple )
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
        ForLoopType< RestoreData >::apply( tuple_, tp );
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
