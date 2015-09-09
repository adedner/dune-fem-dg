#ifndef FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH
#define FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH

#include <string>

namespace Dune
{
namespace Fem
{

  template< class ... StepperArg >
  class SolutionLimiterHandler;


  template< class StepperHead, class... StepperArg >
  class SolutionLimiterHandler< StepperHead, StepperArg ... >
  {
  public:
    typedef std::tuple< typename std::add_pointer< StepperHead >::type,
      typename std::add_pointer< StepperArg >::type... >                               StepperTupleType;

    template< int i >
    struct LimitSolution
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args&& ... a )
      {
        std::get< i >( tuple )->limit( std::forward<Args>(a)... );
      }
    };

    SolutionLimiterHandler( const StepperTupleType& tuple )
      : tuple_( tuple )
    {}

    void step()
    {
      ForLoop< LimitSolution, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
    }

  private:
    const StepperTupleType& tuple_;
  };

  template<>
  class SolutionLimiterHandler<>
  {
  public:
    template< class ... Args >
    SolutionLimiterHandler ( Args && ... ) {}

    void step () {}
  };

}
}

#endif
