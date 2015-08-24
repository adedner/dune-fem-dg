#ifndef FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH
#define FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH

#include <string>

namespace Dune
{
namespace Fem
{
  template< class... StepperArg >
  class CombinedDefaultSolutionLimiterHandler
  {
  public:
    typedef std::tuple< typename std::add_pointer< StepperArg >::type... >                               StepperTupleType;
    typedef typename std::remove_pointer< typename std::tuple_element<0,StepperTupleType>::type >::type  FirstStepperType;


    static_assert( sizeof ... ( StepperArg ) > 0,
                   "CombinedDefaultSolutionLimiterHandler: called with zero steppers" );

    template< int i >
    struct LimitSolution
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get< i >( tuple )->limit( args... );
      }
    };

    CombinedDefaultSolutionLimiterHandler( const StepperTupleType& tuple )
      : tuple_( tuple )
    {}

    void step()
    {
      ForLoop< LimitSolution, 0, sizeof ... ( StepperArg )-1 >::apply( tuple_ );
    }

  private:
    const StepperTupleType& tuple_;

  };


  template< class LimiterOperatorImp >
  class DefaultSolutionLimiterHandler
  {
  public:
    typedef LimiterOperatorImp   LimiterOperatorType;

    DefaultSolutionLimiterHandler( const std::string keyPrefix = "" )
    {}

    template< class SolutionImp >
    void step( SolutionImp& sol )
    {
      limOp_->limit( sol );
    }

    void setLimiter( LimiterOperatorType* limOp )
    {
      limOp_ = limOp;
    }


  private:
    LimiterOperatorType* limOp_;

  };



  class NoSolutionLimiterHandler
  {
  public:
    NoSolutionLimiterHandler( const std::string keyPrefix = "" )
    {}

    template < class ... Args >
    void step( Args&& ... ) {};

    template < class ... Args >
    void setLimiter( Args&& ... ) {};
  };

}
}

#endif
