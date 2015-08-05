#ifndef FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH
#define FEMDG_ALGORITHM_SOLUTIONLIMITERHANDLER_HH

#include <string>

namespace Dune
{
namespace Fem
{

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
