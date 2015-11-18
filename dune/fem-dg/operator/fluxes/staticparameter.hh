#ifndef FEMDG_STATICPARAMETER_HH
#define FEMDG_STATICPARAMETER_HH

#include <string>

namespace Dune
{
namespace Fem
{

  template< class ParametersImp, typename ParametersImp::MethodType::id ident >
  class StaticParameter
    : public ParametersImp
  {
  public:
    StaticParameter( const std::string keyPrefix = "" )
      : ParametersImp( keyPrefix )
    {}

    virtual typename ParametersImp::MethodType::id getMethod() const
    {
      return ident;
    }
  };

}
}
#endif
