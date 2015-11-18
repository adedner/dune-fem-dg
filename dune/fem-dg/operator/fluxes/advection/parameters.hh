#ifndef FEMDG_ADVECTION_FLUX_PARAMETERS_HH
#define FEMDG_ADVECTION_FLUX_PARAMETERS_HH

#include <string>
#include <assert.h>

namespace Dune
{
namespace Fem
{

  struct AdvectionFluxIdentifier
  {
    typedef enum
    {
      none = 0,
      upwind = 1,
      llf = 2,
      general = 3
    } id;
  };


  class AdvectionFluxParameters
    : public Dune::Fem::LocalParameter< AdvectionFluxParameters, AdvectionFluxParameters >
  {
    const std::string keyPrefix_;
  public:
    typedef AdvectionFluxIdentifier MethodType;

    AdvectionFluxParameters( const std::string keyPrefix = "dgadvectionflux." )
      : keyPrefix_( keyPrefix )
    {}

    static std::string methodNames( const MethodType::id mthd )
    {
      const std::string method []
        = { "NONE", "UPWIND" , "LLF", };
      assert( mthd >= MethodType::none && mthd < MethodType::general );
      return method[ mthd ];
    }

    virtual MethodType::id getMethod() const
    {
      const std::string method []
        = { methodNames( MethodType::none ),
            methodNames( MethodType::upwind ),
            methodNames( MethodType::llf )
          };
      return (MethodType::id) Fem::Parameter::getEnum( keyPrefix_ + "method", method );
    }

  };


}
}
#endif
