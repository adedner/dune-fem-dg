#ifndef DUNE_FEM_DG_EULERFLUXES_PARAMETER_HH
#define DUNE_FEM_DG_EULERFLUXES_PARAMETER_HH

// system includes
#include <string>
#include <cmath>

#include "../advection/fluxbase.hh"

// dune-grid includes
#if WELLBALANCE
#include <dune/grid/common/genericreferenceelements.hh>
#endif

#include <dune/fem-dg/operator/fluxes/rotator.hh>

namespace Dune
{
namespace Fem
{
namespace Euler
{

  struct FluxIdentifier
  {
    typedef enum
    {
      llf = 0,
      hll = 1,
      hllc = 2,
      llf2 = 3,
      general = 4
    } id;
  };

  class EulerFluxParameters
    : public Dune::Fem::LocalParameter< EulerFluxParameters, EulerFluxParameters >
  {
    const std::string keyPrefix_;
  public:
    typedef FluxIdentifier MethodType;

    EulerFluxParameters( const std::string keyPrefix = "dgadvectionflux." )
      : keyPrefix_( keyPrefix )
    {}

    static std::string methodNames( const MethodType::id mthd )
    {
      const std::string method []
        = { "LLF", "HLL" , "HLLC", "LLF2" };
      assert( mthd >= MethodType::llf && mthd < MethodType::general );
      std::cout << mthd << std::endl;
      return method[ mthd ];
    }

    virtual MethodType::id getMethod() const
    {
      const std::string method []
        = { methodNames( MethodType::llf ),
            methodNames( MethodType::hll ),
            methodNames( MethodType::hllc ),
            methodNames( MethodType::llf2 )
          };
      return (MethodType::id) Fem::Parameter::getEnum( keyPrefix_ + "method", method );
    }

  };


}
}
}

#endif // file declaration
