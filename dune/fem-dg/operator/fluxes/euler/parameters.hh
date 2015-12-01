#ifndef DUNE_FEM_DG_EULERFLUXES_PARAMETER_HH
#define DUNE_FEM_DG_EULERFLUXES_PARAMETER_HH

// system includes
#include <string>
#include <cmath>

#include <dune/fem/io/parameter.hh>
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

  namespace AdvectionFlux
  {
    enum Enum
    {
      llf = 0,
      hll = 1,
      hllc = 2,
      llf2 = 3,
      general = 4
    };

    //parameters which can be chosen in parameter file
    const Enum        _enums[] = { Enum::llf, Enum::hll, Enum::hllc, Enum::llf2 };
    const std::string _strings[] = { "LLF", "HLL" , "HLLC", "LLF2" };
    static const int  _size = 4;

    //helper class for static parameter selection
    template< Enum id = Enum::general >
    struct Identifier
    {
      static const Enum value = id;
      typedef Enum type;
    };
  }

  template< Euler::AdvectionFlux::Enum id = Euler::AdvectionFlux::Enum::general >
  class AdvectionFluxParameters
    : public Dune::Fem::LocalParameter< AdvectionFluxParameters<id>, AdvectionFluxParameters<id> >
  {
  public:
    typedef typename Euler::AdvectionFlux::Identifier<id> IdType;
  private:
    typedef typename IdType::type                  IdEnum;
  public:

    template< IdEnum ident >
    using GeneralIdType = Euler::AdvectionFlux::Identifier< ident >;

    AdvectionFluxParameters( const std::string keyPrefix = "dgadvectionflux." )
      : keyPrefix_( keyPrefix )
    {}

    static std::string methodNames( const IdEnum mthd )
    {
      for( int i = 0; i < Euler::AdvectionFlux::_size; i++)
        if( Euler::AdvectionFlux::_enums[i] == mthd )
          return Euler::AdvectionFlux::_strings[i];
      assert( false );
      return "invalid advection flux";
    }

    virtual IdEnum getMethod() const
    {
      const int i = Fem::Parameter::getEnum( keyPrefix_ + "method", Euler::AdvectionFlux::_strings );
      return Euler::AdvectionFlux::_enums[i];
    }
  private:
    const std::string keyPrefix_;
  };


}
}
}

#endif // file declaration
