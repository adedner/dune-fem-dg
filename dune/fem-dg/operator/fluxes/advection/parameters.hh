#ifndef FEMDG_ADVECTION_FLUX_PARAMETERS_HH
#define FEMDG_ADVECTION_FLUX_PARAMETERS_HH

#include <string>
#include <assert.h>
#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

  namespace AdvectionFlux
  {
    enum Enum
    {
      none = 0,
      upwind = 1,
      llf = 2,
      general = 3
    };

    //parameters which can be chosen in parameter file
    const Enum        _enums[] = { Enum::none, Enum::upwind, Enum::llf };
    const std::string _strings[] = { "NONE", "UPWIND" , "LLF" };
    static const int  _size = 3;

    //helper class for static parameter selection
    template< Enum id = Enum::general >
    struct Identifier
    {
      static const Enum value = id;
      typedef Enum type;
    };

  }

  template< AdvectionFlux::Enum id = AdvectionFlux::Enum::general >
  class AdvectionFluxParameters
    : public Dune::Fem::LocalParameter< AdvectionFluxParameters<id>, AdvectionFluxParameters<id> >
  {
  public:
    typedef typename AdvectionFlux::Identifier<id> IdType;
  private:
    typedef typename IdType::type                  IdEnum;
  public:

    template< IdEnum ident >
    using GeneralIdType = AdvectionFlux::Identifier< ident >;

    AdvectionFluxParameters( const std::string keyPrefix = "dgadvectionflux." )
      : keyPrefix_( keyPrefix )
    {}

    static std::string methodNames( const IdEnum mthd )
    {
      for( int i = 0; i < AdvectionFlux::_size; i++)
        if( AdvectionFlux::_enums[i] == mthd )
          return AdvectionFlux::_strings[i];
      assert( false );
      return "invalid advection flux";
    }

    virtual IdEnum getMethod() const
    {
      const int i = Fem::Parameter::getEnum( keyPrefix_ + "method", AdvectionFlux::_strings );
      return AdvectionFlux::_enums[i];
    }

  private:
    const std::string keyPrefix_;
  };


}
}
#endif
