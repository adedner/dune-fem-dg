#ifndef FEMDG_ADVECTION_FLUX_PARAMETERS_HH
#define FEMDG_ADVECTION_FLUX_PARAMETERS_HH

#include <string>
#include <assert.h>
#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Namespace containing all parameters to select an advection flux.
   */
  namespace AdvectionFlux
  {
    /**
     * \brief Enum of all known advection flux implementations.
     *
     * \ingroup FemDGParameter
     */
    enum Enum
    {
      //! no flux
      none = 0,
      //! upwind flux
      upwind = 1,
      //! local Lax-Friedrichs flux
      llf = 2,
      //! general flux: parameter selection is done via parameter file!
      general = 3
    };

    //! Contains all known enums for advection fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::none, Enum::upwind, Enum::llf };
    //! Contains all known names of advection fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "NONE", "UPWIND" , "LLF" };
    //! Number of known advection fluxes which can be chosen via parameter file.
    static const int  _size = 3;

    //! Helper class for static parameter selection
    template< Enum id = Enum::general >
    struct Identifier
    {
      static const Enum value = id;
      typedef Enum type;
    };

  }

  /**
   * \brief Parameter class for advection flux parameters.
   *
   * \ingroup ParameterClass
   */
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

    /**
     * \brief Constructor
     *
     * \param[in] keyPrefix key prefix for parameter file.
     */
    AdvectionFluxParameters( const std::string keyPrefix = "dgadvectionflux." )
      : keyPrefix_( keyPrefix )
    {}

    /**
     * \brief returns name of the flux
     *
     * \param[in] mthd enum of Euler flux
     * \returns string which could be used for the selection of a flux in a parameter file.
     */
    static std::string methodNames( const IdEnum mthd )
    {
      for( int i = 0; i < AdvectionFlux::_size; i++)
        if( AdvectionFlux::_enums[i] == mthd )
          return AdvectionFlux::_strings[i];
      assert( false );
      return "invalid advection flux";
    }

    /**
     * \brief returns enum of the flux
     */
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
