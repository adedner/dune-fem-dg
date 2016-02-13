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
      /////////////// standard fluxes ////////////////////////
      //! no flux
      none = 0,
      //! upwind flux
      upwind = 1,
      //! local Lax-Friedrichs flux
      llf = 2,
      //! general flux: parameter selection is done via parameter file!
      general = 3,

      /////////////// euler fluxes //////////////////////////////
      //! the local Lax-Friedrichs flux (with wellbalance option)
      euler_llf = 4,
      //! the Harten, Lax and van Leer (HLL) flux
      euler_hll = 5,
      //! the HLLC flux
      euler_hllc = 6,
      //! general flux: Parameter selection is done via parameter file!
      euler_general = 7

    };

    //! Contains all known enums for advection fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::none, Enum::upwind, Enum::llf, Enum::euler_llf, Enum::euler_hll, Enum::euler_hllc };
    //! Contains all known names of advection fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "NONE", "UPWIND" , "LLF", "EULER-LLF", "EULER-HLL" , "EULER-HLLC" };
    //! Number of known advection fluxes which can be chosen via parameter file.
    static const int  _size = 6;

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
    typedef AdvectionFlux::Enum                  IdEnum;

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
