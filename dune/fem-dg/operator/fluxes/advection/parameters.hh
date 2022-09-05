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
      default_,
      /////////////// standard fluxes ////////////////////////
      //! no flux
      none,
      //! upwind flux
      upwind,
      //! local Lax-Friedrichs flux
      llf,
      //! general flux: parameter selection is done via parameter file!
      general,

      /////////////// euler fluxes //////////////////////////////
      //! the local Lax-Friedrichs flux (with wellbalance option)
      euler_llf,
      //! the Harten, Lax and van Leer (HLL) flux
      euler_hll,
      //! the Harten, Lax and van Leer (HLL) flux with potential temperature formulation
      euler_hll_pt,
      //! the HLLC flux
      euler_hllc,
      //! general flux: Parameter selection is done via parameter file!
      euler_general,

      /////////////// mhd fluxes //////////////////////////////
      mhd_dw,
      mhd_hllem,
      mhd_general,

      //! a flux implemented by the user and provided to the python code
      //! for C++ just add another enum id
      userdefined
    };

    //! Contains all known enums for advection fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::none, Enum::upwind, Enum::llf,
                                   Enum::euler_llf, Enum::euler_hll, Enum::euler_hll_pt, Enum::euler_hllc,
                                   Enum::mhd_dw, Enum::mhd_hllem
                                   };
    //! Contains all known names of advection fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "NONE", "UPWIND" , "LLF",
                                     "EULER-LLF", "EULER-HLL" , "EULER-HLL-PT", "EULER-HLLC",
                                     "MHD-DW", "MHD-HLLEM"
                                     };
    //! Number of known advection fluxes which can be chosen via parameter file.
    static const int  _size = 8;

  }

  /**
   * \brief Parameter class for advection flux parameters.
   *
   * \ingroup ParameterClass
   */
  class AdvectionFluxParameters
    : public Dune::Fem::LocalParameter< AdvectionFluxParameters, AdvectionFluxParameters >
  {
  public:
    typedef AdvectionFlux::Enum                  IdEnum;

    /**
     * \brief Constructor
     *
     * \param[in] keyPrefix key prefix for parameter file.
     */
    AdvectionFluxParameters( const std::string& keyPrefix,
                             const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : keyPrefix_( keyPrefix ),
        parameter_( parameter )
    {}

    AdvectionFluxParameters( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : AdvectionFluxParameters( "dgadvectionflux.", parameter )
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
      const int i = parameter_.getEnum( keyPrefix_ + "method", AdvectionFlux::_strings );
      return AdvectionFlux::_enums[i];
    }

    const Dune::Fem::ParameterReader& parameters() const { return parameter_; }

  private:
    const std::string keyPrefix_;
    const Dune::Fem::ParameterReader parameter_;
  };


}
}
#endif
