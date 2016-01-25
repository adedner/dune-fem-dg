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

  /**
   * \brief Namespace containing all parameters to select an Euler flux.
   */
  namespace AdvectionFlux
  {
    /**
     * \brief Enum of all known Euler flux implementations.
     *
     * \ingroup FemDGParameter
     */
    enum Enum
    {
      //! the local Lax-Friedrichs flux (with wellbalance option)
      llf = 0,
      //! the Harten, Lax and van Leer (HLL) flux
      hll = 1,
      //! the HLLC flux
      hllc = 2,
      //! the local Lax-Friedrichs flux
      llf2 = 3,
      //! general flux: Parameter selection is done via parameter file!
      general = 4
    };

    //! Contains all known enums for Euler fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::llf, Enum::hll, Enum::hllc, Enum::llf2 };
    //! Contains all known names of Euler fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "LLF", "HLL" , "HLLC", "LLF2" };
    //! Number of known Euler fluxes which can be chosen via parameter file.
    static const int  _size = 4;

    //! Helper class for static parameter selection
    template< Enum id = Enum::general >
    struct Identifier
    {
      static const Enum value = id;
      typedef Enum type;
    };
  }

  /**
   * \brief Parameter class for (Euler) advection flux parameters.
   *
   * \ingroup ParameterClass
   */
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
      for( int i = 0; i < Euler::AdvectionFlux::_size; i++)
        if( Euler::AdvectionFlux::_enums[i] == mthd )
          return Euler::AdvectionFlux::_strings[i];
      assert( false );
      return "invalid advection flux";
    }

    /**
     * \brief returns enum of the flux
     */
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
