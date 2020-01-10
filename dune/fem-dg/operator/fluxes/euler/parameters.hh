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
   * \brief Parameter class for (Euler) advection flux parameters.
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
    AdvectionFluxParameters( const std::string keyPrefix = "dgadvectionflux.",
                             const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : keyPrefix_( keyPrefix ),
        parameter_( parameter )
    {}

    /**
     * \brief Constructor
     *
     * \param[in] keyPrefix key prefix for parameter file.
     */
    AdvectionFluxParameters( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : parameter_( parameter )
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
  private:
    const std::string keyPrefix_;
    const Dune::Fem::ParameterReader parameter_;
  };


}
}
}

#endif // file declaration
