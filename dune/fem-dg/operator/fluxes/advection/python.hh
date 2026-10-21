#ifndef FEMDG_ADVECTION_FLUX_PYTHON_USERDEFINED_HH
#define FEMDG_ADVECTION_FLUX_PYTHON_USERDEFINED_HH

#include <string>
#include <assert.h>

#include <dune/fem-dg/operator/fluxes/advection/fluxbase.hh>
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/models/modelwrapper.hh>

namespace Dune
{
 namespace Fem
 {

  /**
   * \brief Defines an interface for advective fluxes passed from Python side.
   *
   * \ingroup AdvectionFluxes
   *
   * \tparam ModelImp type of the analytical model
   * \tparam FluxParameterImp type of the flux parameters
   * \tparam enableRightModel true if two separate models are needed (left/right)
   */
  template <class ModelImp,
            class FluxParameterImp = AdvectionFluxParameters,
            bool enableRightModel  = false >
  class DGAdvectionFluxPythonUserDefined
    : public DGAdvectionFluxBase< ModelImp, FluxParameterImp, enableRightModel >
  {
  public:
    typedef ModelImp ModelType;

    typedef DGAdvectionFluxBase< ModelType, FluxParameterImp, enableRightModel > BaseType;

    //static const int dimRange = ModelType::dimRange;
    //typedef typename ModelType::DomainType         DomainType;
    //typedef typename ModelType::RangeType          RangeType;
    //typedef typename ModelType::JacobianRangeType  JacobianRangeType;
    //typedef typename ModelType::FluxRangeType      FluxRangeType;
    //typedef typename ModelType::FaceDomainType     FaceDomainType;

    typedef FluxParameterImp                       ParameterType;
    //typedef typename ParameterType::IdEnum         IdEnum;

    /**
     * \brief Constructor
     *
     * \param[in] mod analytical model
     * \param[in] parameters  parameter reader
     */
    DGAdvectionFluxPythonUserDefined (const ModelImp& modelImp,
                                     const Dune::Fem::ParameterReader& parameter = Dune::Fem::Parameter::container() )
      : DGAdvectionFluxPythonUserDefined( modelImp, ParameterType( parameter ) )
    {
    }

    /**
     * \brief Constructor
     *
     * \param[in] mod analytical model
     * \param[in] parameters  advection parameters
     */
    DGAdvectionFluxPythonUserDefined (const ModelType& model,
                                      const ParameterType& parameter )
      : BaseType( model, parameter )
    {
    }

    /**
     * \brief Copy Constructor
     *
     * \param[in] other  object to copy
     */
    DGAdvectionFluxPythonUserDefined (const DGAdvectionFluxPythonUserDefined& other )
      : BaseType( other.model_ )
    {
    }
  };

 } // end namespace Fem
} // end namespace Dune
#endif
