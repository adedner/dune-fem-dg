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

  template< class FunctionSpace >
  struct EmptyAdditional
  {
    static const int limitedDimRange = FunctionSpace :: dimRange;
    static const bool hasAdvection = true;
    static const bool hasDiffusion = false;
    static const bool hasStiffSource = false;
    static const bool hasNonStiffSource = false;
    static const bool hasFlux = true;
  };

  /**
   * \brief Defines an interface for advective fluxes.
   *
   * \ingroup AdvectionFluxes
   *
   * \tparam ModelImp type of the analytical model
   * \tparam FluxParameterImp type of the flux parameters
   */
  template <class ModelImp, class Additional=EmptyAdditional<typename ModelImp::DFunctionSpaceType>,
            class FluxParameterImp = AdvectionFluxParameters >
  class DGAdvectionFluxPythonUserDefine
    : public DGAdvectionFluxBase<
           AdvectionModelWrapper< typename ModelImp::GridPartType::GridType,
                                  ModelImp,
                                  Additional,
                                  NoLimiter< typename ModelImp::DFunctionSpaceType::DomainFieldType > >, FluxParameterImp >
  {
    typedef AdvectionModelWrapper< typename ModelImp::GridPartType::GridType,
                                   ModelImp,
                                   Additional,
                                   NoLimiter< typename ModelImp::DFunctionSpaceType::DomainFieldType > > ModelWrapperType;

  public:
    typedef typename std::conditional< std::is_base_of< Dune::Fem::IsFemDGModel, ModelImp >::value,
            ModelImp, ModelWrapperType > :: type ModelType;
                    //typename AdvectionLimiterFunctionSelector< typename ModelImp::DFunctionSpaceType::DomainFieldType, limiterFunctionId > :: type >

    typedef DGAdvectionFluxBase< ModelType, FluxParameterImp  > BaseType;

    static const int dimRange = ModelType::dimRange;
    typedef typename ModelType::DomainType         DomainType;
    typedef typename ModelType::RangeType          RangeType;
    typedef typename ModelType::JacobianRangeType  JacobianRangeType;
    typedef typename ModelType::FluxRangeType      FluxRangeType;
    typedef typename ModelType::FaceDomainType     FaceDomainType;

    typedef FluxParameterImp                       ParameterType;
    typedef typename ParameterType::IdEnum         IdEnum;

    /**
     * \brief Constructor
     *
     * \param[in] mod analytical model
     * \param[in] parameters  parameter reader
     */
    DGAdvectionFluxPythonUserDefine (const ModelImp& modelImp,
                                     const Dune::Fem::ParameterReader& parameter = Dune::Fem::Parameter::container() )
      : DGAdvectionFluxPythonUserDefine( modelImp, ParameterType( parameter ) )
    {
    }

    /**
     * \brief Constructor
     *
     * \param[in] mod analytical model
     * \param[in] parameters  parameter reader
     */
    DGAdvectionFluxPythonUserDefine (const ModelImp& modelImp,
                                     const ParameterType& parameter )
      : BaseType( *(new ModelType(modelImp)), parameter )
    {
      modelPtr_.reset( &this->model_ );
    }

    /**
     * \brief Constructor
     *
     * \param[in] mod analytical model
     * \param[in] parameters  parameter reader
     */
    DGAdvectionFluxPythonUserDefine (const ModelType& model,
                                     const ParameterType& parameter )
      : BaseType( model, parameter )
    {
    }

    /**
     * \brief Copy Constructor
     *
     * \param[in] other  object to copy
     */
    DGAdvectionFluxPythonUserDefine (const DGAdvectionFluxPythonUserDefine& other )
      : BaseType( other.model_ )
    {
    }

  protected:
    std::unique_ptr< const ModelType > modelPtr_;
  };

} // end namespace Fem
} // end namespace Dune
#endif
