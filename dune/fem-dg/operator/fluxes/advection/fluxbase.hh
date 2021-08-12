#ifndef FEMDG_ADVECTION_FLUX_HH
#define FEMDG_ADVECTION_FLUX_HH

#include <string>
#include <assert.h>

#include "parameters.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Defines an interface for advective fluxes.
   *
   * \ingroup AdvectionFluxes
   *
   * \tparam ModelImp type of the analytical model
   * \tparam FluxParameterImp type of the flux parameters
   */
  template <class ModelImp, class FluxParameterImp >
  class DGAdvectionFluxBase
  {
  public:
    typedef ModelImp  ModelType;

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
    DGAdvectionFluxBase (const ModelType& mod,
                         const Dune::Fem::ParameterReader& parameter = Dune::Fem::Parameter::container() )
      : model_(mod),
        param_( parameter )
    {}

    /**
     * \brief Constructor
     *
     * \param[in] mod analytical model
     * \param[in] parameters  parameter reader
     */
    DGAdvectionFluxBase (const ModelType& mod,
                         const ParameterType& parameter )
      : model_(mod),
        param_( parameter )
    {}

    /**
     * \brief Returns the name of the flux.
     */
    static std::string name () { return "AdvectionFlux"; }

    /**
     * \brief Returns the analytical model.
     */
    const ModelType& model () const { return model_; }

    /**
     * \brief Returns the parameters of the flux.
     */
    const ParameterType& parameter () const { return param_; }

    /**
     * \brief Evaluates the numerical flux \f$g(u^+,u^-)\f$.
     *
     * \param[in] left local evaluation of inside entity \f$ E^+ \f$
     * \param[in] right local evaluation of outside entity \f$ E^- \f$
     * \param[in] uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[in] uRight evaluation of the local function, i.e. \f$ u_{E^-}( \hat{x} ) \f$
     * \param[in] jacLeft evaluation of the gradient of the local function, i.e. \f$ \nabla  u_{E^+}( \hat{x} ) \f$
     * \param[in] jacRight evaluation of the gradient of the local function, i.e. \f$ \nabla u_{E^-}( \hat{x} ) \f$
     * \param[out] gLeft the result of the numerical flux \f$ g(u^+,u^-) \f$
     * \param[out] gRight the result of the numerical flux \f$ g(u^+,u^-) \f$
     *
     * \return maximum wavespeed * normal
     */
    template <class LocalEvaluation>
    inline double numericalFlux( const LocalEvaluation& left,
                                 const LocalEvaluation& right,
                                 const RangeType& uLeft,
                                 const RangeType& uRight,
                                 const JacobianRangeType& jacLeft,
                                 const JacobianRangeType& jacRight,
                                 RangeType& gLeft,
                                 RangeType& gRight ) const
    {

      assert( false );
      abort();
      return 0.0;
    }
  protected:
    // NOTE: model is stored as object for various reasons
    // such as thread safety and model states
    ModelType model_;
    const ParameterType param_;
  };


}
}
#endif
