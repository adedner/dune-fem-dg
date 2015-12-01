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
   * \brief defines the advective flux
   *
   * \ingroup AdvectionFluxes
   */
  template <class ModelImp, class FluxParameterImp >
  class DGAdvectionFluxBase
  {
    typedef typename ModelImp::Traits            Traits;

    enum { dimRange = ModelImp::dimRange };
    typedef typename ModelImp::DomainType        DomainType;
    typedef typename ModelImp::RangeType         RangeType;
    typedef typename ModelImp::JacobianRangeType JacobianRangeType;
    typedef typename ModelImp::FluxRangeType     FluxRangeType;
    typedef typename ModelImp::FaceDomainType    FaceDomainType;
    typedef typename ModelImp::EntityType        EntityType;
    typedef typename ModelImp::IntersectionType  IntersectionType;

  public:
    typedef ModelImp                               ModelType;
    typedef FluxParameterImp                       ParameterType;
    typedef typename ParameterType::IdType         IdType;
    typedef typename IdType::type                  IdEnum;


    /**
     * \brief Constructor
     */
    DGAdvectionFluxBase (const ModelType& mod,
                         const ParameterType& parameters = ParameterType() )
      : model_(mod),
        param_( parameters )
    {}

    static std::string name () { return "AdvectionFlux"; }

    const ModelType& model () const { return model_; }

    const ParameterType& parameter () const { return param_; }

    /**
     * \brief evaluates the flux \f$g(u,v)\f$
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
    const ModelType& model_;
    const ParameterType& param_;
  };


}
}
#endif
