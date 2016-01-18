#ifndef FEMDG_UPWIND_FLUX_HH
#define FEMDG_UPWIND_FLUX_HH

#include <cmath>
#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief defines the advective flux using an upwind scheme
   *
   * \ingroup AdvectionFluxes
   */
  template <class ModelImp>
  class UpwindFlux
    : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters< AdvectionFlux::Enum::upwind > >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters< AdvectionFlux::Enum::upwind > >
                                                  BaseType;

    typedef typename ModelImp::Traits             Traits;
    enum { dimRange = ModelImp::dimRange };
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;
    typedef typename ModelImp::EntityType         EntityType;
    typedef typename ModelImp::IntersectionType   IntersectionType;

  public:
    typedef typename BaseType::IdType             IdType;
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    using BaseType::model_;
    /**
     * \brief Constructor
     */
    UpwindFlux( const ModelType& mod, const ParameterType& parameter = ParameterType() )
      : BaseType( mod, parameter )
    {}

    static std::string name () { return "Upwind"; }

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
      const FaceDomainType& x = left.localPosition();

      // get normal from intersection
      const DomainType normal = left.intersection().integrationOuterNormal(x);

      // get velocity
      const DomainType v = model_.velocity( left );
      const double upwind = normal * v;

      if (upwind>0)
        gLeft = uLeft;
      else
        gLeft = uRight;
      gLeft *= upwind;
      gRight = gLeft;
      return std::abs(upwind);
    }
  };



}
}
#endif
