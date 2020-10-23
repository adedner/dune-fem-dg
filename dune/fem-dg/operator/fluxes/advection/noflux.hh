#ifndef DUNE_FEMDG_NOFLUX_HH
#define DUNE_FEMDG_NOFLUX_HH

#include <string>
#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{
  /**
   * \brief Advective flux returning zero flux.
   *
   * \ingroup AdvectionFluxes
   */
  template <class ModelImp>
  class NoFlux
    : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
                                                  BaseType;

    typedef typename ModelImp::Traits             Traits;
    static const int dimRange = ModelImp::dimRange;
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;
  public:
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;
    using BaseType::model_;

    /**
     * \copydoc DGAdvectionFluxBase::DGAdvectionFluxBase()
     */
    NoFlux( const ModelType& mod, const ParameterType& parameter )
      : BaseType( mod, parameter )
    {}

    static std::string name () { return "NoFlux"; }

    /**
     * The numerical upwind flux \f$ g \f$ is defined by
     * \f[ g(u^+,u^-) =  0. \f]
     *
     * This means, we do not have any flux through intersections.
     *
     * \copydoc DGAdvectionFluxBase::numericalFlux()
     */
    template< class LocalEvaluation >
    inline double
    numericalFlux( const LocalEvaluation& left,
                   const LocalEvaluation& right,
                   const RangeType& uLeft,
                   const RangeType& uRight,
                   const JacobianRangeType& jacLeft,
                   const JacobianRangeType& jacRight,
                   RangeType& gLeft,
                   RangeType& gRight) const
    {
      gLeft  = 0;
      gRight = 0;
      return 0;
    }
  };


}
}
#endif
