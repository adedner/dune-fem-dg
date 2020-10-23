#ifndef FEMDG_UPWIND_FLUX_HH
#define FEMDG_UPWIND_FLUX_HH

#include <cmath>
#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Defines an advective flux using an upwind scheme.
   *
   * \ingroup AdvectionFluxes
   *
   * The numerical upwind flux \f$ g \f$ is defined by
   * \f[ g(u^+,u^-) = \begin{cases} F(u^+)\cdot n & F(u^+)\cdot n >0 \\ F(u^-)\cdot n & \text{otherwise} \end{cases}. \f]
   *
   * The function \f$ F \f$ is given by the analytical model.
   */
  template <class ModelImp>
  class UpwindFlux
    : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
                                                  BaseType;

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
    UpwindFlux( const ModelType& mod, const ParameterType& parameter = ParameterType() )
      : BaseType( mod, parameter )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () { return "Upwind"; }

    /**
     * The numerical upwind flux \f$ g \f$ is defined by
     * \f[ g(u^+,u^-) = \begin{cases} F(u^+)\cdot n & F'(u^+)\cdot n >0 \\ F(u^-)\cdot n & \text{otherwise} \end{cases}. \f]
     *
     * \copydoc DGAdvectionFluxBase::numericalFlux()
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
      const DomainType v = model_.velocity( left, uLeft );
      const auto upwind = normal * v;

      if (upwind > 0)
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
