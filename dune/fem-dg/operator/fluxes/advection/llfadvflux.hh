#ifndef FEMDG_LLFADVFLUX_FLUX_HH
#define FEMDG_LLFADVFLUX_FLUX_HH

#include <string>
#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{

  /**
   *  \brief advection flux using local Lax-Friedrichs
   *
   *  \ingroup AdvectionFluxes
   */
  template <class ModelImp>
  class LLFAdvFlux
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
    LLFAdvFlux( const ModelType& mod, const ParameterType& parameter = ParameterType() )
      : BaseType( mod, parameter )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () { return "Local Lax-Friedrichs"; }

    /**
     * The numerical upwind flux \f$ g \f$ is defined by
     * \f[ g(u^+,u^-) =  \frac{1}{2}( F(u^+) + F(u^-) )\cdot n -
     * \frac{1}{2} \max(\lambda_\max(F'(u^+)\cdot n),\lambda_\max(F'(u^-))) (u^- - u^+). \f]
     *
     * where \f$ \lambda_\max \f$ denotes the largest eigenvalue
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
      DomainType normal = left.intersection().integrationOuterNormal(x);
      const auto faceArea = normal.two_norm();
      normal *= 1./faceArea;

      double maxspeedl, maxspeedr;
      double viscparal, viscparar;

      FluxRangeType anaflux;

      // assume that model is set to left entity
      model_.setEntity( left.entity() );

      model_.advection( left, uLeft, jacLeft, anaflux );
      // set gLeft
      anaflux.mv( normal, gLeft );

      model_.maxWaveSpeed( left,  normal, uLeft,  viscparal, maxspeedl );

      model_.setEntity( right.entity() );

      model_.advection( right, uRight, jacRight, anaflux );
      // add to gLeft
      anaflux.umv( normal, gLeft );

      model_.maxWaveSpeed( right, normal, uRight, viscparar, maxspeedr );

      const double maxspeed = std::max( maxspeedl, maxspeedr);
      const double viscpara = std::max( viscparal, viscparar);

      RangeType visc( uRight );
      visc -= uLeft;
      visc *= viscpara;
      gLeft -= visc;

      // multiply with 0.5 for averaging and with face area
      gLeft *= 0.5 * faceArea;
      // conservation property
      gRight = gLeft;

      // model is an internal object so we can avoid this
      //model_.setEntity( left.entity() );

      return maxspeed * faceArea;
    }
  };

}
}
#endif
