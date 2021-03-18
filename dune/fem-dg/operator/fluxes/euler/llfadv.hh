#ifndef DUNE_FEM_DG_EULERLLFADVECTIONFLUXES_HH
#define DUNE_FEM_DG_EULERLLFADVECTIONFLUXES_HH

// system includes
#include <string>
#include <cmath>

#include "parameters.hh"
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

  /**
    *  \brief Local Lax-Friedrichs flux for the euler problem.
    *
    *  \note It is possible to enable a well-belanced scheme via the
    *  preprocessor define `WELLBALANCE`
    *
    *  \ingroups AdvectionFluxes
    */
  template< class ModelImp >
  class EulerLLFFlux
   : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >   BaseType;

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
    EulerLLFFlux(const ModelType& mod, const ParameterType& param = ParameterType() )
      : BaseType( mod, param )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () { return "LLF"; }

   /**
     * The numerical upwind flux \f$ g \f$ is defined by
     * \f[ g(u^+,u^-) =  \frac{1}{2}( F(u^+) + F(u^-) )\cdot n -
     * \frac{1}{2} \max(\lambda_\max(F'(u^+)\cdot n),\lambda_\max(F'(u^-))) (u^- - u^+). \f]
     *
     * where \f$ \lambda_\max \f$ denotes the largest eigenvalue
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
      const FaceDomainType& x = left.localPosition();
      DomainType normal = left.intersection().integrationOuterNormal(x);
      const double len = normal.two_norm();
      normal *= 1./len;

      RangeType visc;
      FluxRangeType anaflux;

      model_.advection( left, uLeft, jacLeft, anaflux );

      // set gLeft
      anaflux.mv( normal, gLeft );

      model_.advection( right, uRight, jacRight, anaflux );
      anaflux.umv( normal, gLeft );

      double maxspeedl, maxspeedr, maxspeed;
      double viscparal, viscparar, viscpara;

      model_.maxWaveSpeed( left,  normal, uLeft,  viscparal, maxspeedl );
      model_.maxWaveSpeed( right, normal, uRight, viscparar, maxspeedr );

      maxspeed = (maxspeedl > maxspeedr) ? maxspeedl : maxspeedr;
      viscpara = (viscparal > viscparar) ? viscparal : viscparar;
      visc = uRight;
      visc -= uLeft;
      visc *= viscpara;
      gLeft -= visc;

      gLeft *= 0.5*len;
      gRight = gLeft;

  #if WELLBALANCE
      const DomainType xGlobal = left.intersection().geometry().global(x);
      const double g = model_.problem().g();

      // calculate geopotential in the grid elements sharing 'it'
      const double z = xGlobal[dimDomain-1];

      // calculate num. flux for pressure
      double pInside, TInside;
      double pOutside, TOutside;
      model_.pressAndTemp( uLeft, pInside, TInside );
      model_.pressAndTemp( uRight, pOutside, TOutside );
      const double pNumFlux = 0.5*(pInside + pOutside);

      const double rhoAvg = 0.5*(uLeft[0] + uRight[0]);
      const double rhoJump = uRight[0] - uLeft[0];

      // add well-balancing terms to vertical momentum flux, scaled with normal
      gLeft[dimDomain]  += 0.5*normal[dimDomain-1]*( -2.*(pNumFlux-pInside))*len;
      gRight[dimDomain] += 0.5*normal[dimDomain-1]*( -2.*(pNumFlux-pOutside))*len;

      gLeft[dimDomain]  -= 0.25*normal[dimDomain-1]*( g*z*rhoJump )*len;
      gRight[dimDomain] += 0.25*normal[dimDomain-1]*( g*z*rhoJump )*len;
  #endif

      return maxspeed * len;
    }
  };



}
}

#endif // file declaration
