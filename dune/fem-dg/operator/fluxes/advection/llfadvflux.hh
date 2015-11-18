#ifndef FEMDG_LLFADVFLUX_FLUX_HH
#define FEMDG_LLFADVFLUX_FLUX_HH

#include <string>
#include "../staticparameter.hh"
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
    : public DGAdvectionFluxBase< ModelImp >
  {
    typedef DGAdvectionFluxBase< ModelImp >       BaseType;

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
    typedef typename BaseType::MethodType         MethodType;
    static const typename MethodType::id method = MethodType::llf;
    typedef typename BaseType::ModelType          ModelType;
    typedef StaticParameter< typename BaseType::ParameterType, method >
                                                  ParameterType;
    using BaseType::model_;

    LLFAdvFlux( const ModelType& mod, const ParameterType& parameter = ParameterType() )
      : BaseType( mod, method, parameter )
    {}

    static std::string name () { return "Local Lax-Friedrichs"; }

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
      const FaceDomainType& x = left.localPoint();
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

      model_.maxSpeed( left,  normal, uLeft,  viscparal, maxspeedl );
      model_.maxSpeed( right, normal, uRight, viscparar, maxspeedr );

      maxspeed = (maxspeedl > maxspeedr) ? maxspeedl : maxspeedr;
      viscpara = (viscparal > viscparar) ? viscparal : viscparar;

      visc = uRight;
      visc -= uLeft;
      visc *= viscpara;
      gLeft -= visc;

      gLeft *= 0.5*len;
      gRight = gLeft;

      return maxspeed;
    }
  };

}
}
#endif
