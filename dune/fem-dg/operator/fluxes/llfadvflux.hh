#ifndef FEMDG_LLFADVFLUX_FLUX_HH
#define FEMDG_LLFADVFLUX_FLUX_HH

#include <dune/fem-dg/operator/fluxes/advectionflux.hh>

namespace Dune
{

  /**
   *  \brief advection flux using local Lax-Friedrichs
   *
   *  \ingroup AdvectionFluxes
   */
  template <class Model>
  class LLFAdvFlux
    : public DGAdvectionFluxBase< Model >
  {
    typedef DGAdvectionFluxBase< Model >          BaseType;
  public:
    typedef Model                                 ModelType;
    typedef typename ModelType::Traits            Traits;
    enum { dimRange = ModelType::dimRange };
    typedef typename ModelType::DomainType        DomainType;
    typedef typename ModelType::RangeType         RangeType;
    typedef typename ModelType::JacobianRangeType JacobianRangeType;
    typedef typename ModelType::FluxRangeType     FluxRangeType;
    typedef typename ModelType::FaceDomainType    FaceDomainType;
    typedef typename ModelType::EntityType        EntityType;
    typedef typename ModelType::IntersectionType  IntersectionType;

    /**
     * \brief Constructor
     */
    LLFAdvFlux(const ModelType& mod) : model_(mod) {}

    static std::string name () { return "LaxFriedrichsFlux"; }

    const ModelType& model() const {return model_;}

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
  protected:
    const ModelType& model_;
  };

}
#endif
