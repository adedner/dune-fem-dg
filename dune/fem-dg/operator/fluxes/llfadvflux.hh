#ifndef FEMDG_LLFADVFLUX_FLUX_HH
#define FEMDG_LLFADVFLUX_FLUX_HH


namespace Dune
{

  /**
   *  \brief advection flux using local Lax-Friedrichs
   *
   *  \ingroup Fluxes
   */
  template <class ModelType>
  class LLFAdvFlux
  {
  public:
    typedef ModelType Model;
    typedef typename Model::Traits Traits;
    enum { dimRange = Model::dimRange };
    typedef typename Model :: DomainType          DomainType;
    typedef typename Model :: RangeType           RangeType;
    typedef typename Model :: JacobianRangeType   JacobianRangeType;
    typedef typename Model :: FluxRangeType       FluxRangeType;
    typedef typename Model :: FaceDomainType      FaceDomainType;
    typedef typename Model :: EntityType          EntityType;
    typedef typename Model :: IntersectionType    IntersectionType;
    /**
     * @brief Constructor
     */
    LLFAdvFlux(const Model& mod) : model_(mod) {}

    static std::string name () { return "LaxFriedrichsFlux"; }

    const Model& model() const {return model_;}

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
    const Model& model_;
  };

}
#endif
