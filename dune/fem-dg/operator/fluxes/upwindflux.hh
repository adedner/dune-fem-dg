#ifndef FEMDG_UPWIND_FLUX_HH
#define FEMDG_UPWIND_FLUX_HH

#include <cmath>
#include <dune/fem-dg/operator/fluxes/advectionflux.hh>

namespace Dune
{
  /**
   * \brief defines the advective flux using an upwind scheme
   *
   * \ingroup AdvectionFluxes
   */
  template <class Model>
  class UpwindFlux
    : public DGAdvectionFluxBase< Model >
  {
    typedef DGAdvectionFluxBase< Model >               BaseType;
  public:
    typedef Model                                  ModelType;
    typedef typename ModelType::Traits             Traits;
    enum { dimRange = ModelType::dimRange };
    typedef typename ModelType::DomainType         DomainType;
    typedef typename ModelType::RangeType          RangeType;
    typedef typename ModelType::JacobianRangeType  JacobianRangeType;
    typedef typename ModelType::FluxRangeType      FluxRangeType;
    typedef typename ModelType::DiffusionRangeType DiffusionRangeType;
    typedef typename ModelType::FaceDomainType     FaceDomainType;
    typedef typename ModelType::EntityType         EntityType;
    typedef typename ModelType::IntersectionType   IntersectionType;

  public:
    /**
     * \brief Constructor
     */
    UpwindFlux( const ModelType& mod )
      : BaseType( mod ),
        model_(mod)
    {}

    static std::string name () { return "UpwindFlux"; }

    const Model& model() const {return model_;}

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
      const FaceDomainType& x = left.localPoint();

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
  protected:
    const Model& model_;
  };
}
#endif
