#ifndef FEMDG_UPWIND_FLUX_HH
#define FEMDG_UPWIND_FLUX_HH

#include <string>
#include <cmath>

namespace Dune
{
  /**
   * @brief defines the advective flux
   */
  template <class ModelType>
  class UpwindFlux {
  public:
    typedef ModelType Model;
    typedef typename Model::Traits Traits;
    enum { dimRange = Model::dimRange };
    typedef typename Model :: DomainType DomainType;
    typedef typename Model :: RangeType RangeType;
    typedef typename Model :: JacobianRangeType JacobianRangeType;
    typedef typename Model :: FluxRangeType FluxRangeType;
    typedef typename Model :: DiffusionRangeType DiffusionRangeType;
    typedef typename Model :: FaceDomainType  FaceDomainType;
    typedef typename Model :: EntityType  EntityType;
    typedef typename Model :: IntersectionType  IntersectionType;

  public:
    /**
     * @brief Constructor
     */
    UpwindFlux(const Model& mod) : model_(mod) {}

    static std::string name () { return "UpwindFlux"; }

    const Model& model() const {return model_;}

    /**
     * @brief evaluates the flux \f$g(u,v)\f$
     *
     * @return maximum wavespeed * normal
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
