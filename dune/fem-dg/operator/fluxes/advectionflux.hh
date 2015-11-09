#ifndef FEMDG_ADVECTION_FLUX_HH
#define FEMDG_ADVECTION_FLUX_HH

#include <string>
#include <assert.h>

namespace Dune
{

  /**
   * \defgroup Fluxes Analytical and Numerical Fluxes
   *
   * A numerical flux is given by...
   *
   * A the moment fluxes can be separated into to categories:
   * An advective flux or a diffusive flux.
   *
   */


  /**
   * \defgroup AdvectionFluxes Advection Fluxes
   *
   * \ingroup Fluxes
   *
   * An advective numerical flux is given by...
   *
   *
   */


  /**
   * \brief defines the advective flux using an upwind scheme
   *
   * \ingroup AdvectionFluxes
   */
  template <class Model>
  class DGAdvectionFluxBase
  {
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

  public:
    /**
     * \brief Constructor
     */
    DGAdvectionFluxBase(const ModelType& mod) : model_(mod) {}

    static std::string name () { return "AdvectionFlux"; }
    const ModelType& model() const {return model_;}

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
      assert( false );
      abort();
      return 0.0;
    }
  protected:
    const ModelType& model_;
  };
}
#endif
