#ifndef DUNE_FEMDG_NOFLUX_HH
#define DUNE_FEMDG_NOFLUX_HH

#include <string>
#include <dune/fem-dg/operator/fluxes/advectionflux.hh>

namespace Dune
{

/**
 * \brief advective flux returning zero
 *
 * \ingroup AdvectionFluxes
 */
template <class Model>
class NoFlux
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
public:
  /**
   * \brief constructor
   */
  NoFlux(const ModelType& mod)
    : BaseType( mod ),
      model_(mod)
  {}

  static std::string name () { return "NoFlux"; }

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

protected:
  const ModelType& model_;
};

}
#endif
