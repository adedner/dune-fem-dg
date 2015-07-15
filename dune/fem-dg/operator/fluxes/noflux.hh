#ifndef DUNE_FEMDG_NOFLUX_HH
#define DUNE_FEMDG_NOFLUX_HH

#include <string>

namespace Dune
{

template <class ModelType>
class NoFlux {
public:
  typedef ModelType Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model :: DomainType DomainType;
  typedef typename Model :: RangeType RangeType;
  typedef typename Model :: FluxRangeType FluxRangeType;
  typedef typename Model :: FaceDomainType  FaceDomainType;
  typedef typename Model :: EntityType  EntityType;
  typedef typename Model :: IntersectionType  IntersectionType;
public:
  /**
   * @brief constructor
   */
  NoFlux(const Model& mod) : model_(mod) {}

  static std::string name () { return "UpwindFlux"; }

  const Model& model() const {return model_;}

  template< class LocalEvaluation >
  inline double
  numericalFlux( const LocalEvaluation& left,
                 const LocalEvaluation& right,
                 const RangeType& uLeft,
                 const RangeType& uRight,
                 RangeType& gLeft,
                 RangeType& gRight) const
  {
    gLeft  = 0;
    gRight = 0;
    return 0;
  }

  /**
   * @brief evaluates the flux \f$g(u,v)\f$
   *
   * @return maximum wavespeed * normal
   */
  template <class QuadratureImp>
  inline double numericalFlux( const IntersectionType& it,
                               const EntityType& inside,
                               const EntityType& outside,
                               const double time,
                               const QuadratureImp& faceQuadInner,
                               const QuadratureImp& faceQuadOuter,
                               const int quadPoint,
                               const RangeType& uLeft,
                               const RangeType& uRight,
                               RangeType& gLeft,
                               RangeType& gRight ) const
  {
    gLeft  = 0;
    gRight = 0;
    return 0;
  }
protected:
  const Model& model_;
};

}
#endif
