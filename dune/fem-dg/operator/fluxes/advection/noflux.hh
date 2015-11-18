#ifndef DUNE_FEMDG_NOFLUX_HH
#define DUNE_FEMDG_NOFLUX_HH

#include <string>
#include "../staticparameter.hh"
#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{
  /**
   * \brief advective flux returning zero
   *
   * \ingroup AdvectionFluxes
   */
  template <class ModelImp>
  class NoFlux
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
    static const typename MethodType::id method = MethodType::none;
    typedef typename BaseType::ModelType          ModelType;
    typedef StaticParameter< typename BaseType::ParameterType, method >
                                                  ParameterType;
    using BaseType::model_;

    /**
     * \brief constructor
     */
    NoFlux( const ModelType& mod, const ParameterType& parameter = ParameterType() )
      : BaseType( mod, method, parameter )
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
  };


}
}
#endif
