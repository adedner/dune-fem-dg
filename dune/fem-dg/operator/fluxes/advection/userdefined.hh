#ifndef FEMDG_ADVECTION_USERDEFINED_HH
#define FEMDG_ADVECTION_USERDEFINED_HH

#include <string>
#include <assert.h>

#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>


namespace Dune
{
namespace Fem
{

  /**
   * \brief class specialization for a general flux chosen by a parameter file.
   *
   * The purpose of this class is to allow the selection of an advection flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::userdefined >
   : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters  >
                                                  BaseType;

    typedef typename ModelImp::Traits             Traits;
    static const int dimRange = ModelImp::dimRange;
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;

  public:
    typedef typename BaseType::IdEnum             IdEnum;
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    /**
     * \brief Constructor
     */
    DGAdvectionFlux (const ModelType& mod,
                     const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () { return "AdvectionFlux (userdefined)"; }

    /**
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
      DUNE_THROW(NotImplemented,"DGAdvactionFlux< userdefined > was not implemented!");
      return 0.0;
    }

  };

} // end namespace Fem
} // end namespace Dune
#endif
