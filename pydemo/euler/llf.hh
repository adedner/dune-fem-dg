#ifndef FEMDG_ADVECTION_USERDEFINED_HH
#define FEMDG_ADVECTION_USERDEFINED_HH

#include <string>

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
    enum { dimRange = ModelImp::dimRange };
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;
  protected:
    using BaseType::model_;

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
    static std::string name () { return "userdefined LLF"; }

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
      // User defined LLF
      const FaceDomainType& x = left.localPosition();
      DomainType normal = left.intersection().integrationOuterNormal(x);
      const auto faceArea = normal.two_norm();
      normal *= 1./faceArea;

      FluxRangeType anaflux;

      model_.advection( left, uLeft, jacLeft, anaflux );
      // set gLeft
      anaflux.mv( normal, gLeft );

      model_.advection( right, uRight, jacRight, anaflux );
      // add to gLeft
      anaflux.umv( normal, gLeft );

      double maxspeedl, maxspeedr;
      double viscparal, viscparar;

      model_.maxSpeed( left,  normal, uLeft,  viscparal, maxspeedl );
      model_.maxSpeed( right, normal, uRight, viscparar, maxspeedr );

      const double maxspeed = std::max( maxspeedl, maxspeedr);
      const double viscpara = std::max( viscparal, viscparar);

      RangeType visc( uRight );
      visc -= uLeft;
      visc *= viscpara;
      gLeft -= visc;

      // multiply with 0.5 for averaging and with face area
      gLeft *= 0.5 * faceArea;
      // conservation property
      gRight = gLeft;

      return maxspeed * faceArea;
    }

  };

} // end namespace Fem
} // end namespace Dune
#endif
