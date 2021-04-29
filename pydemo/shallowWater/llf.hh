#ifndef FEMDG_ADVECTION_USERDEFINED_HH
#define FEMDG_ADVECTION_USERDEFINED_HH

#include <string>

#include <dune/fem-dg/operator/fluxes/advection/python.hh>

namespace Dune
{
  namespace Fem
  {
    template< class ModelImp >
    class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::userdefined >
    : public DGAdvectionFluxPythonUserDefine< ModelImp >
                // EmptyAdditional< typename ModelImp::DFunctionSpaceType> >
    {
      typedef DGAdvectionFluxPythonUserDefine< ModelImp >
        //, EmptyAdditional< typename ModelImp::DFunctionSpaceType> >
        BaseType;

    protected:
      using BaseType::model_;

    public:
      typedef typename BaseType::IdEnum             IdEnum;
      typedef typename BaseType::ModelType          ModelType;
      typedef typename BaseType::ParameterType      ParameterType;

      // ModelType != ModelImp
      enum { dimRange = ModelType::dimRange };
      typedef typename ModelType::DomainType         DomainType;
      typedef typename ModelType::RangeType          RangeType;
      typedef typename ModelType::JacobianRangeType  JacobianRangeType;
      typedef typename ModelType::FluxRangeType      FluxRangeType;
      typedef typename ModelType::FaceDomainType     FaceDomainType;

      /**
       * \brief Constructor
       */
      DGAdvectionFlux (const ModelImp& modelImp,
                       const ParameterType& parameters = ParameterType() )
        : BaseType( modelImp, parameters )
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

        model_.maxWaveSpeed( left,  normal, uLeft,  viscparal, maxspeedl );
        model_.maxWaveSpeed( right, normal, uRight, viscparar, maxspeedr );

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
  }
}
#endif
