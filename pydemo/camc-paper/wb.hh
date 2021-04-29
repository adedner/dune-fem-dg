#pragma once
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem-dg/operator/fluxes/advection/python.hh>

template<class ModelImp>
struct WB : public Dune::Fem::DGAdvectionFluxPythonUserDefine< ModelImp >
{
  typedef Dune::Fem::DGAdvectionFluxPythonUserDefine< ModelImp > BaseType;
  typedef typename BaseType::RangeType          RangeType;
  typedef typename BaseType::JacobianRangeType  JacobianRangeType;
  typedef typename BaseType::FluxRangeType      FluxRangeType;
  typedef typename BaseType::DomainType         DomainType;

  WB (const ModelImp& model, double g)
  : BaseType(model), g_(g) {}

  template< class LocalContext >
  double
  numericalFlux( const LocalContext& left,         const LocalContext& right,
                 const RangeType& uLeft,           const RangeType& uRight,
                 const JacobianRangeType& jacLeft, const JacobianRangeType& jacRight,
                 RangeType& gLeft,                 RangeType& gRight) const
  {
    double speedL, speedR;
    const auto& x = left.localPosition();
    auto normal = left.intersection().integrationOuterNormal(x);
    const double faceArea = normal.two_norm();
    normal *= 1./faceArea;

    double zLeft = uLeft[3];
    double zRight = uRight[3];
    double zStar  = std::max(zLeft,zRight);
    double hLeft  = std::max(0.,uLeft[0]-zStar);
    double hRight = std::max(0.,uRight[0]-zStar);
    {
      model().setEntity(left.entity());
      // modify states (left)
      double z        = zStar - std::max(0.,zStar-uLeft[0]);
      double etaLeft  = hLeft + z;
      double etaRight = hRight + z;
      RangeType vLeft{etaLeft, hLeft/(etaLeft-zLeft)*uLeft[1],
                               hLeft/(etaLeft-zLeft)*uLeft[2],
                       z};
      RangeType vRight{etaRight, hRight/(etaRight-zRight)*uRight[1],
                                 hRight/(etaRight-zRight)*uRight[2],
                       z};
      // compute LLF flux
      speedL = llf(left,vLeft,vRight,jacLeft,jacRight,normal, gLeft);
      // add correction
      gLeft[1] += normal[0] * g_ * etaLeft * (z - zLeft);
      gLeft[2] += normal[1] * g_ * etaLeft * (z - zLeft);
    }
    model().setEntity(right.entity());
    {
      // modify states (right)
      double z        = zStar - std::max(0.,zStar-uRight[0]);
      double etaLeft  = hLeft + z;
      double etaRight = hRight + z;
      RangeType vLeft{etaLeft, hLeft/(etaLeft-zLeft)*uLeft[1],
                               hLeft/(etaLeft-zLeft)*uLeft[2],
                      z};
      RangeType vRight{etaRight, hRight/(etaRight-zRight)*uRight[1],
                                 hRight/(etaRight-zRight)*uRight[2],
                       z};
      // compute LLF flux
      speedR = llf(right,vLeft,vRight,jacLeft,jacRight,normal, gRight);
      // add correction
      gRight[1] += normal[0] * g_ * etaRight * (z - zRight);
      gRight[2] += normal[1] * g_ * etaRight * (z - zRight);
    }
    gLeft  *= faceArea;
    gRight *= faceArea;
    model().setEntity(left.entity());
    return std::max(speedL,speedR)*faceArea;
  }
  private:
  template< class LocalEvaluation >
  double
  llf( const LocalEvaluation& eval,
       const RangeType& uLeft,           const RangeType& uRight,
       const JacobianRangeType& jacLeft, const JacobianRangeType& jacRight,
       const DomainType &normal, RangeType& flux ) const
  {
    // LLF part applied to modified height values
    FluxRangeType anaflux;
    model().advection( eval, uLeft, jacLeft, anaflux );
    anaflux.mv( normal, flux );
    model().advection( eval, uRight, jacRight, anaflux );
    anaflux.umv( normal, flux );

    double maxspeedl, maxspeedr;
    double viscparal, viscparar;
    model().maxWaveSpeed( eval, normal, uLeft,  viscparal, maxspeedl );
    model().maxWaveSpeed( eval, normal, uRight, viscparar, maxspeedr );

    const double maxspeed = std::max( maxspeedl, maxspeedr);
    const double viscpara = std::max( viscparal, viscparar);

    RangeType visc( uRight );
    visc -= uLeft;
    visc *= viscpara;
    flux -= visc;

    flux *= 0.5;
    return maxspeed;
  }

  using BaseType::model;
  double g_;
};
