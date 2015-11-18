#ifndef DUNE_FEM_DG_EULER_FLUXES_HH
#define DUNE_FEM_DG_EULER_FLUXES_HH

// system includes
#include <string>
#include <cmath>

#include "../advection/fluxbase.hh"
#include "eulerflux_impl.hh"
#include "llfadv.hh"

namespace Dune
{
namespace Fem
{
namespace Euler
{


  //template< class Model, int >
  template< class Model, Euler::FluxIdentifier::id >
  class EulerFlux;

  template< class Model >
  class EulerFlux< Model, Euler::FluxIdentifier::llf >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::LLF > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::LLF > > BaseType ;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::MethodType             MethodType;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod,  MethodType::llf, parameters )
    {}
    static std::string name () { return "LLF (Dennis)"; }
  };

  template< class Model >
  class EulerFlux< Model, Euler::FluxIdentifier::hll >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL > > BaseType ;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::MethodType             MethodType;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod, MethodType::hll, parameters )
    {}
    static std::string name () { return "HLL (Dennis)"; }
  };

  template< class Model >
  class EulerFlux< Model, Euler::FluxIdentifier::hllc >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLLC > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLLC > > BaseType ;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::MethodType             MethodType;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod, MethodType::hllc, parameters )
    {}
    static std::string name () { return "HLLC (Dennis)"; }
  };

  template< class Model >
  class EulerFlux< Model, Euler::FluxIdentifier::llf2 >
    : public LLFFlux< Model >
  {
    typedef LLFFlux< Model >  BaseType;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::MethodType             MethodType;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod, const ParameterType& param = ParameterType() )
      : BaseType( mod, param )
    {}
    static std::string name () { return "LLF"; }
  };

  template< class ModelImp >
  class EulerFlux< ModelImp, Euler::FluxIdentifier::general >
    : public DGAdvectionFluxBase< ModelImp, EulerFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, EulerFluxParameters >   BaseType;

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
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    using BaseType::model_;
    using BaseType::method_;

    EulerFlux( const ModelType& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters.getMethod(), parameters ),
        flux_llf_( mod ),
        flux_hll_( mod ),
        flux_hllc_( mod ),
        flux_llf2_( mod )
    {}

    static std::string name () { return "EulerFlux (via parameter file)"; }

    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
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
      switch (method_)
      {
        case MethodType::llf:
          return flux_llf_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case MethodType::hll:
          return flux_hll_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case MethodType::hllc:
          return flux_hllc_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case MethodType::llf2:
          return flux_llf2_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      assert( false );
      return 0.0;
    }

  private:
    EulerFlux< ModelImp, MethodType::llf >  flux_llf_;
    EulerFlux< ModelImp, MethodType::hll >  flux_hll_;
    EulerFlux< ModelImp, MethodType::hllc > flux_hllc_;
    EulerFlux< ModelImp, MethodType::llf2 > flux_llf2_;

  };

}
}
}

#endif // file declaration
