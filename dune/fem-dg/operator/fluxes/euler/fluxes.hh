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
  template< class Model, class AdvectionFluxIdentifierImp >
  class EulerFlux;

  template< class Model >
  class EulerFlux< Model, Euler::AdvectionFlux::Identifier< Euler::AdvectionFlux::Enum::llf > >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::LLF > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::LLF > > BaseType ;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::IdType::type           IdEnum;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
    static std::string name () { return "LLF (Dennis)"; }
  };

  template< class Model >
  class EulerFlux< Model, Euler::AdvectionFlux::Identifier< Euler::AdvectionFlux::Enum::hll > >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL > > BaseType ;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::IdType::type           IdEnum;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
    static std::string name () { return "HLL (Dennis)"; }
  };

  template< class Model >
  class EulerFlux< Model, Euler::AdvectionFlux::Identifier< Euler::AdvectionFlux::Enum::hllc > >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLLC > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLLC > > BaseType ;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::IdType::type           IdEnum;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
    static std::string name () { return "HLLC (Dennis)"; }
  };

  template< class Model >
  class EulerFlux< Model, Euler::AdvectionFlux::Identifier< Euler::AdvectionFlux::Enum::llf2 > >
    : public LLFFlux< Model >
  {
    typedef LLFFlux< Model >  BaseType;
  public:
    typedef typename BaseType::ParameterType          ParameterType;
    typedef typename BaseType::IdType::type           IdEnum;
    typedef typename BaseType::ModelType              ModelType;

    EulerFlux( const Model& mod, const ParameterType& param = ParameterType() )
      : BaseType( mod, param )
    {}
    static std::string name () { return "LLF"; }
  };

  template< class ModelImp >
  class EulerFlux< ModelImp, Euler::AdvectionFlux::Identifier< Euler::AdvectionFlux::Enum::general > >
    : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters< Euler::AdvectionFlux::Enum::general > >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters< Euler::AdvectionFlux::Enum::general > >   BaseType;

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
    typedef typename BaseType::IdType             IdType;
    typedef typename IdType::type                 IdEnum;
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    template< IdEnum ident >
    using GeneralIdType = typename ParameterType::template GeneralIdType< ident >;

    EulerFlux( const ModelType& mod,
               const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters ),
        method_( parameters.getMethod() ),
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
        case IdEnum::llf:
          return flux_llf_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case IdEnum::hll:
          return flux_hll_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case IdEnum::hllc:
          return flux_hllc_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case IdEnum::llf2:
          return flux_llf2_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      assert( false );
      std::cerr << "Error: Advection flux not chosen via parameter file" << std::endl;
      return 0.0;
    }

  private:
    const IdEnum&                                    method_;
    EulerFlux< ModelImp, GeneralIdType< IdEnum::llf > >  flux_llf_;
    EulerFlux< ModelImp, GeneralIdType< IdEnum::hll > > flux_hll_;
    EulerFlux< ModelImp, GeneralIdType< IdEnum::hllc > > flux_hllc_;
    EulerFlux< ModelImp, GeneralIdType< IdEnum::llf2 > > flux_llf2_;

  };

}
  //template< Euler::FluxIdentifier::value id >
  //struct Fluxes< Euler::FluxIdentifier< id > >
  //{
  //  template< class ModelImp >
  //  using Advection=Euler::EulerFlux< ModelImp, id>;
  //};



}
}

#endif // file declaration
