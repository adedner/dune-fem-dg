#ifndef FEMDG_ADVECTION_FLUXES_HH
#define FEMDG_ADVECTION_FLUXES_HH

#include <string>
#include <assert.h>

#include "llfadvflux.hh"
#include "noflux.hh"
#include "upwindflux.hh"


namespace Dune
{
namespace Fem
{

  template< class ModelImp, class AdvectionFluxIdentifierImp >
  class DGAdvectionFlux;


  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Identifier< AdvectionFlux::Enum::llf > >
    : public LLFAdvFlux< ModelImp >
  {
    typedef LLFAdvFlux< ModelImp > BaseType;
  public:
    typedef typename BaseType::ParameterType      ParameterType;
    typedef typename BaseType::IdType::type       IdEnum;
    typedef typename BaseType::ModelType          ModelType;

    DGAdvectionFlux( const ModelType& mod,
                     const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
  };


  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Identifier< AdvectionFlux::Enum::none > >
    : public NoFlux< ModelImp >
  {
    typedef NoFlux< ModelImp > BaseType;
  public:
    typedef typename BaseType::ParameterType      ParameterType;
    typedef typename BaseType::IdType::type       IdEnum;
    typedef typename BaseType::ModelType          ModelType;

    DGAdvectionFlux( const ModelType& mod,
                     const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
  };


  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Identifier< AdvectionFlux::Enum::upwind > >
    : public UpwindFlux< ModelImp >
  {
    typedef UpwindFlux< ModelImp > BaseType;
  public:
    typedef typename BaseType::ParameterType      ParameterType;
    typedef typename BaseType::IdType::type       IdEnum;
    typedef typename BaseType::ModelType          ModelType;

    DGAdvectionFlux( const ModelType& mod,
                     const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
  };

  /**
   * \brief defines the advective flux using an upwind scheme
   *
   * \ingroup AdvectionFluxes
   */
  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Identifier< AdvectionFlux::Enum::general > >
   : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters< AdvectionFlux::Enum::general > >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters< AdvectionFlux::Enum::general >  >
                                                  BaseType;

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

    /**
     * \brief Constructor
     */
    DGAdvectionFlux (const ModelType& mod,
                     const ParameterType& parameters = ParameterType())
      : BaseType( mod, parameters ),
        method_( parameters.getMethod() ),
        flux_none_( mod ),
        flux_llf_( mod ),
        flux_upwind_( mod )
    {}
    static std::string name () { return "AdvectionFlux (via parameter file)"; }

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
        case IdEnum::none:
          return flux_none_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case IdEnum::llf:
          return flux_llf_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case IdEnum::upwind:
          return flux_upwind_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      std::cerr << "Error: Advection flux not chosen via parameter file" << std::endl;
      assert( false );
      return 0.0;
    }

  private:
    const IdEnum&                                    method_;
    DGAdvectionFlux< ModelType, GeneralIdType< IdEnum::none > >   flux_none_;
    DGAdvectionFlux< ModelType, GeneralIdType< IdEnum::llf > >    flux_llf_;
    DGAdvectionFlux< ModelType, GeneralIdType< IdEnum::upwind > > flux_upwind_;

  };


 // template< class FluxIdentifierImp >
 // struct Fluxes;

 // template< AdvectionFluxIdentifier::value id >
 // struct Fluxes< AdvectionFluxIdentifier< id > >
 // {
 //   template< class ModelImp >
 //   using Advection=DGAdvectionFlux< ModelImp, id>;
 // };




}
}
#endif
