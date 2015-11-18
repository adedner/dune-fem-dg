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

  //template< class ModelImp, int >
  template< class ModelImp, AdvectionFluxIdentifier::id >
  class DGAdvectionFlux;


  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFluxIdentifier::llf >
    : public LLFAdvFlux< ModelImp >
  {
    typedef LLFAdvFlux< ModelImp > BaseType;
  public:
    typedef typename BaseType::ParameterType      ParameterType;
    typedef typename ParameterType::MethodType    MethodType;
    typedef typename BaseType::ModelType          ModelType;

    DGAdvectionFlux( const ModelType& mod,
                     const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
  };


  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFluxIdentifier::none >
    : public NoFlux< ModelImp >
  {
    typedef NoFlux< ModelImp > BaseType;
  public:
    typedef typename BaseType::ParameterType      ParameterType;
    typedef typename ParameterType::MethodType    MethodType;
    typedef typename BaseType::ModelType          ModelType;

    DGAdvectionFlux( const ModelType& mod,
                     const ParameterType& parameters = ParameterType() )
      : BaseType( mod, parameters )
    {}
  };


  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFluxIdentifier::upwind >
    : public UpwindFlux< ModelImp >
  {
    typedef UpwindFlux< ModelImp > BaseType;
  public:
    typedef typename BaseType::ParameterType      ParameterType;
    typedef typename ParameterType::MethodType    MethodType;
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
  class DGAdvectionFlux< ModelImp, AdvectionFluxIdentifier::general >
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
    static const typename MethodType::id method = MethodType::general;
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    /**
     * \brief Constructor
     */
    DGAdvectionFlux (const ModelType& mod,
                     const ParameterType& parameters = ParameterType())
      : model_( mod ),
        param_( parameters ),
        method_( param_.getMethod() ),
        flux_none_( mod ),
        flux_llf_( mod ),
        flux_upwind_( mod )
    {}
    static std::string name () { return "AdvectionFlux (via parameter file)"; }

    const ModelType& model() const { return model_; }

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
        case MethodType::none:
          return flux_none_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case MethodType::llf:
          return flux_llf_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
        case MethodType::upwind_:
          return flux_upwind_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      assert( false );
      return 0.0;
    }

  private:
    const ModelType&                                 model_;
    const ParameterType&                             param_;
    typename MethodType::id                          method_;
    DGAdvectionFlux< ModelType, MethodType::none >   flux_none_;
    DGAdvectionFlux< ModelType, MethodType::llf >    flux_llf_;
    DGAdvectionFlux< ModelType, MethodType::upwind > flux_upwind_;

  };


}
}
#endif
