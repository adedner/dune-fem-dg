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

  /**
   * \brief The purpose of this class is to allow the selection of an advection flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp, AdvectionFlux::Enum id >
  class DGAdvectionFlux;


  /**
   * \brief class specialization for the local Lax-Friedrichs flux.
   *
   * The purpose of this class is to allow the selection of an advection flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::llf >
    : public LLFAdvFlux< ModelImp >
  {
    typedef LLFAdvFlux< ModelImp > BaseType;
  public:
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };


  /**
   * \brief class specialization for no advection flux.
   *
   * The purpose of this class is to allow the selection of an advection flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::none >
    : public NoFlux< ModelImp >
  {
    typedef NoFlux< ModelImp > BaseType;
  public:
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };


  /**
   * \brief class specialization for upwind flux.
   *
   * The purpose of this class is to allow the selection of an advection flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::upwind >
    : public UpwindFlux< ModelImp >
  {
    typedef UpwindFlux< ModelImp > BaseType;
  public:
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };


  /**
   * \brief class specialization for a general flux chosen by a parameter file.
   *
   * The purpose of this class is to allow the selection of an advection flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::general >
   : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters  >
                                                  BaseType;

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
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... ),
        method_( this->parameter().getMethod() ),
        flux_none_( this->model(), this->parameter() ),
        flux_llf_( this->model(), this->parameter() ),
        flux_upwind_( this->model(), this->parameter() )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () { return "AdvectionFlux (via parameter file)"; }

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
      if( method_ == IdEnum::upwind )
      {
        return flux_upwind_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      else if( method_ == IdEnum::llf )
      {
        return flux_llf_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      else if( method_ == IdEnum::none )
      {
        return flux_none_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      else
      {
        std::cerr << "Error: Advection flux " << method_ << " not supported!" << std::endl;
        assert( false );
        std::abort();
      }
      return 0.0;
    }

  private:
    const IdEnum                                 method_;
    DGAdvectionFlux< ModelType, IdEnum::none >   flux_none_;
    DGAdvectionFlux< ModelType, IdEnum::llf >    flux_llf_;
    DGAdvectionFlux< ModelType, IdEnum::upwind > flux_upwind_;

  };



}
}
#endif
