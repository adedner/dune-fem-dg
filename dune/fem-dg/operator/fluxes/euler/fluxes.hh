#ifndef DUNE_FEM_DG_EULER_FLUXES_HH
#define DUNE_FEM_DG_EULER_FLUXES_HH

// system includes
#include <string>
#include <cmath>

#include "../advection/fluxbase.hh"
#include "../advection/fluxes.hh"
#include "eulerflux_impl.hh"
#include "llfadv.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief class specialization for the local Lax-Friedrichs flux.
   *
   * The purpose of this class is to allow the selection of an Euler flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class Model >
  class DGAdvectionFlux< Model, AdvectionFlux::Enum::euler_llf >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::LLF > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::LLF > > BaseType ;
  public:
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
    static std::string name () { return "LLF (Dennis)"; }
  };

  /**
   * \brief class specialization for the HLL flux.
   *
   * The purpose of this class is to allow the selection of an Euler flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class Model >
  class DGAdvectionFlux< Model, AdvectionFlux::Enum::euler_hll >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL > > BaseType ;
  public:
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
    static std::string name () { return "HLL (Dennis)"; }
  };

  /**
   * \brief class specialization for the HLLC flux.
   *
   * The purpose of this class is to allow the selection of an Euler flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class Model >
  class DGAdvectionFlux< Model, AdvectionFlux::Enum::euler_hll_pt >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL_PT > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLL_PT > > BaseType ;
  public:
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
    static std::string name () { return "HLL_PT (Dennis)"; }
  };

  /**
   * \brief class specialization for the HLLC flux.
   *
   * The purpose of this class is to allow the selection of an Euler flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class Model >
  class DGAdvectionFlux< Model, AdvectionFlux::Enum::euler_hllc >
    : public EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLLC > >
  {
    typedef EulerFluxImpl< Model, EulerNumFlux::EulerFlux<Model,EulerNumFlux::EulerFluxType::HLLC > > BaseType ;
  public:
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
    static std::string name () { return "HLLC (Dennis)"; }
  };

  /**
   * \brief class specialization for a general flux chosen by a parameter file.
   *
   * The purpose of this class is to allow the selection of an Euler flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::euler_general >
    : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >   BaseType;

    static const int dimRange = ModelImp::dimRange;
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;

  public:
    typedef AdvectionFlux::Enum                   IdEnum;
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    /**
     * \copydoc DGAdvectionFluxBase::DGAdvectionFluxBase()
     */
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... ),
        method_( this->parameter().getMethod() ),
        flux_llf_(this->model() ),
        flux_hll_( this->model() ),
        flux_hllc_( this->model() )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () { return "AdvectionFlux - Euler (via parameter file)"; }

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
      if( IdEnum::euler_llf == method_ )
      {
        return flux_llf_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      else if ( IdEnum::euler_hll == method_ )
      {
        return flux_hll_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      else if ( IdEnum::euler_hllc == method_ )
      {
        return flux_hllc_.numericalFlux( left, right, uLeft, uRight, jacLeft, jacRight, gLeft, gRight );
      }
      else
      {
        std::cerr << "Error: Advection flux not chosen via parameter file" << std::endl;
        assert( false );
        std::abort();
      }
      return 0.0;
    }

  private:
    const IdEnum                                    method_;
    DGAdvectionFlux< ModelImp, IdEnum::euler_llf >  flux_llf_;
    DGAdvectionFlux< ModelImp, IdEnum::euler_hll >  flux_hll_;
    DGAdvectionFlux< ModelImp, IdEnum::euler_hllc > flux_hllc_;

  };

}
}

#endif // file declaration
