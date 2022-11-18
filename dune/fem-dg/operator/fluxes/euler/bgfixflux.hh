#ifndef DUNE_FEM_DG_EULER_BGFIX_FLUX_HH
#define DUNE_FEM_DG_EULER_BGFIX_FLUX_HH

// system includes
#include <string>
#include <cmath>

#include "../advection/fluxbase.hh"
#include "../advection/fluxes.hh"

namespace Dune
{
namespace Fem
{
  /**
   * \brief class specialization for a general flux chosen by a parameter file.
   *
   * The purpose of this class is to allow the selection of an Euler flux
   * via an enum given in AdvectionFlux::Enum.
   */
  template< class ModelImp, AdvectionFlux::Enum id >
  class DGAdvectionFluxBGFix : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >  BaseType;
    typedef DGAdvectionFlux< ModelImp, id >                           NumFluxImplType;

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

    using BaseType::model;

    /**
     * \copydoc DGAdvectionFluxBase::DGAdvectionFluxBase()
     */
    template< class ... Args>
    DGAdvectionFluxBGFix(  Args&&... args )
      : BaseType( std::forward<Args>(args)... ),
        numFlux_( std::forward<Args>(args)... )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () {
      return std::string("(BGFix) + FluxImpl");
    }

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
      // Assumption: background atmosphere is continuous
      auto bgLeft  = model().background( left.entity(), left.quadraturePoint() );
      auto bgRight = model().background( right.entity(), right.quadraturePoint() );

      RangeType uLeftTot( uLeft );
      RangeType uRightTot( uRight );

      uLeftTot  += bgLeft;
      uRightTot += bgRight;

      RangeType gBgLeft, gBgRight;

      // numerical flux
      const double ws = numFlux_.numericalFlux(left, right, uLeftTot, uRightTot,
          jacLeft, jacRight, gLeft, gRight );

      // background flux
      numFlux_.numericalFlux(left, right, bgLeft, bgRight, jacLeft, jacRight, gBgLeft, gBgRight );

      gLeft  -= gBgLeft;
      gRight -= gBgRight;

      return ws;
    }

  protected:
    NumFluxImplType numFlux_;
  };
}
}

#endif // file declaration
