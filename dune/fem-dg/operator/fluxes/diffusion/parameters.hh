#ifndef DUNE_FEM_DG_DIFFUSIONFLUXPARAMETER_HH
#define DUNE_FEM_DG_DIFFUSIONFLUXPARAMETER_HH

#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

  ///////////////////////////////////////////////////////////
  //
  //  Identifier for Diffusion Fluxes for Primal methods
  //
  //////////////////////////////////////////////////////////
  struct DGDiffusionFluxIdentifier
  {
    typedef enum
    {
      cdg2    = 0,  // CDG 2 (Compact Discontinuous Galerkin 2)
      cdg     = 1,  // CDG (Compact Discontinuous Galerkin)
      br2     = 2,  // BR2 (Bassi-Rebay 2)
      ip      = 3,  // IP (Interior Penalty)
      nipg    = 4,  // NIPG (Non-symmetric Interior  Penalty)
      bo      = 5,  // BO (Baumann-Oden)
      general = 6,  // general means all methods chosen via parameter file
      none    = 7   // no diffusion (advection only)
    } id;
  };

  struct DGLiftingFluxIdentifier
  {
    typedef enum
    {
      id_id    = 0,  // int_Omega r([u]).tau  = -int_e [u].{tau}
      id_A     = 1,  // int_Omega r([u]).tau  = -int_e [u].{Atau}
      A_A      = 2   // int_Omega r([u]).Atau = -int_e [u].{Atau}
    } id;
  };

  class DGPrimalFormulationParameters
    : public Fem::LocalParameter< DGPrimalFormulationParameters, DGPrimalFormulationParameters >
  {
    const std::string keyPrefix_;
  public:
    typedef DGDiffusionFluxIdentifier MethodType;
    typedef DGLiftingFluxIdentifier   LiftingType;

    DGPrimalFormulationParameters( const std::string keyPrefix = "dgdiffusionflux." )
      : keyPrefix_( keyPrefix )
    {}

    static std::string methodNames( const MethodType::id mthd )
    {
      const std::string method []
        = { "CDG2", "CDG" , "BR2", "IP" , "NIPG", "BO" };
      assert( mthd >= MethodType::cdg2 && mthd < MethodType::general );
      return method[ mthd ];
    }

    virtual MethodType::id getMethod() const
    {
      const std::string method []
        = { methodNames( MethodType::cdg2 ),
            methodNames( MethodType::cdg ),
            methodNames( MethodType::br2 ),
            methodNames( MethodType::ip ),
            methodNames( MethodType::nipg ),
            methodNames( MethodType::bo )
          };
      return (MethodType::id) Fem::Parameter::getEnum( keyPrefix_ + "method", method );
    }

    static std::string liftingNames( const LiftingType::id mthd )
    {
      const std::string method []
        = { "id_id", "id_A" , "A_A" };
      return method[ mthd ];
    }

    virtual LiftingType::id getLifting() const
    {
      const std::string method []
        = { liftingNames( LiftingType::id_id ),
            liftingNames( LiftingType::id_A ),
            liftingNames( LiftingType::A_A )
          };
      return (LiftingType::id) Fem::Parameter::getEnum( keyPrefix_ + "lifting", method, 0 );
    }

    virtual double penalty() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "penalty" );
    }

    virtual double liftfactor() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "liftfactor" );
    }

    virtual double theoryparameters() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "theoryparameters", 0. );
    }

    template <class DomainType>
    void upwind( DomainType& upwd ) const
    {
      Fem::Parameter::get(keyPrefix_ + "upwind", upwd, upwd);
    }
  };


} // end namespace
} // end namespace
#endif
