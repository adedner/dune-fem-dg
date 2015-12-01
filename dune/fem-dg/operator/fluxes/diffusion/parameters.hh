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

  namespace PrimalDiffusionFlux
  {
    enum class Enum
    {
      cdg2,      // CDG 2 (Compact Discontinuous Galerkin 2)
      cdg,       // CDG (Compact Discontinuous Galerkin)
      br2,       // BR2 (Bassi-Rebay 2)
      ip,        // IP (Interior Penalty)
      nipg,      // NIPG (Non-symmetric Interior  Penalty)
      bo,        // BO (Baumann-Oden)
      general,   // general means all methods chosen via parameter file
      none       // no diffusion (advection only)
    };

    //parameters which can be chosen in parameter file
    const Enum        _enums[] = { Enum::cdg2, Enum::cdg, Enum::br2, Enum::ip, Enum::nipg, Enum::bo };
    const std::string _strings[] = { "CDG2", "CDG" , "BR2", "IP" , "NIPG", "BO" };
    static const int  _size = 6;

    //helper class for static parameter selection
    template< Enum id = Enum::general >
    struct Identifier
    {
      typedef Enum type;
      static const type value = id;
    };
  }

  namespace PrimalDiffusionLifting
  {
    enum class Enum
    {
      id_id,  // int_Omega r([u]).tau  = -int_e [u].{tau}
      id_A,   // int_Omega r([u]).tau  = -int_e [u].{Atau}
      A_A     // int_Omega r([u]).Atau = -int_e [u].{Atau}
    };

    //parameters which can be chosen in parameter file
    const Enum        _enums[] = { Enum::id_id, Enum::id_A, Enum::A_A };
    const std::string _strings[] = { "id_id", "id_A" , "A_A" };
    static const int  _size = 3;

    //helper class for static parameter selection
    struct Identifier
    {
      typedef Enum type;
    };
  }




  template< PrimalDiffusionFlux::Enum id = PrimalDiffusionFlux::Enum::general >
  class DGPrimalDiffusionFluxParameters
    : public Fem::LocalParameter< DGPrimalDiffusionFluxParameters<id>, DGPrimalDiffusionFluxParameters<id> >
  {
  public:
    typedef typename PrimalDiffusionFlux::Identifier<id> IdType;
    typedef typename PrimalDiffusionLifting::Identifier  LiftingType;
  private:
    typedef typename IdType::type                        IdEnum;
    typedef typename LiftingType::type                   LiftingEnum;
  public:

    DGPrimalDiffusionFluxParameters( const std::string keyPrefix = "dgdiffusionflux." )
      : keyPrefix_( keyPrefix )
    {}

    static std::string methodNames( const IdEnum& mthd )
    {
      for( int i = 0; i < PrimalDiffusionFlux::_size; i++)
        if( PrimalDiffusionFlux::_enums[i] == mthd )
          return PrimalDiffusionFlux::_strings[i];
      assert( false );
      return "invalid diffusion flux";
    }

    virtual IdEnum getMethod() const
    {
      const int i = Fem::Parameter::getEnum( keyPrefix_ + "method", PrimalDiffusionFlux::_strings );
      return PrimalDiffusionFlux::_enums[i];
    }

    static std::string liftingNames( const LiftingEnum mthd )
    {
      for( int i = 0; i < PrimalDiffusionLifting::_size; i++)
        if( PrimalDiffusionLifting::_enums[i] == mthd )
          return PrimalDiffusionLifting::_strings[i];
      assert( false );
      return "invalid identifier";
    }

    virtual LiftingEnum getLifting() const
    {
      const int i = Fem::Parameter::getEnum( keyPrefix_ + "lifting", PrimalDiffusionLifting::_strings, 0 );
      return PrimalDiffusionLifting::_enums[i];
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
  private:

    const std::string keyPrefix_;

  };


} // end namespace
} // end namespace
#endif
