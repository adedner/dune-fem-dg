#ifndef DUNE_FEM_DG_DIFFUSIONFLUXPARAMETER_HH
#define DUNE_FEM_DG_DIFFUSIONFLUXPARAMETER_HH

#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Namespace containing all parameters to select a primal diffusion flux.
   */
  namespace PrimalDiffusionFlux
  {
    /**
     * \brief Enum of all known primal diffusion flux implementations.
     *
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! CDG 2 (Compact Discontinuous Galerkin 2) flux.
      cdg2,
      //! CDG (Compact Discontinuous Galerkin) flux.
      cdg,
      //! BR2 (Bassi-Rebay 2) flux.
      br2,
      //! IP (Interior Penalty) flux.
      ip,
      //! NIPG (Non-symmetric Interior  Penalty) flux.
      nipg,
      //! BO (Baumann-Oden) flux.
      bo,
      //! general flux: Parameter selection is done via parameter file!
      general,
      //! no diffusion (advection only) flux.
      none
    };

    //! Contains all known enums for primal diffusion fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::cdg2, Enum::cdg, Enum::br2, Enum::ip, Enum::nipg, Enum::bo };
    //! Contains all known names of primal diffusion fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "CDG2", "CDG" , "BR2", "IP" , "NIPG", "BO" };
    //! Number of known primal diffusion fluxes which can be chosen via parameter file.
    static const int  _size = 6;

    //! Helper class for static parameter selection
    template< Enum id = Enum::general >
    struct Identifier
    {
      typedef Enum type;
      static const type value = id;
    };
  }

  /**
   * \brief Namespace containing all parameters to select a lifting for primal diffusion fluxes.
   */
  namespace PrimalDiffusionLifting
  {
    /**
     * \brief Enum of all known liftings for primal diffusion fluxes.
     *
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      //! \f$ \int_\Omega r([u])\cdot\tau  = -\int_e [u]\cdot\{\tau\} \f$
      id_id,
      //! \f$ \int_\Omega r([u])\cdot\tau  = -\int_e [u]\cdot\{A\tau\} \f$
      id_A,
      //! \f$ \int_\Omega r([u])\cdot A\tau = -\int_e [u]\cdot\{A\tau\} \f$
      A_A
    };

    //! Contains all known enums for liftings of primal diffusion fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::id_id, Enum::id_A, Enum::A_A };
    //! Contains all known names lifting of primal diffusion fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "id_id", "id_A" , "A_A" };
    //! Number of known liftings for primal diffusion fluxes which can be chosen via parameter file.
    static const int  _size = 3;

    //! Helper class for static parameter selection
    struct Identifier
    {
      typedef Enum type;
    };
  }



  /**
   * \brief Parameter class for primal diffusion flux parameters.
   *
   * \ingroup ParameterClass
   */
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

    /**
     * \brief Constructor
     *
     * \param[in] keyPrefix key prefix for parameter file.
     */
    DGPrimalDiffusionFluxParameters( const std::string keyPrefix = "dgdiffusionflux." )
      : keyPrefix_( keyPrefix )
    {}

    /**
     * \brief returns name of the flux
     *
     * \param[in] mthd enum of primal diffusion flux
     * \returns string which could be used for the selection of a flux in a parameter file.
     */
    static std::string methodNames( const IdEnum& mthd )
    {
      for( int i = 0; i < PrimalDiffusionFlux::_size; i++)
        if( PrimalDiffusionFlux::_enums[i] == mthd )
          return PrimalDiffusionFlux::_strings[i];
      assert( false );
      return "invalid diffusion flux";
    }

    /**
     * \brief returns enum of the flux
     */
    virtual IdEnum getMethod() const
    {
      const int i = Fem::Parameter::getEnum( keyPrefix_ + "method", PrimalDiffusionFlux::_strings );
      return PrimalDiffusionFlux::_enums[i];
    }

    /**
     * \brief returns name of the flux
     *
     * \param[in] mthd enum of a lifting of the primal diffusion flux
     * \returns string which could be used for the selection of a lifting in a parameter file.
     */
    static std::string liftingNames( const LiftingEnum mthd )
    {
      for( int i = 0; i < PrimalDiffusionLifting::_size; i++)
        if( PrimalDiffusionLifting::_enums[i] == mthd )
          return PrimalDiffusionLifting::_strings[i];
      assert( false );
      return "invalid identifier";
    }

    /**
     * \brief returns enum of the lifting
     */
    virtual LiftingEnum getLifting() const
    {
      const int i = Fem::Parameter::getEnum( keyPrefix_ + "lifting", PrimalDiffusionLifting::_strings, 0 );
      return PrimalDiffusionLifting::_enums[i];
    }

    //! todo please doc me
    virtual double penalty() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "penalty" );
    }

    //! todo please doc me
    virtual double liftfactor() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "liftfactor" );
    }

    /**
     * \brief Returns whether to use parameters given in the literature.
     */
    virtual double theoryparameters() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "theoryparameters", 0. );
    }

    //! todo please doc me
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
