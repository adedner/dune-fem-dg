#ifndef DUNE_FEM_DG_DIFFUSIONFLUXPARAMETER_HH
#define DUNE_FEM_DG_DIFFUSIONFLUXPARAMETER_HH

#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

  /**
   *  \brief Namespace containing an Enum class to describe the formulation
   *  \note Selection of operators reflecting the formulation see FormulationSelector
   */
  namespace Formulation
  {
    /**
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      primal,
      local
    };
  }

  /**
   * \brief Namespace containing all parameters to select a DG method available in primal formulation
   */
  namespace DiffusionFlux
  {
    /**
     * \brief Enum of all available diffusion flux implementations (primal and local formulation).
     *
     * \ingroup FemDGParameter
     */
    enum class Enum
    {
      default_,
      //! BR1 (Bassi-Rebay 1) flux (local formulation).
      br1,
      //! LDG (Local Discontinuous Galerkin) flux (local formulation).
      ldg,
      //! CDG 2 (Compact Discontinuous Galerkin 2) flux (primal formulation).
      cdg2,
      //! CDG (Compact Discontinuous Galerkin) flux (primal formulation).
      cdg,
      //! BR2 (Bassi-Rebay 2) flux (primal formulation).
      br2,
      //! IP (Interior Penalty) flux (primal formulation).
      ip,
      //! NIPG (Non-symmetric Interior  Penalty) flux (primal formulation).
      nipg,
      //! BO (Baumann-Oden) flux (primal formulation).
      bo,
      //! general flux primal formulation: Parameter selection is done via parameter file!
      primal,
      //! general flux local formulation: Parameter selection is done via parameter file!
      local,
      //! no diffusion (advection only) flux.
      none
    };
  }

  /**
   * \brief Namespace containing all parameters to select a DG method available in primal formulation
   */
  namespace PrimalDiffusionFlux
  {
    using namespace DiffusionFlux;

    //! Contains all known enums for primal diffusion fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::cdg2, Enum::cdg, Enum::br2, Enum::ip, Enum::nipg, Enum::bo };
    //! Contains all known names of primal DG diffusion fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "CDG2", "CDG" , "BR2", "IP" , "NIPG", "BO" };
    //! Number of known primal diffusion fluxes which can be chosen via parameter file.
    static const int  _size = 6;

  }

  namespace LocalDiffusionFlux
  {
    using namespace DiffusionFlux;

    //! Contains all known enums for dual diffusion fluxes which can be chosen via parameter file.
    const Enum        _enums[] = { Enum::br1, Enum::ldg };
    //! Contains all known names of dual diffusion fluxes which can be chosen via parameter file.
    const std::string _strings[] = { "BR1", "LDG" };
    //! Number of known primal diffusion fluxes which can be chosen via parameter file.
    static const int  _size = 2;

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

  }



  /**
   * \brief Parameter class for primal diffusion flux parameters.
   *
   * \ingroup ParameterClass
   */
  class DGPrimalDiffusionFluxParameters
    : public Fem::LocalParameter< DGPrimalDiffusionFluxParameters, DGPrimalDiffusionFluxParameters >
  {
  public:
    typedef PrimalDiffusionFlux::Enum           IdEnum;
    typedef PrimalDiffusionLifting::Enum        LiftingEnum;

    /**
     * \brief Constructor
     *
     * \param[in] keyPrefix key prefix for parameter file.
     */
    DGPrimalDiffusionFluxParameters( const std::string& keyPrefix,
                                     const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : keyPrefix_( keyPrefix ),
        parameter_( parameter )
    {}

    /**
     * \brief Default Constructor
     *
     */
    DGPrimalDiffusionFluxParameters( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : DGPrimalDiffusionFluxParameters( "dgdiffusionflux.", parameter )
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
      const int i = parameter_.getEnum( keyPrefix_ + "method", PrimalDiffusionFlux::_strings, 0 );
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
      const int i = parameter_.getEnum( keyPrefix_ + "lifting", PrimalDiffusionLifting::_strings, 0 );
      return PrimalDiffusionLifting::_enums[i];
    }

    //! todo please doc me
    virtual double penalty() const
    {
      return parameter_.getValue<double>( keyPrefix_ + "penalty", 1.0 );
    }

    //! todo please doc me
    virtual double liftfactor() const
    {
      return parameter_.getValue<double>( keyPrefix_ + "liftfactor", 1.0 );
    }

    /**
     * \brief Returns whether to use parameters given in the literature.
     */
    virtual double theoryparameters() const
    {
      return parameter_.getValue<double>( keyPrefix_ + "theoryparameters", 0. );
    }

    //! todo please doc me
    template <class DomainType>
    void upwind( DomainType& upwd ) const
    {
      parameter_.get(keyPrefix_ + "upwind", upwd, upwd);
    }

    const Dune::Fem::ParameterReader& parameters() const
    {
      return parameter_;
    }
  private:

    const std::string keyPrefix_;
    const Dune::Fem::ParameterReader parameter_;
  };

  /**
   * \brief Parameter class for primal diffusion flux parameters.
   *
   * \ingroup ParameterClass
   */
  class DGLocalDiffusionFluxParameters
    : public Fem::LocalParameter< DGLocalDiffusionFluxParameters, DGLocalDiffusionFluxParameters >
  {
  public:
    typedef LocalDiffusionFlux::Enum           IdEnum;

    /**
     * \brief Constructor
     *
     * \param[in] keyPrefix key prefix for parameter file.
     */
    DGLocalDiffusionFluxParameters( const std::string keyPrefix = "dgdiffusionflux." )
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
      for( int i = 0; i < LocalDiffusionFlux::_size; i++)
        if( LocalDiffusionFlux::_enums[i] == mthd )
          return LocalDiffusionFlux::_strings[i];
      assert( false );
      return "invalid diffusion flux";
    }

    /**
     * \brief returns enum of the flux
     */
    virtual IdEnum getMethod() const
    {
      const int i = Fem::Parameter::getEnum( keyPrefix_ + "method", LocalDiffusionFlux::_strings );
      return LocalDiffusionFlux::_enums[i];
    }

    //! todo please doc me
    virtual double penalty() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "penalty" );
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

  //! formulation selector, either primal or local depending on the flux implementation
  template <DiffusionFlux::Enum flux>
  struct FormulationSelector
  {
    static const Formulation::Enum formId = Formulation::Enum::primal;
  };

  template <>
  struct FormulationSelector< DiffusionFlux::Enum::local >
  {
    static const Formulation::Enum formId = Formulation::Enum::local;
  };

  template <>
  struct FormulationSelector< DiffusionFlux::Enum::ldg >
  {
    static const Formulation::Enum formId = Formulation::Enum::local;
  };

  template <>
  struct FormulationSelector< DiffusionFlux::Enum::br1 >
  {
    static const Formulation::Enum formId = Formulation::Enum::local;
  };


} // end namespace
} // end namespace
#endif
