#ifndef DUNE_FEM_DG_DGPRIMALFLUXES__HH
#define DUNE_FEM_DG_DGPRIMALFLUXES__HH

#include "fluxbase.hh"
#include "dgprimalfluxes.hh"
#include "ldgflux.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief The purpose of this class is to allow the selection of a primal diffusion flux
   * via an enum given in DiffusionFlux::Enum.
   *
   * \warning NIPG and BO are not implemented since these methods are not so
   * interesting, use PrimalDiffusionFlux::Enum::primal for this
   */
  template <class DiscreteFunctionSpaceImp,
            class Model,
            PrimalDiffusionFlux::Enum id >
  class DGPrimalDiffusionFlux;

  /**
   * \brief class specialization for a general primal diffusion flux chosen by a parameter file.
   *
   * The purpose of this class is to allow the selection of an Euler flux
   * via an enum given in DiffusionFlux::Enum.
   */
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Enum::primal >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseType :: ParameterType ParameterType ;

    /**
      * \brief constructor reading parameters
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model,
                           const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, PrimalDiffusionFlux::Enum::primal )
    {
    }
  };

  /**
   * \brief class specialization for the CDG2 diffusion flux.
   *
   * The purpose of this class is to allow the selection of a primal diffusion flux
   * via an enum given in DiffusionFlux::Enum.
   */
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Enum::cdg2 >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseType :: ParameterType  ParameterType;

    /**
      * \brief constructor reading parameters
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model,
                           const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, PrimalDiffusionFlux::Enum::cdg2 )
    {
    }
  };


  /**
   * \brief class specialization for the CDG diffusion flux.
   *
   * The purpose of this class is to allow the selection of a primal diffusion flux
   * via an enum given in DiffusionFlux::Enum.
   */
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Enum::cdg >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseType :: ParameterType  ParameterType;

    /**
      * \brief constructor reading parameters
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model,
                           const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, PrimalDiffusionFlux::Enum::cdg )
    {
    }
  };


  /**
   * \brief class specialization for the BR2 diffusion flux.
   *
   * The purpose of this class is to allow the selection of a primal diffusion flux
   * via an enum given in DiffusionFlux::Enum.
   */
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Enum::br2 >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseType :: ParameterType  ParameterType;

    /**
      * \brief constructor reading parameters
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model,
                           const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, PrimalDiffusionFlux::Enum::br2  )
    {
    }
  };


  /**
   * \brief class specialization for the IP diffusion flux.
   *
   * The purpose of this class is to allow the selection of a primal diffusion flux
   * via an enum given in DiffusionFlux::Enum.
   */
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Enum::ip >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseType :: ParameterType  ParameterType;

    /**
      * \brief constructor reading parameters
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model,
                           const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, PrimalDiffusionFlux::Enum::ip  )
    {
    }
  };

  /**
   * \brief class specialization for no diffusion flux.
   *
   * The purpose of this class is to allow the selection of a primal diffusion flux
   * via an enum given in DiffusionFlux::Enum.
   */
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Enum::none >
  : public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseType :: ParameterType  ParameterType;

  public:
    /**
      * \brief constructor reading parameters
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model,
                           const ParameterType& parameters = ParameterType() )
      : BaseType( model, false, parameters )
    {
    }

    void diffusionFluxName ( std::ostream& out ) const
    {
      out << "none";
    }

    void diffusionFluxLiftFactor ( std::ostream& out ) const {}
    void diffusionFluxPenalty ( std::ostream& out ) const {}
  };


 //======================= DUAL FLUXES ========================


  template <class DiscreteFunctionSpaceImp,
            class Model,
            LocalDiffusionFlux::Enum id >
  class DGLocalDiffusionFlux;

  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGLocalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, LocalDiffusionFlux::Enum::br1 >
    : public LDGDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGLocalDiffusionFluxParameters >
  {
    typedef LDGDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGLocalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp                         DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename BaseType::ParameterType                 ParameterType;

    /**
      * \brief constructor reading parameters
      */
    DGLocalDiffusionFlux( GridPartType& gridPart,
                           const Model& model,
                           const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, LocalDiffusionFlux::Enum::br1 )
    {
    }
  };

  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGLocalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, LocalDiffusionFlux::Enum::ldg >
    : public LDGDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGLocalDiffusionFluxParameters >
  {
    typedef LDGDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGLocalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp                         DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename BaseType::ParameterType                 ParameterType;

    /**
      * \brief constructor reading parameters
      */
    DGLocalDiffusionFlux( GridPartType& gridPart,
                          const Model& model,
                          const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, LocalDiffusionFlux::Enum::ldg )
    {
    }
  };

  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGLocalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, LocalDiffusionFlux::Enum::local >
    : public LDGDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGLocalDiffusionFluxParameters >
  {
    typedef LDGDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGLocalDiffusionFluxParameters >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp                         DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename BaseType::ParameterType                 ParameterType;

    /**
      * \brief constructor reading parameters
      */
    DGLocalDiffusionFlux( GridPartType& gridPart,
                         const Model& model,
                         const ParameterType& parameters = ParameterType() )
      : BaseType( gridPart, model, parameters, LocalDiffusionFlux::Enum::local )
    {
    }
  };

} // end namespace Fem
} // end namespace Dune
#endif
