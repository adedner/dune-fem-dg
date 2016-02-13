#ifndef DUNE_FEM_DG_DGPRIMALFLUXES__HH
#define DUNE_FEM_DG_DGPRIMALFLUXES__HH

#include "fluxbase.hh"
#include "dgprimalfluxes.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief The purpose of this class is to allow the selection of a primal diffusion flux
   * via an enum given in DiffusionFlux::Enum.
   *
   * \warning NIPG and BO are not implemented since these methods are not so
   * interesting, use PrimalDiffusionFlux::Enum::general for this
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
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Enum::general >
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
      : BaseType( gridPart, model, parameters )
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
      : BaseType( gridPart, model, parameters )
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
      : BaseType( gridPart, model, parameters )
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
      : BaseType( gridPart, model, parameters )
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
      : BaseType( gridPart, model, parameters )
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


} // end namespace
} // end namespace
#endif
