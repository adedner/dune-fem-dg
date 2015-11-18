#ifndef DUNE_FEM_DG_DGPRIMALFLUXES__HH
#define DUNE_FEM_DG_DGPRIMALFLUXES__HH

#include "fluxbase.hh"
#include "dgprimalfluxes.hh"

namespace Dune
{
namespace Fem
{

  //! DG primal diffusion flux
  template <class DiscreteFunctionSpaceImp,
            class Model,
            DGDiffusionFluxIdentifier::id >
  class DGPrimalDiffusionFlux;

  //////////////////////////////////////////////////////////
  //
  //  general diffusion flux allows choice of method via Parameter
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, DGDiffusionFluxIdentifier::general >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
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
      : BaseType( gridPart, model, parameters.getMethod(), parameters )
    {
    }
  };

  //////////////////////////////////////////////////////////
  //
  //  specialization for CDG2
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, DGDiffusionFluxIdentifier::cdg2 >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
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
      : BaseType( gridPart, model, DGDiffusionFluxIdentifier::cdg2, parameters )
    {
    }
  };


  //////////////////////////////////////////////////////////
  //
  //  specialization for CDG
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, DGDiffusionFluxIdentifier::cdg >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
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
      : BaseType( gridPart, model, DGDiffusionFluxIdentifier::cdg, parameters )
    {
    }
  };


  //////////////////////////////////////////////////////////
  //
  //  specialization for BR2
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, DGDiffusionFluxIdentifier::br2 >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
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
      : BaseType( gridPart, model, DGDiffusionFluxIdentifier::br2, parameters )
    {
    }
  };


  //////////////////////////////////////////////////////////
  //
  //  specialization for IP
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, DGDiffusionFluxIdentifier::ip >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
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
      : BaseType( gridPart, model, DGDiffusionFluxIdentifier::ip, parameters )
    {
    }
  };

  //////////////////////////////////////////////////////////
  //
  //  specialization for no-diffusion
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, DGDiffusionFluxIdentifier::none >
    : public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model >
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

  //////////////////////////////////////////////////////////
  //
  //  specialization for NIPG and BO are missing since these methods are not so
  //  interesting, use DGDiffusionFluxIdentifier::general for this
  //
  //////////////////////////////////////////////////////////



} // end namespace
} // end namespace
#endif
