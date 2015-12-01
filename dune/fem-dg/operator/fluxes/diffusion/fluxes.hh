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
            class DiffusionFluxIdentifierImp >
  class DGPrimalDiffusionFlux;

  //////////////////////////////////////////////////////////
  //
  //  general diffusion flux allows choice of method via Parameter
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::general > >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::general >  >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::general >  >
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

  //////////////////////////////////////////////////////////
  //
  //  specialization for CDG2
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::cdg2 > >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::cdg2 > >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::cdg2 > >
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


  //////////////////////////////////////////////////////////
  //
  //  specialization for CDG
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::cdg > >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::cdg > >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::cdg > >
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


  //////////////////////////////////////////////////////////
  //
  //  specialization for BR2
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::br2 > >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::br2 > >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::br2 > >
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


  //////////////////////////////////////////////////////////
  //
  //  specialization for IP
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::ip > >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::ip > >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::ip > >
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

  //////////////////////////////////////////////////////////
  //
  //  specialization for no-diffusion
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model>
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, PrimalDiffusionFlux::Identifier< PrimalDiffusionFlux::Enum::none > >
    : public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::none > >
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model, DGPrimalDiffusionFluxParameters< PrimalDiffusionFlux::Enum::none > >
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
