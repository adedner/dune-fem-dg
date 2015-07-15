#ifndef DUNE_FEM_DG_PASSTRAITS_HH
#define DUNE_FEM_DG_PASSTRAITS_HH

// Dune includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

namespace Dune {

  //PassTraits
  //----------

  template <class Traits, int polOrd, int dimR>
  class PassTraits
  {
  public:
    // inherit types from Traits
    typedef typename Traits::InitialDataType InitialDataType;
    typedef typename Traits::ModelType  ModelType ;
    typedef typename Traits::FluxType   FluxType;
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId  = Traits :: PrimalDiffusionFluxId ;

    typedef typename ModelType :: Traits ModelTraits;

    typedef typename ModelTraits  :: GridPartType        GridPartType;
    typedef typename ModelTraits  :: GridType            GridType;
    typedef typename GridType     :: ctype               ctype;

    static const int polynomialOrder = polOrd ;
    static const int dimRange  = dimR ;
    static const int dimDomain = Traits::ModelType::dimDomain ;

    typedef typename ModelTraits::FaceDomainType  FaceDomainType;

    typedef Fem::FunctionSpace< ctype, double, dimDomain, dimRange >      FunctionSpaceType;
    typedef Fem::FunctionSpace< ctype, double, dimDomain, 1 >             PressureFunctionSpaceType;

    // CACHING
    typedef Dune::Fem::CachingQuadrature< GridPartType, 1 >  FaceQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >  VolumeQuadratureType;

    /*
#if DGSCHEME
#if ONB
    #warning using DG space with ONB
    typedef Dune::Fem::DiscontinuousGalerkinSpace
#elif LAG
    #warning using DG space with Lagrange base functions
    typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace
#elif LEG
    #warning using DG space with hierarich Legendre base functions
    typedef Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace
#else
#define USE_STRONG_BC
    #warning using p-adaptive DG space
    typedef Dune::Fem::PAdaptiveDGSpace
#define PADAPTSPACE
#endif
#else
#define USE_STRONG_BC
#if LAGRANGESPACE
    #warning using Lagrange space
    typedef Dune::Fem::LagrangeDiscreteFunctionSpace
#else
#define USE_STRONG_BC
    #warning using p-adaptive Lagrange space
    typedef Dune::Fem::PAdaptiveLagrangeSpace
#define PADAPTSPACE
#endif
#endif
    */
    typedef Dune::Fem::DiscontinuousGalerkinSpace
       < FunctionSpaceType, GridPartType, polOrd, Dune::Fem::CachingStorage > DiscreteFunctionSpaceType;
    typedef Dune::Fem::DiscontinuousGalerkinSpace < PressureFunctionSpaceType, GridPartType,
            polOrd, Dune::Fem::CachingStorage > PressureDiscreteFunctionSpaceType;

    static const bool symmetricSolver = true ;
//#if WANT_ISTL
    typedef Solvers<DiscreteFunctionSpaceType, istl,  symmetricSolver> SolversType;
/*
#elif WANT_PETSC
    typedef Solvers<DiscreteFunctionSpaceType, petsc, symmetricSolver> SolversType;
#else
    typedef Solvers<DiscreteFunctionSpaceType, fem,   symmetricSolver> SolversType;
#endif
*/

    typedef typename SolversType :: DiscreteFunctionType       DestinationType;
    typedef typename SolversType :: LinearOperatorType         LinearOperatorType;
    typedef typename SolversType :: LinearInverseOperatorType  LinearInverseOperatorType;
  };

} // end namespace Dune
#endif
