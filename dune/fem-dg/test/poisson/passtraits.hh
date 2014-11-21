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

namespace Dune {

  //PassTraits
  //----------

  template <class Model,int dimRange,int polOrd>
  class PassTraits
  {
  public:
    typedef typename Model :: Traits                                 ModelTraits;
    typedef typename ModelTraits :: GridPartType                     GridPartType;
    typedef typename GridPartType :: GridType                        GridType;
    typedef typename GridType :: ctype                               ctype;
    static const int dimDomain = Model :: Traits :: dimDomain;

    template <class Grid, int dim >
    struct FaceQuadChooser
    {
      typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > Type;
    };

#if HAVE_UG
    template <int dim>
    struct FaceQuadChooser< const Dune::UGGrid< dim >, dim >
    {
      typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > Type;
    };
#endif

    // CACHING
    typedef typename FaceQuadChooser< const GridType, GridType::dimension > :: Type  FaceQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >     VolumeQuadratureType;

    // Allow generalization to systems
    typedef Dune::Fem::FunctionSpace< ctype, double, dimDomain, dimRange >      FunctionSpaceType;
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
       < FunctionSpaceType, GridPartType, polOrd, Dune::Fem::CachingStorage > DiscreteFunctionSpaceType;

    static const bool symmetricSolver = true ;
#if WANT_ISTL
    typedef Solvers<DiscreteFunctionSpaceType, istl,  symmetricSolver> SolversType;
#elif WANT_PETSC
    typedef Solvers<DiscreteFunctionSpaceType, petsc, symmetricSolver> SolversType;
#else
    typedef Solvers<DiscreteFunctionSpaceType, fem,   symmetricSolver> SolversType;
#endif

    typedef typename SolversType :: DiscreteFunctionType       DestinationType;
    typedef typename SolversType :: LinearOperatorType         LinearOperatorType;
    typedef typename SolversType :: LinearInverseOperatorType  LinearInverseOperatorType;

  };

} // end namespace Dune
#endif
