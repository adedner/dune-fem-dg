#ifndef DUNE_FEM_DG_PASSTRAITS_HH
#define DUNE_FEM_DG_PASSTRAITS_HH

// Dune includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/fvspace.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>


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

    //typedef Fem::ElementQuadrature< GridPartType, 0 >                     VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 0 >                     VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 >                     FaceQuadratureType;
    //typedef ElementQuadrature< GridPartType, 1 >                     FaceQuadratureType;

    // Allow generalization to systems
    typedef Fem::FunctionSpace< ctype, double, dimDomain, dimRange >      FunctionSpaceType;
    typedef Fem::DiscontinuousGalerkinSpace< 
                                        FunctionSpaceType,
                                        GridPartType, polOrd,
                                        Fem::CachingStorage >             DiscreteFunctionSpaceType;
    //typedef Fem::LegendreDiscontinuousGalerkinSpace< FunctionSpaceType,
    //                                    GridPartType, polOrd,
    //                                    Fem::CachingStorage >             DiscreteFunctionSpaceType;
    typedef Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >    DestinationType;

    // Indicator for Limiter
    typedef Fem::FunctionSpace< ctype, double, dimDomain, 3> FVFunctionSpaceType;
    typedef Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Fem::SimpleStorage> IndicatorSpaceType;
    typedef Fem::AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;

    typedef AdaptationHandler< GridType, FunctionSpaceType >  AdaptationHandlerType ;
  };

} // end namespace Dune 
#endif
