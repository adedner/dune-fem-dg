#ifndef DUNE_FEM_DG_PASSTRAITS_HH
#define DUNE_FEM_DG_PASSTRAITS_HH

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>

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

    typedef typename FunctionSpaceType :: DomainType         DomainType;
    typedef typename FunctionSpaceType :: RangeType          RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType  JacobianRangeType;
    typedef typename FunctionSpaceType :: RangeFieldType     RangeFieldType ;
    typedef typename FunctionSpaceType :: DomainFieldType    DomainFieldType ;

    //static const int dimRange  = ModelTraits::dimRange;
    //static const int dimDomain = ModelTraits::dimDomain;
    //typedef Fem::ElementQuadrature< GridPartType, 0 >                     VolumeQuadratureType;
    //typedef ElementQuadrature< GridPartType, 1 >                     FaceQuadratureType;

    typedef Fem::CachingQuadrature< GridPartType, 0 >    VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 >    FaceQuadratureType;

    // Allow generalization to systems
    typedef Fem::DiscontinuousGalerkinSpace<
                                        FunctionSpaceType,
                                        GridPartType, polynomialOrder,
                                        Fem::CachingStorage >             DiscreteFunctionSpaceType;
    //typedef Fem::LegendreDiscontinuousGalerkinSpace< FunctionSpaceType,
    //                                    GridPartType, polOrd,
    //                                    Fem::CachingStorage >             DiscreteFunctionSpaceType;
    typedef Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >    DestinationType;

    // Indicator for Limiter
    typedef Fem::FunctionSpace< ctype, double, ModelTraits::dimDomain, 3> FVFunctionSpaceType;
    typedef Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Fem::SimpleStorage> IndicatorSpaceType;
    typedef Fem::AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;

    typedef AdaptationHandler< GridType, FunctionSpaceType >  AdaptationHandlerType ;
  };

} // end namespace Dune
#endif
