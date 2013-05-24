#ifndef DUNE_FEM_DG_OPERATORBASE_HH
#define DUNE_FEM_DG_OPERATORBASE_HH

#include <string>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

// dune-fem-dg includes 
#include <dune/fem-dg/pass/dgpass.hh>
#include <dune/fem-dg/operator/dg/passtraits.hh>

#ifdef USE_SMP_PARALLEL
#include <dune/fem-dg/pass/threadpass.hh>
#endif

namespace Dune {  
    
  // DGAdvectionDiffusionOperatorBase
  //---------------------------------

  template< class Traits > 
  class DGAdvectionDiffusionOperatorBase : 
    public Fem::SpaceOperatorInterface 
      < typename PassTraits< 
          typename Traits :: Model, Traits::dimRange, Traits :: polOrd > :: DestinationType >
  {
    enum { u = Traits :: u, 
           cdgpass  = Traits :: cdgpass };

    enum { polOrd = Traits :: polOrd };

    typedef Fem::SpaceOperatorInterface < typename PassTraits<
                      typename Traits :: Model, Traits::dimRange, Traits :: polOrd > ::
                      DestinationType > BaseType ;
    
  public:
    using BaseType :: operator () ;

    typedef typename Traits :: NumFluxType NumFluxType;
    typedef typename Traits :: Model Model;

    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef typename Model::Traits::GridType GridType;

    typedef typename Traits :: DiscreteModelType DiscreteModelType;
    typedef typename DiscreteModelType :: DiffusionFluxType DiffusionFluxType;

    typedef typename DiscreteModelType::Traits AdvTraits;
    
    typedef typename AdvTraits::DiscreteFunctionType AdvDFunctionType;
    typedef typename AdvTraits::IndicatorType        IndicatorType ;
    
    typedef Fem::StartPass< AdvDFunctionType, u
#ifdef USE_SMP_PARALLEL
         , NonBlockingCommHandle< AdvDFunctionType >   
#endif
      > Pass0Type;

    typedef 
#ifdef USE_SMP_PARALLEL
      ThreadPass < 
#endif
      LocalCDGPass< DiscreteModelType, Pass0Type, cdgpass >
#ifdef USE_SMP_PARALLEL
      , true // non-blocking communication 
        > 
#endif
    Pass1Type;

    typedef typename AdvTraits::DomainType AdvDomainType;
    typedef typename AdvTraits::DiscreteFunctionSpaceType AdvDFunctionSpaceType;
    typedef typename AdvTraits::DestinationType AdvDestinationType;

    typedef AdvDomainType DomainType;
    typedef typename AdvTraits::GridPartType GridPartType;
    typedef AdvDFunctionSpaceType DiscreteFunctionSpaceType;
    typedef AdvDestinationType DestinationType;

    typedef typename DiscreteModelType :: AdaptationHandlerType  AdaptationHandlerType;

  public:
    DGAdvectionDiffusionOperatorBase( GridPartType& gridPart , const NumFluxType& numf ) 
      : model_( numf.model() )
      , numflux_( numf )
      , gridPart_( gridPart )
      , space_( gridPart_ )
      , discreteModel_( model_, numflux_, DiffusionFluxType( gridPart_, model_ ) )
      , startPass_()
      , pass1_( discreteModel_, startPass_, space_ )
    {}

    IndicatorType* indicator() { return 0; }

    void setAdaptationHandler( AdaptationHandlerType& adHandle, double weight = 1 ) 
    {
#ifdef USE_SMP_PARALLEL
      // also set adaptation handler to the discrete models in the thread pass 
      {
        pass1_.setAdaptationHandler( adHandle, weight );
      }
#else
      {
        // set adaptation handle to discrete model 
        discreteModel_.setAdaptationHandler( adHandle, weight );
      }
#endif
    }

    void setTime(const double time){
	    pass1_.setTime( time );
    }

    double timeStepEstimate() const {
      return pass1_.timeStepEstimate();
    }

    //! evaluate the spatial operator 
    void operator()( const DestinationType& arg, DestinationType& dest ) const 
    {
      pass1_( arg, dest );
    }

    //! only evaluate fluxes of operator 
    void evaluateOnly( const DestinationType& arg ) const {
      // only apply operator without storing result, for evalaution 
      // of the aposteriori error estimator mainly 
      DestinationType* emptyPtr = 0 ;
	    pass1_( arg, *emptyPtr );
    }

    inline const DiscreteFunctionSpaceType& space() const {
	    return space_;
    } 
    inline DiscreteFunctionSpaceType& space() {
	    return space_;
    } 

    inline void switchupwind() 
    { 
      // call upwind switcher on pass (in case its a thread pass)
      pass1_.switchUpwind();
    }

    inline double maxAdvectionTimeStep() const 
    {
      return discreteModel_.maxAdvectionTimeStep();
    } 
    inline double maxDiffusionTimeStep() const 
    {
      return discreteModel_.maxDiffusionTimeStep();
    } 

    template <class Matrix> 
    inline void operator2Matrix( Matrix& matrix, DestinationType& rhs ) const 
    {
      pass1_.operator2Matrix( matrix , rhs );
    }

    template <class Entity, class Intersection, class Quadrature>
    inline void flux(const DestinationType &u, 
                     const Entity &entity, const Entity &nb,
                     const Intersection &intersection, 
                     const Quadrature &faceQuadInner, const Quadrature &faceQuadOuter,
                     const int l,
                     typename DestinationType::RangeType &fluxEn,
                     typename DestinationType::RangeType &fluxNb) const
    {
      pass1_.flux(u,entity,nb,intersection,faceQuadInner,faceQuadOuter,l,fluxEn,fluxNb);
    }

    inline void limit( DestinationType& U ) const {}

    inline double computeTime() const 
    {
      return pass1_.computeTime();
    }

    inline size_t numberOfElements () const 
    {
      return pass1_.numberOfElements();
    }
    
    void printmyInfo(std::string filename) const {}

    virtual std::string description() const = 0;

  protected:
    const Model& model_;
    const NumFluxType& numflux_;
    GridPartType& gridPart_;
    AdvDFunctionSpaceType space_;
    DiscreteModelType discreteModel_;
    Pass0Type startPass_;
    Pass1Type pass1_;
  };

}
#endif
