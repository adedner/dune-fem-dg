#ifndef DUNE_DGOPERATORBASE_HH
#define DUNE_DGOPERATORBASE_HH

#include <string>

#include "../pass/ippass.hh"
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#ifdef NSMOD_USE_SMP_PARALLEL
#include <dune/fem-dg/pass/threadpass.hh>
#endif

namespace Dune {  


  // DGAdvectionDiffusionOperatorBase
  //---------------------------------

  template< class Traits > 
  class DGAdvectionDiffusionOperatorBase : 
    public SpaceOperatorInterface 
      < typename PassTraits< 
          typename Traits :: Model, Traits::dimRange, Traits :: polOrd > :: DestinationType >
  {
    enum { u = Traits :: u, 
           cdgpass  = Traits :: cdgpass };

    enum { polOrd = Traits :: polOrd };
    
  public:
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
    
    typedef StartPass< AdvDFunctionType, u > Pass0Type;
    typedef 
#ifdef NSMOD_USE_SMP_PARALLEL
      ThreadPass<
#endif
        LocalCDGPass< DiscreteModelType, Pass0Type, cdgpass >
#ifdef NSMOD_USE_SMP_PARALLEL 
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
    DGAdvectionDiffusionOperatorBase( GridType& grid , const NumFluxType& numf ) 
      : grid_( grid )
      , model_( numf.model() )
      , numflux_( numf )
      , gridPart_( grid_ )
      , space_( gridPart_ )
      , discreteModel_( model_, numflux_, DiffusionFluxType( gridPart_, model_ ) )
      , startPass_()
      , pass1_( discreteModel_, startPass_, space_ )
    {}

    IndicatorType* indicator() { return 0; }

    void setAdaptationHandler( AdaptationHandlerType& adHandle, double weight = 1 ) 
    {
      discreteModel_.setAdaptationHandler( adHandle, weight );
    }

    void setTime(const double time) {
	    pass1_.setTime( time );
    }

    double timeStepEstimate() const {
	    return pass1_.timeStepEstimate();
    }

    void operator()( const DestinationType& arg, DestinationType& dest ) const {
	    pass1_( arg, dest );
    }

    inline const DiscreteFunctionSpaceType& space() const {
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

    inline void limit(const DestinationType& arg,DestinationType& dest) const
    {}

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
    GridType& grid_;
    const Model& model_;
    const NumFluxType& numflux_;
    GridPartType gridPart_;
    AdvDFunctionSpaceType space_;
    DiscreteModelType discreteModel_;
    Pass0Type startPass_;
    Pass1Type pass1_;
  };

}
#endif
