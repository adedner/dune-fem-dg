#ifndef DUNE_FEM_DG_OPERATORBASE_HH
#define DUNE_FEM_DG_OPERATORBASE_HH

#include <string>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/pass/insertfunction.hh>

// dune-fem-dg includes
#include <dune/fem-dg/pass/dgpass.hh>
#include <dune/fem-dg/operator/dg/passtraits.hh>
#include <dune/fem-dg/misc/parameterkey.hh>

#ifdef USE_SMP_PARALLEL
#include <dune/fem/misc/threads/domainthreaditerator.hh>
#include <dune/fem/misc/threads/threaditerator.hh>
#include <dune/fem-dg/pass/threadpass.hh>
#endif

namespace Dune {


  // DGAdvectionDiffusionOperatorBase
  //---------------------------------

  template< class Traits >
  class DGAdvectionDiffusionOperatorBase :
    public Fem::SpaceOperatorInterface< typename Traits :: DestinationType >
  {
    enum { u = Traits :: u,
           cdgpass  = Traits :: cdgpass };

    enum { polynomialOrder = Traits :: polynomialOrder };

    typedef Fem::SpaceOperatorInterface< typename Traits :: DestinationType >  BaseType ;

  public:
    using BaseType :: operator () ;

    typedef typename Traits :: FluxType           AdvectionFluxType;
    typedef typename Traits :: ModelType          ModelType;
    typedef typename ModelType :: ProblemType     ProblemType ;

    enum { dimRange  = Traits::dimRange };
    enum { dimDomain = Traits::dimDomain };

    typedef typename Traits :: GridType GridType;
    typedef typename Traits :: DiscreteModelType DiscreteModelType;

    typedef typename DiscreteModelType :: DiffusionFluxType DiffusionFluxType;

    typedef typename DiscreteModelType::Traits AdvTraits;

    typedef typename AdvTraits::DestinationType        AdvDFunctionType;
    // for convenience (not used here)
    typedef AdvDFunctionType                           IndicatorType;
    typedef typename AdvTraits::GridPartType GridPartType;

    typedef Fem::StartPass< AdvDFunctionType, u
#ifdef USE_SMP_PARALLEL
         , NonBlockingCommHandle< AdvDFunctionType >
#endif
      > Pass0Type;

    typedef typename ModelType :: ModelParameter ModelParameter;
    typedef typename Traits :: ExtraParameterTupleType ExtraParameterTupleType;

    template <class Tuple, int i>
    struct InsertFunctions
    {
      typedef InsertFunctions< Tuple, i-1 > PreviousInsertFunctions;
      typedef typename PreviousInsertFunctions :: PassType PreviousPass ;
      typedef typename std::remove_pointer< typename std::tuple_element< i-1, Tuple > ::type > :: type  DiscreteFunction;
      static const int passId = std::tuple_element< i-1, ModelParameter >::type::value;
      typedef Dune::Fem::InsertFunctionPass< DiscreteFunction, PreviousPass, passId > PassType;

      static std::shared_ptr< PassType > createPass( Tuple& tuple )
      {
        std::shared_ptr< PreviousPass > previousPass = PreviousInsertFunctions::createPass( tuple );
        const DiscreteFunction* df = std::get< i-1 >( tuple );
        return std::shared_ptr< PassType > ( new PassType( df, previousPass ) );
      }
    };

    template <class Tuple>
    struct InsertFunctions< Tuple, 0 >
    {
      typedef Pass0Type PassType;
      static std::shared_ptr< PassType > createPass( Tuple& tuple )
      {
        return std::shared_ptr< PassType > ( new PassType() );
      }
    };

    typedef InsertFunctions< ExtraParameterTupleType, std::tuple_size< ExtraParameterTupleType >::value > InsertFunctionsType;
    typedef typename InsertFunctionsType :: PassType InsertFunctionPassType;

    typedef
#ifdef USE_SMP_PARALLEL
      ThreadPass <
#endif
      LocalCDGPass< DiscreteModelType, InsertFunctionPassType, cdgpass >
#ifdef USE_SMP_PARALLEL
      , Fem::DomainDecomposedIteratorStorage< GridPartType >
    //, Fem::ThreadIterator< GridPartType >
      , true // non-blocking communication
        >
#endif
    Pass1Type;

    typedef typename AdvTraits::DiscreteFunctionSpaceType AdvDFunctionSpaceType;
    typedef typename AdvTraits::DestinationType AdvDestinationType;

    typedef AdvDFunctionSpaceType DiscreteFunctionSpaceType;
    typedef AdvDestinationType DestinationType;

    typedef typename DiscreteModelType :: AdaptationType  AdaptationType;

  public:
    DGAdvectionDiffusionOperatorBase( GridPartType& gridPart, ProblemType& problem,
                                      ExtraParameterTupleType& tuple,
                                      const std::string name = "" )
      : model_( problem )
      , numflux_( model_ )
      , gridPart_( gridPart )
      , space_( gridPart_ )
      , discreteModel_( model_, numflux_,
                        DiffusionFluxType( gridPart_, model_, DGPrimalFormulationParameters( ParameterKey::generate( name, "dgdiffusionflux." ) ) ) )
      , previousPass_( InsertFunctionsType::createPass( tuple ) )
      , pass1_( discreteModel_, *previousPass_, space_ )
    {}

    IndicatorType* indicator() { return 0; }

    void setAdaptation( AdaptationType& adHandle, double weight = 1 )
    {
#ifdef USE_SMP_PARALLEL
      // also set adaptation handler to the discrete models in the thread pass
      {
        pass1_.setAdaptation( adHandle, weight );
      }
#else
      {
        // set adaptation handle to discrete model
        discreteModel_.setAdaptation( adHandle, weight );
      }
#endif
    }

    void setTime(const double time) {
	    pass1_.setTime( time );
    }

    double timeStepEstimate() const {
	    return pass1_.timeStepEstimate();
    }

    //! evaluate the spatial operator
    void operator()( const DestinationType& arg, DestinationType& dest ) const {
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
    inline void limit( const DestinationType& arg, DestinationType& U ) const {}

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

    const ModelType& model() const { return model_; }

  protected:
    ModelType           model_;
    AdvectionFluxType   numflux_;
    GridPartType& gridPart_;

    AdvDFunctionSpaceType space_;
    DiscreteModelType discreteModel_;
    std::shared_ptr< InsertFunctionPassType > previousPass_;
    Pass1Type pass1_;
  };

}
#endif
