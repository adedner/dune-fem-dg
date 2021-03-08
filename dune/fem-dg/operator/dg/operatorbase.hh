#ifndef DUNE_FEM_DG_OPERATORBASE_HH
#define DUNE_FEM_DG_OPERATORBASE_HH

#include <string>
#include <memory>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

// dune-fem-dg includes
#include <dune/fem-dg/pass/insertfunction.hh>
#include <dune/fem-dg/pass/dgpass.hh>
#include <dune/fem-dg/misc/parameterkey.hh>

#include <dune/fem/misc/threads/domainthreaditerator.hh>
#include <dune/fem/misc/threads/threaditerator.hh>
#include <dune/fem-dg/pass/threadpass.hh>

namespace Dune
{
namespace Fem
{

  // DGAdvectionDiffusionOperatorBase
  //---------------------------------

  /**
   * \brief advection diffusion operator
   *
   * \note This operator is based on the Pass-Concept
   *
   * \ingroup PassBased
   * \ingroup PassOperator
   */
  template< class Traits >
  class DGAdvectionDiffusionOperatorBase :
    public Fem::SpaceOperatorInterface< typename Traits::DestinationType >
  {
    enum { u = Traits::u,
           cdgpass  = Traits::cdgpass };

    enum { polynomialOrder = Traits::polynomialOrder };

    typedef Fem::SpaceOperatorInterface< typename Traits::DestinationType >  BaseType;

    static const bool threading = Traits :: threading ;

    //note: ExtraParameterTuple contains non pointer types from now on
    template <class Tuple, class StartPassImp, unsigned long int i >
    struct InsertFunctions
    {
      static_assert( i <= std::tuple_size< Tuple >::value,
                     "Not enough extra discrete functions specified (see: ExtraParameterTuple) requested by model." );

      typedef InsertFunctions< Tuple, StartPassImp, i-1 >                        PreviousInsertFunctions;
      typedef typename PreviousInsertFunctions::PassType                         PreviousPass;

      typedef std::tuple_element_t< i-1, Tuple >                                 DiscreteFunction;

      typedef Dune::Fem::InsertFunctionPass< DiscreteFunction, PreviousPass, i-1 > PassType;

      template< class ExtraArgImp >
      static std::shared_ptr< PassType > createPass( const std::shared_ptr< ExtraArgImp >& tuple )
      {
        static const int size = ExtraArgImp::size;
        //this has to fit!
        static_assert( i <= size,
                       "Not enough Discrete Function Type in container!");

        //static_assert( std::is_same<std::tuple_element_t<i-1,ExtraArgImp>,std::tuple_element_t<i-1,Tuple> >::value,
        //               "Discrete Function Type has to have the same type!");

        auto previousPass = PreviousInsertFunctions::createPass( tuple );
        //maybe empty!
        //const std::shared_ptr<DiscreteFunction> df = std::get< i-1 >( tuple );
        const auto df = (*tuple)( _index<i-1>() ).get();
        return std::make_shared<PassType>( df, previousPass );
      }
    };

    template <class Tuple, class StartPassImp >
    struct InsertFunctions< Tuple, StartPassImp, 0 >
    {
      typedef StartPassImp            PassType;

      template< class ExtraArgImp >
      static decltype(auto) createPass( ExtraArgImp& tuple )
      {
        return std::make_shared<PassType>();
      }
    };
  public:
    using BaseType::operator () ;

    // dummy method for a troubled cell indicator to be passed to the
    // limited advection operator
    typedef void* TroubledCellIndicatorType;
    void setTroubledCellIndicator(TroubledCellIndicatorType indicator) {}

    typedef typename Traits::ModelType                    ModelType;
    typedef typename ModelType::ProblemType               ProblemType ;

    typedef typename Traits::GridType                     GridType;
    typedef typename Traits::DiscreteModelType            DiscreteModelType;

    typedef typename DiscreteModelType::AdvectionFluxType AdvectionFluxType;
    typedef typename DiscreteModelType::DiffusionFluxType DiffusionFluxType;

    typedef typename DiscreteModelType::Traits            AdvTraits;

    typedef typename AdvTraits::DestinationType           AdvDFunctionType;
    // for convenience (not used here)
    typedef AdvDFunctionType                              IndicatorType;
    typedef typename AdvTraits::GridPartType              GridPartType;

    // select non-blocking communication handle
    typedef typename
      std::conditional< threading,
          NonBlockingCommHandle< AdvDFunctionType >,
          EmptyNonBlockingComm > :: type NonBlockingCommHandleType;

    typedef Fem::StartPass< AdvDFunctionType, u, NonBlockingCommHandleType >  Pass0Type;

    typedef typename Traits::ExtraParameterTupleType      ExtraParameterTupleType;

    typedef InsertFunctions< ExtraParameterTupleType, Pass0Type, ModelType::modelParameterSize >
                                                          InsertFunctionsType;

    typedef typename InsertFunctionsType::PassType        InsertFunctionPassType;

    typedef Fem::ThreadIterator< GridPartType >           ThreadIteratorType;

    typedef LocalCDGPass< DiscreteModelType, InsertFunctionPassType, cdgpass >   InnerPassType;
    typedef typename std::conditional< threading,
         ThreadPass< InnerPassType, ThreadIteratorType, true /* nonblockingcomm */ >,
         InnerPassType > :: type                                                 Pass1Type;

    typedef typename AdvTraits::DiscreteFunctionSpaceType AdvDFunctionSpaceType;
    typedef typename AdvTraits::DestinationType AdvDestinationType;

    typedef AdvDFunctionSpaceType DiscreteFunctionSpaceType;
    typedef AdvDestinationType DestinationType;

    typedef typename DiscreteModelType::AdaptationType  AdaptationType;

  public:
    template< class ExtraParameterTupleImp >
    DGAdvectionDiffusionOperatorBase( GridPartType& gridPart,
                                      const ModelType& model,
                                      const AdvectionFluxType& numFlux,
                                      ExtraParameterTupleImp& tuple,
                                      const std::string name = "",
                                      const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : gridPart_( gridPart )
      , model_( model )
      , numFlux_( numFlux )
      , space_( gridPart_ )
      , discreteModel_( model_, numFlux_, DiffusionFluxType( gridPart_, model_, parameter ) )
      , previousPass_( InsertFunctionsType::createPass( tuple ) )
      , pass1_( discreteModel_, *previousPass_, space_ )
      , counter_(0)
    {}

    IndicatorType* indicator() { return 0; }

    void setAdaptation( AdaptationType& adHandle, double weight = 1 )
    {
      // also set adaptation handler to the discrete models in the thread pass
      pass1_.setAdaptation( adHandle, weight );
    }

    void setTime(const double time) {
	    pass1_.setTime( time );
    }

    double timeStepEstimate() const {
	    return pass1_.timeStepEstimate();
    }

    //! evaluate the spatial operator
    void operator()( const DestinationType& arg, DestinationType& dest ) const {
      called();
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
    int counter() const {return counter_;}
    void called() const {counter_++;}

    inline void switchupwind()
    {
      // call upwind switcher on pass (in case its a thread pass)
      pass1_.switchUpwind();
    }

    template <class Entity, class Intersection, class Quadrature>
    inline void numericalFlux(const DestinationType &u,
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

    const DiscreteModelType& discreteModel() const { return discreteModel_; }

  protected:
    GridPartType&                              gridPart_;
    const ModelType&                           model_;
    const AdvectionFluxType&                   numFlux_;

    AdvDFunctionSpaceType                      space_;
    DiscreteModelType                          discreteModel_;
    std::shared_ptr< InsertFunctionPassType >  previousPass_;
    Pass1Type                                  pass1_;
    mutable int                                counter_;
  };

}
}
#endif
