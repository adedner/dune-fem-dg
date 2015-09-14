#ifndef FEMDG_DEFAULTADAPTHANDLER_HH
#define FEMDG_DEFAULTADAPTHANDLER_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem/space/common/adaptmanager.hh>

namespace Dune
{
namespace Fem
{

  template< class IndicatorImp, class GradientIndicatorImp >
  class AdaptIndicator
  {
  public:
    typedef uint64_t                                                                           UInt64Type;

    typedef IndicatorImp                                                                       IndicatorType;
    typedef GradientIndicatorImp                                                               GradientIndicatorType;
    typedef typename IndicatorType::DestinationType                                            DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType                           DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType                                   GridPartType;
    typedef typename GridPartType::GridType                                                    GridType;
    typedef Dune::AdaptationParameters                                                         AdaptationParametersType;

    typedef Dune::AdaptationHandler< GridType, typename DiscreteFunctionSpaceType::FunctionSpaceType >
                                                                                               AdaptationHandlerType;

    template< class Problem, class ExtraTupleParameter >
    AdaptIndicator( DiscreteFunctionType& sol, Problem& problem, const ExtraTupleParameter& tuple, const std::string keyPrefix = "" )
    : sol_( sol ),
      adaptationHandler_( nullptr ),
      keyPrefix_( keyPrefix ),
      adaptParam_( AdaptationParametersType( Dune::ParameterKey::generate( keyPrefix, "fem.adaptation." ) ) ),
      indicator_( const_cast<GridPartType&>(sol_.gridPart()), problem, tuple, keyPrefix_ ),
      gradientIndicator_( sol_.space(), problem, adaptParam_ )
    {}

    bool adaptive() const
    {
      return adaptParam_.adaptive();
    }

    size_t numberOfElements() const
    {
      return gradientIndicator_.numberOfElements();
    }

    UInt64Type globalNumberOfElements() const
    {
      if( adaptive() && adaptationHandler_ )
      {
        UInt64Type globalElements = adaptationHandler_->globalNumberOfElements() ;
        return globalElements;
      }
      return 0;
    }

    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp )
    {
      // create adaptation handler in case of apost indicator
      if( adaptive() )
      {
        if( !adaptationHandler_ && adaptParam_.aposterioriIndicator() )
        {
          adaptationHandler_.reset( new AdaptationHandlerType( const_cast<GridType&>(sol_.gridPart().grid()), tp ) );
          indicator_.setAdaptation( *adaptationHandler_ );
        }
      }
    }

    void finalize()
    {
      adaptationHandler_.reset( nullptr );
    }

    void estimateMark( const bool initialAdapt = false )
    {
      if( adaptive() )
      {
        if( useIndicator() )
        {
          //clear old values of indicator
          adaptationHandler_->resetStatus();

          indicator_.evaluateOnly( sol_ );
          // mark all entities depending on error
          adaptationHandler_->markEntities( initialAdapt );
        }
        else if( useGradientIndicator() )
        {
          gradientIndicator_.estimateAndMark( sol_ );
        }
      }
    }

    void postAdapt()
    {
      if( adaptive() )
        if( useIndicator() )
          adaptationHandler_->resetStatus();

    }

    void preAdapt()
    {}

    int minNumberOfElements() const
    {
      assert( adaptationHandler_ );
      return adaptationHandler_->minNumberOfElements();
    }

    int maxNumberOfElements() const
    {
      assert( adaptationHandler_ );
      return adaptationHandler_->maxNumberOfElements();
    }

    const int finestLevel() const
    {
      return adaptParam_.finestLevel();
    }

  private:
    AdaptationHandlerType* adaptationHandler()
    {
      assert( adaptationHandler_ );
      return adaptationHandler_;
    }

    bool useIndicator() const
    {
      return (adaptationHandler_.get()!=nullptr);
    }

    bool useGradientIndicator() const
    {
      return adaptParam_.gradientBasedIndicator();
    }

    DiscreteFunctionType&                         sol_;
    std::unique_ptr< AdaptationHandlerType >      adaptationHandler_;
    const std::string                             keyPrefix_;
    const AdaptationParametersType                adaptParam_;
    IndicatorType                                 indicator_;
    GradientIndicatorType                         gradientIndicator_;
  };


  template< class ... StepperArg >
  class AdaptHandler;


  template< class StepperHead, class... StepperArg >
  class AdaptHandler< StepperHead, StepperArg... >
  {
    typedef uint64_t                                                                           UInt64Type;

    typedef std::tuple< typename std::add_pointer< StepperHead >::type, typename std::add_pointer< StepperArg >::type... > StepperTupleType;
    typedef typename StepperHead::GridType GridType;

    typedef Dune::Fem::RestrictProlongDefaultTuple< typename StepperHead::DiscreteFunctionType, typename StepperArg::DiscreteFunctionType ... >      RestrictionProlongationType;

    typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType >              AdaptationManagerType;

    typedef Dune::AdaptationParameters                                                         AdaptationParametersType;

    struct EstimateMark {
      template<class T, class... Args > static void applyImpl( T e, Args&& ... a )
      { e->estimateMark( std::forward<Args>(a)... ); }
    };
    struct SetAdaptation {
      template<class T, class... Args > static void applyImpl( T e, Args&& ... a )
      { e->setAdaptation( std::forward<Args>(a)... ); }
    };
    struct PreAdapt {
      template<class T, class... Args > static void applyImpl( T e, Args&& ... a )
      { e->preAdapt( std::forward<Args>(a)... ); }
    };
    struct PostAdapt {
      template<class T, class... Args > static void applyImpl( T e, Args&& ... a )
      { e->postAdapt( std::forward<Args>(a)... ); }
    };
    struct Finalize {
      template<class T, class... Args > static void applyImpl( T e, Args&& ... a )
      { e->finalize( std::forward<Args>(a)... ); }
    };
    struct MinMaxNumElements {
      template<class T, class... Args > static void applyImpl( T e, int& min, int& max, Args&& ... a )
      {
        min = std::min( min, e->minNumberOfElements( std::forward<Args>(a)... ) );
        max = std::max( max, e->maxNumberOfElements( std::forward<Args>(a)... ) );
      }
    };
    struct NumberOfElements {
      template<class T, class... Args > static void applyImpl( T e, int& max, Args&& ... a )
      { max = std::max( max, e->numberOfElements( std::forward<Args>(a)... ) ); }
    };
    struct GlobalNumberOfElements {
      template<class T, class... Args > static void applyImpl( T e, int& max, Args&& ... a )
      { max = std::max( max, e->globalNumberOfElements( std::forward<Args>(a)... ) ); }
    };
    struct FinestLevel {
      template<class T, class... Args > static void applyImpl( T e, int& max, Args&& ... a )
      { max = std::max( max, e->finestLevel( std::forward<Args>(a)... ) ); }
    };
    struct Adaptive {
      template<class T, class... Args > static void applyImpl( T e, bool& adaptive, Args&& ... a )
      { adaptive |= e->adaptive( std::forward<Args>(a)... ); }
    };

    template< class Caller >
    class LoopCallee
    {
      template<class C, class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      getAdaptIndicator( T, Args&& ... ){}
      template<class C, class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      getAdaptIndicator( T elem, Args &&... a )
      {
        if( elem->adaptIndicator() )
          C::applyImpl(elem->adaptIndicator(), std::forward<Args>(a)... );
      }
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... Args >
        static void apply ( Tuple &tuple, Args&& ... a )
        {
          getAdaptIndicator< Caller >( std::get<i>( tuple ), std::forward<Args>(a)... );
        }
      };
    };



  public:

    AdaptHandler( StepperTupleType& tuple )
    : tuple_( tuple ),
      rp_( nullptr ),
      adaptationManager_(),
      keyPrefix_( "" ),
      adaptParam_( AdaptationParametersType( Dune::ParameterKey::generate( keyPrefix_, "fem.adaptation." ) ) )
    {
      setRestrProlong( Std::index_sequence_for< StepperHead, StepperArg ... >() );
      if( adaptive() )
        rp_->setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
    }

    template< std::size_t ... i >
    void setRestrProlong( Std::index_sequence< i ... > )
    {
      rp_.reset( new RestrictionProlongationType( *std::get< i >( tuple_ )->adaptationSolution()... ) );
    }

    bool adaptive () const
    {
      bool adaptive = false;
      ForLoop< LoopCallee<Adaptive>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, adaptive );
      return adaptive;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      if( tp.timeStep() % adaptParam_.adaptCount() == 0 )
      {
        estimateMark( false );
        adapt();
      }
    }

    void init()
    {
      estimateMark( true );
      adapt();
    }

    size_t numberOfElements() const
    {
      int numElements = 0;
      ForLoop< LoopCallee<NumberOfElements>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, numElements );
      return numElements;
    }

    UInt64Type globalNumberOfElements() const
    {
      if( adaptive() )
      {
        UInt64Type globalElements = 0;
        ForLoop< LoopCallee<GlobalNumberOfElements>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, globalElements );
        if( Dune::Fem::Parameter::verbose () )
        {
          double min = std::numeric_limits< double >::max;
          double max = 0.0;
          ForLoop< LoopCallee<MinMaxNumElements>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, min, max );
           std::cout << "grid size (sum,min,max) = ( "
            << globalElements << " , " << min << " , " << max << ")" << std::endl;
        }
        return globalElements;
      }
      return 0;
    }

    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp )
    {
      ForLoop< LoopCallee<SetAdaptation>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
    }

    void finalize()
    {
      ForLoop< LoopCallee<Finalize>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
    }

    double& adaptationTime()
    {
      adaptationTime_ = adaptive() ? adaptationManager().adaptationTime() : 0.0;
      return adaptationTime_;
    }

    double& loadBalanceTime()
    {
      loadBalanceTime_ = adaptive() ? adaptationManager().loadBalanceTime() : 0.0;
      return loadBalanceTime_;
    }

    const int finestLevel() const
    {
      int finestLevel = 0;
      ForLoop< LoopCallee<FinestLevel>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, finestLevel );
      return finestLevel;
    }

  protected:
    GridType &grid () { return std::get< 0 >( tuple_ )->grid(); }

    AdaptationManagerType& adaptationManager()
    {
      if( !adaptationManager_ )
        adaptationManager_.reset( new AdaptationManagerType( grid(), *rp_ ) );
      return *adaptationManager_;
    }

    void estimateMark( const bool initialAdaptation = false )
    {
      if( adaptive() )
        ForLoop< LoopCallee<EstimateMark>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_, initialAdaptation );
    }

    void adapt()
    {
      if( adaptive() )
      {
        //int sequence = getSequence( get<0>( tuple_ ) );

        ForLoop< LoopCallee<PreAdapt>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
        adaptationManager().adapt();
        ForLoop< LoopCallee<PostAdapt>::template Apply, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );

        //TODO include limiterHandler
        //if( sequence !=  getSequence( get<0>( tuple_ ) ) )
          //limiterHandler_.step( get<0>( tuple_ )->adaptationSolution() );
      }
    }

  private:

    //template< class T >
    //static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value, int >::type
    //getSequence( T ){}
    //template< class T >
    //static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value, int >::type
    //getSequence( T elem )
    //{
    //  if( elem->adaptationSolution() )
    //    return elem->adaptationSolution()->space().sequence();
    //  return 0;
    //}

    StepperTupleType&                         tuple_;
    std::unique_ptr< RestrictionProlongationType > rp_;
    std::unique_ptr< AdaptationManagerType >  adaptationManager_;
    const std::string                         keyPrefix_;
    const AdaptationParametersType            adaptParam_;
    double                                    adaptationTime_;
    double                                    loadBalanceTime_;
  };


  template<>
  class AdaptHandler<>
  {
    typedef uint64_t                                                                           UInt64Type;
  public:

    template< class ... Args >
    AdaptHandler ( Args && ... ) {}

    template< class ... Args >
    bool adaptive( Args&& ... ) const { return false; }

    template< class ... Args >
    void setIndicator( Args&& ... ) {}

    template< class ... Args >
    void step( Args&& ... ) {}

    template< class ... Args >
    void init( Args&& ... ) {}

    template< class ... Args >
    size_t numberOfElements( Args&& ... ) const { return 0; }

    template< class ... Args >
    UInt64Type globalNumberOfElements( Args&& ... ) const { return 0; }

    template< class ... Args >
    void setAdaptation( Args&& ... ){}

    template< class ... Args >
    void finalize( Args&& ... ) {}

    template< class ... Args >
    const double adaptationTime( Args&& ... ) const { return 0.0; }

    template< class ... Args >
    const double loadBalanceTime( Args&& ... ) const { return 0.0; }

    template< class ... Args >
    const double finestLevel( Args&& ... ) const { return 0.0; }
  };


}
}
#endif
