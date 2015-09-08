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

    template< int i >
    struct EstimateMark
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, Args && ... args )
      {
        if( elem->adaptIndicator() )
          elem->adaptIndicator()->estimateMark( args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), args... );
      }
    };
    template< int i >
    struct SetAdaptation
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, Args && ... args )
      {
        if( elem->adaptIndicator() )
          elem->adaptIndicator()->setAdaptation( args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), args... );
      }
    };
    template< int i >
    struct PreAdapt
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, Args && ... args )
      {
        if( elem->adaptIndicator() )
          elem->adaptIndicator()->preAdapt( args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), args... );
      }
    };
    template< int i >
    struct PostAdapt
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, Args && ... args )
      {
        if( elem->adaptIndicator() )
          elem->adaptIndicator()->postAdapt( args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), args... );
      }
    };
    template< int i >
    struct Finalize
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, Args && ... args )
      {
        if( elem->adaptIndicator() )
          elem->adaptIndicator()->finalize( args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), args... );
      }
    };
    template< int i >
    struct MinMaxNumElements
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, int& min, int& max, Args && ... args )
      {
        if( elem->adaptIndicator() )
        {
          min = std::min( min, elem->adaptIndicator()->minNumberOfElements( args... ) );
          max = std::max( max, elem->adaptIndicator()->maxNumberOfElements( args... ) );
        }
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, double& min, double& max, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), min, max, args... );
      }
    };
    template< int i >
    struct NumberOfElements
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, int& max, Args && ... args )
      {
        if( elem->adaptIndicator() )
          max = std::max( max, elem->adaptIndicator()->numberOfElements( args... ) );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, int& max, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), max, args... );
      }
    };
    template< int i >
    struct GlobalNumberOfElements
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, UInt64Type& max, Args && ... args )
      {
        if( elem->adaptIndicator() )
          max = std::max( max, elem->adaptIndicator()->globalNumberOfElements( args... ) );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, UInt64Type& max, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), max, args... );
      }
    };
    template< int i >
    struct FinestLevel
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, int& max, Args && ... args )
      {
        if( elem->adaptIndicator() )
          max = std::max( max, elem->adaptIndicator()->finestLevel( args... ) );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, int& max, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), max, args... );
      }
    };
    template< int i >
    struct Adaptive
    {
      template<class T, class... Args >
      static typename enable_if< std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T, Args&& ... ){}
      template<class T, class... Args >
      static typename enable_if< !std::is_void< typename std::remove_pointer<T>::type::AdaptIndicatorType >::value >::type
      adaptIndicator( T elem, bool& adaptive, Args && ... args )
      {
        if( elem->adaptIndicator() )
          adaptive |= elem->adaptIndicator()->adaptive( args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, bool& adaptive, Args && ... args )
      {
        adaptIndicator( std::get<i>( tuple ), adaptive, args... );
      }
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
      ForLoop< Adaptive, 0, sizeof ... ( StepperArg ) >::apply( tuple_, adaptive );
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
      ForLoop< NumberOfElements, 0, sizeof ... ( StepperArg ) >::apply( tuple_, numElements );
      return numElements;
    }

    UInt64Type globalNumberOfElements() const
    {
      if( adaptive() )
      {
        UInt64Type globalElements = 0;
        ForLoop< GlobalNumberOfElements, 0, sizeof ... ( StepperArg ) >::apply( tuple_, globalElements );
        if( Dune::Fem::Parameter::verbose () )
        {
          double min = std::numeric_limits< double >::max;
          double max = 0.0;
          ForLoop< MinMaxNumElements, 0, sizeof ... ( StepperArg ) >::apply( tuple_, min, max );
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
      ForLoop< SetAdaptation, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
    }

    void finalize()
    {
      ForLoop< Finalize, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
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
      ForLoop< FinestLevel, 0, sizeof ... ( StepperArg ) >::apply( tuple_, finestLevel );
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
        ForLoop< EstimateMark, 0, sizeof ... ( StepperArg ) >::apply( tuple_, initialAdaptation );
    }

    void adapt()
    {
      if( adaptive() )
      {
        //int sequence = get<0>( tuple_ )->adaptationSolution()->space().sequence();

        ForLoop< PreAdapt, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
        adaptationManager().adapt();
        ForLoop< PostAdapt, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );

        //TODO include limiterHandler
        //if( sequence !=  get<0>( tuple_ )->adaptationSolution()->space().sequence() )
          //limiterHandler_.step( get<0>( tuple_ )->adaptationSolution() );
      }
    }

  private:

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
