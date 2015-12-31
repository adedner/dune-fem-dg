#ifndef FEMDG_DEFAULTADAPTHANDLER_HH
#define FEMDG_DEFAULTADAPTHANDLER_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/common/forloop.hh>
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
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>

namespace Dune
{
namespace Fem
{

  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class AdaptHandler;


  template< class AlgTupleImp, std::size_t... Ints >
  class AdaptHandler< AlgTupleImp, Std::index_sequence< Ints... > >
  {
    template< class TupleType > struct RPDefaultTupleExtractor;
    template< class ... Args > struct RPDefaultTupleExtractor< std::tuple< Args... > >
    { typedef Dune::Fem::RestrictProlongDefaultTuple< typename std::remove_pointer< Args >::type::DiscreteFunctionType... > type; };

    typedef AlgTupleImp                                                                        AlgTupleType;

    typedef Std::index_sequence< Ints... >                                                     IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                                    TupleReducerType;
    typedef typename TupleReducerType::type                                                    TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef uint64_t                                                                           UInt64Type;

    typedef typename std::remove_pointer< typename std::tuple_element< 0, TupleType >::type >::type::GridType
                                                                                               GridType;


    typedef typename RPDefaultTupleExtractor< TupleType >::type                                RestrictionProlongationType;

    typedef AdaptationManager< GridType, RestrictionProlongationType >                         AdaptationManagerType;

    typedef AdaptationParameters                                                               AdaptationParametersType;



    struct EstimateMark {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->estimateMark( std::forward<Args>(a)... ); }
    };
    struct SetAdaptation {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->setAdaptation( std::forward<Args>(a)... ); }
    };
    struct PreAdapt {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->preAdapt( std::forward<Args>(a)... ); }
    };
    struct PostAdapt {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->postAdapt( std::forward<Args>(a)... ); }
    };
    struct Finalize {
      template<class T, class... Args > static void apply( T e, Args&& ... a )
      { e->finalize( std::forward<Args>(a)... ); }
    };
    struct MinMaxNumElements {
      template<class T, class... Args > static void apply( T e, int& min, int& max, Args&& ... a )
      {
        min = std::min( min, e->minNumberOfElements( std::forward<Args>(a)... ) );
        max = std::max( max, e->maxNumberOfElements( std::forward<Args>(a)... ) );
      }
    };
    struct NumberOfElements {
      template<class T, class... Args > static void apply( T e, int& max, Args&& ... a )
      { max = std::max( max, e->numberOfElements( std::forward<Args>(a)... ) ); }
    };
    struct GlobalNumberOfElements {
      template<class T, class... Args > static void apply( T e, int& max, Args&& ... a )
      { max = std::max( max, e->globalNumberOfElements( std::forward<Args>(a)... ) ); }
    };
    struct FinestLevel {
      template<class T, class... Args > static void apply( T e, int& max, Args&& ... a )
      { max = std::max( max, e->finestLevel( std::forward<Args>(a)... ) ); }
    };
    struct Adaptive {
      template<class T, class... Args > static void apply( T e, bool& adaptive, Args&& ... a )
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
          C::apply(elem->adaptIndicator(), std::forward<Args>(a)... );
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




    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0, numAlgs - 1 >;



  public:

    AdaptHandler( AlgTupleType& tuple )
    : tuple_( TupleReducerType::apply( tuple ) ),
      rp_( nullptr ),
      adaptationManager_(),
      keyPrefix_( "" ),
      adaptParam_( AdaptationParametersType( ParameterKey::generate( keyPrefix_, "fem.adaptation." ) ) )
    {

      setRestrProlong( IndexSequenceType() );
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
      ForLoopType< Adaptive >::apply( tuple_, adaptive );
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
      ForLoopType< NumberOfElements >::apply( tuple_, numElements );
      return numElements;
    }

    UInt64Type globalNumberOfElements() const
    {
      if( adaptive() )
      {
        UInt64Type globalElements = 0;
        ForLoopType< GlobalNumberOfElements >::apply( tuple_, globalElements );
        if( Dune::Fem::Parameter::verbose () )
        {
          double min = std::numeric_limits< double >::max;
          double max = 0.0;
          ForLoopType< MinMaxNumElements >::apply( tuple_, min, max );
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
      ForLoopType< SetAdaptation >::apply( tuple_, tp );
    }

    void finalize()
    {
      ForLoopType< Finalize >::apply( tuple_ );
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
      ForLoopType< FinestLevel >::apply( tuple_, finestLevel );
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
        ForLoopType< EstimateMark >::apply( tuple_, initialAdaptation );
    }

    void adapt()
    {
      if( adaptive() )
      {
        //int sequence = getSequence( get<0>( tuple_ ) );

        ForLoopType< PreAdapt >::apply( tuple_ );
        adaptationManager().adapt();
        ForLoopType< PostAdapt >::apply( tuple_ );

        //TODO include limiterHandler
        //if( sequence !=  getSequence( get<0>( tuple_ ) ) )
          //limiterHandler_.step( get<0>( tuple_ )->adaptationSolution() );
      }
    }

  private:
    TupleType                                 tuple_;
    std::unique_ptr< RestrictionProlongationType > rp_;
    std::unique_ptr< AdaptationManagerType >  adaptationManager_;
    const std::string                         keyPrefix_;
    const AdaptationParametersType            adaptParam_;
    double                                    adaptationTime_;
    double                                    loadBalanceTime_;
  };


  template< class TupleImp >
  class AdaptHandler< TupleImp, Std::index_sequence<> >
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
