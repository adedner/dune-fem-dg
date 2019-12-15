#ifndef FEMDG_DEFAULTADAPTCALLER_HH
#define FEMDG_DEFAULTADAPTCALLER_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/algorithm/caller/postprocessing.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/misc/integral_constant.hh>
#include "interface.hh"


namespace Dune
{
namespace Fem
{

  /**
   * \brief Caller class managing the adaptation process.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename std::make_index_sequence< std::tuple_size< AlgTupleImp >::value > >
  class AdaptCaller;


  /**
   * \brief Specialization of a caller class managing the adaptation process.
   *
   * \ingroup Callers
   *
   * This class manages the adaptation process for a tuple of sub-algorithms.
   * For each sub-algorithm adaptation can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef AdaptCaller< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                       std::index_sequence< 0, 2 > >
   *                                           MyCaller;
   * \endcode
   * This would enable adaptation for `Alg1` and `Alg3`.
   *
   * \warning Since all sub-algorithms depends on the same (refined) grid, a
   * partial selection of adaptive algorithms would not work.
   *
   * \todo Improve the situation for partially selected adaptive algorithm tuples.
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam std::index_sequence< Ints... > Index sequence for enabling the checkpointing feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class AdaptCaller< AlgTupleImp, std::index_sequence< Ints... > >
    : public CallerInterface
  {
    template< class TupleType > struct RPDefaultTupleExtractor;
    template< class ... Args > struct RPDefaultTupleExtractor< std::tuple< Args... > >
    { typedef Dune::Fem::RestrictProlongDefaultTuple< typename Args::element_type::DiscreteFunctionType... > type; };

    typedef AlgTupleImp                                                                        AlgTupleType;

    typedef std::index_sequence< Ints... >                                                     IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                                    TupleReducerType;
    typedef typename TupleReducerType::type                                                    TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef uint64_t                                                                           UInt64Type;

    typedef typename std::tuple_element_t< 0, TupleType >::element_type::GridType              GridType;


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
      static typename std::enable_if< std::is_void< typename T::element_type::AdaptIndicatorType >::value >::type
      getAdaptIndicator( T&, Args&& ... ){}
      template<class C, class T, class... Args >
      static typename std::enable_if< !std::is_void< typename T::element_type::AdaptIndicatorType >::value >::type
      getAdaptIndicator( T& elem, Args &&... a )
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

    /**
     * \brief Constructor.
     *
     * \param tuple Tuple of all sub-algorithms.
     */
    AdaptCaller( AlgTupleType& tuple )
    : tuple_( TupleReducerType::apply( tuple ) ),
      rp_( nullptr ),
      adaptationManager_(),
      keyPrefix_( "" ),
      adaptParam_( AdaptationParametersType( ParameterKey::generate( keyPrefix_, "fem.adaptation." ) ) ),
      adaptationTime_(0),
      loadBalanceTime_(0)
    {
      const auto& elem = std::get<0>( tuple_ );
      std::cout << elem.use_count() << std::endl;
      //const auto& elemRay = *elem;
      const auto& sol = elem->solution();
      std::cout << "new container: name: " << sol.name() << ", size: " << sol.size() << std::endl;

      setRestrProlong( IndexSequenceType() );
      if( adaptive() )
        rp_->setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
    }


    /**
     * \brief Prepare an initial refined grid.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void initializeEnd( AlgImp* alg, int loop, TimeProviderImp& tp)
    {
      if( adaptive() )
      {
        // adapt the grid to the initial data
        for( int startCount = 0; startCount < finestLevel(); ++ startCount )
        {
          // call initial adaptation
          estimateMark( true );
          adapt( alg, loop, tp );

          // setup problem again
          alg->initialize( loop, tp );

          // some info in verbose mode
          if( Fem::Parameter::verbose() )
          {
            if( tp.timeStepValid() )
            {
              std::cout << "Start adaptation: step " << startCount << ",  dt = " << tp.deltaT() << ",  grid size: " << alg->gridSize()
                        << std::endl;
            }
          }
        }
      }
    }

    /**
     * \brief Prepare an initial refined grid.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     */
    template< class AlgImp >
    void initializeEnd( AlgImp* alg, int loop )
    {
      if( adaptive() )
      {
        // some info in verbose mode
        if( Fem::Parameter::verbose() )
        {
          std::cout << "Start adaptation: grid size: " << alg->gridSize()
                    << std::endl;
        }
      }
    }

    /**
     * \brief Calls the estimate, mark and adaptation routines to refine the grid.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void solveStart( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( needsAdaptation( alg, loop, tp ) )
      {
        ForLoopType< PreAdapt >::apply( tuple_ );
        estimateMark( false );
        adapt( alg, loop, tp );
        ForLoopType< PostAdapt >::apply( tuple_ );
      }
    }

    /**
     * \brief Calls the estimate, mark and adaptation routines to refine the grid.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     */
    template< class AlgImp >
    void solveStart( AlgImp* alg, int loop )
    {
      if( adaptive() )
      {
        ForLoopType< PreAdapt >::apply( tuple_ );
        estimateMark();
        adaptationManager().adapt();
        ForLoopType< PostAdapt >::apply( tuple_ );
      }
    }

    /**
     * \brief finalize all indicators
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void finalizeStart( AlgImp* alg, int loop, TimeProviderImp& tp)
    {
      ForLoopType< Finalize >::apply( tuple_ );
    }

    /**
     * \brief finalize all indicators
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     */
    template< class AlgImp >
    void finalizeStart( AlgImp* alg, int loop )
    {
      ForLoopType< Finalize >::apply( tuple_ );
    }

    /**
     * \brief Returns true, if all sub-algorithms are adaptive.
     */
    bool adaptive () const
    {
      bool adaptive = false;
      ForLoopType< Adaptive >::apply( tuple_, adaptive );
      return adaptive;
    }


    template< class AlgImp, class TimeProviderImp >
    bool needsAdaptation( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( !tp.timeStepValid() )
        return false;
      return( tp.timeStep() % adaptParam_.adaptCount() == 0 );
    }

    template< class AlgImp >
    bool needsAdaptation( AlgImp* alg, int loop )
    {
      //TODO set a better condition, here (from padaptindicator?)...
      static int maxIteration = 0;
      maxIteration++;
      return( maxIteration < adaptParam_.coarsestLevel() );
    }

    /**
     * \brief Returns the number of elements.
     */
    size_t numberOfElements() const
    {
      int numElements = 0;
      ForLoopType< NumberOfElements >::apply( tuple_, numElements );
      return numElements;
    }

    /**
     * \brief Returns the number of elements.
     */
    UInt64Type globalNumberOfElements() const
    {
      if( adaptive() )
      {
        UInt64Type globalElements = 0;
        ForLoopType< GlobalNumberOfElements >::apply( tuple_, globalElements );
        if( Dune::Fem::Parameter::verbose () )
        {
          double min = std::numeric_limits< double >::max();
          double max = 0.0;
          ForLoopType< MinMaxNumElements >::apply( tuple_, min, max );
           std::cout << "grid size (sum,min,max) = ( "
            << globalElements << " , " << min << " , " << max << ")" << std::endl;
        }
        return globalElements;
      }
      return 0;
    }

    /**
     * \brief Set adaptation manager of sub-algorithms.
     *
     * \param[in] tp time provider
     */
    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp )
    {
      ForLoopType< SetAdaptation >::apply( tuple_, tp );
    }

    /**
     * \brief Set adaptation manager of sub-algorithms.
     */
    void setAdaptation()
    {
      ForLoopType< SetAdaptation >::apply( tuple_ );
    }

    /**
     * \brief Returns the adaptation time used by adaptation manager.
     */
    double& adaptationTime()
    {
      adaptationTime_ = adaptive() ? adaptationManager().adaptationTime() : 0.0;
      return adaptationTime_;
    }

    /**
     * \brief Returns the load balancing time.
     */
    double& loadBalanceTime()
    {
      loadBalanceTime_ = adaptive() ? adaptationManager().loadBalanceTime() : 0.0;
      return loadBalanceTime_;
    }

  protected:
    template< std::size_t ... i >
    void setRestrProlong( std::index_sequence< i ... > )
    {
      const auto& elem = std::get<0>( tuple_ );
      std::cout << elem.use_count() << std::endl;
      auto& elemRay = *elem;
      const auto& sol = elemRay.solution();
      std::cout << "new container: name: " << sol.name() << ", size: " << sol.size() << std::endl;
      std::cout << "new container: name: " << elemRay.adaptationSolution()->name() << ", size: " << elemRay.adaptationSolution()->size() << std::endl;
      rp_.reset( new RestrictionProlongationType( *(std::get< i >( tuple_ )->adaptationSolution() )... ) );
    }

    int finestLevel() const
    {
      //int res = 0;
      //Dune::Hybrid::forEach(tuple_,
      //[&](auto i)
      //{
      //  typedef typename std::is_void< typename std::tuple_element_t<i,TupleType>::element_type::AdaptIndicatorType >::type hasAdaptIndicatorType;
      //  Dune::Hybrid::ifElse( hasAdaptIndicatorType(), [&](auto i){ res=std::max(res,std::get<i>(tuple_)->finestLevel()); } );
      //} );
      //return res;
      int finestLevel = 0;
      ForLoopType< FinestLevel >::apply( tuple_, finestLevel );
      return finestLevel;
    }

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

    template< class AlgImp, class TimeProviderImp >
    void adapt( AlgImp* alg, int loop, TimeProviderImp& tp)
    {
      if( adaptive() )
      {
        int sequence = getSequence( std::get<0>( tuple_ ) );

        adaptationManager().adapt();

        if( sequence !=  getSequence( std::get<0>( tuple_ ) ) )
          alg->postProcessing().solveEnd( alg, loop, tp );
      }
    }


    template< class T >
    static typename std::enable_if< std::is_void< typename T::element_type::AdaptIndicatorType >::value, int >::type
    getSequence( T& ){}
    template< class T >
    static typename std::enable_if< !std::is_void< typename T::element_type::AdaptIndicatorType >::value, int >::type
    getSequence( T& elem )
    {
      if( elem->adaptationSolution() )
        return elem->adaptationSolution()->space().sequence();
      return 0;
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


  /**
   * \brief Specialization of a caller class without adaptation.
   *
   * \ingroup Callers
   */
  template< class TupleImp >
  class AdaptCaller< TupleImp, std::index_sequence<> >
    : public CallerInterface
  {
    typedef uint64_t                          UInt64Type;
  public:

    template< class ... Args >
    AdaptCaller ( Args && ... ) {}

    template< class ... Args >
    bool adaptive( Args&& ... ) const { return false; }

    template< class ... Args >
    size_t numberOfElements( Args&& ... ) const { return 0; }

    template< class ... Args >
    UInt64Type globalNumberOfElements( Args&& ... ) const { return 0; }

    template< class ... Args >
    void setAdaptation( Args&& ... ){}

    template< class ... Args >
    double adaptationTime( Args&& ... ) const { return 0.0; }

    template< class ... Args >
    double loadBalanceTime( Args&& ... ) const { return 0.0; }
  };

}
}
#endif
