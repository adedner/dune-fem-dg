#ifndef FEMDG_DEFAULTADAPTHANDLER_HH
#define FEMDG_DEFAULTADAPTHANDLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem/space/common/adaptmanager.hh>

namespace Dune
{
namespace Fem
{


  template< class GridPartImp, class DiscreteFunctionImp, class RestrictionProlongationImp, class IndicatorImp,
            class GradientIndicatorImp, class SolutionLimiterHandlerImp = NoSolutionLimiterHandler >
  class DefaultAdaptHandler
  {
    typedef GridPartImp                                                                        GridPartType;
    typedef typename GridPartImp::GridType                                                     GridType;
    typedef Dune::AdaptationHandler< GridType,
              typename DiscreteFunctionImp::DiscreteFunctionSpaceType::FunctionSpaceType >     AdaptationHandlerType;
    typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationImp >               AdaptationManagerType;
    typedef uint64_t                                                                           UInt64Type;
    typedef Dune::AdaptationParameters                                                         AdaptationParametersType;
    typedef IndicatorImp                                                                       IndicatorType;
    typedef GradientIndicatorImp                                                               GradientIndicatorType;
  public:

    DefaultAdaptHandler( GridPartType& gridPart, DiscreteFunctionImp& sol,
                         SolutionLimiterHandlerImp& limiterHandler, const std::string keyPrefix = "" )
    : gridPart_( gridPart ),
      sol_( sol ),
      adaptationHandler_(),
      rp_( sol ),
      adaptationManager_(),
      indicator_( nullptr ),
      gradientIndicator_( nullptr ),
      keyPrefix_( keyPrefix ),
      adaptParam_( AdaptationParametersType( Dune::ParameterKey::generate( keyPrefix, "fem.adaptation." ) ) ),
      limiterHandler_( limiterHandler )
    {
      if( adaptive() )
        rp_.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
    }

    bool adaptive() const
    {
      return true;
    }

    template< class Problem, class ExtraTupleParameter >
    void setIndicator( Problem& problem, const ExtraTupleParameter& tuple )
    {
      indicator_.reset( new IndicatorType( gridPart_, problem, tuple, keyPrefix_ ) );
      gradientIndicator_.reset( new GradientIndicatorType( sol_.space(), problem, adaptParam_ ) );
    }


    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      if( tp.timeStep() % adaptParam_.adaptCount() == 0 )
        doEstimateMarkAdapt( false );
    }

    void init()
    {
      doEstimateMarkAdapt( true );
    }

    size_t numberOfElements() const
    {
      if( gradientIndicator_ )
        return gradientIndicator_->numberOfElements();
      return 0;
    }

    UInt64Type globalNumberOfElements() const
    {
      if( adaptive() )
      {
        UInt64Type globalElements = adaptationHandler_->globalNumberOfElements() ;
        if( Dune::Fem::Parameter::verbose () )
        {
           std::cout << "grid size (sum,min,max) = ( "
            << globalElements << " , "
            << adaptationHandler_->minNumberOfElements() << " , "
            << adaptationHandler_->maxNumberOfElements() << ")" << std::endl;
        }
        return globalElements;
      }
      return 0;
    }

    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp )
    {
      // create adaptation handler in case of apost indicator
      if( adaptive() && indicator_ )
      {
        if( !adaptationHandler_ && adaptParam_.aposterioriIndicator() )
        {
          adaptationHandler_.reset( new AdaptationHandlerType( gridPart_.grid(), tp ) );
          indicator_->setAdaptation( *adaptationHandler_ );
        }
      }
    }

    void finalize()
    {
      adaptationHandler_.reset();
    }

    const double adaptationTime()
    {
      return adaptive() ? adaptationManager().adaptationTime() : 0.0;
    }

    const double loadBalanceTime()
    {
      return adaptive() ? adaptationManager().loadBalanceTime() : 0.0;
    }

    const double finestLevel() const
    {
      return adaptParam_.finestLevel();
    }

  protected:
    AdaptationManagerType& adaptationManager()
    {
      if( !adaptationManager_ )
        adaptationManager_.reset( new AdaptationManagerType( gridPart_.grid(), rp_ ) );
      return *adaptationManager_;
    }

    void doEstimateMarkAdapt( const bool initialAdaptation = false )
    {
      if( adaptive() )
      {
        // get grid sequence before adaptation
        const int sequence = sol_.space().sequence();

        if( adaptationHandler_ && indicator_ )
        {
          // call operator once to calculate indicator
          indicator_->evaluateOnly( sol_ );

          // do marking and adaptation
          adaptationHandler_->adapt( adaptationManager(), initialAdaptation );
        }
        else if( adaptParam_.gradientBasedIndicator() && gradientIndicator_ )
        {
          gradientIndicator_->estimateAndMark( sol_ );
          adaptationManager().adapt();
        }
        else if( adaptParam_.shockIndicator() )
        {
          // marking has been done by limiter
          adaptationManager().adapt();
        }

        // if grid has changed then limit solution again
        if( sequence != sol_.space().sequence() )
        {
          limiterHandler_.step( sol_ );
        }
      }
    }

  private:

    GridPartImp&                              gridPart_;
    DiscreteFunctionImp&                      sol_;
    std::unique_ptr< AdaptationHandlerType >  adaptationHandler_;
    RestrictionProlongationImp                rp_;
    std::unique_ptr< AdaptationManagerType >  adaptationManager_;
    std::unique_ptr< IndicatorType >          indicator_;
    std::unique_ptr< GradientIndicatorType >  gradientIndicator_;
    const std::string                         keyPrefix_;
    const AdaptationParametersType            adaptParam_;
    SolutionLimiterHandlerImp&                limiterHandler_;
  };

  class NoAdaptHandler
  {
    typedef uint64_t                                                                           UInt64Type;
  public:
    template< class ... Args >
    NoAdaptHandler( Args&& ... )
    {}

    template< class ... Args >
    bool adaptive( Args&& ... ) const { return false; }

    template< class ... Args >
    void setIndicator( Args&& ... ) {}

    template< class ... Args >
    void step( Args&& ... ) {}

    template< class ... Args >
    void init( Args&& ... ) {}

    template< class ... Args >
    size_t numberOfElements( Args&& ... ) const
    { return 0; }

    template< class ... Args >
    UInt64Type globalNumberOfElements( Args&& ... ) const
    { return 0; }

    template< class ... Args >
    void setAdaptation( Args&& ... ){}

    template< class ... Args >
    void finalize( Args&& ... ) {}


    template< class ... Args >
    const double adaptationTime() const { return 0.0; }

    template< class ... Args >
    const double loadBalanceTime() const { return 0.0; }

    template< class ... Args >
    const double finestLevel() const { return 0.0; }
  };

}
}
#endif
