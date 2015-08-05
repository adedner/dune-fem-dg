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


  template< class GridImp, class DiscreteFunctionImp, class RestrictionProlongationImp, class IndicatorTupleImp, class SolutionLimiterHandlerImp = NoSolutionLimiterHandler >
  class DefaultAdaptHandler
  {
  public:
    typedef Dune::AdaptationHandler< GridImp,
              typename DiscreteFunctionImp::DiscreteFunctionSpaceType::FunctionSpaceType >     AdaptationHandlerType;
    typedef Dune::Fem::AdaptationManager< GridImp, RestrictionProlongationImp >                AdaptationManagerType;
    typedef uint64_t                                                                           UInt64Type;
    typedef Dune::AdaptationParameters                                                         AdaptationParametersType;

    DefaultAdaptHandler( GridImp& grid, DiscreteFunctionImp& sol,
                         SolutionLimiterHandlerImp& limiterHandler, const std::string keyPrefix = "" )
    : grid_( grid ),
      sol_( sol ),
      adaptationHandler_(),
      rp_( sol ),
      adaptationManager_(),
      indicatorTuple_( std::make_tuple( nullptr, nullptr ) ),
      keyPrefix_( keyPrefix ),
      adaptParam_( AdaptationParametersType( Dune::ParameterKey::generate( keyPrefix, "fem.adaptation." ) ) ),
      limiterHandler_( limiterHandler )
    {
      if( adaptive() )
        rp_.setFatherChildWeight( Dune::DGFGridInfo<GridImp> :: refineWeight() );
    }

    AdaptationManagerType& adaptationManager()
    {
      if( !adaptationManager_ )
        adaptationManager_.reset( new AdaptationManagerType( grid_, rp_ ) );
      return *adaptationManager_;
    }

    bool adaptive() const
    {
      return true;
    }

    template< class IndicatorImp, class GradientIndicatorImp >
    void setIndicator( IndicatorImp* ind, GradientIndicatorImp* gradInd )
    {
      indicatorTuple_ = std::make_tuple( ind, gradInd );
    }

    void step()
    {
      if( adaptive() )
        doEstimateMarkAdapt( *get<0>(indicatorTuple_), *get<1>(indicatorTuple_), false );
    }

    void init()
    {
      if( adaptive() )
        doEstimateMarkAdapt( *get<0>(indicatorTuple_), *get<1>(indicatorTuple_), true );
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

    template< class TimeProviderImp, class IndicatorImp >
    void setAdaptation( TimeProviderImp& tp, IndicatorImp& indicator )
    {
      // create adaptation handler in case of apost indicator
      if( adaptive() )
      {
        if( ! adaptationHandler_ && adaptParam_.aposterioriIndicator() )
        {
          adaptationHandler_.reset( new AdaptationHandlerType( grid_, tp ) );
          indicator.setAdaptation( *adaptationHandler_ );
        }
      }
    }

    void finalize()
    {
      adaptationHandler_.reset();
    }

    const AdaptationParametersType& params() const
    { return adaptParam_; }

  protected:

    template <class IndicatorOperator, class GradientIndicator>
    void doEstimateMarkAdapt( const IndicatorOperator& dgIndicator,
                              GradientIndicator& gradientIndicator,
                              const bool initialAdaptation = false )
    {
      if( adaptive() )
      {
        // get grid sequence before adaptation
        const int sequence = sol_.space().sequence();

        if( adaptationHandler_ )
        {
          // call operator once to calculate indicator
          dgIndicator.evaluateOnly( sol_ );

          // do marking and adaptation
          adaptationHandler_->adapt( adaptationManager(), initialAdaptation );
        }
        else if( adaptParam_.gradientBasedIndicator() )
        {
          gradientIndicator.estimateAndMark( sol_ );
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

    GridImp&                                  grid_;
    DiscreteFunctionImp&                      sol_;
    std::unique_ptr< AdaptationHandlerType >  adaptationHandler_;
    RestrictionProlongationImp                rp_;
    std::unique_ptr< AdaptationManagerType >  adaptationManager_;
    IndicatorTupleImp                         indicatorTuple_;
    const std::string                         keyPrefix_;
    const AdaptationParametersType            adaptParam_;
    SolutionLimiterHandlerImp&                limiterHandler_;
  };

  template< class GridImp, class RestrictionProlongationImp >
  class NoAdaptHandler
  {
    public:
    typedef Dune::Fem::AdaptationManager< GridImp, RestrictionProlongationImp >                          AdaptationManagerType;
    typedef uint64_t UInt64Type;
    typedef Dune::AdaptationParameters                                                         AdaptationParametersType;

    template< class Arg1, class Arg2, class Arg3, class ... Args >
    NoAdaptHandler( Arg1&, Arg2&, Arg3, const std::string keyPrefix = "",Args&& ... )
    : adaptationManager_(nullptr),
      adaptParam_( AdaptationParametersType( Dune::ParameterKey::generate( keyPrefix, "fem.adaptation." ) ) )
    {}

    template< class ... Args >
    NoAdaptHandler( Args&& ... )
    : adaptationManager_(nullptr),
      adaptParam_( AdaptationParametersType( Dune::ParameterKey::generate( "", "fem.adaptation." ) ) )
    {}

    AdaptationManagerType& adaptationManager() const
    {
      return *adaptationManager_;
    }

    template< class ... Args >
    bool adaptive( Args&& ... ) const { return false; }

    template< class ... Args >
    void setIndicator( Args&& ... ) {}

    template< class ... Args >
    void step( Args&& ... ) {}

    template< class ... Args >
    void init( Args&& ... ) {}

    template< class ... Args >
    UInt64Type globalNumberOfElements( Args&& ... ) const
    { return 0; }

    template< class ... Args >
    void setAdaptation( Args&& ... ){}

    template< class ... Args >
    void finalize( Args&& ... ) {}

    template< class ... Args >
    const AdaptationParametersType& params() const { return adaptParam_; }


    private:
    mutable std::unique_ptr< AdaptationManagerType >  adaptationManager_;
    const AdaptationParametersType            adaptParam_;
  };

}
}
#endif
