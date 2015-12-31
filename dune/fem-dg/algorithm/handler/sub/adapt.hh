#ifndef FEMDG_SUBDEFAULTADAPTHANDLER_HH
#define FEMDG_SUBDEFAULTADAPTHANDLER_HH

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

  template< class... IndicatorArgs >
  class AdaptIndicator;

  template< class IndicatorImp, class EstimatorImp >
  class AdaptIndicator< IndicatorImp, EstimatorImp >
  {
  public:
    typedef uint64_t                                                                           UInt64Type;

    typedef IndicatorImp                                                                       IndicatorType; //using adaptationhandler
    typedef EstimatorImp                                                                       EstimatorType; //estimator interface
    typedef typename IndicatorType::DestinationType                                            DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType                           DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType                                   GridPartType;
    typedef typename GridPartType::GridType                                                    GridType;
    typedef AdaptationParameters                                                               AdaptationParametersType;

    typedef AdaptationHandler< GridType, typename DiscreteFunctionSpaceType::FunctionSpaceType >
                                                                                               AdaptationHandlerType;

    template< class Problem, class ExtraTupleParameter >
    AdaptIndicator( DiscreteFunctionType& sol, Problem& problem, const ExtraTupleParameter& tuple, const std::string keyPrefix = "" )
    : sol_( sol ),
      adaptationHandler_( nullptr ),
      keyPrefix_( keyPrefix ),
      adaptParam_( AdaptationParametersType( ParameterKey::generate( keyPrefix, "fem.adaptation." ) ) ),
      indicator_( const_cast<GridPartType&>(sol_.gridPart()), problem, tuple, keyPrefix_ ),
      estimator_( sol_.space(), problem, adaptParam_ )
    {}

    bool adaptive() const
    {
      return adaptParam_.adaptive();
    }

    size_t numberOfElements() const
    {
      return estimator_.numberOfElements();
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
          estimator_.estimateAndMark( sol_ );
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
    EstimatorType                                 estimator_;
  };


  template<>
  class AdaptIndicator<>
  {
  public:
   typedef uint64_t                          UInt64Type;

    template< class... Args >
    AdaptIndicator( Args&& ... args )
    {}

    bool adaptive() const { return false; }

    size_t numberOfElements() const { return 0; }

    UInt64Type globalNumberOfElements() const { return 0; }

    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp ){}

    void finalize() {}

    void estimateMark( const bool initialAdapt = false ) {}

    void postAdapt() {}

    void preAdapt() {}

    int minNumberOfElements() const { return 0; }

    int maxNumberOfElements() const { return 0; }

    const int finestLevel() const { return 0; }

  };


  template< class Obj >
  class AdaptIndicatorOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    AdaptIndicatorOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

  template<>
  class AdaptIndicatorOptional< void >
    : public OptionalNullPtr< AdaptIndicator<> >
  {
    typedef OptionalNullPtr< AdaptIndicator<> >    BaseType;
  public:
    template< class... Args >
    AdaptIndicatorOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };



}
}
#endif
