#ifndef FEMDG_SUBDEFAULTADAPTCALLER_HH
#define FEMDG_SUBDEFAULTADAPTCALLER_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>

namespace Dune
{
namespace Fem
{
  /**
   * \brief Adaptation indicator for the estimation and marking of the entities.
   */
  template< class... IndicatorArgs >
  class AdaptIndicator;


  /**
   * \brief Adaptation indicator doing no indication and marking of the entities.
   */
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

    void preAdapt() {}

    void estimateMark( const bool initialAdapt = false ) {}

    void postAdapt() {}

    void finalize() {}

    int minNumberOfElements() const { return 0; }

    int maxNumberOfElements() const { return 0; }

    int finestLevel() const { return 0; }
  };


  /**
   * \brief Specialization doing estimation and marking of the entities for adaptation.
   *
   * \tparam IndicatorImp an adaptation handler
   * \tparam EstimatorImp an indicator using the Estimator interface
   */
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

    typedef typename IndicatorType :: AdvectionFluxType                                        AdvectionFluxType;

    // time provider
    typedef Fem::TimeProviderBase TimeProviderType;

    typedef AdaptationHandler< GridType, typename DiscreteFunctionSpaceType::FunctionSpaceType >
                                                                                               AdaptationHandlerType;
    typedef typename IndicatorType::ExtraParameterTupleType                                    ExtraParameterTupleType;

    template< class Model, class ExtraParameterTupleImp >
    AdaptIndicator( DiscreteFunctionType& sol, Model& model, const AdvectionFluxType& numFlux,
                    const ExtraParameterTupleImp& tuple,
                    const std::string& name = "AdaptIndDG",
                    const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : AdaptIndicator( sol.space(), model, numFlux, tuple, name, parameter )
    {
      sol_ = &sol;
    }

    template< class Model, class ExtraParameterTupleImp >
    AdaptIndicator( const DiscreteFunctionSpaceType& space, Model& model, const AdvectionFluxType& numFlux,
                    const ExtraParameterTupleImp& tuple,
                    const std::string& name = "AdaptIndDG",
                    const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : sol_( nullptr ),
      space_( space ),
      adaptationHandler_( nullptr ),
      adaptParam_( parameter ),
      // this is a DG operator (see operator/dg)
      indicator_( const_cast<GridPartType&>(space.gridPart()), model, numFlux, tuple, name, parameter ),
      estimator_( space, model.problem(), adaptParam_ ),
      timeProvider_( nullptr )
    {}

    /**
     * \brief Returns true if adaptive.
     */
    bool adaptive() const
    {
      return adaptParam_.adaptive();
    }

    /**
     * \brief Returns the number of elements.
     */
    size_t numberOfElements() const
    {
      return estimator_.numberOfElements();
    }

    /**
     * \brief Returns global number of elements.
     */
    UInt64Type globalNumberOfElements() const
    {
      if( adaptive() && adaptationHandler_ )
      {
        UInt64Type globalElements = adaptationHandler_->globalNumberOfElements() ;
        return globalElements;
      }
      return 0;
    }

    /**
     * \brief Set adaptation handler.
     *
     * \param[in] tp time provider
     */
    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp )
    {
      // create adaptation handler in case of apost indicator
      if( adaptive() )
      {
        timeProvider_ = &tp;
        if( !adaptationHandler_ && adaptParam_.aposterioriIndicator() )
        {
          adaptationHandler_.reset( new AdaptationHandlerType( const_cast<GridType&>(space_.gridPart().grid()) ) );
          indicator_.setAdaptation( *adaptationHandler_ );
        }
      }
    }

    /**
     * \brief Clean up indicator.
     */
    void finalize()
    {
      adaptationHandler_.reset( nullptr );
    }

    /**
     * \brief estimate the local error and mark all entities
     *
     * \param[in] initialAdapt true if this is the initial adaptation (needed for adaptation handler)
     */
    void estimateMark( const bool initialAdapt = false )
    {
      if( sol_ )
      {
        const double dt = timeProvider_ ? timeProvider_->deltaT() : 0;
        estimateMark( *sol_, dt, initialAdapt );
      }
      else
        DUNE_THROW(InvalidStateException,"AdaptIndicator: sol_ pointer is not set");
    }

    /**
     * \brief estimate the local error and mark all entities
     *
     * \param[in] initialAdapt true if this is the initial adaptation (needed for adaptation handler)
     */
    void estimateMark( const DiscreteFunctionType& solution, const double dt, const bool initialAdapt = false )
    {
      if( adaptive() )
      {
        if( useIndicator() )
        {
          //clear old values of indicator
          adaptationHandler_->resetStatus( dt );

          indicator_.evaluateOnly( solution );
          // mark all entities depending on error
          adaptationHandler_->markEntities( initialAdapt );
        }
        else if( useGradientIndicator() )
        {
          estimator_.estimateAndMark( solution );
        }
      }
    }

    /**
     * \brief Reset all indicators.
     */
    void postAdapt()
    {
      if( adaptive() )
        if( useIndicator() )
          adaptationHandler_->resetStatus();

    }

    /**
     * \brief Prepare estimator.
     *
     * Not used.
     */
    void preAdapt()
    {}

    /**
     * \brief Returns the minmal number of element.
     */
    int minNumberOfElements() const
    {
      assert( adaptationHandler_ );
      return adaptationHandler_->minNumberOfElements();
    }

    /**
     * \brief Returns the maximal number of element.
     */
    int maxNumberOfElements() const
    {
      assert( adaptationHandler_ );
      return adaptationHandler_->maxNumberOfElements();
    }

    /**
     * \brief Returns the maximum finest level of the grid.
     */
    int finestLevel() const
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

    DiscreteFunctionType*                         sol_;
    const DiscreteFunctionSpaceType&              space_;
    std::unique_ptr< AdaptationHandlerType >      adaptationHandler_;
    const std::string                             keyPrefix_;
    const AdaptationParametersType                adaptParam_;
    IndicatorType                                 indicator_;
    EstimatorType                                 estimator_;

    TimeProviderType*                             timeProvider_;
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
