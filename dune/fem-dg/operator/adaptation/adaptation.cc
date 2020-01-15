#ifndef DUNE_ADAPTATIONOBJECT_CC
#define DUNE_ADAPTATIONOBJECT_CC

#include "adaptation.hh"
#include <dune/fem/misc/gridwidth.hh>

namespace Dune
{
namespace Fem
{

  //! constructor
  template< class GridImp, class FunctionSpace >
  AdaptationHandler< GridImp, FunctionSpace >::
  AdaptationHandler ( GridType &grid,
                      const AdaptationParameters &param )
    : ComputeMinMaxVolume( GridPartType( grid ),
                           param.coarsestLevel( DGFGridInfo< GridType >::refineStepsForHalf() ),
                           param.finestLevel  ( DGFGridInfo< GridType >::refineStepsForHalf() ) )
    , grid_( grid )
    , gridPart_( grid_ )
#ifdef USE_ALUGRID_MPACCESS
    , mpAccess_( MPIHelper::getCommunicator() )
#endif
    , indicator_( grid_, 0, 0.0 )
    , globalTolerance_( param.refinementTolerance() )
    , coarsenTheta_( param.coarsenPercentage() )
    // make intial error count as 2.5 percent of the total error
    , initialTheta_( 0.025 )
    , globalNumElements_( 0 )
    , localNumElements_( 0 )
    , endTime_( param.endTime() )
    , deltaT_( 0 )
    , verbose_( Fem::Parameter::verbose() && param.verbose() )
  {
    const bool verboseOutput = Fem::Parameter::verbose();

    //std::cout << finestVolume() << "  " << coarsestVolume() << std::endl;

    resetStatus();
    if( verboseOutput )
      std::cout << "AdaptationHandler created! \n";

    /*
       // calculate global min of grid width to scale tolerance
       double gridWidth = GridWidth::calcGridWidth( gridPart_ );

       // get global minimum of macro grid width
       double macroGridWidth = grid_.comm().min( gridWidth );

       // globalTolerance_ *= macroGridWidth;
     */

    // scale tolerance with domain volume
    // globalTolerance_ *= volumeOfDomain();
  }

  //! constructor
  template< class GridImp, class FunctionSpace >
  AdaptationHandler< GridImp, FunctionSpace >::
  AdaptationHandler ( const AdaptationHandler &other )
    : grid_( other.grid_ )
    , gridPart_( grid_ )
    , indicator_( other.indicator_ )
    , globalTolerance_( other.globalTolerance_ )
    , coarsenTheta_( other.coarsenTheta_ )
    , initialTheta_( other.initialTheta_ )
    , globalNumElements_( other.globalNumElements_ )
    , localNumElements_( other.localNumElements_ )
    , endTime_( other.endTime_ )
    , deltaT_( other.deltaT_ )
    , verbose_( other.verbose_ )
  {}

  //! clear indicator
  template< class GridImp, class FunctionSpace >
  double
  AdaptationHandler< GridImp, FunctionSpace >::
  volumeOfDomain () const
  {
    double volume = 0;
    // type of iterator, i.e. leaf iterator

    for( const auto& entity : elements( gridPart_, Dune::Partitions::interior ) )
    {
      // sum up the volume
      volume += entity.geometry().volume();
    }

    // return volume of computational domain
    return grid_.comm().sum( volume );
  }

  //! clear indicator
  template< class GridImp, class FunctionSpace >
  void
  AdaptationHandler< GridImp, FunctionSpace >::
  clearIndicator ()
  {
    indicator_.resize();
    indicator_.fill( 0.0 );
  }

  //! add another adaptation handlers indicator container
  template< class GridImp, class FunctionSpace >
  AdaptationHandler< GridImp, FunctionSpace > &
  AdaptationHandler< GridImp, FunctionSpace >::
  operator+= ( const ThisType &other )
  {
    // add all indicator entries
    {
      typedef typename IndicatorType::Iterator IteratorType;
      typedef typename IndicatorType::ConstIterator ConstIteratorType;
      const IteratorType endit = indicator_.end();
      IteratorType it = indicator_.begin();

      ConstIteratorType oIt = other.indicator_.begin();
      for( IteratorType it = indicator_.begin(); it != endit; ++it, ++oIt )
        (*it) += *oIt;
    }
    return *this;
  }

  //! initialize localIndicator with en
  template< class GridImp, class FunctionSpace >
  template< class Entity >
  typename AdaptationHandler< GridImp, FunctionSpace >::LocalIndicatorType
  AdaptationHandler< GridImp, FunctionSpace >::
  localIndicator ( const Entity &entity )
  {
    // convert the given entity to an entity of the grid
    // for wrapped entities the cast to the host entity is necessary
    const GridEntityType &gridEntity = Fem::gridEntity( entity );
    return LocalIndicatorType( this, &indicator_[ gridEntity ] );
  }


  //! add value to local indicator, use setEntity before
  template< class GridImp, class FunctionSpace >
  void
  AdaptationHandler< GridImp, FunctionSpace >::
  addToLocalIndicator ( LocalIndicatorDataType &indicator,
                        const FullRangeType &error, const double h ) const
  {
    //const double dt = deltaT();
    //const double factor = ( h + dt  ) * dt;
    //indicator += (factor * error.two_norm());
    indicator += error.two_norm();
  }

  template< class GridImp, class FunctionSpace >
  void
  AdaptationHandler< GridImp, FunctionSpace >::
  addToLocalIndicator ( const GridEntityType &en, const FullRangeType &error, const double h )
  {
    addToLocalIndicator( indicator_[ en ], error, h );
  }

  template< class GridImp, class FunctionSpace >
  void
  AdaptationHandler< GridImp, FunctionSpace >::
  setLocalIndicator ( const GridEntityType &en, const FullRangeType &error )
  {
    assert( singleThreadMode() );
    indicator_[ en ] = error[ 0 ];
  }

  template< class GridImp, class FunctionSpace >
  double
  AdaptationHandler< GridImp, FunctionSpace >::
  getLocalIndicator ( const GridEntityType &en ) const
  {
    return indicator_[ en ];
  }

  //! calculate sum of local errors
  template< class GridImp, class FunctionSpace >
  double
  AdaptationHandler< GridImp, FunctionSpace >::
  getSumEstimator () const
  {
    assert( singleThreadMode() );
    double sum = 0.0;

    typedef typename IndicatorType::ConstIterator IteratorType;
    const IteratorType endit = indicator_.end();
    for( IteratorType it = indicator_.begin(); it != endit; ++it )
      sum += *it;

    // global sum of estimator
    sum = grid_.comm().sum( sum );
    return sum;
  }

  //! calculate max of local errors
  template< class GridImp, class FunctionSpace >
  double
  AdaptationHandler< GridImp, FunctionSpace >::
  getMaxEstimator () const
  {
    assert( singleThreadMode() );
    double max = 0.0;

    {
      typedef typename IndicatorType::ConstIterator IteratorType;
      const IteratorType endit = indicator_.end();
      IteratorType it = indicator_.begin();

      // initialzie with first entry
      if( it != endit )
        max = *it;

      for(; it != endit; ++it )
        max = std::max( *it, max );
    }

    // global sum of estimator
    return grid_.comm().max( max );
  }

  template< class GridImp, class FunctionSpace >
  int
  AdaptationHandler< GridImp, FunctionSpace >::
  localNumberOfElements () const
  {
    assert( localNumElements_ > 0 );
    return localNumElements_;
  }

  template< class GridImp, class FunctionSpace >
  typename AdaptationHandler< GridImp, FunctionSpace >::UInt64Type
  AdaptationHandler< GridImp, FunctionSpace >::
  globalNumberOfElements () const
  {
    assert( globalNumElements_ > 0 );
    return globalNumElements_;
  }

  template< class GridImp, class FunctionSpace >
  size_t
  AdaptationHandler< GridImp, FunctionSpace >::
  minNumberOfElements () const
  {
    return minNumElements_;
  }

  template< class GridImp, class FunctionSpace >
  size_t
  AdaptationHandler< GridImp, FunctionSpace >::
  maxNumberOfElements () const
  {
    return maxNumElements_;
  }

  template< class GridImp, class FunctionSpace >
  double
  AdaptationHandler< GridImp, FunctionSpace >::
  getLocalInTimeTolerance () const
  {
    double dt = deltaT();
    return (1. - initialTheta_) * globalTolerance_ * globalTolerance_* (dt / endTime_);
  }

  template< class GridImp, class FunctionSpace >
  double
  AdaptationHandler< GridImp, FunctionSpace >::
  getInitialTolerance () const
  {
    assert( singleThreadMode() );
    const double globalNumberElements = globalNumberOfElements();
    const double dt = deltaT();
    const double localInTimeTol = initialTheta_ * globalTolerance_ * globalTolerance_* (dt / endTime_);
    const double localTol = localInTimeTol / globalNumberElements;

    //std::cout << "Return intial tol " << localTol << std::endl;
    if( verbose() )
    {
      // this requires global communication
      const double globalErr = getMaxEstimator();

      std::cout << std::endl;
      std::cout << "   LocalEst_max = " <<  globalErr
                << "   Tol_local = " << localTol
                << "   Tol = " << localInTimeTol
                << "   GlobalTol = " << globalTolerance_
                << "   Num El: " <<  globalNumberElements << "\n";
    }

    return localTol;
  }

  template< class GridImp, class FunctionSpace >
  double
  AdaptationHandler< GridImp, FunctionSpace >::
  getLocalTolerance () const
  {
    assert( singleThreadMode() );
    const double localInTimeTol = getLocalInTimeTolerance();
    const double globalNumberElements = globalNumberOfElements();

    // apply equi distribution strategy
    double localTol = localInTimeTol / globalNumberElements;

    if( verbose() )
    {
      // this requires global communication
      const double globalErr = getMaxEstimator();

      std::cout << std::endl;
      std::cout << "   LocalEst_max = " <<  globalErr
                << "   Tol_local = " << localTol
                << "   Tol = " << localInTimeTol
                << "   GlobalTol = " << globalTolerance_
                << "   Num El: " <<  globalNumberElements << "\n";
    }

    return localTol;
  }

  // --markEntities
  template< class GridImp, class FunctionSpace >
  void
  AdaptationHandler< GridImp, FunctionSpace >::
  markEntities ( const bool initialAdapt )
  {
    assert( singleThreadMode() );
    // get local refine tolerance
    //const double refineTol = ( initialAdapt ) ? getInitialTolerance() : getLocalTolerance();
    const double refineTol = getMaxEstimator() * globalTolerance_ ;
    // get local coarsen tolerance
    const double coarsenTol = refineTol * coarsenTheta_;

    for( const auto& entity : elements( gridPart_ ) )
    {
      // get local error indicator
      const double localIndicator = getLocalIndicator( entity );
      // get entity level
      const double volume = entity.geometry().volume();

      //std::cout << localIndicator << "   " << refineTol << std::endl;
      //std::cout << finestVolume() << "  " << coarsestVolume() << std::endl;

      // if indicator larger than localTol mark for refinement
      if( (localIndicator > refineTol) && (volume > finestVolume()) )
        // mark for refinement
        grid_.mark( REFINE, entity );
      else if( (localIndicator < coarsenTol) && (volume < coarsestVolume()) )
        // mark for coarsening
        grid_.mark( COARSEN, entity );
      else
        // for for nothing
        grid_.mark( NONE, entity );
    }
  }


  //- --adapt
  template< class GridImp, class FunctionSpace >
  template< class AdaptationManagerType >
  void AdaptationHandler< GridImp, FunctionSpace >::
  adapt ( AdaptationManagerType &am, const bool initialAdapt )
  {
    assert( singleThreadMode() );
    // if adaptation is enabled
    if( am.adaptive() )
    {
      // mark all entities depending on error
      markEntities( initialAdapt );

      // do adaptation
      am.adapt();

      // clear indicator and calc number of elements
      resetStatus();
    }
  }

  //! clear indicator and caculate new number of elements
  template< class GridImp, class FunctionSpace >
  void
  AdaptationHandler< GridImp, FunctionSpace >::
  resetStatus ( const double dt )
  {
    assert( singleThreadMode() );

    // clear error indicator
    clearIndicator();

    // re-calculate number of leaf elements
    // and number of level elements
    countElements();

    // output new number of elements
    if( verbose() )
      std::cout << "   Adaptation: Num El = " << globalNumElements_ << "\n";

    // set new dt if provided (default = 0)
    setDeltaT( dt );
  }

  //! count number of overall leaf entities
  template< class GridImp, class FunctionSpace >
  void
  AdaptationHandler< GridImp, FunctionSpace >::
  countElements ()
  {
    assert( singleThreadMode() );

    // count elements
    size_t count = 0;

    // count interior elements, sum will yield global number of elements after
    // communication
    for( const auto& entity : elements( gridPart_, Dune::Partitions::interior ) )
    {
      DUNE_UNUSED_PARAMETER( entity );
      ++count;
    }

    // number of elements that I have
    localNumElements_ = count;

#ifdef USE_ALUGRID_MPACCESS
    typename ALU3DSPACE MpAccessLocal::minmaxsum_t minmaxsum = mpAccess_.minmaxsum( double(count) );
    globalNumElements_ = (UInt64Type) minmaxsum.sum;
    minNumElements_ = (size_t) minmaxsum.min;
    maxNumElements_ = (size_t) minmaxsum.max;
#else
    double minMax[ 2 ] = { double(count), count > 0 ? 1.0/double(count) : 0.0 } ;
    grid_.comm().max( &minMax[ 0 ], 2 );
    maxNumElements_ = (size_t) minMax[ 0 ];
    minNumElements_ = (size_t) (minMax[ 1 ] > 0) ? (1.0/minMax[ 1 ]) : 0;
    globalNumElements_ = count ;
    globalNumElements_ = grid_.comm().sum( globalNumElements_ );
#endif
  }

} // end namespace
} // end namespace Dune
#endif
