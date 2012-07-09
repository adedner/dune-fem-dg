#ifndef DUNE_ADAPTATIONOBJECT_CC
#define DUNE_ADAPTATIONOBJECT_CC

#include "adaptation.hh"
#include <dune/fem/misc/gridwidth.hh>

namespace Dune 
{

//! constructor
template <class GridImp, class FunctionSpace>
AdaptationHandler<GridImp, FunctionSpace> ::  
AdaptationHandler (GridType &grid, 
                   TimeProviderType &timeProvider,
                   const AdaptationParameters& param ) 
  : grid_(grid)
  , gridPart_(grid_)
  , indicator_( grid_, 0 ) // grid , codimension 
  , enIndicator_( 0 )
  , nbIndicator_( 0 )
  , timeProvider_(timeProvider)
  , globalTolerance_ ( param.refinementTolerance() )
  , coarsenTheta_( param.coarsenPercentage() )
  , finestLevel_( param.finestLevel( DGFGridInfo<GridType>::refineStepsForHalf() ) )
  , coarsestLevel_( param.coarsestLevel( DGFGridInfo<GridType>::refineStepsForHalf() ) )
  , globalNumElements_ (0)
  , maxNumberOfElementsAllowed_( getMaxNumberOfElements() )  
  , localNumElements_(0)
  , endTime_( param.endTime() )  
  , maxLevelCounter_()
  , verbose_( Parameter :: verbose() && param.verbose() )
{
  const bool verboseOutput = Parameter :: verbose() ;

  // set default values
  initialTheta_ = 0.995;
  
  maxLevFlag_ = 1;

  resetStatus();
  if( verboseOutput )
  {
    std::cout << "AdaptationHandler created! \n";
  }

  /*
  // calculate global min of grid width to scale tolerance  
  double gridWidth = GridWidth::calcGridWidth( gridPart_ );
  
  // get global minimum of macro grid width 
  double macroGridWidth = grid_.comm().min( gridWidth );

  // globalTolerance_ *= macroGridWidth;
  */

  // scale tolerance with domain volume 
  globalTolerance_ *= volumeOfDomain();
}

//! clear indicator 
template <class GridImp, class FunctionSpace>
double 
AdaptationHandler<GridImp, FunctionSpace> ::  
volumeOfDomain() const
{
  double volume = 0;
  // type of iterator, i.e. leaf iterator 
  typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
  
  const IteratorType endit = gridPart_.template end< 0 > ();
  for (IteratorType it = gridPart_.template begin< 0 > (); 
       it != endit; ++it)
  {
    // entity 
    const GridEntityType& entity = *it;

    // sum up the volume 
    volume += entity.geometry().volume();
  }

  // return volume of computational domain 
  return grid_.comm().sum( volume );
}

//! clear indicator 
template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
clearIndicator()
{
  // set all entries to zero
  indicator_.clear();
}

//! initialize localIndicator with en 
template <class GridImp, class FunctionSpace>
template <class Entity>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
setEntity( const Entity& entity )
{
  // convert the given entity to an entity of the grid 
  // for wrapped entities the cast to the host entity is necessary 
  const GridEntityType& gridEntity = Fem :: gridEntity( entity );
  enIndicator_ = & indicator_[ gridEntity ];
}
  
//! initialize localIndicator with en 
template <class GridImp, class FunctionSpace>
template <class Entity>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
setNeighbor( const Entity& neighbor )
{
  // convert the given entity to an entity of the grid 
  // for wrapped entities the cast to the host entity is necessary 
  const GridEntityType& gridNeighbor = Fem :: gridEntity( neighbor );
  nbIndicator_ = & indicator_[ gridNeighbor ];
}
  
//! add value to local indicator, use setEntity before 
template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
addToLocalIndicator( LocalIndicatorType& indicator, 
                     const FullRangeType& error, const double h )
{
  const double dt = timeProvider_.deltaT();
  const double factor = ( h + dt  ) * dt ;
  indicator += (factor * error.two_norm());  
}

//! add value to local indicator, use setEntity before 
template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
addToEntityIndicator(const FullRangeType& error, const double h )
{
  assert( enIndicator_ );
  addToLocalIndicator( *enIndicator_, error, h );
}

//! add value to local indicator, use setEntity before 
template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
addToNeighborIndicator(const FullRangeType& error, const double h )
{
  assert( nbIndicator_ );
  addToLocalIndicator( *nbIndicator_, error, h );
}

template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
addToLocalIndicator(const GridEntityType &en, const FullRangeType& error, const double h )
{
  addToLocalIndicator( indicator_[ en ], error, h );
  return;
}

template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
setLocalIndicator(const GridEntityType &en, const FullRangeType& error)
{
  indicator_[ en ] = error[ 0 ] ;
  return;
}

template <class GridImp, class FunctionSpace>
double  
AdaptationHandler<GridImp, FunctionSpace> ::  
getLocalIndicator(const GridEntityType &en) const
{
  return indicator_[ en ].value();
}

  //! calculate sum of local errors 
template <class GridImp, class FunctionSpace>
double  
AdaptationHandler<GridImp, FunctionSpace> ::  
getSumEstimator() const 
{    
  double sum = 0.0;

  typedef typename IndicatorType :: ConstIterator   IteratorType ;
  const IteratorType endit = indicator_.end();
  for(IteratorType it = indicator_.begin(); it != endit; ++it)
  {
    sum += (*it).value();
  }

  // global sum of estimator 
  sum = grid_.comm().sum( sum );
  return sum;
}

  //! calculate max of local errors 
template <class GridImp, class FunctionSpace>
double  
AdaptationHandler<GridImp, FunctionSpace> ::  
getMaxEstimator() const 
{    
  double max = 0.0;

  {
    typedef typename IndicatorType :: ConstIterator   IteratorType ;
    const IteratorType endit = indicator_.end();
    IteratorType it = indicator_.begin(); 

    // initialzie with first entry 
    if ( it != endit ) max = (*it).value();
    
    for( ; it != endit; ++it)
    {
      max = std::max((*it).value(), max);
    }
  }

  // global sum of estimator 
  max = grid_.comm().max(max);
  return max;
}

template <class GridImp, class FunctionSpace>
int  
AdaptationHandler<GridImp, FunctionSpace> ::  
localNumberOfElements () const 
{
  assert( localNumElements_ > 0 );
  return localNumElements_;
}

template <class GridImp, class FunctionSpace>
int   
AdaptationHandler<GridImp, FunctionSpace> ::  
globalNumberOfElements () const 
{
  assert( globalNumElements_ > 0 );
  return globalNumElements_;
}

template <class GridImp, class FunctionSpace>
double  
AdaptationHandler<GridImp, FunctionSpace> ::  
getLocalInTimeTolerance () const 
{
  double dt = timeProvider_.deltaT();
  return (1. - initialTheta_) * globalTolerance_ * globalTolerance_* (dt / endTime_);
  //return globalTolerance_;
}

template <class GridImp, class FunctionSpace>
double  
AdaptationHandler<GridImp, FunctionSpace> ::  
getLocalTolerance () const
{
  const double localInTimeTol = getLocalInTimeTolerance ();
  const double globalNumberElements = globalNumberOfElements();
  const double elementsOnMaxLevel = maxLevelCounter_[finestLevel_]; 
  
  /*
  if( elementsOnMaxLevel < 
        0.1 * globalNumberElements ) 
  {
    globalTolerance_ *= 0.01;
  }
  else if( elementsOnMaxLevel > 
           0.5 * globalNumberElements )
  {
    globalTolerance_ *= 1.05;
  }
  */

  /*
  double factor = 1.0 ;

  if( globalNumberElements > 1.25 * maxNumberOfElementsAllowed_ ) 
  {
    // (tanh( ( 2.0*pi * x - pi))+1)*0.5
    double x = globalNumberElements / maxNumberOfElementsAllowed_ ;
    //factor = 1.0 + 0.5 * (std::sin(M_PI*x + 0.5*M_PI) + 1.0);
    factor = 1.0 + 0.5 * (std::tanh(M_PI*x - M_PI) + 1.0);
    std::cout << factor << " factor|percent " << x << std::endl;
  }
  else if ( globalNumberElements < 0.25 * maxNumberOfElementsAllowed_ )
  {
    double x = globalNumberElements / maxNumberOfElementsAllowed_ ;
    //factor = 0.5 * (std::sin(M_PI*x + 0.5*M_PI) + 1.0);
    factor = 0.5 * (std::tanh(M_PI*x - M_PI) + 1.0);
    std::cout << factor << " factor|percent " << x << std::endl;
  }

  // adjust global tolerance such that maximal number of elements 
  // is kept more or less constant 
  globalTolerance_ *= factor ;
  */

  double localTol = localInTimeTol / globalNumberElements ;

  if( verbose() )
  {
    // this requires global communication 
    const double globalErr = getMaxEstimator(); 

    std::cout << "Level counters: ";
    for(size_t i=0; i<maxLevelCounter_.size(); ++i ) 
    {
      double percent = maxLevelCounter_[ i ]/ globalNumberElements ;
      std::cout << "Level[ " << i << " ] = " << percent << "  ";
    }
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
template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
markEntities ()
{
  // get local refine tolerance 
  const double refineTol = getLocalTolerance();
  // get local coarsen tolerance 
  const double coarsenTol = refineTol * coarsenTheta_;

  // type of iterator, i.e. leaf iterator 
  typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
  
  const IteratorType endit = gridPart_.template end< 0 > ();
  for (IteratorType it = gridPart_.template begin< 0 > (); 
       it != endit; ++it)
  {
    // entity 
    const GridEntityType& entity = *it;

    // get local error indicator 
    const double localIndicator = getLocalIndicator(entity);
    // get entity level 
    const int level = entity.level() ;
    
    // if indicator larger than localTol mark for refinement 
    if( (localIndicator > refineTol) && (level < finestLevel_) )
    {
      // mark for refinement 
      grid_.mark(REFINE, entity);
    }
    else if ( (localIndicator < coarsenTol) && (level > coarsestLevel_) )
    {
      // mark for coarsening 
      grid_.mark(COARSEN, entity);
    }
    else
    {
      // for for nothing 
      grid_.mark(NONE, entity);
    }
  }
  return;
}


//- --adapt 
template <class GridImp, class FunctionSpace>
template <class AdaptationManagerType>
void AdaptationHandler<GridImp, FunctionSpace> ::  
adapt( AdaptationManagerType& am )
{
  // if adaptation is enabled 
  if( am.adaptive() )
  {
    // mark all entities depending on error
    markEntities(); 

    // do adaptation 
    am.adapt();

    // clear indicator and calc number of elements 
    resetStatus();
  }
}

//! clear indicator and caculate new number of elements 
template <class GridImp, class FunctionSpace>
void 
AdaptationHandler<GridImp, FunctionSpace> ::  
resetStatus() 
{
  // clear error indicator 
  clearIndicator();

  // re-calculate number of leaf elements 
  // and number of level elements 
  globalNumElements_ = countElements();

  // output new number of elements 
  if( verbose() )
  {
    std::cout << "   Adaptation: Num El = " << globalNumElements_ << "\n";
  }
}

  //! count number of overall leaf entities 
template <class GridImp, class FunctionSpace>
int 
AdaptationHandler<GridImp, FunctionSpace> ::  
countElements() const 
{
  if(  maxLevelCounter_.size() > 0)
    maxLevelCounter_.clear();

  maxLevelCounter_.resize(finestLevel_ + 1);
  for(size_t i=0; i<maxLevelCounter_.size(); ++i)
    maxLevelCounter_[i] = 0;
  
  // count elements 
  size_t count = 0;
  // type of iterator, i.e. leaf iterator 
  typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
  
  const IteratorType endit = gridPart_.template end< 0 > ();
  for (IteratorType it = gridPart_.template begin< 0 > (); 
       it != endit; ++it)
  {
    ++count;
    const int level = it->level();
    assert( level < int(maxLevelCounter_.size()) );
    ++maxLevelCounter_[ level ];
  }

  // number of elements that I have
  localNumElements_ = count;

  {
    size_t commSize = maxLevelCounter_.size();
    int* commBuff = new int [commSize+1]; 
    
    for(size_t i=0; i<commSize; ++i)
    {
      commBuff[i] = maxLevelCounter_[i];
    }
    commBuff[commSize] = count;
      
    grid_.comm().sum(commBuff,commSize+1);
    
    for(size_t i=0; i<commSize; ++i)
    {
      maxLevelCounter_[i] = commBuff[i]; 
    }

    count = commBuff[commSize];
    delete [] commBuff;
  }
  // return element count 
  return count;
}
  
} // end namespace Dune 
#endif
