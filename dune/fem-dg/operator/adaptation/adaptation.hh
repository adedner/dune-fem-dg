/***********************************************************************************************
 
   Sourcefile:  adaptation.cc

   Titel:       grid and function adaptation due to error indicator

   Decription:  classes: Adaptation


***********************************************************************************************/
#ifndef DUNE_ADAPTATIONOBJECT_HH
#define DUNE_ADAPTATIONOBJECT_HH

// include restricion, prolongation and adaptation operator classes for discrete functions
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/space/fvspace.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/solver/timeprovider.hh>

namespace Dune 
{

// class for the organization of the adaptation prozess
template <class GridImp, class ProblemFunctionSpace >
class AdaptationHandler
{
public:
  enum { COARSEN = -1, NONE = 0, REFINE = 1 };

  typedef GridImp GridType ;
  typedef DGAdaptiveLeafGridPart< GridType > GridPartType ;

  // useful enums and typedefs
  enum { dim = GridType :: dimension };
  enum { dimworld = GridType::dimensionworld };
  
  typedef typename ProblemFunctionSpace :: RangeType FullRangeType ;

  // initialize functionspace, etc., for the indicator function
  typedef typename ToNewDimRangeFunctionSpace< 
      ProblemFunctionSpace, 1 > :: Type IndicatorFuncSpaceType;

  //! type of discrete function space for indicator  
  typedef FiniteVolumeSpace<IndicatorFuncSpaceType, 
                   GridPartType, 0, CachingStorage> DiscreteFunctionSpaceType;

  typedef DiscreteFunctionSpaceType IndicatorDiscreteFunctionSpaceType;

  // discrete function type of adaptive functions
  typedef typename DiscreteFunctionSpaceType::RangeType           RangeType;
  typedef typename GridType :: template Codim<0> :: Entity        EntityType; 
  typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType; 

  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType ;

  typedef AdaptiveDiscreteFunction< IndicatorDiscreteFunctionSpaceType >            IndicatorDiscreteFunctionType; 

  typedef typename IndicatorDiscreteFunctionType::LocalFunctionType IndicatorLocalFuncType;

  //typedef GridTimeProviderType< GridType > TimeProviderType ;
  typedef TimeProviderBase TimeProviderType ;

  // interface for adaptation operator 
  typedef AdaptationManagerInterface   AdaptInterfaceType;

private:  
  // no copying 
  AdaptationHandler (const AdaptationHandler&); 
public:
  typedef IndicatorLocalFuncType LocalIndicatorType;
  //! constructor
  AdaptationHandler (GridType &grid, 
                     TimeProviderType &timeProvider ) ;

  //const TimeProviderType &timeProvider() const;
  TimeProviderType* timeProvider () { return &timeProvider_; } 
public: 
  //! clear indicator 
  void clearIndicator();

  //! initialize enIndicator with en 
  void setEntity(const EntityType & en);
  
  //! initialize nbIndicator with en 
  void setNeighbor(const EntityType & en);
  
  //! add value to local indicator, use setEntity before 
  void addToLocalIndicator(const FullRangeType& error, const double h );

  //! add value to local indicator, use setNeighbor before 
  void addToNeighborIndicator( const FullRangeType& error, const double h );

  //! add value to local indicator, use setNeighbor before 
  void addToLocalIndicator( LocalIndicatorType& indicator,
        const FullRangeType& error, const double h );

  void addToLocalIndicator(const EntityType &en, const FullRangeType& error );

  void setLocalIndicator(const EntityType &en, const FullRangeType& error);

  double getLocalIndicator(const EntityType &en) const;

  //! calculate sum of local errors 
  double getSumEstimator() const ;

  //! calculate max of local errors 
  double getMaxEstimator() const ;

  // overall number of leaf elements 
  int globalNumberOfElements () const ;

  // number of local leaf elements 
  int localNumberOfElements () const ;

  double getLocalInTimeTolerance () const ;

  double getLocalTolerance () const;

  // --markEntities
  void markEntities ();

  //- --adapt 
  template <class AdaptationManagerType> 
  void adapt( AdaptationManagerType& );

  //! reset status of indicator and count elements 
  void resetStatus ();

  //! count number of overall leaf entities 
  int countElements() const ;

  //! export indicator function
  IndicatorDiscreteFunctionType& indicator () ;

  //! module interface for intialize 
  void initialize() 
  {
    clearIndicator();
  }
  
  //! module interface for one time step  
  void solveTimeStep() 
  {
    adapt();
  }

private:
  //! grid part, has grid and ind set 
  GridType&  grid_;
  GridPartType gridPart_;

  //! function space for discrete solution
  IndicatorDiscreteFunctionSpaceType indicatorSpace_;

  //! indicator function
  IndicatorDiscreteFunctionType indicator_;

  //! local function of indicator_ for entity 
  LocalIndicatorType enIndicator_;
  //! local function of indicator_ for neighbor 
  LocalIndicatorType nbIndicator_;

  //! timestep size in time discretization parameters und endTime 
  TimeProviderType & timeProvider_;

  //! parameters for adaptation
  mutable double globalTolerance_;
  //double localInTimeTolerance_;
  //double localToleranceOrig_;
  //double localTolerance_;
  double coarsenTheta_;
  double initialTheta_;

  double alphaSigSet_;
  double tolSigSet_;
  double tolSigSetInv_;
  double tolMaxLevSet_;
  int    numSigSet_;
  int    numMaxLev_;
  int    maxLevFlag_;
  double maxLevAlpha_;
  double maxLevBeta_;

  int finestLevel_; 
  int coarsestLevel_; 

  int globalNumElements_;
  mutable int localNumElements_;

  double endTime_;

  mutable std::vector<int> maxLevelCounter_;
};

} // end namespace Dune 

#include "adaptation.cc"
#endif
