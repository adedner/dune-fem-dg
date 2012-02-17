/***********************************************************************************************
 
   Sourcefile:  adaptation.cc

   Titel:       grid and function adaptation due to error indicator

   Decription:  classes: Adaptation


***********************************************************************************************/
#ifndef DUNE_ADAPTATIONOBJECT_HH
#define DUNE_ADAPTATIONOBJECT_HH

// include restricion, prolongation and adaptation operator classes for discrete functions
#include <dune/grid/utility/persistentcontainer.hh>
 
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/solver/timeprovider.hh>

namespace Dune 
{

struct AdaptationParameters
: public LocalParameter< AdaptationParameters, AdaptationParameters >
{
  int markStrategy_ ; 

  int getStrategy() const 
  {
    const std::string names[] = { "shockind", "apost" , "grad" };
    // default value is gradient 
    return Parameter :: getEnum("fem.adaptation.markingStrategy", names, 2 );
  }

  AdaptationParameters() 
    : markStrategy_( getStrategy() )
  {}

  //! simulation end time 
  virtual double endTime() const 
  {
    return Parameter :: getValue< double >("femhowto.endTime" );
  }

  //! retujrn refinement tolerance 
  virtual double refinementTolerance() const 
  {
    return Parameter :: getValue<double> ("fem.adaptation.refineTolerance");
  }

  //! return percentage of refinement tolerance used for coarsening tolerance
  virtual double coarsenPercentage() const 
  { 
    return Parameter :: getValue<double> ("fem.adaptation.coarsenPercent", 0.1 );
  }

  //! return product of refinementTolerance and coarsenPercentage 
  virtual double coarsenTolerance () const 
  {
    return refinementTolerance() * coarsenPercentage();
  }

  //! return maximal level achieved by refinement  
  virtual int finestLevel( const int refineStepsForHalf ) const 
  { 
    return refineStepsForHalf * 
      Parameter :: getValue<int>("fem.adaptation.finestLevel" );
  }

  //! return minimal level achieved by refinement  
  virtual int coarsestLevel( const int refineStepsForHalf ) const 
  { 
    return refineStepsForHalf * 
      Parameter :: getValue<int>("fem.adaptation.coarsestLevel", 0 );
  }

  //! return depth for refining neighbors of a cell marked for refinement
  virtual int neighborRefLevel() const 
  { 
    return Parameter :: getValue<int>("fem.adaptation.grad.neighborRefLevel", 1 );
  }

  //! return true if marking strategy is based on shock indicator 
  virtual bool shockIndicator() const 
  {
    return markStrategy_ == 0;
  }

  //! return true if marking strategy is based on shock indicator 
  virtual bool gradientBasedIndicator() const 
  {
    return markStrategy_ == 2;
  }

  virtual bool aposterioriIndicator() const 
  {
    return markStrategy_ == 1;
  }
};

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
      ProblemFunctionSpace, 1 > :: Type  FunctionSpaceType;

  // discrete function type of adaptive functions
  typedef typename FunctionSpaceType :: RangeType           RangeType;
  typedef typename GridType :: template Codim<0> :: Entity        GridEntityType; 
  typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType; 

  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType ;

  // local indicator 
  struct LocalIndicator 
  {
    LocalIndicator() : value_( -1.0 ) {} 
    double value_;
    LocalIndicator& operator += ( const double& value ) 
    {
      value_ += value; 
      return *this;
    }

    LocalIndicator& operator = ( const double& value ) 
    {
      value_ = value; 
      return *this;
    }
    double value() const { return value_; }
  };

  // type of indicator stored by for each entity 
  typedef LocalIndicator LocalIndicatorType;

  typedef PersistentContainer< GridType, LocalIndicatorType > IndicatorType ;

  // time provider 
  typedef TimeProviderBase TimeProviderType ;

  // interface for adaptation operator 
  typedef AdaptationManagerInterface   AdaptInterfaceType;

private:  
  // no copying 
  AdaptationHandler (const AdaptationHandler&); 
public:
  //! constructor
  AdaptationHandler (GridType &grid, 
                     TimeProviderType &timeProvider,
                     const AdaptationParameters& param = AdaptationParameters() ) ;

  //const TimeProviderType &timeProvider() const;
  TimeProviderType* timeProvider () { return &timeProvider_; } 
public: 
  //! clear indicator 
  void clearIndicator();

  //! initialize enIndicator with en 
  template <class Entity> 
  void setEntity(const Entity& en);
  
  //! initialize nbIndicator with en 
  template <class Entity>
  void setNeighbor(const Entity& en);
  
  //! add value to local indicator, use setEntity before 
  void addToEntityIndicator(const FullRangeType& error, const double h );

  //! add value to local indicator, use setNeighbor before 
  void addToNeighborIndicator( const FullRangeType& error, const double h );

  //! add to local indicator for given entity 
  void addToLocalIndicator( LocalIndicatorType& indicator, const FullRangeType& error, const double h );

  //! add to local indicator for given entity 
  void addToLocalIndicator(const GridEntityType &en, const FullRangeType& error, const double h );

  //! det local indicator for given entity 
  void setLocalIndicator(const GridEntityType &en, const FullRangeType& error);

  //! return local indicator for given entity 
  double getLocalIndicator(const GridEntityType &en) const;

  //! calculate sum of local errors 
  double getSumEstimator() const ;

  //! calculate max of local errors 
  double getMaxEstimator() const ;

  //! overall number of leaf elements 
  int globalNumberOfElements () const ;

  //! number of local leaf elements 
  int localNumberOfElements () const ;

  //! get local in time tolerance  
  double getLocalInTimeTolerance () const ;

  //! get local tolerance 
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

  //! indicator function
  IndicatorType indicator_;

  //! local function of indicator_ for entity 
  LocalIndicatorType* enIndicator_;
  //! local function of indicator_ for neighbor 
  LocalIndicatorType* nbIndicator_;

  //! timestep size in time discretization parameters und endTime 
  TimeProviderType & timeProvider_;

  //! parameters for adaptation
  mutable double globalTolerance_;
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
