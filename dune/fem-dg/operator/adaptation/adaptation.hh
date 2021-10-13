/***********************************************************************************************

   Sourcefile:  adaptation.cc

   Titel:       grid and function adaptation due to error indicator

   Decription:  classes: Adaptation


***********************************************************************************************/
#ifndef DUNE_ADAPTATIONOBJECT_HH
#define DUNE_ADAPTATIONOBJECT_HH

// include restricion, prolongation and adaptation operator classes for discrete functions
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/adaptationmanager.hh>

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/solver/timeprovider.hh>

#include <dune/fem-dg/operator/adaptation/utility.hh>

#if HAVE_DUNE_ALUGRID && HAVE_MPI
#define USE_ALUGRID_MPACCESS
#include <dune/alugrid/3d/alugrid.hh>
#endif

namespace Dune
{
namespace Fem
{

  // class for the organization of the adaptation prozess
  template< class GridImp, class ProblemFunctionSpace >
  class AdaptationHandler
    : public ComputeMinMaxVolume
  {
    typedef AdaptationHandler< GridImp, ProblemFunctionSpace > ThisType;

  public:
    using ComputeMinMaxVolume :: coarsestVolume ;
    using ComputeMinMaxVolume :: finestVolume ;
    using ComputeMinMaxVolume :: computeGlobalMinMax ;

    enum { COARSEN = -1, NONE = 0, REFINE = 1 };

    typedef GridImp GridType;
    typedef Fem::DGAdaptiveLeafGridPart< GridType > GridPartType;

    // useful enums and typedefs
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename ProblemFunctionSpace::RangeType FullRangeType;

    // initialize functionspace, etc., for the indicator function
    typedef typename Fem::ToNewDimRangeFunctionSpace<
      ProblemFunctionSpace, 1 >::Type FunctionSpaceType;

    // discrete function type of adaptive functions
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename GridType::template Codim< 0 >::Entity GridEntityType;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

    // type of indicator stored by for each entity
    typedef double LocalIndicatorDataType;

    typedef PersistentContainer< GridType, LocalIndicatorDataType > IndicatorType;

    class LocalIndicator
    {
      const ThisType *adaptation_;
      LocalIndicatorDataType *localIndicator_;

    public:
      LocalIndicator ()
        : adaptation_( nullptr ),
          localIndicator_( nullptr )
      {}

      LocalIndicator ( const ThisType *adaptation, LocalIndicatorDataType *indicator )
        : adaptation_( adaptation ),
          localIndicator_( indicator )
      {}

      LocalIndicator ( const LocalIndicator &other )
        : adaptation_( other.adaptation_ ),
          localIndicator_( other.localIndicator_ )
      {}

      LocalIndicator &operator= ( const LocalIndicator &other )
      {
        // make sure we are considering the same adaptation object
        adaptation_     = other.adaptation_;
        localIndicator_ = other.localIndicator_;
        return *this;
      }

      //! reset local indicator
      void reset ()
      {
        localIndicator_ = nullptr;
      }

      //! add to local indicator
      void add ( const FullRangeType &error, const double h )
      {
        assert( localIndicator_ );
        assert( adaptation_ );
        adaptation_->addToLocalIndicator( *localIndicator_, error, h );
      }

      //! add to local indicator if localIndicator is valid
      void addChecked ( const FullRangeType &error, const double h )
      {
        if( localIndicator_ )
          add( error, h );
      }
    };

    // type of local indicator
    typedef LocalIndicator LocalIndicatorType;

    // time provider
    typedef Fem::TimeProviderBase TimeProviderType;

    // interface for adaptation operator
    typedef Fem::AdaptationManagerInterface AdaptInterfaceType;

    // type of 64 bit unsigned integer
    typedef uint64_t UInt64Type;

  public:
    //! constructor
    AdaptationHandler ( GridType &grid,
                        const AdaptationParameters &param = AdaptationParameters() );

    // copy constructor
    AdaptationHandler ( const AdaptationHandler & );

    TimeProviderType *timeProvider () {
      std::abort();
      return nullptr;
    }

    double deltaT () const { return deltaT_; }
    void setDeltaT( const double dt ) { deltaT_ = dt; }

    //! clear indicator
    void clearIndicator();

    //! return local indicator object
    template< class Entity >
    LocalIndicatorType localIndicator( const Entity &entity );

    //! add another AdaptationHandlers indicator
    ThisType &operator+= ( const ThisType &other );

    //! reset nbIndicator to null pointer
    void resetNeighbor();

    //! add value to local indicator, use setEntity before
    void addToEntityIndicator( const FullRangeType &error, const double h );

    //! add value to local indicator, use setNeighbor before
    void addToNeighborIndicator( const FullRangeType &error, const double h );

    //! add to local indicator for given entity
    void addToLocalIndicator( LocalIndicatorDataType &indicator, const FullRangeType &error, const double h ) const;

    //! add to local indicator for given entity
    void addToLocalIndicator( const GridEntityType &en, const FullRangeType &error, const double h );

    //! det local indicator for given entity
    void setLocalIndicator( const GridEntityType &en, const FullRangeType &error );

    //! return local indicator for given entity
    double getLocalIndicator( const GridEntityType &en ) const;

    //! calculate sum of local errors
    double getSumEstimator() const;

    //! calculate max of local errors
    double getMaxEstimator() const;

    //! overall number of leaf elements
    UInt64Type globalNumberOfElements () const;

    //! min number of leaf element that one process has
    size_t minNumberOfElements () const;

    //! max number of leaf element that one process has
    size_t maxNumberOfElements () const;

    //! number of local leaf elements
    int localNumberOfElements () const;

    //! get local in time tolerance
    double getLocalInTimeTolerance () const;

    //! get initial tolerance
    double getInitialTolerance () const;

    //! get local tolerance
    double getLocalTolerance () const;

    // --markEntities
    void markEntities ( const bool initialAdapt );

    //- --adapt
    template< class AdaptationManagerType >
    void adapt( AdaptationManagerType &, const bool initialAdapt = false );

    //! reset status of indicator and count elements
    void resetStatus ( const double dt = 0.0 );

    //! module interface for intialize
    void initialize ()
    {
      clearIndicator();
    }

    //! module interface for one time step
    void solveTimeStep ()
    {
      adapt();
    }

    //! return true if verbosity mode is enabled
    bool verbose () const { return verbose_; }

  protected:
    //! count number of overall leaf entities
    void countElements();

    int thread () const { return Dune::Fem::MPIManager::thread(); }
    bool singleThreadMode () const { return Dune::Fem::MPIManager::singleThreadMode(); }

    // return volume of computational domain
    double volumeOfDomain () const;

    //! grid part, has grid and ind set
    GridType &grid_;
    GridPartType gridPart_;

#ifdef USE_ALUGRID_MPACCESS
    ALU3DSPACE MpAccessMPI mpAccess_;
#endif

    //! persistent container holding local indicators
    IndicatorType indicator_;

    //! parameters for adaptation
    mutable double globalTolerance_;
    const double coarsenTheta_;
    const double initialTheta_;

    UInt64Type globalNumElements_;
    size_t minNumElements_;
    size_t maxNumElements_;

    mutable int localNumElements_;

    double endTime_;
    double deltaT_;

    const bool verbose_;
  };

} // end namespace
} // end namespace Dune

#include "adaptation.cc"
#undef USE_ALUGRID_MPACCESS
#endif
