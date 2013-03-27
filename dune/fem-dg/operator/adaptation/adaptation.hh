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
#include <dune/fem/space/common/adaptmanager.hh>

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/solver/timeprovider.hh>


namespace Dune
{

  struct AdaptationParameters
    : public Fem::LocalParameter< AdaptationParameters, AdaptationParameters >
  {
    int markStrategy_;

    int getStrategy () const
    {
      const std::string names[] = { "shockind", "apost", "grad" };
      // default value is gradient
      return Fem::Parameter::getEnum( "fem.adaptation.markingStrategy", names, 2 );
    }

    AdaptationParameters ()
      : markStrategy_( getStrategy() )
    {}

    //! simulation end time
    virtual double endTime () const
    {
      return Fem::Parameter::getValue< double >( "femhowto.endTime" );
    }

    //! retujrn refinement tolerance
    virtual double refinementTolerance () const
    {
      return Fem::Parameter::getValue< double >( "fem.adaptation.refineTolerance" );
    }

    //! return percentage of refinement tolerance used for coarsening tolerance
    virtual double coarsenPercentage () const
    {
      return Fem::Parameter::getValue< double >( "fem.adaptation.coarsenPercent", 0.1 );
    }

    //! return product of refinementTolerance and coarsenPercentage
    virtual double coarsenTolerance () const
    {
      return refinementTolerance() * coarsenPercentage();
    }

    //! return maximal level achieved by refinement
    virtual int finestLevel ( const int refineStepsForHalf ) const
    {
      return refineStepsForHalf *
             Fem::Parameter::getValue< int >( "fem.adaptation.finestLevel" );
    }

    //! return minimal level achieved by refinement
    virtual int coarsestLevel ( const int refineStepsForHalf ) const
    {
      return refineStepsForHalf *
             Fem::Parameter::getValue< int >( "fem.adaptation.coarsestLevel", 0 );
    }

    //! return depth for refining neighbors of a cell marked for refinement
    virtual int neighborRefLevel () const
    {
      return Fem::Parameter::getValue< int >( "fem.adaptation.grad.neighborRefLevel", 1 );
    }

    //! return true if marking strategy is based on shock indicator
    virtual bool shockIndicator () const
    {
      return markStrategy_ == 0;
    }

    //! return true if marking strategy is based on shock indicator
    virtual bool gradientBasedIndicator () const
    {
      return markStrategy_ == 2;
    }

    //! return true if aposteriori indicator is enabled
    virtual bool aposterioriIndicator () const
    {
      return markStrategy_ == 1;
    }

    //! return true if verbosity mode is enabled
    virtual bool verbose () const { return Fem::Parameter::getValue< bool >( "fem.adaptation.verbose", false ); }
  };

// class for the organization of the adaptation prozess
  template< class GridImp, class ProblemFunctionSpace >
  class AdaptationHandler
  {
    typedef AdaptationHandler< GridImp, ProblemFunctionSpace > ThisType;

  public:
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
    typedef typename GridType::template Codim< 0 >::EntityPointer EntityPointerType;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

    // local indicator
    struct LocalIndicatorData
    {
      LocalIndicatorData () : value_( 0.0 ) {}
      double value_;
      LocalIndicatorData &operator+= ( const double &value )
      {
        value_ += value;
        return *this;
      }

      LocalIndicatorData &operator= ( const double &value )
      {
        value_ = value;
        return *this;
      }
      double value () const { return value_; }
    };


    // type of indicator stored by for each entity
    typedef LocalIndicatorData LocalIndicatorDataType;

    typedef PersistentContainer< GridType, LocalIndicatorDataType > IndicatorType;

    class LocalIndicator
    {
      const ThisType *adaptation_;
      LocalIndicatorDataType *localIndicator_;

    public:
      LocalIndicator ()
        : adaptation_( 0 ),
        localIndicator_( 0 ) {}

      LocalIndicator ( const ThisType *adaptation, LocalIndicatorDataType *indicator )
        : adaptation_( adaptation ),
        localIndicator_( indicator ) {}

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
        localIndicator_ = 0;
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
                        TimeProviderType &timeProvider,
                        const AdaptationParameters &param = AdaptationParameters() );

    // copy constructor
    AdaptationHandler ( const AdaptationHandler & );

    //const TimeProviderType &timeProvider() const;
    TimeProviderType *timeProvider () { return &timeProvider_; }

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
    void resetStatus ();

    //! count number of overall leaf entities
    UInt64Type countElements() const;

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
    int thread () const { return Dune::Fem::ThreadManager::thread(); }
    bool singleThreadMode () const { return Dune::Fem::ThreadManager::singleThreadMode(); }

    // return volume of computational domain
    double volumeOfDomain () const;

    //! grid part, has grid and ind set
    GridType &grid_;
    GridPartType gridPart_;

    //! persistent container holding local indicators
    IndicatorType indicator_;

    //! timestep size in time discretization parameters und endTime
    TimeProviderType &timeProvider_;

    //! parameters for adaptation
    mutable double globalTolerance_;
    const double coarsenTheta_;
    const double initialTheta_;

    const int finestLevel_;
    const int coarsestLevel_;

    UInt64Type globalNumElements_;

    mutable int localNumElements_;

    double endTime_;

    mutable std::vector< UInt64Type > maxLevelCounter_;

    const bool verbose_;
  };

} // end namespace Dune

#include "adaptation.cc"
#endif
