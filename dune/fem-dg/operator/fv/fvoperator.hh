#ifndef DUNE_FEMDG_FVOPERATOR_HH
#define DUNE_FEMDG_FVOPERATOR_HH

#include <limits>
#include <dune/common/version.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/misc/threads/threaditerator.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem-dg/pass/context.hh>

namespace Dune
{

namespace Fem
{

namespace detail
{
  template <class DiscreteFunction>
  struct AddLocalEvaluate
  {
    typedef typename DiscreteFunction::DiscreteFunctionType::EntityType  EntityType;
    AddLocalEvaluate(DiscreteFunction &w)
    : w_(w) {}

    template <class LocalDofs>
    void operator () (const EntityType& entity, const LocalDofs& wLocal ) const
    {
      w_.addLocalDofs( entity, wLocal );
    }
    DiscreteFunction &w_;
  };

  // FVOperatorImpl
  template< class DiscreteFunction, class Model, class NumFlux >
  class FVOperatorImpl
  {

  public:
    typedef Model     ModelType;
    typedef NumFlux   NumFluxType;

    // first we extract some types
    typedef DiscreteFunction             DiscreteFunctionType;
    typedef DiscreteFunction             DestinationType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

    static_assert( DiscreteFunctionSpaceType::polynomialOrder == 0, "This only works for FiniteVolumeSpace" );

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef typename GridPartType::Grid Grid;
    static const int dim = GridPartType::dimension;
    static const int dimworld = GridPartType::dimensionworld;
    static const int dimRange = Model::dimRange;
    typedef typename Grid::ctype ctype;
    static const bool isCartesian = Fem::GridPartCapabilities::isCartesian< GridPartType >::v;


    // only apply the scheme to interior elements
    static const Dune :: PartitionIteratorType ptype = Dune :: InteriorBorder_Partition ;

    // types of codim zero entity iterator and geometry
    typedef typename GridPartType::template Codim< 0 >:: template Partition< ptype > :: Iterator  Iterator;
    typedef typename Iterator::Entity                         Entity;
    typedef typename Entity::Geometry                         Geometry;

    // type of intersections and corresponding geometries
    typedef typename GridPartType::Intersection   Intersection;

    //typedef ElementQuadrature< GridPartType, 0 >  VolumeQuadratureType;
    typedef CachingQuadrature< GridPartType, 1 >  FaceQuadratureType;

    // types of vectors
    typedef typename ModelType::DomainType          DomainType;
    typedef typename ModelType::FaceDomainType      FaceDomainType;
    typedef typename ModelType::RangeType           RangeType;
    typedef typename ModelType::JacobianRangeType   JacobianRangeType;

    typedef PointContext< Entity >                LocalEvalEntityType;
    typedef PointContext< Entity, Intersection >  LocalEvalIntersectionType;

  public:
    /** \brief constructor
     *
     *  \param[in]  gridPart  gridPart to operate on
     *  \param[in]  model       discretization of the Model
     */
    FVOperatorImpl ( const DiscreteFunctionSpaceType& space, const Model &model, const NumFluxType& numFlux )
    : space_( space ),
      model_( model ),
      numFlux_( numFlux ),
      time_( 0 ),
      dtEst_( 0 ),
      jacLeft_( 0 ),
      jacRight_( 0 ),
      center_( 0 ),
      faceCenters_(),
      localFaceCenter_( 0.5 )
    {
      typedef typename GridPartType::ctype ctype;
      {
        const auto& geomTypes = space_.geomTypes( 0 );
        assert( geomTypes.size() == 1 );
        const auto& refElem = Dune::referenceElement<ctype, GridPartType::dimension>( geomTypes[ 0 ] );
        center_ = refElem.position( 0, 0 );

        const int faces = refElem.size( 1 );
        faceCenters_.resize( faces );
        for( int i=0; i<faces; ++i )
        {
          faceCenters_[ i ] = refElem.position( i, 1 );
          //std::cout << "Face " << i << " center = " << faceCenters_[ i] <<  std::endl;
        }
      }
      {
        const auto& geomTypes = space_.geomTypes( 1 );
        assert( geomTypes.size() == 1 );
        localFaceCenter_ = Dune::referenceElement<ctype, GridPartType::dimension-1>( geomTypes[ 0 ] ).position( 0, 0 );
      }
    }

    template <class GridFunction, class Iterators>
    void evaluate( const GridFunction& u, DestinationType& w, const Iterators& iterators ) const
    {
      AddLocalEvaluate< DestinationType > addLocalEvaluate(w);
      evaluate( u, w, iterators, addLocalEvaluate );
    }

    /** \brief obtain the grid view for this scheme
     *
     *  \returns the grid view
     */
    const GridPartType &gridPart () const
    {
      return space().gridPart();
    }

    const DiscreteFunctionSpaceType& space() const
    {
      return space_;
    }

    void setTime( const double time )
    {
      time_ = time;
    }

    double timeStepEstimate() const
    {
      return dtEst_;
    }

    size_t numberOfElements () const { return numberOfElements_; }

  protected:
    template <class GridFunction, class Iterators, class Functor>
    void evaluate( const GridFunction& u, DestinationType& w, const Iterators& iterators, Functor& addLocalDofs ) const
    {
      if( u.space().order() > 0 )
      {
        // higher order needs more sophisticated eval
        applyImpl( u, iterators, addLocalDofs, std::true_type() );
      }
      else
      {
        applyImpl( u, iterators, addLocalDofs, std::false_type() );
      }
    }

    template <class GridFunction, class Iterators, class Functor, bool value>
    void applyImpl( const GridFunction& u, const Iterators& iterators, Functor& addLocalDofs, std::integral_constant<bool, value> ) const;

    size_t index( const Entity& entity ) const
    {
      return gridPart().indexSet().index( entity );
    }

    template <class GridFunction, class Functor, class LocalFunction, bool value>
    double applyLocal (const GridFunction& u,
                       const Entity &entiy,
                       Functor& addLocalDofs,
                       LocalFunction& uEn, LocalFunction& uNb,
                       const double dt,
                       std::integral_constant<bool, value> higherOrder ) const;

    const DiscreteFunctionSpaceType& space_;
    // copy model and numFlux for thread safety
    ModelType    model_;
    NumFluxType  numFlux_;
    double time_;
    mutable double dtEst_;

    JacobianRangeType jacLeft_, jacRight_;
    DomainType center_;
    std::vector< DomainType > faceCenters_;
    FaceDomainType localFaceCenter_;

    mutable size_t numberOfElements_;

  }; // end FVOperatorImpl

  template< class DiscreteFunction, class Model, class NumFlux >
  template< class GridFunction, class Iterators, class Functor, bool value >
  inline void FVOperatorImpl< DiscreteFunction, Model, NumFlux >
    ::applyImpl( const GridFunction& u, const Iterators& iterators,
                 Functor& addLocalDofs, std::integral_constant< bool, value > higherOrder ) const
  {
    numberOfElements_ = 0;

    // time step size (using std:min(.,dt) so set to maximum)
    double ws = std::numeric_limits<double>::infinity();

    typedef ConstLocalFunction< GridFunction >  LocalFunctionType ;

    LocalFunctionType uEn( u );
    LocalFunctionType uNb( u );

    // compute update vector and optimum dt in one grid traversal
    const auto endit = iterators.end();
    for( auto it = iterators.begin(); it != endit; ++it )
    {
      ws = applyLocal( u, *it, addLocalDofs, uEn, uNb, ws, higherOrder );
      ++ numberOfElements_;
    }

    dtEst_ = ws;

    // return time step estimate
    // return  ws;
  }

  template< class DiscreteFunction, class Model, class NumFlux >
  template< class GridFunction, class Functor, class LocalFunction, bool value >
  inline double FVOperatorImpl< DiscreteFunction, Model, NumFlux >
    ::applyLocal ( const GridFunction& u,
                   const Entity &entity,
                   Functor& addLocalDofs,
                   LocalFunction& uEn, LocalFunction& uNb,
                   const double prevDt,
                   std::integral_constant< bool, value >
                 ) const
  {
    // initialize model
    model_.setEntity( entity );

    double dt = prevDt;

    auto enIndex = index( entity );

    const Geometry &geo = entity.geometry();

    static const bool higherOrder = value ;

    // cell volume
    const ctype enVolume = geo.volume();

    // 1 over cell volume
    const ctype enVolume_1 = 1.0/enVolume;

    RangeType uLeft;
    RangeType uRight;

    // local entity update
    RangeType enUpdate( 0 );

    // for first order already evaluate here
    if constexpr ( ! higherOrder )
    {
      u.getLocalDofs( entity, uLeft );
    }
    else
    {
      uEn.init( entity );
    }

    if( model_.hasNonStiffSource() )
    {
      LocalEvalEntityType left( entity, center_, center_, enVolume );
      if constexpr ( higherOrder )
      {
        uEn.evaluate( center_, uLeft );
      }
      model_.nonStiffSource( left, uLeft, jacLeft_, enUpdate );
    }

    // the following only makes sense if the model has a flux implemented
    if ( model_.hasFlux() )
    {
      // run through all intersections with neighbors and boundary
      const auto iitend = gridPart().iend( entity );
      for( auto iit = gridPart().ibegin( entity ); iit != iitend; ++iit )
      {
        const Intersection& intersection = *iit;

        // local context interior
        LocalEvalIntersectionType left( entity, intersection, faceCenters_[ intersection.indexInInside() ], localFaceCenter_, enVolume );

        // handle interior face
        if( intersection.neighbor() )
        {
          // access neighbor
          const Entity &neighbor = intersection.outside();

          auto nbIndex = index( neighbor );
          // compute flux from one side only
          if( neighbor.partitionType() != InteriorEntity ||
              enIndex < nbIndex )
          {
            const ctype nbVolume   = ( isCartesian ) ? enVolume   : neighbor.geometry().volume();
            const ctype nbVolume_1 = ( isCartesian ) ? enVolume_1 : 1.0/nbVolume;

            // local context neighbor
            LocalEvalIntersectionType right( neighbor, intersection, faceCenters_[ intersection.indexInOutside() ], localFaceCenter_, nbVolume );

            // evaluate data for higher order
            if constexpr ( higherOrder )
            {
              uNb.init( neighbor );
              std::abort();
              // TODO: need face center here
              uEn.evaluate( center_, uLeft );
              uNb.evaluate( center_, uRight );
            }
            else
            {
              u.getLocalDofs( neighbor, uRight );
            }

            // apply numerical flux
            RangeType fluxLeft, fluxRight;
            const double waveSpeed = numFlux_.numericalFlux( left, right, uLeft, uRight, jacLeft_, jacRight_, fluxLeft, fluxRight );

            // calc update of entity
            enUpdate.axpy( -enVolume_1, fluxLeft );

            // calc update of neighbor
            fluxRight *=  nbVolume_1;
            addLocalDofs( neighbor, fluxRight );

            // compute dt restriction
            dt = std::min( dt, std::min( enVolume, nbVolume ) / waveSpeed );
          }
        }
        // handle boundary face
        else
        {
          // evaluate data for higher order
          if constexpr ( higherOrder )
          {
            // evaluate data
            uEn.evaluate( center_, uLeft );
          }

          const bool hasBndValue = model_.hasBoundaryValue( left );
          RangeType fluxLeft;
          double waveSpeed ;
          if( hasBndValue )
          {
            RangeType fluxRight;
            model_.boundaryValue( left, uLeft, uRight );
            // apply numerical flux
            waveSpeed = numFlux_.numericalFlux( left, left, uLeft, uRight, jacLeft_, jacRight_, fluxLeft, fluxRight );
          }
          else
          {
            // use boundary flux from model
            waveSpeed = model_.boundaryFlux( left, uLeft, jacLeft_, fluxLeft );
          }

          // apply update
          enUpdate.axpy( -enVolume_1, fluxLeft );

          // compute dt restriction
          dt = std::min( dt, enVolume / waveSpeed );
        }
      } // end all intersections
    } // end hasFlux

    // finally update entity value
    addLocalDofs( entity, enUpdate );
    return dt;
  }
} // end namespace detail

// DGAdvectionDiffusionOperatorBase
//---------------------------------

/**
 * \brief advection diffusion operator
 *
 * \note This operator is based on the Pass-Concept
 *
 * \ingroup PassBased
 * \ingroup PassOperator
 */
template< class Traits >
class FVOperator :
  public Fem::SpaceOperatorInterface< typename Traits::DestinationType >
{
  //enum { u = Traits::u,
  //       cdgpass  = Traits::cdgpass };

  enum { polynomialOrder = Traits::polynomialOrder };
  static_assert( polynomialOrder == 0, "FVOperator only works for FiniteVolumeSpace" );

  typedef Fem::SpaceOperatorInterface< typename Traits::DestinationType >  BaseType;

  static const bool threading = Traits :: threading ;

public:
  using BaseType::operator () ;

  // dummy method for a troubled cell indicator to be passed to the
  // limited advection operator
  typedef void* TroubledCellIndicatorType;
  void setTroubledCellIndicator(TroubledCellIndicatorType indicator) {}

  typedef typename Traits::ModelType                    ModelType;
  typedef typename ModelType::ProblemType               ProblemType ;

  typedef typename Traits::GridType                     GridType;

  typedef typename Traits::AdvectionFluxType AdvectionFluxType;

  // for convenience (not used here)
  typedef typename Traits::GridPartType                 GridPartType;
  typedef typename Traits::LimiterIndicatorType         LimiterIndicatorType;
  typedef typename Traits::AdaptationHandlerType        AdaptationType;

  typedef typename Traits::DestinationType DestinationType;
  typedef typename DestinationType::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

  typedef detail::FVOperatorImpl< DestinationType, ModelType, AdvectionFluxType >    FVOperatorType;

  // threading types
  typedef Fem::ThreadIterator< GridPartType >                ThreadIteratorType;

public:
  template< class ExtraParameterTupleImp >
  FVOperator( GridPartType& gridPart,
              const ModelType& model,
              const AdvectionFluxType& numFlux,
              ExtraParameterTupleImp& tuple,
              const std::string name = "",
              const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : FVOperator( gridPart, model, numFlux, name, parameter )
  {}

  FVOperator( GridPartType& gridPart,
              const ModelType& model,
              const AdvectionFluxType& numFlux,
              const std::string name = "",
              const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : gridPart_( gridPart )
    , space_( gridPart_ )
    , iterators_( gridPart_ )
    , impl_( space_, model, numFlux )
    , wTmp_()
    , numberOfElements_( 0 )
    , counter_(0)
    , communicate_( true )
  {
    const size_t threadM1 = MPIManager::numThreads()-1;
    if( wTmp_.size() < threadM1 )
    {
      wTmp_.clear();
      for( size_t i=0; i<threadM1; ++i )
        wTmp_.emplace_back( DestinationType( "wTmp", space_ ) );
    }
  }

  void setAdaptation( AdaptationType& adHandle, double weight = 1 )
  {
    // TODO: implement
  }

  // no indicator available for FV schemes
  const LimiterIndicatorType* indicator () const { return nullptr; }

  void enableIndicator() {}
  void disableIndicator() {}

  void setTime(const double time)
  {
    // use current thread number to obtain correct sizes
    const size_t size = MPIManager::numThreads();
    for( size_t i=0; i<size; ++i )
    {
      impl_[i].setTime( time );
    }
  }

  double timeStepEstimate() const { return timeStepEstimate_; }

  //! evaluate the spatial operator
  void operator()( const DestinationType& u, DestinationType& w ) const
  {
    evaluate( u, w );
  }

  //! evaluate the spatial operator
  template <class GridFunction>
  void operator()( const GridFunction& u, DestinationType& w ) const
  {
    evaluate( u, w );
  }

  //! only evaluate fluxes of operator
  void evaluateOnly( const DestinationType& arg ) const {}

  inline const DiscreteFunctionSpaceType& space() const {
    return space_;
  }

  inline DiscreteFunctionSpaceType& space() {
    return space_;
  }

  int counter() const { return counter_; }
  void called() const { ++counter_; }

  inline void switchupwind()
  {
  }

  template <class Entity, class Intersection, class Quadrature>
  inline void numericalFlux(const DestinationType &u,
                            const Entity &entity, const Entity &nb,
                            const Intersection &intersection,
                            const Quadrature &faceQuadInner, const Quadrature &faceQuadOuter,
                            const int l,
                            typename DestinationType::RangeType &fluxEn,
                            typename DestinationType::RangeType &fluxNb) const
  {
  }


  inline bool hasLimiter () const { return false; }

  inline void limit( DestinationType& U ) const {}
  inline void limit( const DestinationType& arg, DestinationType& U ) const {}

  inline double computeTime() const
  {
    return 0.0;
  }

  inline size_t numberOfElements () const
  {
    return numberOfElements_;
  }

  void printmyInfo(std::string filename) const {}

  virtual std::string description() const { return std::string(""); }

protected:
  const FVOperatorType& impl() const { return *impl_; }

  // update number of interior elements as sum over threads
  double gatherTimeStepEstimate () const
  {
    double dt = std::numeric_limits<double>::infinity();
    numberOfElements_ = 0;
    // use current thread number to obtain correct sizes
    const size_t size = MPIManager::numThreads();
    for( size_t i=0; i<size; ++i )
    {
      dt = std::min( dt, impl_[ i ].timeStepEstimate() );
      numberOfElements_ += impl_[ i ].numberOfElements();
    }
    return dt;
  }

  void gather( DestinationType& w ) const
  {
    // use current thread number to obtain correct sizes
    // thread 0 is using w directly, therefore numThreads-1
    const size_t size = MPIManager::numThreads()-1;
    for( size_t i=0; i<size; ++i )
    {
      w += wTmp_[ i ];
    }
  }

  template <class GridFunction>
  void evaluate( const GridFunction& u, DestinationType& w ) const
  {
    // update counter
    called();

    iterators_.update();
    auto doEval = [this, &u, &w] ()
    {
      int thread = MPIManager::thread();
      DestinationType& wTmp = (thread == 0) ? w : this->wTmp_[ thread-1 ];
      wTmp.clear();

      this->impl().evaluate( u, wTmp, this->iterators_ );
    };

    try {
      // execute in parallel
      MPIManager :: run ( doEval );

      // collect all values from wTmp
      gather( w );

      // update time step estimate and numberOfElements
      timeStepEstimate_ = gatherTimeStepEstimate();
    }
    catch ( const SingleThreadModeError& e )
    {
      // reset w from previous entries
      w.clear();
      // re-run in single thread mode if previous attempt failed
      impl().evaluate( u, w, space_ );

      // update time step estimate
      timeStepEstimate_ = impl().timeStepEstimate();
      numberOfElements_ = impl().numberOfElements();
    }

    // synchronize result
    if( communicate_ )
      w.communicate();
  }

  GridPartType&                       gridPart_;
  DiscreteFunctionSpaceType           space_;
  mutable ThreadIteratorType          iterators_;
  ThreadSafeValue< FVOperatorType >   impl_;
  mutable std::vector< DestinationType >  wTmp_;
  mutable double                      timeStepEstimate_;
  mutable size_t                      numberOfElements_;
  mutable int                         counter_;
  bool                                communicate_;
};

} // end namespace Fem
} // end namespace Dune

#endif // #ifndef FVSCHEME_HH
