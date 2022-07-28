#ifndef DUNE_FEM_DG_THREADPASS_HH
#define DUNE_FEM_DG_THREADPASS_HH

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/quadrature/cornerpointset.hh>

#include <dune/fem-dg/pass/pass.hh>

#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{
namespace Fem
{

  struct NonBlockingCommParameter
  {
    static bool nonBlockingCommunication()
    {
      // non-blocking communication is only avaiable in smp mode
      // because the implementation is done in ThreadPass
      return Fem :: Parameter :: getValue< bool > ("femdg.nonblockingcomm", false );
    }
  };

  template < class DiscreteFunction >
  class DeleteCommunicatedDofs : public Fem::ParallelScalarProduct < DiscreteFunction >
  {
    typedef Fem::ParallelScalarProduct < DiscreteFunction >  BaseType;

  public:
    typedef typename BaseType :: DiscreteFunctionType       DiscreteFunctionType;
    typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;


    //! constructor taking space
    explicit DeleteCommunicatedDofs( const DiscreteFunctionSpaceType &space )
      : BaseType( space )
    {
    }

    //! delete ghost values again, otherwise the Newton solver
    //! of the implicit ODE solvers wont converge
    void deleteCommunicatedDofs( DiscreteFunctionType& df ) const
    {
#if HAVE_MPI
      const auto& auxiliaryDofs = this->auxiliaryDofs();

      // don't delete the last since this is the overall size
      const int auxiliarySize = auxiliaryDofs.size() - 1;
      for(int auxiliary = 0; auxiliary<auxiliarySize; ++auxiliary)
      {
        typedef typename DiscreteFunctionType :: DofBlockPtrType DofBlockPtrType;
        DofBlockPtrType block = df.block( auxiliaryDofs[ auxiliary ] );
        const int blockSize = DiscreteFunctionType :: DiscreteFunctionSpaceType :: localBlockSize ;
        for(int l = 0; l<blockSize; ++l )
          (*block)[ l ] = 0;
      }
#endif
    }
  };

  template < class DestinationType >
  class NonBlockingCommHandle
  {
    typedef typename DestinationType :: DiscreteFunctionSpaceType :: CommunicationManagerType
          :: NonBlockingCommunicationType  NonBlockingCommunicationType;

    mutable std::unique_ptr< NonBlockingCommunicationType > nonBlockingComm_;
    const bool useNonBlockingComm_ ;
  public:
    NonBlockingCommHandle()
      : nonBlockingComm_(),
        useNonBlockingComm_( NonBlockingCommParameter :: nonBlockingCommunication() )
      {}

    NonBlockingCommHandle( const NonBlockingCommHandle& other )
      : nonBlockingComm_(),
        useNonBlockingComm_( other.useNonBlockingComm_ )
    {}

    bool nonBlockingCommunication () const { return useNonBlockingComm_; }

    ~NonBlockingCommHandle()
    {
      // make sure all communications have been finished
      assert( ! nonBlockingComm_ );
    }

    // send data
    void initComm( const DestinationType& dest ) const
    {
      if( nonBlockingCommunication() && ! nonBlockingComm_ )
      {
        nonBlockingComm_.reset( new NonBlockingCommunicationType(
            dest.space().communicator().nonBlockingCommunication() ) );

        // perform send operation
        nonBlockingComm_->send( dest );
      }
    }

    // receive data
    void receiveComm( const DestinationType& destination ) const
    {
      if( nonBlockingCommunication() && nonBlockingComm_ )
      {
        DestinationType& dest = const_cast< DestinationType& > ( destination );
        nonBlockingComm_->receive( dest );
        nonBlockingComm_.reset();
      }
    }

    // cleanup possibly overwritten ghost values
    void finalizeComm( const DestinationType& dest ) const
    {
      // only delete non-interior dofs in non-blocking mode
      if( nonBlockingCommunication() )
      {
        // make sure communication has been finished
        assert( ! nonBlockingComm_ );
        DeleteCommunicatedDofs< DestinationType > delDofs( dest.space() );
        delDofs.deleteCommunicatedDofs( const_cast< DestinationType& > ( dest ) );
      }
    }
  };

  /**
   * \brief Pass which turns a simple pass into a threading pass.
   *
   * \ingroup Pass
   *
   * \tparam InnerPass The pass which should be turned into a threading pass.
   * \tparam ThreadIterator An iterator.
   * \tparam nonblockingcomm Boolean indicating whether blocking or non blocking communication is used.
   */
  template < class InnerPass,
             class ThreadIterator,
             bool nonblockingcomm = true >
  class ThreadPass :
    public LocalPass< typename InnerPass :: DiscreteModelType,
                      typename InnerPass :: PreviousPassType,
                      InnerPass :: passId >
  {
    typedef ThreadPass< InnerPass, ThreadIterator, nonblockingcomm > ThisType;
  public:
    typedef InnerPass                                                     InnerPassType;
    typedef typename InnerPass::DiscreteModelType                         DiscreteModelType;
    typedef typename InnerPass::PreviousPassType                          PreviousPassType;

    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelType, PreviousPassType, InnerPass :: passId>  BaseType;

    // Types from the base class
    typedef typename BaseType::EntityType                                 EntityType;
    typedef typename BaseType::ArgumentType                               ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType           DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType      VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType        FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Iterator over the space
    typedef typename DiscreteFunctionSpaceType::IteratorType              IteratorType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType                  GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType              GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType                DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType                 RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType         JacobianRangeType;

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType               IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType                       IntersectionType;
    typedef typename GridType::template Codim<0>::Geometry                Geometry;

    // Various other types
    typedef typename DestinationType::LocalFunctionType                   LocalFunctionType;

    typedef NonBlockingCommHandle< DestinationType >                      NonBlockingCommHandleType ;

    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };

    // type of local id set
    typedef typename GridPartType::IndexSetType                           IndexSetType;

    // type of thread iterators (e.g. Fem::DomainDecomposedIteratorStorage or Fem::ThreadIterator)
    typedef ThreadIterator                                                ThreadIteratorType;

  protected:
    using BaseType::spc_;
    using BaseType::previousPass_ ;

  public:
    using BaseType::pass;

    /** \brief Constructor
     * \param discreteModel Actual discrete model
     * \param pass Previous pass
     * \param spc Space belonging to the discrete function local to this pass
     * \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order
     * \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order
     */
    ThreadPass(const DiscreteModelType& discreteModel,
               PreviousPassType& pass,
               const DiscreteFunctionSpaceType& spc,
               const int volumeQuadOrd = -1,
               const int faceQuadOrd = -1)
      : ThreadPass( discreteModel, pass, spc, Dune::Fem::Parameter::container(), volumeQuadOrd, faceQuadOrd )
    {}

    /** \brief Constructor
     * \param discreteModel Actual discrete model
     * \param pass Previous pass
     * \param spc Space belonging to the discrete function local to this pass
     * \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order
     * \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order
     */
    ThreadPass(const DiscreteModelType& discreteModel,
               PreviousPassType& pass,
               const DiscreteFunctionSpaceType& spc,
               const Dune::Fem::ParameterReader &parameter,
               const int volumeQuadOrd = -1,
               const int faceQuadOrd = -1) :
      BaseType(pass, spc),
      delDofs_( spc ),
      iterators_( spc.gridPart() ),
      singleDiscreteModel_( discreteModel ),
      discreteModels_( Fem::MPIManager::maxThreads() ),
      passes_( Fem::MPIManager::maxThreads() ),
      passComputeTime_( Fem::MPIManager::maxThreads(), 0.0 ),
      firstStage_( false ),
      arg_(0), dest_(0),
      nonBlockingComm_(),
      volumeQuadOrd_( volumeQuadOrd ),
      faceQuadOrd_( faceQuadOrd ),
      firstCall_( true ),
      requireCommunication_( true ),
      parameter_( &parameter ),
      sumComputeTime_( parameter.getValue<bool>("fem.parallel.sumcomputetime", false ) )
    {
      // initialize thread pass here since it otherwise fails when parameters
      // are passed from the Python side
#if 0
      // initialize each thread pass by the thread itself to avoid NUMA effects
      {
        // lambda calling parallel code
        auto threadedRun = [this] () { this->runThread(); };

        // see threadpool.hh
        Fem :: MPIManager :: runLocked( threadedRun );
      }
#else
      {
        // fall back if the above does not work
        const int maxThreads = Fem::MPIManager::maxThreads();
        for(int i=0; i<maxThreads; ++i)
        {
          createInnerPass( i, i == 0 );
        }
      }
#endif

#ifndef NDEBUG
      {
        // check that all objects have been created
        const int maxThreads = Fem::MPIManager::maxThreads();
        for(int i=0; i<maxThreads; ++i)
        {
          assert( discreteModels_[ i ] );
          assert( passes_[ i ] );
        }
      }

      if( Fem :: Parameter :: verbose() )
        std::cout << "Thread Pass initialized\n";
#endif
      // get information about communication
      requireCommunication_ = passes_[ 0 ]->requireCommunication();
      parameter_ = nullptr;
    }

    template <class TroubledCellIndicatorType>
    void setTroubledCellIndicator( TroubledCellIndicatorType indicator )
    {
      const int maxThreads = Fem::MPIManager::maxThreads();
      for(int i=0; i<maxThreads; ++i)
      {
        passes_[ i ]->setTroubledCellIndicator( indicator );
      }
    }

    virtual ~ThreadPass () {}

    template <class AdaptationType>
    void setAdaptation( AdaptationType& adHandle, double weight )
    {
      const int maxThreads = Fem::MPIManager::maxThreads();
      for(int thread=0; thread<maxThreads; ++thread)
      {
        discreteModels_[ thread ]->setAdaptation(
            adHandle, weight, &iterators_.filter( thread ) );
        // add filter in thread parallel versions
      }
    }

    //! call appropriate method on all internal passes
    void enable() const
    {
      const int maxThreads = Fem::MPIManager::maxThreads();
      for(int thread=0; thread<maxThreads; ++thread)
      {
        pass( thread ).enable();
      }
    }

    //! call appropriate method on all internal passes
    void disable() const
    {
      const int maxThreads = Fem::MPIManager::maxThreads();
      for(int thread=0; thread<maxThreads; ++thread)
      {
        pass( thread ).disable();
      }
    }

    //! print tex info
    void printTexInfo(std::ostream& out) const
    {
      BaseType::printTexInfo(out);
      pass( 0 ).printTexInfo(out);
    }

    //! Estimate for the timestep size
    double timeStepEstimateImpl() const
    {
      double dtMin = pass( 0 ).timeStepEstimateImpl();
      const int numThreads = Fem::MPIManager::numThreads();
      for( int i = 1; i < numThreads ; ++i)
      {
        dtMin = std::min( dtMin, pass( i ).timeStepEstimateImpl() );
      }
      return dtMin;
    }

  protected:
    enum SkipMode { skipNone, skipInterior, skipNonInterior };
    //! returns true for flux evaluation if neighbor
    //! is on same thread as entity
    template <SkipMode mode>
    struct NBChecker
    {
      const ThreadIteratorType& storage_;
      const int myThread_;

#ifndef NDEBUG
      mutable int counter_;
      mutable int nonEqual_;
#endif
      NBChecker( const ThreadIteratorType& st,
                 const int myThread )
        : storage_( st ),
          myThread_( myThread )
#ifndef NDEBUG
          , counter_( 0 )
          , nonEqual_( 0 )
#endif
      {}

      // returns true if niehhbor can be updated
      bool operator () ( const EntityType& en, const EntityType& nb ) const
      {
#ifndef NDEBUG
        ++counter_;
        if( myThread_ != storage_.thread( nb ) )
          ++nonEqual_;
#endif
        // storage_.thread can also return negative values in which case the
        // update of the neighbor is skipped, e.g. for ghost elements
        return myThread_ == storage_.thread( nb );
      }

      //! return true if actually some intersections would be skipped
      bool isActive () const { return mode != skipNone ; }

      //! returns true if the intersection with neighbor nb should be skipped
      bool skipIntersection( const EntityType& nb ) const
      {
        // noskip means all intersections are considered
        switch( mode )
        {
          case skipNone:         return false;
          case skipInterior:     return nb.partitionType() == InteriorEntity;
          case skipNonInterior:  return nb.partitionType() != InteriorEntity;
          default: return false ;
        }
        return false ;
      }
    };

    InnerPass& pass( const int thread ) const
    {
      assert( (int) passes_.size() > thread );
      return *( passes_[ thread ] );
    }

    using BaseType::time ;
    using BaseType::computeTime_ ;
    using BaseType::destination_ ;
    using BaseType::destination;

  public:
    using BaseType::receiveCommunication;
    using BaseType::space;

    //! switch upwind direction
    void switchUpwind()
    {
      const int maxThreads = Fem::MPIManager::maxThreads();
      for(int i=0; i<maxThreads; ++i )
        discreteModels_[ i ]->switchUpwind();
    }

    //! overload compute method to use thread iterators
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      // reset number of elements (see LocalPass)
      this->numberOfElements_ = 0;

      // set time for all passes, this is used in prepare of pass
      // and therefore has to be done before prepare is called
      const int numThreads = Fem::MPIManager::numThreads();
      for(int i=0; i<numThreads; ++i )
      {
        // set time to each pass
        pass( i ).setTime( time() );
      }

      try
      {
        // update thread iterators in case grid changed
        iterators_.update();

        // call prepare before parallel area
        const int numThreads = Fem::MPIManager::numThreads();
        pass( 0 ).prepare( arg, dest, true );
        passComputeTime_[ 0 ] = 0.0 ;
        for(int i=1; i<numThreads; ++i )
        {
          // prepare pass (make sure pass doesn't clear dest, this will conflict)
          pass( i ).prepare( arg, dest, false );
          passComputeTime_[ i ] = 0.0 ;
        }
        firstStage_ = true ;

        arg_  = &arg ;
        dest_ = &dest ;

        // lambda calling parallel code
        auto threadedRun = [this] () { this->runThread(); };

        ////////////////////////////////////////////////////////////
        // BEGIN PARALLEL REGION, first stage, element integrals
        ////////////////////////////////////////////////////////////
        {
          // see threadpool.hh
          Fem :: MPIManager :: run( threadedRun );
        }
        /////////////////////////////////////////////////
        // END PARALLEL REGION
        /////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // BEGIN PARALLEL REGION, second stage, surface integrals
        // only for non-blocking communication
        ////////////////////////////////////////////////////////////
        if( useNonBlockingCommunication() )
        {
          // mark second stage
          firstStage_ = false ;

          // see threadpool.hh
          Fem :: MPIManager :: run( threadedRun );
        }
        /////////////////////////////////////////////////
        // END PARALLEL REGION
        /////////////////////////////////////////////////

        double accCompTime = 0.0;
        double ratioMaster = 1.0;
        for(int i=0; i<numThreads; ++i )
        {
          // get number of elements
          this->numberOfElements_ += pass( i ).numberOfElements();

          if( sumComputeTime_ )
          {
            accCompTime += passComputeTime_[ i ];
          }
          else
          {
            // accumulate time
            accCompTime = std::max( passComputeTime_[ i ], accCompTime );
          }

          // thread 0 should have longer compute time since also communication has to be done
          if( passComputeTime_[ 0 ] > 0 )
            ratioMaster = std::min( ratioMaster, double(passComputeTime_[ i ] / passComputeTime_[ 0 ] ) );
        }

        //std::cout << "ratio = " << ratioMaster << std::endl;
        // store ration information for next partitioning
        //iterators_.setMasterRatio( ratioMaster );

        // increase compute time
        computeTime_ += accCompTime ;

      } // end if first call
      catch (const Dune::Fem::SingleThreadModeError& e )
      {
        // reset all passes
        for(int i=0; i<numThreads; ++i )
        {
          pass( i ).finalize( arg, dest, false );
        }

        // for the first call we need to receive data already here,
        // since the flux calculation is done at once
        if( useNonBlockingCommunication() )
        {
          // RECEIVE DATA, send was done on call of operator() (see pass.hh)
          receiveCommunication( arg );
        }

        // Compute pass in single thread mode
        pass( 0 ).compute( arg, dest );
        // get number of elements
        this->numberOfElements_ = pass( 0 ).numberOfElements();
      }

      // if useNonBlockingComm_ is disabled then communicate here if communication is required
      if( requireCommunication_ && ! nonBlockingComm_.nonBlockingCommunication() )
      {
        assert( dest_ );

        if( dest_ )
        {
          // communicate calculated function
          dest_->communicate();
        }
      }

      // remove pointers
      arg_  = 0;
      dest_ = 0;
    }

    //! return true if communication is necessary and non-blocking should be used
    bool useNonBlockingCommunication() const
    {
      return requireCommunication_ && nonBlockingComm_.nonBlockingCommunication();
    }

    void initComm() const
    {
      if( useNonBlockingCommunication() && destination_ )
        nonBlockingComm_.initComm( destination() );
    }

    void receiveComm() const
    {
      if( useNonBlockingCommunication() && destination_ )
        nonBlockingComm_.receiveComm( destination() );
    }

    //! parallel section of compute
    void createInnerPass(const int thread, const bool isMainThread ) const
    {
      // initialization (called from constructor of this class)
      if( ! passes_[ thread ] )
      {
        // use separate discrete discrete model for each thread
        discreteModels_[ thread ].reset( new DiscreteModelType( singleDiscreteModel_ ));
        const bool verbose = Dune::Fem::Parameter::verbose() && isMainThread;
        if( parameter_ )
        {
          // create dg passes, the last bool disables communication in the pass itself
          passes_[ thread ].reset(new InnerPassType( *discreteModels_[ thread ], previousPass_, space(), *parameter_, volumeQuadOrd_, faceQuadOrd_, verbose ));
        }
        else
        {
          // create dg passes, the last bool disables communication in the pass itself
          passes_[ thread ].reset(new InnerPassType( *discreteModels_[ thread ], previousPass_, space(), volumeQuadOrd_, faceQuadOrd_, verbose ));
        }

      }
    }

    //! parallel section of compute
    void runThread() const
    {
      const int thread = Fem::MPIManager::thread() ;
      // make sure thread 0 is master thread
      assert( (thread == 0) == Fem::MPIManager::isMainThread() );

      // initialization (called from constructor of this class)
      if( ! passes_[ thread ] )
      {
        createInnerPass( thread, Fem::MPIManager::isMainThread() );
        return ;
      }

      //! get pass for my thread
      InnerPassType& myPass = pass( thread );

      // stop time
      Dune::Timer timer ;

      const bool computeInteriorIntegrals = firstStage_;

      // Iterator is of same type as the space iterator
      typedef typename ThreadIteratorType :: IteratorType Iterator;

      if( useNonBlockingCommunication() )
      {
        if ( computeInteriorIntegrals )
        {
          // create NB checker, skip if neighbor is NOT interior
          NBChecker< skipNonInterior > nbChecker( iterators_, thread );

          const Iterator endit = iterators_.end();
          for (Iterator it = iterators_.begin(); it != endit; ++it)
          {
            const auto entity = *it;
            assert( iterators_.thread( entity ) == thread );
            myPass.applyLocalInterior( entity, nbChecker );
          }

          // receive ghost data (only master thread)
          if( thread == 0 && requireCommunication_ )
          {
            // RECEIVE DATA, send was done on call of operator() (see pass.hh)
            receiveCommunication( *arg_ );
          }
        }
        else
        {
          // create NB checker, skip if neighbors is interior
          NBChecker< skipInterior > nbChecker( iterators_, thread );

          const Iterator endit = iterators_.end();
          for (Iterator it = iterators_.begin(); it != endit; ++it)
          {
            const auto entity = *it;
            assert( iterators_.thread( entity ) == thread );
            myPass.applyLocalProcessBoundary( entity, nbChecker );
          }

          assert( arg_ );
          // dest can also be null pointer
          // when the operator is evaluated only
          // for evaluation of the estimators

          // finalize pass (make sure communication is done in case of thread parallel
          // program, this would give conflicts)
          myPass.finalize(*arg_, *dest_, false );
        }
      }
      else
      {
        // create NB checker, noSkip of intersections
        NBChecker< skipNone > nbChecker( iterators_, thread );

        const Iterator endit = iterators_.end();
        for (Iterator it = iterators_.begin(); it != endit; ++it)
        {
          const auto entity = *it;
          assert( iterators_.thread( entity ) == thread );
          myPass.applyLocal( entity, nbChecker );
        }

        assert( arg_ );
        // dest can also be null pointer
        // when the operator is evaluated only
        // for evaluation of the estimators

        // finalize pass (make sure communication is not done in case of thread parallel
        // program, this would give conflicts)
        myPass.finalize(*arg_, *dest_, false );
      }

      // accumulate compute time for this thread
      passComputeTime_[ thread ] += timer.elapsed();
    }


    //! In the preparations, store pointers to the actual arguments and
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // we use prepare of internal operator
      std::abort();
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      // we use finalize of internal operator
      std::abort();
    }

    void applyLocal( const EntityType& en) const
    {
      // we use applyLocal of internal operator
      std::abort();
    }

  private:
    ThreadPass();
    ThreadPass(const ThreadPass&);
    ThreadPass& operator=(const ThreadPass&);

  protected:
    // create an instance of the parallel scalarproduct here to avoid
    // deleting on every call of finalizeComm
    DeleteCommunicatedDofs< DestinationType > delDofs_;

    mutable ThreadIteratorType iterators_;
    const DiscreteModelType& singleDiscreteModel_;
    mutable std::vector< std::unique_ptr< DiscreteModelType > > discreteModels_;
    mutable std::vector< std::unique_ptr< InnerPassType > > passes_;
    mutable std::vector< double > passComputeTime_;
    mutable bool firstStage_;

    // temporary variables
    mutable const ArgumentType* arg_;
    mutable DestinationType* dest_;

    // non-blocking communication handler
    NonBlockingCommHandleType nonBlockingComm_;

    const int volumeQuadOrd_;
    const int faceQuadOrd_;
    mutable bool firstCall_;
    bool requireCommunication_;

    const Dune::Fem::ParameterReader* parameter_;

    const bool sumComputeTime_;
  };
//! @}

} // end namespace
} // end namespace Dune

#endif
