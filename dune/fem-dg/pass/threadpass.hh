#ifndef DUNE_THREADPASS_HH
#define DUNE_THREADPASS_HH

#include <dune/common/fvector.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/misc/threaditerator.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/dgmodelcaller.hh>
#include <dune/fem/pass/pass.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/allgeomtypes.hh> 
#include <dune/fem/space/common/arrays.hh> 

#include "threadhandle.hh"
#include "domaindecomposed.hh"

namespace Dune {

  struct NonBlockingCommParameter 
  {
    static bool nonBlockingCommunication() 
    {
      // non-blocking communication is only avaiable in smp mode
#ifdef USE_SMP_PARALLEL
      return Fem :: Parameter :: getValue< bool > ("femdg.nonblockingcomm", false );
#else 
      return false;
#endif
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
      typedef typename BaseType :: SlaveDofsType  SlaveDofsType;
      SlaveDofsType &slaves = this->slaveDofs();

      // don't delete the last since this is the overall Size 
      const int slaveSize = slaves.size() - 1;
      for(int slave = 0; slave<slaveSize; ++slave)
      {
        typedef typename DiscreteFunctionType :: DofBlockPtrType DofBlockPtrType;
        DofBlockPtrType block = df.block( slaves[ slave ] );
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

    mutable NonBlockingCommunicationType* nonBlockingComm_;
    const bool useNonBlockingComm_ ;
  public:
    NonBlockingCommHandle() 
      : nonBlockingComm_( 0 ),
        useNonBlockingComm_( NonBlockingCommParameter :: nonBlockingCommunication() )
      {}

    NonBlockingCommHandle( const NonBlockingCommHandle& other ) 
      : nonBlockingComm_( 0 ),
        useNonBlockingComm_( other.useNonBlockingComm_ ) 
    {}

    bool nonBlockingCommunication () const { return useNonBlockingComm_; }

    ~NonBlockingCommHandle() 
    {
      // make sure all communications have been finished
      assert( nonBlockingComm_ == 0 );
    }

    // send data 
    void initComm( const DestinationType& dest ) const 
    {
      if( nonBlockingCommunication() && nonBlockingComm_ == 0 )
      {
        nonBlockingComm_ = new NonBlockingCommunicationType(
            dest.space().communicator().nonBlockingCommunication() );

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
        delete nonBlockingComm_;
        nonBlockingComm_ = 0;
      }
    }

    // cleanup possibly overwritten ghost values 
    void finalizeComm( const DestinationType& dest ) const 
    {
      // only delete non-interior dofs in non-blocking mode 
      if( nonBlockingCommunication() ) 
      {
        // make sure communication has been finished
        assert( nonBlockingComm_ == 0 );
        DeleteCommunicatedDofs< DestinationType > delDofs( dest.space() );
        delDofs.deleteCommunicatedDofs( const_cast< DestinationType& > ( dest ) );
      }
    }
  };

  template < class InnerPass, bool nonblockingcomm = true > 
  class ThreadPass :
    public Fem::LocalPass< typename InnerPass :: DiscreteModelType,
                           typename InnerPass :: PreviousPassType, 
                           InnerPass :: passId > 
  {
    typedef ThreadPass< InnerPass > ThisType;
  public:
    typedef InnerPass InnerPassType;
    typedef typename InnerPass :: DiscreteModelType  DiscreteModelType;
    typedef typename InnerPass :: PreviousPassType   PreviousPassType;

    //- Typedefs and enums
    //! Base class
    typedef Fem::LocalPass< DiscreteModelType, PreviousPassType, InnerPass :: passId>  BaseType;

    // Types from the base class
    typedef typename BaseType::EntityType  EntityType;
    typedef typename EntityType :: EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Iterator over the space
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef Fem::CombinedSelector< ThisType , SelectorType > CombinedSelectorType;
    typedef Fem::DGDiscreteModelCaller< DiscreteModelType 
                                   , ArgumentType 
                                   , CombinedSelectorType
                                 > DiscreteModelCallerType;

    typedef NonBlockingCommHandle< DestinationType > NonBlockingCommHandleType ;

    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };

    // type of local id set 
    typedef typename GridPartType::IndexSetType IndexSetType; 

    //typedef Fem::ThreadIterator< DiscreteFunctionSpaceType > ThreadIteratorType;
    typedef Fem::DomainDecomposedIteratorStorage< DiscreteFunctionSpaceType > ThreadIteratorType;

    // type of adaptation handler 
    typedef typename DiscreteModelType :: AdaptationHandlerType AdaptationHandlerType ;

  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order 
    ThreadPass(DiscreteModelType& problem, 
               PreviousPassType& pass, 
               const DiscreteFunctionSpaceType& spc,
               const int volumeQuadOrd = -1,
               const int faceQuadOrd = -1) :
      BaseType(pass, spc),
      delDofs_( spc ),
      iterators_( spc ),
      singleProblem_( problem ),
      problems_( Fem::ThreadManager::maxThreads() ),
      passes_( Fem::ThreadManager::maxThreads() ),
      passComputeTime_( Fem::ThreadManager::maxThreads(), 0.0 ),
      passStage_( Fem::ThreadManager::maxThreads(), false ),
      arg_(0), dest_(0),
      nonBlockingComm_(),
      numberOfElements_( 0 ),
      firstCall_( true ),
      sumComputeTime_( Fem :: Parameter :: getValue<bool>("fem.parallel.sumcomputetime", false ) )
    {
      const int maxThreads = Fem::ThreadManager::maxThreads();
      for(int i=0; i<maxThreads; ++i)
      {
        // use serparate discrete problem for each thread 
        problems_[ i ] = new DiscreteModelType( problem );
        // create dg passes, the last bool disables communication in the pass itself
        passes_[ i ]   = new InnerPassType( *problems_[ i ], pass, spc, volumeQuadOrd, faceQuadOrd, false );
      }
#ifndef NDEBUG
      if( Fem :: Parameter :: verbose() )
        std::cout << "Thread Pass initialized\n";
#endif
    }

    virtual ~ThreadPass () 
    {
      for(int i=0; i<Fem::ThreadManager::maxThreads(); ++i)
      {
        delete passes_[ i ];
        delete problems_[ i ];
      }
    }

    void setAdaptationHandler( AdaptationHandlerType& adHandle, double weight ) 
    {
      const int maxThreads = Fem::ThreadManager::maxThreads();
      for(int thread=0; thread<maxThreads; ++thread)
      {
        problems_[ thread ]->setAdaptationHandler( adHandle, 
#ifdef NSMOD_USE_SMP_PARALLEL
            iterators_.filter( thread ), // add filter in thread parallel versions 
#endif
            weight );
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
      const int maxThreads = Fem::ThreadManager::maxThreads();
      for( int i = 1; i < maxThreads ; ++i)
      {
        dtMin = std::min( dtMin, pass( i ).timeStepEstimateImpl() );
      }
      return dtMin;
    }

  protected:  
    //! returns true for flux evaluation if neighbor 
    //! is on same thread as entity  
    struct NBChecker 
    {
      const ThreadIteratorType& storage_;
      const int myThread_; 

#ifndef NDEBUG
      mutable int counter_; 
      mutable int nonEqual_;
#endif
      NBChecker( const ThreadIteratorType& st, 
                 const int myThread = Fem::ThreadManager::thread() ) 
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
    };

    InnerPass& pass( const int thread ) const 
    {
      assert( (int) passes_.size() > thread );
      return *( passes_[ thread ] );
    }

    using BaseType :: time ;
    using BaseType :: computeTime_ ;
    using BaseType :: destination_ ;
    using BaseType :: destination ;
    using BaseType :: receiveCommunication ;

  public:  
    //! return number of elements visited on last application
    size_t numberOfElements () const 
    {
      return numberOfElements_; 
    }

    //! switch upwind direction
    void switchUpwind() 
    {
      const int maxThreads = Fem::ThreadManager::maxThreads();
      for(int i=0; i<maxThreads; ++i ) 
        problems_[ i ]->switchUpwind(); 
    }

    //! overload compute method to use thread iterators 
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      const bool updateAlso = (& dest != 0);
      if( updateAlso ) 
      {
        // clear destination 
        dest.clear();
      }

      // reset number of elements 
      numberOfElements_ = 0;

      // set time for all passes, this is used in prepare of pass 
      // and therefore has to be done before prepare is called
      const int maxThreads = Fem::ThreadManager::maxThreads();
      for(int i=0; i<maxThreads; ++i ) 
      {
        // set time to each pass 
        pass( i ).setTime( time() );
      }

      // for the first call we only run on one thread to avoid 
      // clashes with the singleton storages for quadratures 
      // and base function caches etc.
      // after one grid traversal everything should be set up 
      if( firstCall_ ) 
      {
        //! get pass for my thread  
        InnerPassType& myPass = pass( 0 );
        //std::cout << "Thread Pass :: First Call !!! " << this << std::endl;

        // stop time 
        Timer timer ;

        // pepare 
        myPass.prepare( arg, dest );

        // for the first call we need to receive data already here,
        // since the flux calculation is done at once
        if( nonBlockingComm_.nonBlockingCommunication() ) 
        {
          // RECEIVE DATA, send was done on call of operator() (see pass.hh)
          receiveCommunication( arg );
        }

        // Iterator is of same type as the space iterator 
        typedef typename DiscreteFunctionSpaceType :: IteratorType Iterator;
        const Iterator endit = iterators_.space().end();
        for (Iterator it = iterators_.space().begin(); it != endit; ++it)
        {
          myPass.applyLocal( *it );
        }

        // finalize pass 
        myPass.finalize(arg, dest);

        // get number of elements 
        numberOfElements_ = myPass.numberOfElements(); 
        // store time 
        computeTime_ += timer.elapsed();

        // set tot false since first call has been done
        firstCall_ = false ;
      }
      else 
      {
        // update thread iterators in case grid changed 
        iterators_.update();

        // call prepare before parallel area 
        const int maxThreads = Fem::ThreadManager::maxThreads();
        for(int i=0; i<maxThreads; ++i ) 
        {
          // prepare pass (make sure pass doesn't clear dest, this will conflict)
          pass( i ).prepare( arg, dest );
          passComputeTime_[ i ] = 0.0 ;
          passStage_[ i ] = true ;
        }
        
        arg_  = &arg ; 
        dest_ = &dest ;

        ////////////////////////////////////////////////////////////
        // BEGIN PARALLEL REGION, first stage, element integrals  
        ////////////////////////////////////////////////////////////
        {
          // see threadhandle.hh 
          Fem :: ThreadHandle :: run( *this ); 
        }
        /////////////////////////////////////////////////
        // END PARALLEL REGION 
        /////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // BEGIN PARALLEL REGION, second stage, surface integrals 
        // only for non-blocking communication 
        ////////////////////////////////////////////////////////////
        if( nonBlockingComm_.nonBlockingCommunication() ) 
        {
          // mark second stage 
          for(int i=0; i<maxThreads; ++i ) 
            passStage_[ i ] = false ;

          // RECEIVE DATA, send was done on call of operator() (see pass.hh)
          receiveCommunication( arg );

          // see threadhandle.hh 
          Fem :: ThreadHandle :: run( *this ); 
        }
        /////////////////////////////////////////////////
        // END PARALLEL REGION 
        /////////////////////////////////////////////////

        arg_  = 0;
        dest_ = 0;

        double accCompTime = 0.0;
        for(int i=0; i<maxThreads; ++i ) 
        {
          // get number of elements 
          numberOfElements_ += pass( i ).numberOfElements(); 

          if( sumComputeTime_ ) 
          {
            accCompTime += passComputeTime_[ i ];
          }
          else 
          {
            // accumulate time 
            accCompTime = std::max( passComputeTime_[ i ], accCompTime );
          }
        }

        // increase compute time 
        computeTime_ += accCompTime ;

      } // end if first call 

      // set max time steps 
      setMaxTimeSteps();

      // if useNonBlockingComm_ is disabled then communicate here
      if( ! nonBlockingComm_.nonBlockingCommunication() && updateAlso ) 
      {
        // communicate calculated function 
        dest.communicate();
      }
    }

    void initComm() const 
    {
      if( nonBlockingComm_.nonBlockingCommunication() && destination_ ) 
        nonBlockingComm_.initComm( destination() );
    }

    void receiveComm() const
    {
      if( nonBlockingComm_.nonBlockingCommunication() && destination_ ) 
        nonBlockingComm_.receiveComm( destination() );
    }

    //! parallel section of compute 
    void runThread() const
    {
      const int thread = Fem::ThreadManager::thread() ;
      //! get pass for my thread  
      InnerPassType& myPass = pass( thread );

      // stop time 
      Timer timer ;

      const bool computeElementIntegral = passStage_[ thread ];

      // create NB checker 
      NBChecker nbChecker( iterators_, thread );

      // Iterator is of same type as the space iterator 
      typedef typename ThreadIteratorType :: IteratorType Iterator;

      if( nonBlockingComm_.nonBlockingCommunication() ) 
      {
        if ( computeElementIntegral ) 
        {
          const Iterator endit = iterators_.end();
          for (Iterator it = iterators_.begin(); it != endit; ++it)
          {
            assert( iterators_.thread( *it ) == thread );
            myPass.elementIntegral( *it );
          }
        }
        else 
        {
          const Iterator endit = iterators_.end();
          for (Iterator it = iterators_.begin(); it != endit; ++it)
          {
            assert( iterators_.thread( *it ) == thread );
            myPass.surfaceIntegral( *it, nbChecker );
          }

          assert( arg_ );
          // dest can also be null pointer 
          // when the operator is evaluated only 
          // for evaluation of the estimators 

          // finalize pass (make sure communication is done in case of thread parallel
          // program, this would give conflicts)
          myPass.finalize(*arg_, *dest_);
        }
      }
      else 
      {
        const Iterator endit = iterators_.end();
        for (Iterator it = iterators_.begin(); it != endit; ++it)
        {
          assert( iterators_.thread( *it ) == thread );
          myPass.applyLocal( *it, nbChecker );
        }

        assert( arg_ );
        // dest can also be null pointer 
        // when the operator is evaluated only 
        // for evaluation of the estimators 

        // finalize pass (make sure communication is not done in case of thread parallel
        // program, this would give conflicts)
        myPass.finalize(*arg_, *dest_);
      }

      passComputeTime_[ thread ] += timer.elapsed();
    }


    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // we use prepare of internal operator
      abort();
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      // we use finalize of internal operator
      abort();
    }

    void applyLocal( const EntityType& en) const 
    {
      // we use applyLocal of internal operator
      abort();
    }
  protected:
    void setMaxTimeSteps() const
    {
      const int maxThreads = Fem::ThreadManager::maxThreads();
      double maxAdvStep = 0;
      double maxDiffStep = 0;
      for(int i=0; i<maxThreads; ++i ) 
      {
        maxAdvStep  = std::max( maxAdvStep,  problems_[ i ]->maxAdvectionTimeStep() );
        maxDiffStep = std::max( maxDiffStep, problems_[ i ]->maxDiffusionTimeStep() );
      }

      // set time steps to single problem 
      singleProblem_.setMaxTimeSteps( maxAdvStep, maxDiffStep );
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
    DiscreteModelType& singleProblem_;
    std::vector< DiscreteModelType* > problems_; 
    std::vector< InnerPassType* > passes_;
    mutable std::vector< double > passComputeTime_;
    mutable std::vector< bool   > passStage_;

    // temporary variables 
    mutable const ArgumentType* arg_; 
    mutable DestinationType* dest_; 

    // non-blocking communication handler 
    NonBlockingCommHandleType nonBlockingComm_;

    mutable size_t numberOfElements_;
    mutable bool firstCall_;
    const bool sumComputeTime_;
  };
//! @}  
} // end namespace Dune

#endif
