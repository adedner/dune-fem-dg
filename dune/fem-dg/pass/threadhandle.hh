#ifndef DUNE_FEM_PTHREADCLASS_HH
#define DUNE_FEM_PTHREADCLASS_HH

#include <cassert>
#include <mutex>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

namespace Dune
{
namespace Fem
{

  class ThreadHandle
  {
    class ObjectIF
    {
    protected:
      ObjectIF() {}
    public:
      virtual ~ObjectIF() {}
      virtual void run() = 0;
    };

    template <class Object>
    class ObjectWrapper : public ObjectIF
    {
      Object& obj_;
      std::mutex* mutex_;
    public:
      ObjectWrapper( Object& obj, std::mutex* mtx = nullptr )
        : obj_( obj ), mutex_( mtx ) {}

      void run ()
      {
        // if mutex available run in locked mode
        if( mutex_ )
          mutex_->lock();

        obj_.runThread();

        // unlock if locked before
        if( mutex_ )
          mutex_->unlock();
      }
    };

#ifdef USE_PTHREADS
    ////////////////////////////////////////////
    // class ThreadHandleObject
    ////////////////////////////////////////////
    class ThreadHandleObject
    {
      ObjectIF* objPtr_;
      pthread_barrier_t* barrierBegin_ ;
      pthread_barrier_t* barrierEnd_ ;
      pthread_t threadId_ ;
      int maxThreads_ ;
      int threadNumber_ ;

      bool isSlave () const { return threadNumber_ > 0; }

    public:
      // constructor creating thread with given thread number
      ThreadHandleObject(pthread_barrier_t* barrierBegin,
                         pthread_barrier_t* barrierEnd,
                         const int maxThreads,
                         const int threadNumber )
        : objPtr_( 0 ),
          barrierBegin_ ( barrierBegin ),
          barrierEnd_ ( barrierEnd ),
          threadId_( 0 ),
          maxThreads_( maxThreads ),
          threadNumber_( threadNumber )
      {
        assert( threadNumber > 0 );
      }

      // constructor creating master thread
      explicit ThreadHandleObject(pthread_barrier_t* barrierBegin,
                                  pthread_barrier_t* barrierEnd,
                                  const int maxThreads)
        : objPtr_( 0 ),
          barrierBegin_ ( barrierBegin ),
          barrierEnd_ ( barrierEnd ),
          threadId_( pthread_self() ),
          maxThreads_( maxThreads ),
          threadNumber_( 0 )
      {
      }

      // copy constructor
      ThreadHandleObject(const ThreadHandleObject& other)
        : objPtr_( other.objPtr_ ),
          barrierBegin_( other.barrierBegin_ ),
          barrierEnd_( other.barrierEnd_ ),
          threadId_( other.threadId_ ),
          maxThreads_( other.maxThreads_ ),
          threadNumber_( other.threadNumber_ )
      {}

      // assigment operator
      ThreadHandleObject& operator = ( const ThreadHandleObject& other)
      {
        objPtr_       = other.objPtr_ ;
        barrierBegin_ = other.barrierBegin_ ;
        barrierEnd_   = other.barrierEnd_ ;
        threadId_     = other.threadId_;
        maxThreads_   = other.maxThreads_;
        threadNumber_ = other.threadNumber_;
        return *this;
      }

      // Create the thread and start work
      void start( ObjectIF* obj )
      {
        // init object
        objPtr_ = obj;

        if( isSlave() )
        {
          // if thread has not been initialized
          if( threadId_ == 0 )
          {
            // create a joinable thread
            pthread_create(&threadId_, 0, &ThreadHandleObject::startThread, (void *) this);
          }
        }
        else
        {
          // on master thread there is no need to start an extra thread
          run();
        }
      }

      //! return 1 of thread is stoped, 0 otherwise
      int stoped() const
      {
        return ( objPtr_ == 0 ) ? 1 : 0 ;
      }

      // do the work
      void run()
      {
        assert( barrierBegin_ );
        assert( barrierEnd_ );

        // wait for all threads
        pthread_barrier_wait( barrierBegin_ );

        // when object pointer is set call run, else terminate
        if( objPtr_ )
        {
          objPtr_->run();
        }
        else
        {
          // this is terminating the threads
          return ;
        }

        // work finished, set objPtr to zero
        objPtr_ = 0 ;

        // wait for all threads
        pthread_barrier_wait( barrierEnd_ );

        // when thread is not master then
        // just call run and wait at barrier
        if( isSlave() )
        {
          run();
        }
      }

      //! destroy thread by calling pthread_join
      void destroy()
      {
        if( isSlave() )
          pthread_join(threadId_, 0);
      }

    private:
      // This is the static class function that serves as a
      // C style function pointer for the pthread_create call
      static void* startThread(void *obj)
      {
        // set maxThreads and threadNumber for slave thread
        ThreadManager :: initThread( ((ThreadHandleObject *) obj)->maxThreads_, ((ThreadHandleObject *) obj)->threadNumber_ );

        // do the work
        ((ThreadHandleObject *) obj)->run();

        return 0;
      }
    }; // end ThreadHandleObject
    ////////////////////////////////////////////////////
    //  end ThreadHandleObject
    ////////////////////////////////////////////////////

    std::vector< ThreadHandleObject > threads_;
    pthread_barrier_t waitBegin_ ;
    pthread_barrier_t waitEnd_ ;
    const int maxThreads_ ;

  private:
    // prohibit copying
    ThreadHandle( const ThreadHandle& );
    // default constructor
    ThreadHandle()
      : threads_()
      , waitBegin_()
      , waitEnd_()
      , maxThreads_( ThreadManager :: maxThreads() )
    {
      // initialize barrier
      pthread_barrier_init( &waitBegin_, 0, maxThreads_ );

      // initialize barrier
      pthread_barrier_init( &waitEnd_, 0, maxThreads_ );

      // initialize slave threads
      for(int i=1; i<maxThreads_; ++i)
      {
        // create thread handles for pthreads
        threads_.push_back( ThreadHandleObject( &waitBegin_, &waitEnd_, maxThreads_, i ) );
      }

      // insert master thread at last because this thread creates
      // all other threads before it start its calculations
      threads_.push_back( ThreadHandleObject( &waitBegin_, &waitEnd_, maxThreads_ ) );
    } // end constructor

    //! start all threads to do the job
    void startThreads( ObjectIF* obj = 0 )
    {
      // set number of active threads
      ThreadManager :: initMultiThreadMode( maxThreads_ );

      // start threads, this will call the runThread method
      for(int i=0; i<maxThreads_; ++i)
      {
        threads_[ i ].start( obj );
      }

      // wait until all threads are done
      int count = 0;
      while( count < maxThreads_ )
      {
        count = 0;
        // join threads
        for(int i=0; i<maxThreads_; ++i)
        {
          count += threads_[ i ].stoped() ;
        }
      }

      // activate initSingleThreadMode again
      Fem :: ThreadManager :: initSingleThreadMode();
    }

    //! run all threads
    template <class Object>
    void runThreads( Object& obj, std::mutex* mtx = nullptr )
    {
      // create object wrapper
      ObjectWrapper< Object > objPtr( obj, mtx );

      // start parallel execution
      startThreads( & objPtr ) ;
    }

    // return instance of ThreadHandle
    static ThreadHandle& instance()
    {
      static std::unique_ptr< ThreadHandle > handle( new ThreadHandle() );
      return *handle;
    }

  public:
    //! destructor deleting threads
    ~ThreadHandle()
    {
      // start threads with null object which will terminate each sub thread
      startThreads () ;

      // call thread join
      for(int i=0; i<maxThreads_; ++i)
      {
        threads_[ i ].destroy();
      }

      // destroy barrier
      pthread_barrier_destroy( &waitEnd_ );
      // destroy barrier
      pthread_barrier_destroy( &waitBegin_ );
    }

#endif // end HAVE_PTHREAD

  public:
    template <class Object>
    static void run ( Object& obj )
    {
      // this routine should not be called in multiThreadMode, since
      // this routine is actually starting the multiThreadMode
      if( ! ThreadManager :: singleThreadMode() )
        DUNE_THROW(InvalidStateException,"ThreadHandle :: run called from thread parallel region!");

#ifdef USE_PTHREADS
      if( ThreadManager :: pthreads )
      {
        // pthread version
        instance().runThreads( obj );
      }
      else
#endif
      {
        // OpenMP parallel region
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          obj.runThread();
        }
      }
    }

    template <class Object>
    static void runLocked ( Object& obj )
    {
      // this routine should not be called in multiThreadMode, since
      // this routine is actually starting the multiThreadMode
      if( ! ThreadManager :: singleThreadMode() )
        DUNE_THROW(InvalidStateException,"ThreadHandle :: run called from thread parallel region!");

      // run threads in blocking mode
      std::mutex mtx;

#ifdef USE_PTHREADS
      if( ThreadManager :: pthreads )
      {
        // pthread version
        instance().runThreads( obj, &mtx );
      }
      else
#endif
      {
        // use mutex to lock call
        typename ThreadHandle::template ObjectWrapper< Object > objCaller( obj, &mtx );

        // OpenMP parallel region
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          objCaller.run();
        }
      }
    }

  };

} // end namespace Fem

} // end namespace Dune
#endif
