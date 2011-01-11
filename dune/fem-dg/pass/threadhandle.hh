#ifndef DUNE_FEM_PTHREADCLASS_HH
#define DUNE_FEM_PTHREADCLASS_HH

#if defined USE_PTHREADS && ! HAVE_PTHREAD  
#error "pthreads not found, reconfigure!"
#endif

#include <cassert> 
#include <dune/fem/misc/threadmanager.hh>

namespace Dune { 

namespace Fem {

class ThreadHandle 
{
#ifdef USE_PTHREADS 
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
  public:  
    ObjectWrapper( Object& obj ) : obj_( obj ) {}
    void run () { obj_.runThread(); }
  };

  ////////////////////////////////////////////
  // class ThreadHandleObject 
  ////////////////////////////////////////////
  class ThreadHandleObject
  {
    ObjectIF* objPtr_;
    pthread_barrier_t* barrier_ ;
    pthread_t threadId_ ;
    int threadNumber_ ;

    bool isMaster() const { return threadNumber_ == 0; }

  public:
    // constructor creating thread with given thread number 
    explicit ThreadHandleObject(pthread_barrier_t* barrier,
                                const int threadNumber ) 
      : objPtr_( 0 ), 
        barrier_ ( barrier ),
        threadId_( 0 ),
        threadNumber_( threadNumber )
    {
      assert( threadNumber > 0 );
    }

    // constructor creating master thread 
    explicit ThreadHandleObject( pthread_barrier_t* barrier )
      : objPtr_( 0 ),
        barrier_( barrier ),
        threadId_( pthread_self() ),
        threadNumber_( 0 )
    {
    }

    // copy constructor 
    ThreadHandleObject(const ThreadHandleObject& other) 
      : objPtr_( other.objPtr_ ),
        barrier_( other.barrier_ ),
        threadId_( other.threadId_ ),
        threadNumber_( other.threadNumber_ )
    {}

    // assigment operator 
    ThreadHandleObject& operator = ( const ThreadHandleObject& other) 
    {
      objPtr_       = other.objPtr_ ;
      barrier_      = other.barrier_ ;
      threadId_     = other.threadId_;
      threadNumber_ = other.threadNumber_;
      return *this;
    }

    // Create the thread and start work
    void start( ObjectIF& obj ) 
    {
      // init object 
      objPtr_ = & obj ;

      // on master thread there is no need to start an extra thread 
      if( isMaster() ) 
      {
        run();
      }
      else 
      {
        // if thread has not been initialized 
        if( threadId_ == 0 )
        {
          // create a joinable thread 
          pthread_create(&threadId_, 0, &ThreadHandleObject::startThread, (void *) this);
          // the master thread is also adding the threadnumber for the given ids 
          ThreadManager :: setThreadNumber( threadId_, threadNumber_ );
        }
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
      assert( barrier_ );
      // wait for all threads 
      pthread_barrier_wait( barrier_ );

      // when object pointer is set call run 
      //if( objPtr_ ) 
      {
        assert( objPtr_ );
        objPtr_->run();
      }
      
      // work finished, set objPtr to zero 
      objPtr_ = 0 ;

      // when thread is not master then 
      // just call run and wait at barrier 
      if( ! isMaster() ) 
      {
        run();
      }
    }

    /*
    //! destroy thread by calling pthread_join  
    void finish() 
    {
      // work finished 
      objPtr_ = 0;
      return ;

      // run on master to reach barrier 
      if( isMaster() ) 
      {
        run();
      }
    }

    //! destroy thread by calling pthread_join  
    void destroy() 
    {
      //if( ! isMaster() )
      //  pthread_join(threadId_, 0);
    }
    */
  private:
    // This is the static class function that serves as a 
    // C style function pointer for the pthread_create call
    static void* startThread(void *obj)
    {
      // do the work
      ((ThreadHandleObject *) obj)->run();
      return 0;
    }
  }; // end ThreadHandleObject 
  ////////////////////////////////////////////////////
  //  end ThreadHandleObject
  ////////////////////////////////////////////////////

  std::vector< ThreadHandleObject > threads_;
  pthread_barrier_t waitAll_ ;
  const int maxThreads_ ;

private:  
  // prohibit copying 
  ThreadHandle( const ThreadHandle& );
  // default constructor 
  ThreadHandle() 
    : threads_()
    , waitAll_()
    , maxThreads_( ThreadManager :: maxThreads() )
  {
    // initialize barrier 
    pthread_barrier_init( &waitAll_, 0, maxThreads_ );

    // initialize slave threads 
    for(int i=1; i<maxThreads_; ++i)
    {
      // create thread handles for pthreads 
      threads_.push_back( ThreadHandleObject( &waitAll_, i ) );
    }

    // insert master thread at last because this thread creates 
    // all other threads before it start its calculations 
    threads_.push_back( ThreadHandleObject( &waitAll_ ) );
  } // end constructor 

  /*
  // destructor deleting threads 
  ~ThreadHandle() 
  {
    // set number of active threads 
    ThreadManager :: initMultiThreadMode( maxThreads_ );

    for(int i=0; i<maxThreads_; ++i)
    {
      threads_[ i ].finish();
    }

    // wait for all to finish  
    wait();

    // activate initSingleThreadMode again 
    Fem :: ThreadManager :: initSingleThreadMode();

    // call thread join 
    for(int i=0; i<maxThreads_; ++i)
    {
      threads_[ i ].destroy();
    }

    // destroy barrier 
    pthread_barrier_destroy( &waitAll_ );
  }
  */

  //! wait until all threads are stoped 
  void wait() const 
  {
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
  }

  template <class Object>
  void runThreads( Object& obj ) 
  {
    // create object wrapper 
    ObjectWrapper< Object > objPtr( obj );

    // set number of active threads 
    ThreadManager :: initMultiThreadMode( maxThreads_ );

    // start threads, this will call the runThread method  
    for(int i=0; i<maxThreads_; ++i)
    {
      threads_[ i ].start( objPtr );
    }

    // wait until all threads are done 
    wait();

    // activate initSingleThreadMode again 
    Fem :: ThreadManager :: initSingleThreadMode();
  }

  // return instance of ThreadHandle
  static ThreadHandle& instance() 
  {
    static ThreadHandle handle;
    return handle;
  }

#endif

public:  
  template <class Object> 
  static void run ( Object& obj ) 
  {
#ifdef USE_PTHREADS
    {
      // pthread version 
      instance().runThreads( obj );
    }
#else 
    // OpenMP parallel region 
#ifdef _OPENMP 
#pragma omp parallel
#endif
    {
      obj.runThread();
    }
#endif
  }

};

} // end namespace Fem  

} // end namespace Dune 
#endif
