#ifndef DUNE_FEM_PTHREADCLASS_HH
#define DUNE_FEM_PTHREADCLASS_HH

#if defined USE_PTHREADS && ! HAVE_PTHREAD  
#error "pthreads not found, reconfigure!"
#endif

#include <cassert> 
#include <dune/fem/misc/threadmanager.hh>

namespace Dune { 

namespace Fem {

template <class Object> 
class ThreadHandle 
{
#ifdef USE_PTHREADS 
  ////////////////////////////////////////////
  // class ThreadHandleObject 
  ////////////////////////////////////////////
  class ThreadHandleObject
  {
    Object* objPtr_;
    pthread_barrier_t* barrier_ ;
    pthread_t threadId_ ;
    int threadNumber_ ;
    bool finished_ ;

    bool isMaster() const { return threadNumber_ == 0; }

  public:
    // constructor creating thread with given thread number 
    explicit ThreadHandleObject(Object& obj, const int threadNumber ) 
      : objPtr_( &obj ), 
        barrier_ ( 0 ),
        threadId_( threadNumber ),
        threadNumber_( threadNumber ),
        finished_( false )
    {
      assert( threadNumber > 0 );
    }

    // constructor creating master thread 
    explicit ThreadHandleObject( Object& obj ) 
      : objPtr_( &obj ),
        barrier_( 0 ),
        threadId_( pthread_self() ),
        threadNumber_( 0 ),
        finished_( false )
    {
    }

    // copy constructor 
    ThreadHandleObject(const ThreadHandleObject& other) 
      : objPtr_( other.objPtr_ ),
        barrier_( other.barrier_ ),
        threadId_( other.threadId_ ),
        threadNumber_( other.threadNumber_ ),
        finished_( other.finished_ ) 
    {}

    // assigment operator 
    ThreadHandleObject& operator = ( const ThreadHandleObject& other) 
    {
      objPtr_       = other.objPtr_ ;
      barrier_      = other.barrier_ ;
      threadId_     = other.threadId_;
      threadNumber_ = other.threadNumber_;
      finished_     = other.finished_;
      return *this;
    }

    // Create the thread and start work
    void start( pthread_barrier_t& barrier ) 
    {
      // init object 
      init( barrier );

      // on master thread there is no need to start an extra thread 
      if( isMaster() ) 
      {
        run();
      }
      else 
      {
        // create a joinable thread 
        pthread_create(&threadId_, 0, &ThreadHandleObject::startThread, (void *) this);
        // the master thread is also adding the threadnumber for the given ids 
        ThreadManager :: setThreadNumber( threadId_, threadNumber_ );
      }
    }

    // stop thread and join, return 1 when finished  
    int stoped() 
    {
      if( finished_ ) 
      {
        if( ! isMaster() ) 
          pthread_join(threadId_, 0);

        return 1;
      }
      return 0;
    }

    // do the work 
    void run() 
    {
      // wait for all threads 
      assert( barrier_ );
      pthread_barrier_wait( barrier_ );

      // run thread work 
      objPtr_->runThread();
      
      // work finished 
      finished_ = true ;
      barrier_ = 0;
    }

  private:
    void init( pthread_barrier_t& barrier ) 
    {
      barrier_ = & barrier ;
      finished_ = false ;
    }

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
#else // else ifdef USE_PTHREAD

  // object to call runThread 
  Object& obj_;
#endif // end ifdef USE_PTHREAD

public:
  // constructor 
  ThreadHandle( Object& obj ) 
#ifdef USE_PTHREADS
    : threads_()
  {
    const int maxThreads = ThreadManager :: maxThreads() ;
    for(int i=1; i<maxThreads; ++i)
    {
      // create thread handles for pthreads 
      threads_.push_back( ThreadHandleObject( obj, i ) );
    }

    // insert master thread at last because this thread creates 
    // all other threads before it start its calculations 
    threads_.push_back( ThreadHandleObject( obj ) );
#else 
    : obj_( obj )
  {
#endif
  } // end constructor 

  void run () 
  {
#ifdef USE_PTHREADS
    // get number of threads 
    const int maxThreads = ThreadManager :: maxThreads() ;
    
    // initialize barrier 
    pthread_barrier_t barrier ;
    pthread_barrier_init( &barrier , 0, maxThreads );

    // set number of active threads 
    ThreadManager :: initMultiThreadMode( maxThreads );

    // start threads, this will call the runThread method  
    for(int i=0; i<maxThreads; ++i)
    {
      threads_[ i ].start( barrier );
    }

    // wait until all threads are done 
    int count = 0; 
    while( count < maxThreads )
    {
      count = 0;
      // join threads 
      for(int i=0; i<maxThreads; ++i)
      {
        count += threads_[ i ].stoped() ;
      }
    }

    // activate initSingleThreadMode again 
    Fem :: ThreadManager :: initSingleThreadMode();

#else 
    // OpenMP parallel region 
#ifdef _OPENMP 
#pragma omp parallel
#endif
    {
      obj_.runThread();
    }
#endif
  }

};

} // end namespace Fem  

} // end namespace Dune 
#endif
