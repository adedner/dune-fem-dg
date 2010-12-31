#ifndef DUNE_FEM_PTHREADCLASS_HH
#define DUNE_FEM_PTHREADCLASS_HH

#if HAVE_PTHREAD 
#include <pthread.h>
#endif

#include <cassert> 
#include <dune/fem/misc/threadmanager.hh>

namespace Dune { 

namespace Fem {

template <class Object> 
class ThreadHandle 
{
#if HAVE_PTHREAD 
  ////////////////////////////////////////////
  // class ThreadHandleObject 
  ////////////////////////////////////////////
  class ThreadHandleObject
  {
    struct ObjectStorage
    {
      Object* objPtr_;
      pthread_barrier_t* barrier_ ;
      int threadNumber_ ;
      bool finished_ ;

      ObjectStorage(Object* obj, const int tn )
        : objPtr_( obj ), 
          barrier_( 0 ),
          threadNumber_( tn ),
          finished_( false ) 
      {
      }

      ObjectStorage( const ObjectStorage& other )
        : objPtr_( other.objPtr_ ),
          barrier_( other.barrier_ ),
          threadNumber_( other.threadNumber_ ),
          finished_( other.finished_ ) 
      {
      }

      ObjectStorage& operator = ( const ObjectStorage& other ) 
      {
        objPtr_       = other.objPtr_;
        barrier_      = other.barrier_ ;
        threadNumber_ = other.threadNumber_;
        finished_     = other.finished_;
        return *this;
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

      void init( pthread_barrier_t& barrier ) 
      {
        barrier_ = & barrier ;
        finished_ = false ;
      }

      bool finished() const { return finished_; }
      bool isMaster() const { return threadNumber_ == 0; }
      int threadNumber() const { return threadNumber_ ; }
    };

    ObjectStorage obj_;
    pthread_t threadId_ ;
  public:
    explicit ThreadHandleObject(Object& obj, const int threadNumber ) 
      : obj_( &obj, threadNumber ),
        threadId_( threadNumber )
    {
      assert( threadNumber > 0 );
    }

    explicit ThreadHandleObject( Object& obj ) 
      : obj_( &obj, 0 ),
        threadId_( pthread_self() )
    {
    }

    ThreadHandleObject(const ThreadHandleObject& other) 
      : obj_( other.obj_ ),
        threadId_( other.threadId_ )
    {}

    ThreadHandleObject& operator = ( const ThreadHandleObject& other) 
    {
      obj_ = other.obj_ ;
      threadId_ = other.threadId_;
      return *this;
    }

    // Create the thread and start work
    void start( pthread_barrier_t& barrier ) 
    {
      // init object 
      obj_.init( barrier );

      // on master thread there is no need to start an extra thread 
      if( obj_.isMaster() ) 
      {
        obj_.run();
      }
      else 
      {
        // create a joinable thread 
        pthread_create(&threadId_, 0, &ThreadHandleObject::startThread, (void *) &obj_);
        // the master thread is also adding the threadnumber for the given ids 
        ThreadManager :: setThreadNumber( threadId_, obj_.threadNumber() );
      }
    }

    // stop thread and join, return 1 when finished  
    int stoped() 
    {
      if( obj_.finished() ) 
      {
        if( ! obj_.isMaster() ) 
          pthread_join(threadId_, 0);

        return 1;
      }
      return 0;
    }

  private:
    // This is the static class function that serves as a 
    // C style function pointer for the pthread_create call
    static void* startThread(void *obj)
    {
      // do the work
      ((ObjectStorage *) obj)->run();
      return 0;
    }
  }; // end ThreadHandleObject 
  ////////////////////////////////////////////////////
  //  end ThreadHandleObject
  ////////////////////////////////////////////////////

  std::vector< ThreadHandleObject > threads_;
#endif // end if HAVE_PTHREAD

public:
  // constructor 
  ThreadHandle( Object& obj ) 
#if HAVE_PTHREAD 
    : threads_()
  {
    const int maxThreads = ThreadManager :: maxThreads() ;
    for(int i=1; i<maxThreads; ++i)
    {
      // create thread handles for pthreads 
      threads_.push_back( ThreadHandleObject( obj, i ) );
    }

    // insert master thread at last 
    threads_.push_back( ThreadHandleObject( obj ) );
#else 
  {
#endif
  }

  void run () 
  {
#if HAVE_PTHREAD 
    // get number of threads 
    const int maxThreads = ThreadManager :: maxThreads() ;
    
    // initialize barrier 
    pthread_barrier_t barrier ;
    pthread_barrier_init( &barrier , NULL, maxThreads );

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
    {
      std::cerr << "ERROR: pthreads not available! " << std::endl;
      abort();
    }
#endif
  }

};

} // end namespace Fem  

} // end namespace Dune 
#endif
