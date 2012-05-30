#ifndef DUNE_FEM_DG_RUNFILE_HH
#define DUNE_FEM_DG_RUNFILE_HH

#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/threadmanager.hh>
#include <dune/fem-dg/pass/threadpass.hh>

namespace Dune {

  template <class GridType>
  class RunFile : public AutoPersistentObject 
  {
    typedef typename GridType :: Traits :: CollectiveCommunication CommunicatorType;
    const CommunicatorType& comm_; 
    const std::string runFileName_;
    const int writeRunFile_; // 0 don't, 1 only speedup file, 2 write all runfiles 
                             // 3 only write 0, others at end, 4 all files at end 
    std::ostream* runfile_;

    std::vector< double > times_ ;
    double elements_;
    double maxDofs_;
    size_t timesteps_;

    // write in milli seconds
    inline size_t inMS(const double t) 
    {
      return (size_t (t * 1e3));
    }

    void writeHeader(std::ostream& runfile) 
    {
      // write header 
      runfile << "# Time          ";
      runfile << "   dt         ";
      runfile << "  Elements   ";
      runfile << "        dg   ";
      runfile << "       ode      ";
      runfile << "  adapt      ";
      runfile << "     lb   ";
      runfile << "       all  ";
      runfile << "      indi   ";
      runfile << "   limfunc   ";
      runfile << "     limit (in ms)  " << std::endl;
      runfile.flush();
    }

    std::string runFileName(const int rank) const 
    {
      std::stringstream runfile;
      runfile << Parameter :: commonOutputPath() << "/run." << rank; 
      return runfile.str();
    }

    std::ostream* createRunFile( const int rank, 
                                 const int writeId, 
                                 const bool newStart ) 
    {
      // in case of no writing or only speedup table don't create runfile
      if( writeId <= 1 ) return 0;

      bool writeAtOnce = ( writeId > 2 );
      // when writeId == 2 then only for rank 0 write file every time step 
      // this is for monitoring issues
      if( rank == 0 && writeId == 3 ) writeAtOnce = false ;

      if( writeAtOnce ) 
      {
        return new std::stringstream();
      }
      else 
      {
        std::ofstream* file = new std::ofstream( runFileName_.c_str(), ( newStart ) ? std::ios::out : std::ios::app );
        if( ! file ) 
        {
          std::cerr << "Couldn't open run file <"<<runFileName_<<">, ciao!" << std::endl;
          abort();
        }
        return file;
      }
    }
  public:  
    RunFile( const CommunicatorType& comm, const bool newStart )
      : comm_( comm )
      , runFileName_( runFileName( comm_.rank() ) )
      , writeRunFile_( Parameter :: getValue< int > ("fem.parallel.runfile", 0 ) )
      , runfile_( createRunFile( comm_.rank(), writeRunFile_, newStart ) ) 
      , times_() 
      , elements_( 0.0 )
      , maxDofs_( 0.0 )
      , timesteps_( 0 )
    {
      if( runfile_ && newStart ) 
      {
        writeHeader( *runfile_ );
      }
    }

    //! destructor 
    ~RunFile() 
    {
      delete runfile_;
    }


  protected:
    template <class T>
    void writeVectors(std::ostream& file, 
                      const std::string& descr, 
                      const std::vector< T >& sumTimes,
                      const std::vector< T >& maxTimes,
                      const std::vector< T >& minTimes ) const 
    {
      const size_t size = sumTimes.size();
      file << "#########################################" << std::endl ;
      file << "# Sum " << descr << std::endl ;
      for(size_t i=0; i<size-1; ++i)
      {
        file << sumTimes[ i ] << "  ";
      }
      file << std::endl;
      file << "# Max " << descr << std::endl ;
      for(size_t i=0; i<size-1; ++i)
      {
        file << maxTimes[ i ] << "  ";
      }
      file << std::endl;
      file << "# Min " << descr << std::endl ;
      for(size_t i=0; i<size-1; ++i)
      {
        file << minTimes[ i ] << "  ";
      }
      file << std::endl;
    }

  public:  
    void flush() 
    {
      // if write is > 0 then create speedup file 
      if( writeRunFile_ )
      {
        times_.push_back( elements_ );
        const size_t size = times_.size();
        std::vector< double > sumTimes ( times_ );
        std::vector< double > maxTimes ( times_ );
        std::vector< double > minTimes ( times_ );

        // sum, max, and min for all procs 
        comm_.sum( &sumTimes[ 0 ], size );
        comm_.min( &minTimes[ 0 ], size );

        maxTimes.push_back( maxDofs_ );
        comm_.max( &maxTimes[ 0 ], maxTimes.size() );

        maxDofs_ = maxTimes.back();
        maxTimes.pop_back();

        if( comm_.rank() == 0 && timesteps_ > 0 ) 
        {
          const int maxThreads = Fem :: ThreadManager :: maxThreads ();
          const double tasks = comm_.size() * maxThreads ;

          { // adjust elements to be the average element number  
            size_t i = size - 1 ;
            sumTimes[ i ] /= (double) timesteps_;
            maxTimes[ i ] /= (double) timesteps_;
            minTimes[ i ] /= (double) timesteps_;
          }

          std::stringstream runfile;
          runfile << Parameter :: commonOutputPath() << "/speedup." << comm_.size(); 
          std::ofstream file ( runfile.str().c_str() );
          if( file ) 
          {
            const double averageElements = sumTimes[ size - 1 ] / tasks ;

            const bool nonBlocking = 
#ifdef NSMOD_USE_SMP_PARALLEL
              NonBlockingCommHelper :: nonBlockingCommunication() ;
#else 
              false ;
#endif

            file << "# Procs = " << comm_.size() << " * " << maxThreads << " (MPI * threads)" << std::endl ;
            const char* commType = nonBlocking ? "asynchron" : "standard";
            file << "# Comm: " << commType << std::endl;
            file << "# Timesteps = " << timesteps_ << std::endl ;
            file << "# Max DoFs (per element): " << maxDofs_ << std::endl;
            file << "# Elements / timestep: sum    max    min    average  " << std::endl;
            file << sumTimes[ size-1 ] << "  " << maxTimes[ size-1 ] << "  " << minTimes[ size-1 ] << "  " << ((size_t)averageElements) << std::endl;
            file << "# DG       ODE     ADAPT    LB      TIMESTEP    " << std::endl ;

            // multiply sumTimes with maxThhreads since the sum would be to small otherwise 
            for(size_t i=0; i<size; ++i)
            {
              sumTimes[ i ] *= maxThreads ;
#if HAVE_BLUEGENE_P_ARCH
              // for some reason the time on bluegene is 
              // devided by number of threads 
              sumTimes[ i ] *= maxThreads ;
              maxTimes[ i ] *= maxThreads ;
              minTimes[ i ] *= maxThreads ;
#endif
            }
            { 
              std::vector<size_t> sumTimesElem(size);
              std::vector<size_t> maxTimesElem(size);
              std::vector<size_t> minTimesElem(size);

              for(size_t i=0; i<size; ++i)
              {
                sumTimesElem[ i ] = inMS( sumTimes[ i ] * averageElements );
                maxTimesElem[ i ] = inMS( maxTimes[ i ] * averageElements );
                minTimesElem[ i ] = inMS( minTimes[ i ] * averageElements );
              }
              {
                std::string descr("(time of all timesteps in ms)");
                writeVectors( file, descr, sumTimesElem, maxTimesElem, minTimesElem );
              }
              for(size_t i=0; i<size; ++i)
              {
                sumTimesElem[ i ] /= timesteps_;
                maxTimesElem[ i ] /= timesteps_;
                minTimesElem[ i ] /= timesteps_;
              }
              {
                std::string descr("(average time / timestep in ms)");
                writeVectors( file, descr, sumTimesElem, maxTimesElem, minTimesElem );
              }
            }

            // devide by timesteps 
            for(size_t i=0; i<size; ++i)
            {
              sumTimes[ i ] /= (double) timesteps_;
              maxTimes[ i ] /= (double) timesteps_;
              minTimes[ i ] /= (double) timesteps_;
            }

            {
              std::string descr( "( average time / timestep / element in sec )" );
              writeVectors( file, descr, sumTimes, maxTimes, minTimes );
            }
          }
        } // end speedup file 

        if( runfile_ ) 
        {
          std::stringstream* str = dynamic_cast< std::stringstream* > (runfile_); 
          if( str ) 
          {
            std::ofstream file( runFileName_.c_str() );

            if( ! file ) 
            {
              std::cerr << "Couldn't open run file <"<<runFileName_<<">, ciao!" << std::endl;
              abort();
            }

            file << str->str();
            file.flush();
            file.close();
          }
        }
      }
    }

    //! write timestep data 
    inline void write( const double t, 
                       const double ldt, 
                       const size_t nElements,
                       const size_t maxDofs,
                       const double dgOperatorTime,
                       const double odeSolve,
                       const double adaptTime,
                       const double lbTime,
                       const double timeStepTime,
                       const std::vector<double>& limitSteps = std::vector<double>() )
    {
      std::vector< double > times( 5 + limitSteps.size(), 0.0 );
      times[ 0 ] = dgOperatorTime ; 
      times[ 1 ] = odeSolve ;
      times[ 2 ] = adaptTime ;
      times[ 3 ] = lbTime ;
      times[ 4 ] = timeStepTime ;
      for(size_t i=5; i<limitSteps.size(); ++i)
        times[ i ] = limitSteps[ i ];

      maxDofs_ = std::max( double(maxDofs), maxDofs_ );

      write( t, ldt, nElements, times );
    }

    //! clone of write method 
    inline void write( const double t, 
                       const double ldt, 
                       const size_t nElements,
                       const std::vector<double>& times) 
    {
      if( writeRunFile_ ) 
      {
        const size_t size = times.size() ;
        const size_t oldsize = times_.size();
        if( oldsize < size  )
        {
          times_.resize( size ); 
          for( size_t i=oldsize; i<size; ++i) 
            times_[ i ] = 0;
        }

        const double elems = nElements ;
        elements_ += elems ;
        for(size_t i=0; i<size; ++i ) 
          times_[ i ] += times[ i ] / elems ; 

        ++timesteps_ ;

        if( runfile_ ) 
        {
          std::ostream& runfile = (*runfile_);
          const int space = 12;
          runfile << std::scientific << t  << "  ";
          runfile << std::setw(space) << ldt << "  ";
          runfile << std::setw(space) << nElements << " ";
          for(size_t i=0; i<size; ++i) 
            runfile << std::setw(space) << inMS( times[ i ] ) << " ";
          runfile << std::endl;

          runfile.flush();
        }
      }
    }

    //! backup routine 
    void backup() const 
    {
      typedef PersistenceManager :: BackupStreamType  BackupStreamType ;
      BackupStreamType& stream = PersistenceManager :: backupStream();

      stream << elements_ ;
      stream << maxDofs_ ;
      stream << timesteps_ ;
      const size_t tSize = times_.size();
      stream << tSize ;
      for( size_t i=0; i<tSize; ++i ) 
      {
        stream << times_[ i ];
      }
    }

    //! restore routine 
    void restore () 
    {
      typedef PersistenceManager :: RestoreStreamType  RestoreStreamType ;
      RestoreStreamType& stream = PersistenceManager :: restoreStream();

      stream >> elements_ ;
      stream >> maxDofs_ ;
      stream >> timesteps_ ;

      size_t tSize;
      stream >> tSize ;

      times_.resize( tSize );
      for( size_t i=0; i<tSize; ++i ) 
      {
        stream >> times_[ i ];
      }
    }
  }; // end class runfile

} // end namespace Dune 
#endif
