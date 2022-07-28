#ifndef DUNE_FEM_DG_DIAGNOSTICS_HH
#define DUNE_FEM_DG_DIAGNOSTICS_HH

#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem-dg/pass/threadpass.hh>

namespace Dune
{
namespace Fem
{

  class Diagnostics
    : public Fem::AutoPersistentObject
  {
    typedef Fem :: MPIManager :: Communication CommunicatorType;
    const CommunicatorType& comm_;
    const std::string runFileName_;
    const int writeDiagnostics_; // 0 don't, 1 only speedup file, 2 write all diagnosticss
                             // 3 only write 0, others at end, 4 all files at end
    std::ostream* diagnostics_;

    std::vector< double > times_ ;
    std::vector< double > timesPerElem_ ;
    double elements_;
    double maxDofs_;
    size_t timesteps_;

    static const size_t width = 12;

    // write in milli seconds
    inline size_t inMS(const double t) const
    {
      return (size_t (t * 1e3));
    }

    void writeHeader(std::ostream& diagnostics)
    {
      // write header
      diagnostics << "# Time          ";
      diagnostics << "   dt         ";
      diagnostics << "  Elements   ";
      diagnostics << "        dg   ";
      diagnostics << "       ode      ";
      diagnostics << "  adapt      ";
      diagnostics << "     lb   ";
      diagnostics << "       all  ";
      diagnostics << "      indi   ";
      diagnostics << "   limfunc   ";
      diagnostics << "     limit (in ms)  " << std::endl;
      diagnostics.flush();
    }

    std::string runFileName(const int rank, const std::string keyPrefix = "" ) const
    {
      std::stringstream diagnostics;
      diagnostics << Fem :: Parameter :: commonOutputPath() << "/run" << keyPrefix << "." << rank;
      return diagnostics.str();
    }

    std::ostream* createDiagnostics( const int rank,
                                     const int writeId,
                                     const bool newStart )
    {
      // in case of no writing or only speedup table don't create diagnostics
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
    // no copy and assignment
    Diagnostics( const Diagnostics& ) = delete;
    Diagnostics& operator=( const Diagnostics& ) = delete;

    explicit Diagnostics( const bool newStart, const std::string keyPrefix = "" )
      : comm_( Fem :: MPIManager :: comm() )
      , runFileName_( runFileName( comm_.rank(), keyPrefix ) )
      , writeDiagnostics_( Fem :: Parameter :: getValue< int > ("fem.parallel.diagnostics", 0 ) )
      , diagnostics_( createDiagnostics( comm_.rank(), writeDiagnostics_, newStart ) )
      , times_()
      , timesPerElem_()
      , elements_( 0.0 )
      , maxDofs_( 0.0 )
      , timesteps_( 0 )
    {
      if( diagnostics_ && newStart )
      {
        writeHeader( *diagnostics_ );
      }
    }

    //! destructor
    ~Diagnostics()
    {
      delete diagnostics_;
    }


  protected:
    template <class T>
    void writeVectors(std::ostream& file,
                      const std::string& descr,
                      const std::vector< T >& sumTimes,
                      const std::vector< T >& maxTimes,
                      const std::vector< T >& minTimes,
                      const std::string& sumDescr,
                      const std::string& maxDescr ) const
    {
      const size_t size = sumTimes.size();
      file.precision(6);
      file << std::scientific ;
      file << "########################################################################################" << std::endl ;
      file << "# Sum " << descr << sumDescr << std::endl;
      for(size_t i=0; i<size; ++i)
      {
        file << std::setw(width) << sumTimes[ i ] << "   ";
      }
      file << std::endl;
      file << "# Max " << descr << maxDescr << std::endl;
      for(size_t i=0; i<size; ++i)
      {
        file << std::setw(width) << maxTimes[ i ] << "   ";
      }
      file << std::endl;
      file << "# Min " << descr << maxDescr << std::endl;
      for(size_t i=0; i<size; ++i)
      {
        file << std::setw(width) << minTimes[ i ] << "   ";
      }
      file << std::endl;
    }

  public:
    void flush() const
    {
      // if write is > 0 then create speedup file
      if( writeDiagnostics_ )
      {
        const size_t size = times_.size();

        std::vector< double > sumTimes ( times_ );
        std::vector< double > maxTimes ( times_ );
        std::vector< double > minTimes ( times_ );

        sumTimes.reserve( 2*size+1 );
        maxTimes.reserve( 2*size+1 );
        minTimes.reserve( 2*size+1 );

        for( size_t i=0; i<size; ++ i )
        {
          sumTimes.push_back( timesPerElem_[i] );
          maxTimes.push_back( timesPerElem_[i] );
          minTimes.push_back( timesPerElem_[i] );
        }

        // add number of elements
        sumTimes.push_back( elements_ );
        maxTimes.push_back( elements_ );
        minTimes.push_back( elements_ );

        // sum, max, and min for all procs
        comm_.sum( &sumTimes[ 0 ], sumTimes.size() );
        comm_.min( &minTimes[ 0 ], minTimes.size() );

        maxTimes.push_back( maxDofs_ );
        comm_.max( &maxTimes[ 0 ], maxTimes.size() );

        const double maxDofs = maxTimes.back();
        maxTimes.pop_back();

        if( comm_.rank() == 0 && timesteps_ > 0 )
        {
          const int maxThreads = Fem :: MPIManager :: maxThreads ();
          const double timesteps = double(timesteps_);

          const size_t bufferSize = sumTimes.size();
          const double sumElements = sumTimes[ bufferSize - 1 ] / timesteps;
          const double maxElements = maxTimes[ bufferSize - 1 ] / timesteps;
          const double minElements = minTimes[ bufferSize - 1 ] / timesteps;

          std::vector< double > sumTimesPerElem( size, 0.0 );
          std::vector< double > maxTimesPerElem( size, 0.0 );
          std::vector< double > minTimesPerElem( size, 0.0 );

          for( size_t i=0; i<size; ++i )
          {
            sumTimesPerElem[ i ] = sumTimes[ i + size ];
            maxTimesPerElem[ i ] = maxTimes[ i + size ];
            minTimesPerElem[ i ] = minTimes[ i + size ];
          }
          sumTimes.resize( size );
          maxTimes.resize( size );
          minTimes.resize( size );

          std::stringstream diagnostics;
          diagnostics << Fem :: Parameter :: commonOutputPath() << "/speedup." << comm_.size();
          std::ofstream file ( diagnostics.str().c_str() );
          if( file )
          {
            const double tasks = comm_.size();// * maxThreads ;
            const double averageElements = sumElements / tasks ;

            // get information about communication type
            const bool nonBlocking = NonBlockingCommParameter :: nonBlockingCommunication() ;

            file << "# Procs = " << comm_.size() << " * " << maxThreads << " (MPI * threads)" << std::endl ;
            const char* commType = nonBlocking ? "asynchronous" : "standard";
            file << "# Comm: " << commType << std::endl;
            file << "# Timesteps = " << timesteps_ << std::endl ;
            file << "# Max DoFs (per element): " << maxDofs << std::endl;
            file << "# Elements / timestep / MPI tasks:" << std::endl;
            file << "#" << std::setw(width-1) << "sum";
            file << std::setw(width) << "max";
            file << std::setw(width) << "min";
            file << std::setw(width) << "average" << std::endl;
            file << std::setw(width) << sumElements;
            file << std::setw(width) << maxElements;
            file << std::setw(width) << minElements;
            file << std::setw(width) << averageElements << std::endl;
            file << "#########################################################################################" << std::endl;
            file << "#";
            file << std::setw(width) << "DG    ";
            file << std::setw(width+2) << "ODE   ";
            file << std::setw(width+1) << "ADAPT ";
            file << std::setw(width+1) << "   LB ";
            file << std::setw(width+5) << "TIMESTEP" << std::endl;

            // multiply sumTimes with maxThhreads since the sum would be to small otherwise
            for(size_t i=0; i<size; ++i)
            {
              sumTimes[ i ] *= maxThreads ;
              sumTimesPerElem[ i ] *= maxThreads;
            }
            {
              std::string sumDescr("(opt: stays constant over #core increase)");
              std::string maxDescr("(opt: inversely proportional to #core increase)");
              {
                std::string descr("(time of all timesteps in sec)");
                writeVectors( file, descr, sumTimes, maxTimes, minTimes, sumDescr, maxDescr );
              }
              for(size_t i=0; i<size; ++i)
              {
                sumTimes[ i ] /= timesteps;
                maxTimes[ i ] /= timesteps;
                minTimes[ i ] /= timesteps;
              }
              {
                std::string descr("(time / timesteps in sec)");
                writeVectors( file, descr, sumTimes, maxTimes, minTimes, sumDescr, maxDescr );
              }
            }

            std::string sumDescrPerElem("(opt: grows with #core)");
            std::string maxDescrPerElem("(opt: stays constant over #core increase)");
            {
              std::string descr( "(time / elements in sec )" );
              writeVectors( file, descr, sumTimesPerElem, maxTimesPerElem, minTimesPerElem, sumDescrPerElem, maxDescrPerElem );
            }

            // devide per elem times by timesteps
            for(size_t i=0; i<size; ++i)
            {
              sumTimesPerElem[ i ] /= timesteps;
              maxTimesPerElem[ i ] /= timesteps;
              minTimesPerElem[ i ] /= timesteps;
            }

            {
              std::string descr( "(time / (timesteps * elements) in sec )" );
              writeVectors( file, descr, sumTimesPerElem, maxTimesPerElem, minTimesPerElem, sumDescrPerElem, maxDescrPerElem );
            }
          }
        } // end speedup file

        if( diagnostics_ )
        {
          std::stringstream* str = dynamic_cast< std::stringstream* > (diagnostics_);
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
      times[ 0 ] = dgOperatorTime;
      times[ 1 ] = odeSolve;
      times[ 2 ] = adaptTime;
      times[ 3 ] = lbTime;
      times[ 4 ] = timeStepTime ;
      for(size_t i=0; i<limitSteps.size(); ++i)
        times[ i+5 ] = limitSteps[ i ];

      maxDofs_ = std::max( double(maxDofs), maxDofs_ );

      write( t, ldt, nElements, times );
    }

    //! clone of write method
    inline void write( const double t,
                       const double ldt,
                       const size_t nElements,
                       const std::vector<double>& times )
    {
      if( writeDiagnostics_ )
      {
        const size_t size = times.size() ;
        const size_t oldsize = times_.size();
        if( oldsize < size  )
        {
          times_.resize( size );
          timesPerElem_.resize( size );
          for( size_t i=oldsize; i<size; ++i)
          {
            times_[ i ] = 0;
            timesPerElem_[ i ] = 0;
          }
        }

        elements_ += double( nElements );
        for(size_t i=0; i<size; ++i )
        {
          times_[ i ]        += times[ i ] ;
          timesPerElem_[ i ] += times[ i ] / double( nElements );
        }

        ++timesteps_ ;

        if( diagnostics_ )
        {
          std::ostream& diagnostics = (*diagnostics_);
          const int space = 12;
          diagnostics << std::scientific << t  << "  ";
          diagnostics << std::setw(space) << ldt << "  ";
          diagnostics << std::setw(space) << nElements << " ";
          for(size_t i=0; i<size; ++i)
            diagnostics << std::setw(space) << inMS( times[ i ] ) << " ";
          diagnostics << std::endl;

          diagnostics.flush();
        }
      }
    }

    //! backup routine
    void backup() const
    {
      // flush run file
      flush();

      typedef Fem :: PersistenceManager :: BackupStreamType  BackupStreamType ;
      BackupStreamType& stream = Fem :: PersistenceManager :: backupStream();

      stream << elements_ ;
      stream << maxDofs_ ;
      stream << timesteps_ ;
      const size_t tSize = times_.size();
      stream << tSize ;
      for( size_t i=0; i<tSize; ++i )
      {
        stream << times_[ i ];
      }
      for( size_t i=0; i<tSize; ++i )
      {
        stream << timesPerElem_[ i ];
      }
    }

    //! restore routine
    void restore ()
    {
      typedef Fem :: PersistenceManager :: RestoreStreamType  RestoreStreamType ;
      RestoreStreamType& stream = Fem :: PersistenceManager :: restoreStream();

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
      timesPerElem_.resize( tSize );
      for( size_t i=0; i<tSize; ++i )
      {
        stream >> timesPerElem_[ i ];
      }
    }
  }; // end class diagnostics

} // end namespace
} // end namespace Dune
#endif
