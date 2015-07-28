#ifndef FEMHOWTO_BASE_HH
#define FEMHOWTO_BASE_HH

#include <sstream>
#include <sys/time.h>
#include <sys/resource.h>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/common/timer.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>

#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>

typedef Dune::Fem::Parameter  ParameterType ;

struct EocDataOutputParameters :
       public Dune::Fem::LocalParameter<Dune::Fem::DataWriterParameters,EocDataOutputParameters>
{
  protected:
  std::string loop_;
  public:
  EocDataOutputParameters(int loop, const std::string& name) {
    std::stringstream ss;
    ss << name << loop;
    loop_ = ss.str();
  }
  EocDataOutputParameters(const EocDataOutputParameters& other)
  : loop_(other.loop_) {}

  std::string path() const {
    return loop_;
  }
};



struct EocParameters :
       public Dune::Fem::LocalParameter<EocParameters,EocParameters>
{
  protected:
  std::string keyPrefix_;

  public:
  EocParameters( const std::string keyPrefix = "fem.eoc.")
    : keyPrefix_( keyPrefix )
  {}

  const EocDataOutputParameters dataOutputParameters( int loop, const std::string& name ) const
  {
    return EocDataOutputParameters( loop, name );
  }

  virtual int steps() const
  {
    checkOldParameterUsed( "femhowto.eocSteps" );
    return Dune::Fem::Parameter::getValue<int>( keyPrefix_ + "steps", 1);
  }

  virtual std::string outputPath() const
  {
    checkOldParameterUsed( "femhowto.eocOutputPath" );
    return Dune::Fem::Parameter::getValue<std::string>( keyPrefix_ + "outputpath", "./");
  }

  virtual std::string fileName() const
  {
    checkOldParameterUsed( "femhowto.eocFileName" );
    return Dune::Fem::Parameter::getValue<std::string>( keyPrefix_ + "filename", "eoc" );
  }

  virtual int quadOrder() const
  {
    return Dune::Fem::Parameter::getValue< int >( keyPrefix_ + "quadorder", -1 );
  }

  private:

  //helper function for deprecation warnings
  void checkOldParameterUsed( const std::string old ) const
  {
    if( Dune::Fem::Parameter::exists( old ) )
    {
      std::cerr << "WARNING: deprecated key, please update your parameter '" << old << "'!" << std::endl;
      //enforce parameter update by throwing an abort();
      abort();
    }
  }

};


static inline std::string checkPointRestartFileName ()
{
  static bool initialized = false ;
  static std::string checkFileName ;
  if( ! initialized )
  {
    const std::string key( "fem.io.checkpointrestartfile" );
    if( ParameterType :: exists( key ) )
    {
      checkFileName = ParameterType::getValue<std::string> ( key );
    }
    else
      checkFileName = "" ;
    initialized = true ;
  }

  return checkFileName ;
}

template< class HGridType >
Dune::GridPtr< HGridType > initialize( const std::string& problemDescription, const std::string name = "" )
{
#ifdef ALUGRID_CONSTRUCTION_WITH_STREAMS
  // precision of out streams (here ALUGrid backup streams)
  const int precision = ParameterType :: getValue< int > ("fem.io.precision", 16);
  ALU3DSPACE ALUGridExternalParameters :: precision () = precision;
#endif

  // grid pointer object
  Dune::GridPtr< HGridType > gridptr ;

  const std::string checkPointRestartFile = checkPointRestartFileName();

  // if checkpoint restart file has been given
  if( checkPointRestartFile.size() > 0 )
  {
    if( 0 == Dune::Fem::MPIManager :: rank () )
    {
      std::cout << std::endl;
      std::cout << "********************************************************" << std::endl;
      std::cout << "**   Restart from checkpoint `" << checkPointRestartFile << "' " << std::endl;
      std::cout << "********************************************************" << std::endl;
      std::cout << std::endl;
    }

    // restore grid from checkpoint
    gridptr = Dune::Fem::CheckPointer< HGridType > :: restoreGrid( checkPointRestartFile, -1,
                                                                   Dune::Fem::CheckPointerParameters( Dune::ParameterKey::generate( "", "fem.io." ) ) );
  }
  else  // normal new start
  {
    // ----- read in runtime parameters ------
    const std::string filekey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
    const std::string filename = ParameterType::commonInputPath() + "/" + ParameterType::getValue< std::string >(  filekey );

    // initialize grid with given macro file name
    gridptr = Dune::GridPtr< HGridType >( filename );
    ParameterType::appendDGF( filename );

    // load balance grid in case of parallel runs
    gridptr->loadBalance();

    // and refine the grid globally
    Dune::AdaptationParameters adaptParam( Dune::ParameterKey::generate( "", "fem.adaptation." ) );
    const int startLevel = adaptParam.coarsestLevel();
    for(int level=0; level < startLevel ; ++level)
      Dune::Fem::GlobalRefine::apply(*gridptr, 1 );
  }

  // initialize FemEoc if eocSteps > 1
  EocParameters eocParam( Dune::ParameterKey::generate( "", "fem.eoc." ) );
  Dune::Fem::FemEoc::initialize( eocParam.outputPath(), eocParam.fileName(), problemDescription );

  return gridptr;
}

//! get memory in MB
inline double getMemoryUsage()
{
  struct rusage info;
  getrusage( RUSAGE_SELF, &info );
  return (info.ru_maxrss / 1024.0);
}

/////////////////////////////////////////////////////////////////////////////
//
//  compute algorithm
//
/////////////////////////////////////////////////////////////////////////////
template <class Algorithm>
void compute(Algorithm& algorithm)
{
  typedef typename Algorithm::GridType GridType;
  GridType& grid = algorithm.grid();

  typename Algorithm::IOTupleType dataTup = algorithm.dataTuple() ;

  typedef Dune::Fem::DataOutput<GridType, typename Algorithm::IOTupleType> DataOutputType;
  DataOutputType dataOutput( grid, dataTup );

  typedef typename Algorithm::SolverMonitorType SolverMonitorType;

  // initialize FemEoc if eocSteps > 1
  EocParameters eocParam( Dune::ParameterKey::generate( "", "fem.eoc." ) );
  Dune::Fem::FemEoc::initialize( eocParam.outputPath(), eocParam.fileName(), problemDescription );

  Dune::Fem::FemEoc::clear();


  const unsigned int femTimerId = Dune::FemTimer::addTo("timestep");
  for(int eocloop=0; eocloop < eocParam.eocSteps(); ++eocloop )
  {
    // start fem timer (enable with -DFEMTIMER)
    Dune::FemTimer :: start(femTimerId);

    // additional timer
    Dune::Timer timer ;

    // call algorithm and return solver statistics and some info
    SolverMonitorType monitor = algorithm.solve( eocloop );

    // get run time
    const double runTime = timer.elapsed();

    // also get times for FemTimer if enabled
    Dune::FemTimer::stop(femTimerId);
    Dune::FemTimer::printFile("./timer.out");
    Dune::FemTimer::reset(femTimerId);

    // finalize the algorithm
    algorithm.finalize( eocloop );

    std::stringstream eocInfo ;
    // generate EOC information
    Dune::Fem::FemEoc::write( monitor.gridWidth, monitor.elements, runTime, monitor.timeSteps,
                              monitor.doubleValues(), monitor.intValues(), eocInfo );

    // in verbose mode write EOC info to std::cout
    if( ParameterType :: verbose() )
    {
      std::cout << std::endl << "EOC info: " << std::endl << eocInfo.str() << std::endl;
    }

    // write eoc step
    dataOutput.writeData( eocloop );

    // Refine the grid for the next EOC Step. If the scheme uses adaptation,
    // the refinement level needs to be set in the algorithms' initialize method.
    if( eocParam.eocSteps() > 1 )
      Dune::Fem::GlobalRefine::apply(grid,Dune::DGFGridInfo<GridType>::refineStepsForHalf());

  } /***** END of EOC Loop *****/

  // FemTimer cleanup
  Dune::FemTimer::removeAll();
}

#endif
