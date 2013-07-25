#ifndef FEMHOWTO_BASE_HH
#define FEMHOWTO_BASE_HH

#include <sstream>

#include <dune/common/misc.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/common/timer.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>

typedef Dune::Fem::Parameter  ParameterType ;

//! FemEoc initialization helper
struct InitFemEoc 
{
  //! return number of eoc steps 
  static inline int eocSteps() 
  {
    return ParameterType::getValue<int>("femhowto.eocSteps", 1);
  }

  //! initialize FemEoc if eocSteps is > 1 
  static inline void initializeFemEoc( const std::string& problemDescription )
  {
    if( eocSteps() > 1 ) 
    {
      // output of error and eoc information
      std::string eocOutPath = 
        ParameterType::getValue<std::string>("femhowto.eocOutputPath", std::string("./"));

      std::string eocFileName = 
        ParameterType::getValue<std::string>("femhowto.eocFileName", std::string("eoc"));

      Dune::Fem::FemEoc::initialize(eocOutPath, eocFileName, problemDescription);
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
Dune::GridPtr< HGridType > initialize( const std::string& problemDescription )
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
    gridptr = Dune::Fem::CheckPointer< HGridType > :: restoreGrid( checkPointRestartFile );
  }
  else  // normal new start 
  {
    // ----- read in runtime parameters ------
    const std::string filekey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
    const std::string filename = ParameterType::getValue< std::string >( filekey ); /*@\label{base:param0}@*/

    // initialize grid with given macro file name 
    gridptr = Dune::GridPtr< HGridType >( filename );
    ParameterType::appendDGF( filename );

    // load balance grid in case of parallel runs 
    gridptr->loadBalance();

    // and refine the grid globally
    const int startLevel = ParameterType::getValue<int>("fem.adaptation.coarsestLevel", 0);
    for(int level=0; level < startLevel ; ++level)
      Dune::Fem::GlobalRefine::apply(*gridptr, 1 ); /*@\label{fv:globalRefine1}@*/
  }

  // initialize FemEoc if eocSteps > 1 
  InitFemEoc :: initializeFemEoc( problemDescription );

  return gridptr;
} /*@LST0E@*/

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
  const typename Algorithm::DiscreteSpaceType& space = algorithm.space();
  const typename Algorithm::GridPartType& gridPart = space.gridPart();
  typedef typename Algorithm::GridPartType::GridType GridType;
  GridType& grid = algorithm.grid();

  // get number of eoc steps 
  const int eocSteps = InitFemEoc :: eocSteps ();

  typename Algorithm::IOTupleType dataTup = algorithm.dataTuple() ; 

  typedef Dune::Fem::DataOutput<GridType, typename Algorithm::IOTupleType> DataOutputType;
  DataOutputType dataOutput( grid, dataTup );

  typedef typename Algorithm::SolverMonitorType SolverMonitorType;

  const unsigned int femTimerId = Dune::FemTimer::addTo("timestep");
  for(int eocloop=0; eocloop < eocSteps; ++eocloop )
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

    // calculate grid width
    const double h = Dune::Fem::GridWidth::calcGridWidth(gridPart);

    // finalize the algorithm 
    algorithm.finalize( eocloop );

    // only do this if we have more than 1 eoc step
    if( eocSteps > 1 ) 
    {
      std::stringstream eocInfo ;
      // generate EOC information
      Dune::Fem::FemEoc::write(h,grid.size(0), runTime, 
                      monitor.timeSteps, monitor.avgTimeStep, monitor.minTimeStep,
                      monitor.maxTimeStep, monitor.total_newton_iterations, monitor.total_ils_iterations, 
                      monitor.max_newton_iterations, monitor.max_ils_iterations, eocInfo );
       
      // in verbose mode write EOC info to std::cout  
      if( ParameterType :: verbose() )
      {
        std::cout << std::endl << "EOC info: " << std::endl << eocInfo.str() << std::endl;
      }

      // write eoc step 
      dataOutput.writeData( eocloop );

      // Refine the grid for the next EOC Step. If the scheme uses adaptation,
      // the refinement level needs to be set in the algorithms' initialize method.
      Dune::Fem::GlobalRefine::apply(grid,Dune::DGFGridInfo<GridType>::refineStepsForHalf());
    }

  } /***** END of EOC Loop *****/

  // FemTimer cleanup
  Dune::FemTimer::removeAll();
}

#endif
