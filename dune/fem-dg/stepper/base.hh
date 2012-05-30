#ifndef FEMHOWTO_BASE_HH
#define FEMHOWTO_BASE_HH

#include <sstream>

#include <dune/common/misc.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>

//! FemEoc initialization helper
struct InitFemEoc 
{
  //! return number of eoc steps 
  static inline int eocSteps() 
  {
    return Dune::Parameter::getValue<int>("femhowto.eocSteps", 1);
  }

  //! initialize FemEoc if eocSteps is > 1 
  static inline void initializeFemEoc( const std::string& problemDescription )
  {
    if( eocSteps() > 1 ) 
    {
      // output of error and eoc information
      std::string eocOutPath = 
        Dune::Parameter::getValue<std::string>("femhowto.eocOutputPath", std::string("./"));

      Dune::FemEoc::initialize(eocOutPath, "eoc", problemDescription);
    }
  }
};

inline std::string checkPointRestartFileName () 
{
  const std::string key( "fem.io.checkpointrestartfile" );
  if( Dune :: Parameter :: exists( key ) )
  {
    return Dune :: Parameter::getValue<std::string> ( key );
  }
  else 
    return std::string( "" );
}

template< class HGridType >
Dune::GridPtr< HGridType > initialize( const std::string& problemDescription )
{ 
#ifdef ALUGRID_CONSTRUCTION_WITH_STREAMS 
  // precision of out streams (here ALUGrid backup streams)
  const int precision = Dune :: Parameter :: getValue< int > ("fem.io.precision", 16);
  ALUGridSpace :: ALUGridExternalParameters :: precision () = precision;
#endif

  // grid pointer object 
  Dune::GridPtr< HGridType > gridptr ; 

  const std::string checkPointRestartFile = checkPointRestartFileName();

  // if checkpoint restart file has been given
  if( checkPointRestartFile.size() > 0 ) 
  {
     // restore grid from checkpoint
     gridptr = Dune::CheckPointer< HGridType > :: restoreGrid( checkPointRestartFile );
  }
  else  // normal new start 
  {
    // ----- read in runtime parameters ------
    const std::string filekey = Dune::IOInterface::defaultGridKey( HGridType::dimension );
    const std::string filename = Dune::Parameter::getValue< std::string >( filekey ); /*@\label{base:param0}@*/

    // initialize grid with given macro file name 
    gridptr = Dune::GridPtr< HGridType >( filename );
    Dune::Parameter::appendDGF( filename );

    // load balance grid in case of parallel runs 
    gridptr->loadBalance();

    // and refine the grid globally
    const int startLevel = Dune::Parameter::getValue<int>("fem.adaptation.coarsestLevel", 0);
    for(int level=0; level < startLevel ; ++level)
      Dune::GlobalRefine::apply(*gridptr, 1 ); /*@\label{fv:globalRefine1}@*/
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
  typedef typename Algorithm::DiscreteFunctionType DiscreteFunctionType;
  typename Algorithm::DiscreteSpaceType& space = algorithm.space();
  typename Algorithm::GridPartType& gridPart = space.gridPart();
  typedef typename Algorithm::GridPartType::GridType GridType;
  GridType& grid = gridPart.grid();

  // get number of eoc steps 
  const int eocSteps = InitFemEoc :: eocSteps ();

  typename Algorithm::IOTupleType dataTup ( &algorithm.solution() );
  typedef Dune::DataOutput<GridType, typename Algorithm::IOTupleType> DataOutputType;
  DataOutputType dataOutput( grid, dataTup );

  const unsigned int femTimerId = Dune::FemTimer::addTo("timestep");
  for(int eocloop=0; eocloop < eocSteps; ++eocloop)
  {
    Dune::FemTimer :: start(femTimerId);

    // do one step
    double avgTimeStep = 0;
    double minTimeStep = 0;
    double maxTimeStep = 0;
    size_t counter = 0;
    int total_newton_iterations = 0;
    int total_ils_iterations = 0;
    int max_newton_iterations = 0;
    int max_ils_iterations = 0;

    algorithm( avgTimeStep, minTimeStep, maxTimeStep,
               counter, total_newton_iterations, total_ils_iterations,
               max_newton_iterations, max_ils_iterations );

    double runTime = Dune::FemTimer::stop(femTimerId);

    Dune::FemTimer::printFile("./timer.out");
    Dune::FemTimer::reset(femTimerId);

    // calculate grid width
    const double h = Dune::GridWidth::calcGridWidth(gridPart);

    // finalize the algorithm 
    algorithm.finalize( eocloop );

    // only do this if we have more than 1 eoc step
    if( eocSteps > 1 ) 
    {
      if( Dune::Parameter :: verbose() )
        Dune::FemEoc::write(h,grid.size(0), runTime, counter,avgTimeStep,minTimeStep,
                      maxTimeStep, total_newton_iterations, total_ils_iterations, 
                      max_newton_iterations, max_ils_iterations, std::cout);
      else
        Dune::FemEoc::write(h,grid.size(0),runTime,counter,avgTimeStep,minTimeStep,
                      maxTimeStep,total_newton_iterations,total_ils_iterations,
                      max_newton_iterations, max_ils_iterations);

      // write eoc step 
      dataOutput.writeData(eocloop);

      // Refine the grid for the next EOC Step. If the scheme uses adaptation,
      // the refinement level needs to be set in the algorithms' initialize method.
      Dune::GlobalRefine::apply(grid,Dune::DGFGridInfo<GridType>::refineStepsForHalf());
    }

  } /***** END of EOC Loop *****/

  // FemTimer cleanup
  Dune::FemTimer::removeAll();
}

#endif
