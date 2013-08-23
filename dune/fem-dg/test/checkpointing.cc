//#undef ENABLE_MPI
// include host specific macros set by configure script   /*@LST0S@*/
#include <config.h>

#include "checkpointing.hh"

/**
 * @brief main function for the LocalDG Advection-Diffusion application
 *
 * \c main starts the Simulation of an advection-diffusion pde with
 * the localdg method with EOC analysis and visual output to grape, paraview or
 * gnuplot.
 * \attention The localdg implementation uses the \c Dune::Pass
 * concept.
 *
 * @param argc number of arguments from command line
 * @param argv array of arguments from command line
 * @param envp array of environmental variables
 * @return 0 we don't program bugs. :)
 */
int main(int argc, char ** argv, char ** envp) {          /*@LST0S@*/
  typedef Dune::GridSelector :: GridType GridType;

  /* Initialize MPI (always do this even if you are not using MPI) */
  Dune::Fem::MPIManager :: initialize( argc, argv );

#if HAVE_ALUGRID 
  ALUGridSpace :: ALUGridExternalParameters :: precision () = 4;
#endif

  try 
  {

    // *** Initialization
    Dune::Fem::Parameter::append(argc,argv);                           /*@\label{dg:param0}@*/
    if (argc >= 2) {
      // the first parameter is the parameter file 
      Dune::Fem::Parameter::append( argv[ 1 ] );
    } 
    else 
    {
      Dune::Fem::Parameter::append("parameter");                       /*@\label{dg:paramfile}@*/
    }                                                       /*@\label{dg:param1}@*/

  // ProblemType is a Dune::Function that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  typedef Dune::U0< GridType > ProblemType;
  ProblemType problem;

  const bool newStart = (argc < 3);
  Dune::GridPtr<GridType> gridptr;
  const char* checkPointFile = 0; 
  std::string problemDescr( std::string("CheckPointing") + problem.description() );  
  if( newStart ) 
  {
    // Note to me: problem description is for FemEOC
    gridptr = initialize< GridType >( problemDescr ); /*@\label{dg:gridinit}@*/
  }
  else 
  {
    // assume that second parameter is the checkoint file  
    checkPointFile = argv[ 2 ];
    if( 0 == Dune::Fem::MPIManager :: rank () ) 
    {
      std::cout << std::endl;
      std::cout << "********************************************************" << std::endl;
      std::cout << "**   Restart from checkpoint `" << checkPointFile << "' " << std::endl;
      std::cout << "********************************************************" << std::endl;
      std::cout << std::endl;
    }

    // restore grid from checkpoint
    gridptr = Dune::Fem::CheckPointer< GridType > :: restoreGrid( checkPointFile );

    // initialize FemEoc 
    initializeFemEoc( problemDescr );
  }

  // get grid reference 
  GridType & grid = *gridptr;

  /*********************************************
   * EOC Loop                                  *
   ********************************************/
  Stepper<GridType, ProblemType :: BaseType> stepper(grid, problem, checkPointFile);  /*@\label{dg:stepper}@*/

  compute(stepper);                                                   /*@\label{dg:evolve}@*/

  Dune::Fem::Parameter::write("parameter.log");

  }
  catch (Dune::Exception &e) {                            /*@\label{dg:catch0}@*/
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }                                                      /*@\label{dg:catch1}@*/

  return 0;
}

