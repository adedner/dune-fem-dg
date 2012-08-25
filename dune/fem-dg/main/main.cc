#if defined GRIDDIM
#ifndef CODEDIM
#define CODEDIM GRIDDIM
#endif
#endif

// configure macros 
#include <config.h>

// streams for backup
#include <dune/fem-dg/misc/streams.hh>

// local includes
#include <dune/fem-dg/main/main.hh>
#include <dune/fem-dg/main/codegen.hh>

#include <dune/fem/misc/threadmanager.hh>
#include <dune/fem-dg/pass/threadpass.hh>

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/solver/odesolver.hh>

#if HAVE_PETSC
#include <petsc.h>
#endif

#if HAVE_SLEPC
#include <slepc.h>
#endif


#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
#include <dune/fem/space/basefunctions/codegen.hh>

std::string autoFilename(const int dim, const int polynomialOrder ) 
{
  std::stringstream name; 
  name << "autogeneratedcode_" << dim << "d_k" << polynomialOrder << ".hh";
  return name.str();
}

void cleanup(const int k, const std::string& filename ) 
{
  std::cerr << "Code for k="<< k << " generated!! " << std::endl;
  Dune::Fem::CodegenInfo :: instance().dumpInfo();
  Dune::Fem::CodegenInfo :: instance().clear();
  std::remove( filename.c_str() );
}

void codegen()
{
  // only one thread for codegen 
  Dune :: Fem :: ThreadManager :: setMaxNumberThreads( 1 );

  try 
  {
    const int dim = CODEDIM;
    std::cout << "Generating Code \n";
    std::string filename;
#ifdef ONLY_ONE_P
    try {
      filename = autoFilename( dim, POLORDER ) ;
      Dune::Fem::CodegenInfo :: instance().setFileName( filename );
      DG_ONE_P :: simulate();  
    }
    catch (Dune::Fem::CodegenInfoFinished) 
    {
      cleanup( POLORDER , filename );
    }
#else 
#ifndef ONLY_P1_P2
    try {
      filename = autoFilename( dim, 0 ) ;
      Dune::Fem::CodegenInfo :: instance().setFileName( filename );
      DG_P0 :: simulate();  
    }
    catch (Dune::Fem::CodegenInfoFinished) 
    {
      cleanup( 0 , filename );
    }
#endif
    try {
      filename = autoFilename( dim, 1 ) ;
      Dune::Fem::CodegenInfo :: instance().setFileName( filename );
      DG_P1 :: simulate();  
    }
    catch (Dune::Fem::CodegenInfoFinished) 
    {
      cleanup( 1 , filename );
    }
    try 
    {
      filename = autoFilename( dim, 2 ) ;
      Dune::Fem::CodegenInfo :: instance().setFileName( filename );
      DG_P2 :: simulate();  
    }
    catch (Dune::Fem::CodegenInfoFinished) 
    {
      cleanup( 2 , filename );
    }
#ifndef ONLY_P1_P2
    try 
    {
      filename = autoFilename( dim, 3 ) ;
      Dune::Fem::CodegenInfo :: instance().setFileName( filename );
      DG_P3 :: simulate();  
    }
    catch (Dune::Fem::CodegenInfoFinished) 
    {
      cleanup( 3 , filename );
    }
#endif
#ifdef ALSO_P4
    try 
    {
      filename = autoFilename( dim, 4 ) ;
      Dune::Fem::CodegenInfo :: instance().setFileName( filename );
      DG_P4 :: simulate();  
    }
    catch (Dune::Fem::CodegenInfoFinished) 
    {
      cleanup( 4 , filename );
    }
#endif
#endif
  }
  catch (Dune::Fem::CodegenInfoFinished) {}


  //////////////////////////////////////////////////
  //  write include header 
  //////////////////////////////////////////////////

  std::ofstream file( "autogeneratedcode.hh" );

  if( file ) 
  {
    for( int dim=1; dim<4; ++dim ) 
    {
      std::stringstream dimfilename;
      dimfilename << "./autogeneratedcode/autogeneratedcode_" << dim << "d.hh";
      file << "#if CODEDIM == " << dim << std::endl;
      file << "#include \"" << dimfilename.str() << "\"" << std::endl;
      file << "#endif" << std::endl;
      
      std::ofstream dimfile( dimfilename.str().c_str() );
      if( dimfile )  
      {
        const int maxPolOrder = ( dim == 3 ) ? 4 : 8 ;
        // max polorder is 4 in 3d and 8 in 2d at the moment 
        for( int i=0; i <= maxPolOrder; ++ i ) 
        {
          dimfile << "#if POLORDER == " << i << std::endl;
          dimfile << "#include \"" << autoFilename( dim, i ) << "\"" << std::endl;
          dimfile << "#endif" << std::endl;
        }
      }
    }
  }

  // dump all information 
  std::cerr << "All automated code generated, bye, bye !! " << std::endl;
}
#endif

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
 * @return 0 we don't program bugs. :)
 */
int main(int argc, char ** argv) 
{     

  /* Initialize MPI (always do this even if you are not using MPI) */
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try {

#if HAVE_PETSC
    static char help[] = "Petsc-Slepc init";
#endif
#if HAVE_SLEPC
    SlepcInitialize(&argc,&argv,(char*)0,help);
#elif HAVE_PETSC 
    PetscInitialize(&argc,&argv,(char *)0,help);
#endif

  // *** Initialization
  Dune::Fem::Parameter::append(argc,argv);                           /*@\label{dg:param0}@*/
  if (argc >= 2) 
  {
    Dune::Fem::Parameter::append(argv[1]);
  } 
  else 
  {
    Dune::Fem::Parameter::append("parameter");                       /*@\label{dg:paramfile}@*/
  }                                                       /*@\label{dg:param1}@*/

  // get number of desired threads (default is 1)
  int numThreads = Dune::Fem::Parameter::getValue< int >("fem.parallel.numberofthreads", 1);
  Dune :: Fem :: ThreadManager :: setMaxNumberThreads( numThreads );

  int polynomialOrder = 1;
  polynomialOrder = Dune::Fem::Parameter :: getValue("femhowto.polynomialOrder", polynomialOrder );

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
  codegen();
#else 

  // code for only one polynomial degree 
#ifdef ONLY_ONE_P 
  {
    DG_ONE_P :: simulate();  
  }
#else 

  // code for several polynomial degrees 
#ifndef ONLY_P1_P2
  if( polynomialOrder == 0 ) 
  {
    DG_P0 :: simulate();  
  }
#endif
  if( polynomialOrder == 1 ) 
  {
    DG_P1 :: simulate();  
  }
  else if( polynomialOrder == 2 )
  {
    DG_P2 :: simulate();  
  }
#ifndef ONLY_P1_P2
  else if( polynomialOrder == 3 )
  {
    DG_P3 :: simulate();  
  }
#endif
#ifdef ALSO_P4
  else if( polynomialOrder == 4 )
  {
    DG_P4 :: simulate();  
  }
#endif
  else 
  {
    std::cerr << "Got wrong polynomialOrder " << polynomialOrder << ", valid are 1,2 ! " << std::endl;
    abort();
  }
#endif

  // write parameters used 
  Dune::Fem::Parameter::write("parameter.log");
#endif
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
