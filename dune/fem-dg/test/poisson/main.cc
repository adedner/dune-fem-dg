// configure macros
#include <config.h>

#include "problemcreator.hh"
#include <dune/fem-dg/misc/simulator.hh>


int main(int argc, char ** argv)
{
  /* Initialize MPI (always do this even if you are not using MPI) */
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try {
    // *** Initialization
    Dune::Fem::Parameter::append(argc,argv);
    if (argc >= 2)
      Dune::Fem::Parameter::append(argv[1]);
    else
      Dune::Fem::Parameter::append("parameter");

    // write parameters used (before simulation starts)
    Dune::Fem::Parameter::write("parameter.log");

    typedef Dune::GridSelector :: GridType GridType;
    PoissonProblemCreator< GridType > problem;

    // run simulation
    Simulator::run( problem );
  }
  catch (const Dune::Exception &e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }
  return 0;
}