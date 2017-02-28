// configure macros
#include <config.h>

#include <dune/fem-dg/misc/simulator.hh>
#include "algorithmcreator.hh"

int main(int argc, char ** argv)
{
  /* Initialize MPI (always do this even if you are not using MPI) */
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try {
    // read Parameters
    if( !readParameters( argc, argv ) )
      return 1;

    Dune::Fem::PoissonAlgorithmCreator<Dune::GridSelector::GridType> algorithmCreator;

    // run simulation
    Dune::Fem::Simulator::run( algorithmCreator );
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
