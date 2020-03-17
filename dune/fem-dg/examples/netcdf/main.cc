// configure macros
#include <config.h>

#include <dune/fem-dg/misc/simulator.hh>

#if HAVE_NETCDF
#include <netcdf.h>

int main(int argc, char ** argv)
{
  /* Initialize MPI (always do this even if you are not using MPI) */
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try {
    // read Parameters
    if( !readParameters( argc, argv ) )
      return 1;

    int ncid = -1;
    int err = nc_create("./ncfile", NC_NOCLOBBER, &ncid );

    // inquire dimensions etc.
    int dims  = -1;
    int nvars = -1;
    int natts = -1;
    int unlimdimid = -1;
    err = nc_inq(ncid, &dims, &nvars, &natts, &unlimdimid);

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
#else
int main(int argc, char ** argv)
{
  return 0;
}
#endif
