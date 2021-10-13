#ifndef DUNE_FEM_PTHREADCLASS_HH
#define DUNE_FEM_PTHREADCLASS_HH
#warning "Deprecated header, use #include <dune/fem/misc/mpimanager.hh> instead!"
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{
  namespace Fem
  {
    //! deprecated typedef, use ThreadPool
    typedef MPIManager ThreadPool;
  }
}
#endif
