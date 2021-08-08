#ifndef DUNE_FEM_PTHREADCLASS_HH
#define DUNE_FEM_PTHREADCLASS_HH
#warning "Deprecated header, use #include <dune/fem/misc/threads/threadpool.hh> instead!"
#include <dune/fem/misc/threads/threadpool.hh>

namespace Dune
{
  namespace Fem
  {
    //! deprecated typedef, use ThreadPool
    typedef ThreadPool ThreadPool;
  }
}
#endif
