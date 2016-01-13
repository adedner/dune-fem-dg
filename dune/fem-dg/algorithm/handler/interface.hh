#ifndef FEMDG_ALGORITHM_HANDLERINTERFACE_HH
#define FEMDG_ALGORITHM_HANDLERINTERFACE_HH

namespace Dune
{
namespace Fem
{
  struct HandlerInterface
  {
    template< class... A > void initialize_pre( A&&... ) {}
    template< class... A > void initialize_post( A&&... ) {}

    template< class... A > void preSolve_pre( A&&... ) {}
    template< class... A > void preSolve_post( A&&... ) {}

    template< class... A > void solve_pre( A&&... ) {}
    template< class... A > void solve_post( A&&... ) {}

    template< class... A > void postSolve_pre( A&&... ) {}
    template< class... A > void postSolve_post( A&&... ) {}

    template< class... A > void finalize_pre( A&&... ) {}
    template< class... A > void finalize_post( A&&... ) {}
  };

}
}

#endif
