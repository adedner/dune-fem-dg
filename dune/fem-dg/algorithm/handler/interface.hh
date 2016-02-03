#ifndef FEMDG_ALGORITHM_HANDLERINTERFACE_HH
#define FEMDG_ALGORITHM_HANDLERINTERFACE_HH

namespace Dune
{
namespace Fem
{
  /**
   * \brief Interface class for handlers.
   *
   * A handler is a class which provides predefined methods
   * which are called by an algorithm in a predefined order.
   *
   * An algorithm requires at least five basic methods from an sub-algorithm
   * * initialize()
   * * preSolve()
   * * solve()
   * * postSolve()
   * * finalize()
   *
   * To modify the behaviour of subalgorithms additional classes (called handlers) are
   * used.
   *
   * The next figure shows how it should work.
   *
   * ![Example for some handlers which are clipped into the algorithm](algorithms_scheme.png)
   */
  struct HandlerInterface
  {
    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.initialize()`
     */
    template< class... A > void initialize_pre( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.initialize()`
     */
    template< class... A > void initialize_post( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.preSolve()`
     */
    template< class... A > void preSolve_pre( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.preSolve()`
     */
    template< class... A > void preSolve_post( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.solve()`
     */
    template< class... A > void solve_pre( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.solve()`
     */
    template< class... A > void solve_post( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.postSolve()`
     */
    template< class... A > void postSolve_pre( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.postSolve()`
     */
    template< class... A > void postSolve_post( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.finalize()`
     */
    template< class... A > void finalize_pre( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.finalize()`
     */
    template< class... A > void finalize_post( A&&... ) {}
  };

}
}

#endif
