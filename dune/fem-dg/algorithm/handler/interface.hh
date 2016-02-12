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
    template< class... A > void initializeStart( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.initialize()`
     */
    template< class... A > void initializeEnd( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.preSolve()`
     */
    template< class... A > void preSolveStart( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.preSolve()`
     */
    template< class... A > void preSolveEnd( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.solve()`
     */
    template< class... A > void solveStart( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.solve()`
     */
    template< class... A > void solveEnd( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.postSolve()`
     */
    template< class... A > void postSolveStart( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.postSolve()`
     */
    template< class... A > void postSolveEnd( A&&... ) {}

    /**
     * \brief Handler interface method which is called shortly _before_ `algorithm.finalize()`
     */
    template< class... A > void finalizeStart( A&&... ) {}
    /**
     * \brief Handler interface method which is called shortly _after_ `algorithm.finalize()`
     */
    template< class... A > void finalizeEnd( A&&... ) {}
  };

}
}

#endif
