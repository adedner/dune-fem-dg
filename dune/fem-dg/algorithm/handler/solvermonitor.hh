#ifndef FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH
#define FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH

#include <string>
#include <dune/fem-dg/algorithm/monitor.hh>

namespace Dune
{
namespace Fem
{

  class DefaultSolverMonitorHandler
  {
  public:
    typedef SolverMonitor<1>     SolverMonitorType;

    DefaultSolverMonitorHandler( const std::string keyPrefix = "" )
      : solverMonitor_()
    {}


    void stepPrint()
    {
      std::cout << ",  Newton: " << *solverMonitor_.newton_iterations
                << "  ILS: " << *solverMonitor_.ils_iterations
                << "  OC: " << *solverMonitor_.operator_calls << std::endl;
    }

    template< class OdeSolverMonitorImp, class TimeProviderImp >
    void step( OdeSolverMonitorImp& odeSolverMonitor, TimeProviderImp& tp )
    {
      *solverMonitor_.newton_iterations     = odeSolverMonitor.newtonIterations_;
      *solverMonitor_.ils_iterations        = odeSolverMonitor.linearSolverIterations_;
      *solverMonitor_.max_newton_iterations = odeSolverMonitor.maxNewtonIterations_ ;
      *solverMonitor_.max_ils_iterations    = odeSolverMonitor.maxLinearSolverIterations_;
      *solverMonitor_.operator_calls        = odeSolverMonitor.spaceOperatorCalls_;

      //set time step size to monitor
      solverMonitor_.setTimeStepInfo( tp );
    }

    void finalize( const double gridWidth, const double gridSize )
    {
      solverMonitor_.finalize( gridWidth, gridSize );
    }

    SolverMonitorType& monitor()
    {
      return solverMonitor_;
    }

  private:
    SolverMonitorType solverMonitor_;
  };


  class DefaultSteadyStateSolverMonitorHandler
  {
  public:
    typedef SolverMonitor<1>     SolverMonitorType;

    DefaultSteadyStateSolverMonitorHandler( const std::string keyPrefix = "" )
      : solverMonitor_()
    {}


    void stepPrint()
    {
      std::cout << ",  Newton: " << *solverMonitor_.newton_iterations
                << "  ILS: " << *solverMonitor_.ils_iterations << std::endl;
    }

    template< class OdeSolverMonitorImp, class TimeProviderImp >
    void step( OdeSolverMonitorImp& odeSolverMonitor, TimeProviderImp& tp )
    {
    }

    template< class LinearSolverMonitorImp >
    void finalize( const double gridWidth, const double gridSize, const LinearSolverMonitorImp& solver )
    {
      *solverMonitor_.gridWidth = gridWidth;
      *solverMonitor_.elements  = gridSize;

      //*solverMonitor_.ils_iterations = invDgOperator_->iterations();
      //*solverMonitor_.totalNewtonIterations_ = solver.iterations();
      //*solverMonitor_.totalLinearSolverIterations_ = solver.linearIterations();

      *solverMonitor_.newton_iterations     = 0;
      *solverMonitor_.ils_iterations        = solver.iterations();
      *solverMonitor_.max_newton_iterations = 0;
      *solverMonitor_.max_ils_iterations    = 0;
      *solverMonitor_.operator_calls        = 0;
    }

    SolverMonitorType& monitor()
    {
      return solverMonitor_;
    }

  private:
    SolverMonitorType solverMonitor_;
  };

  class NoSolverMonitorHandler
  {
  public:
    NoSolverMonitorHandler( const std::string keyPrefix = "" )
    {}

    template <class ... Args>
    void stepPrint(Args&& ... ) const {};

    template <class ... Args>
    void step(Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };

}
}

#endif
