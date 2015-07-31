#ifndef FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH
#define FEMDG_ALGORITHM_DIAGNOSTICSHANDLER_HH

#include <dune/fem-dg/misc/diagnostics.hh>

namespace Dune
{
namespace Fem
{

  class DefaultDiagnosticsHandler
  {
  public:

    DefaultDiagnosticsHandler()
      : diagnostics_( true )
    {}

    template< class TimeProviderImp, class DiscreteSolutionImp, class OdeSolverMonitorImp, class TimerImp >
    void write( TimeProviderImp& tp, const DiscreteSolutionImp& solution,
                const OdeSolverMonitorImp& odeSolverMonitor,
                const TimerImp& timer )
    {
      double maxNumDofs = solution.space().blockMapper().maxNumDofs() * solution.space().localBlockSize;
      const double nElements = odeSolverMonitor.numberOfElements_;
      const double ldt = tp.deltaT();

      std::vector<double> times(4);
      times[0] = maxNumDofs;
      times[1] = odeSolverMonitor.operatorTime_;
      times[2] = odeSolverMonitor.odeSolveTime_;
      times[3] = timer.elapsed();

      diagnostics_.write( tp.time() + ldt, ldt, nElements, times );
    }

    template< class TimeProviderImp, class DiscreteSolutionImp, class OdeSolverMonitorImp, class TimerImp, class AdaptationManagerImp >
    void write( TimeProviderImp& tp, const DiscreteSolutionImp& solution,
                const OdeSolverMonitorImp& odeSolverMonitor,
                const TimerImp& timer,
                const AdaptationManagerImp& adaptManager )
    {
      double maxNumDofs = solution.space().blockMapper().maxNumDofs() * solution.space().localBlockSize;
      const double nElements = odeSolverMonitor.numberOfElements_;
      const double ldt = tp.deltaT();

      std::vector<double> times(6);
      times[0] = maxNumDofs;
      times[1] = odeSolverMonitor.operatorTime_;
      times[2] = odeSolverMonitor.odeSolveTime_;
      times[3] = adaptManager.adaptationTime();
      times[4] = adaptManager.loadBalanceTime();
      times[5] = timer.elapsed();

      diagnostics_.write( tp.time() + ldt, ldt, nElements, times );
    }


    void finalize() const
    {
      diagnostics_.flush();
    }


  private:
    Diagnostics diagnostics_;
  };



  class NoDiagnosticsHandler
  {
  public:
    template <class ... Args>
    void write( Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };

}
}

#endif
