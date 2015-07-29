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

    template< class TimeProviderImp, class OdeSolverMonitorImp >
    void write( TimeProviderImp& tp, const OdeSolverMonitorImp& odeSolverMonitor )
    {
      const double nElements = odeSolverMonitor.numberOfElements_;
      const double ldt = tp.deltaT();

      std::vector<double> times(3);
      times[0] = odeSolverMonitor.operatorTime_;
      times[1] = odeSolverMonitor.odeSolveTime_;

      diagnostics_.write( tp.time() + ldt, ldt, nElements, times );
    }


    template< class TimeProviderImp, class OdeSolverMonitorImp >
    void write( TimeProviderImp& tp, const double maxNumDofs,
                const OdeSolverMonitorImp& odeSolverMonitor )
    {
      const double nElements = odeSolverMonitor.numberOfElements_;
      const double ldt = tp.deltaT();

      std::vector<double> times(3);
      times[0] = maxNumDofs;
      times[1] = odeSolverMonitor.operatorTime_;
      times[2] = odeSolverMonitor.odeSolveTime_;

      diagnostics_.write( tp.time() + ldt, ldt, nElements, times );
    }


    template< class TimeProviderImp, class OdeSolverMonitorImp >
    void write( TimeProviderImp& tp, const double maxNumDofs,
                const OdeSolverMonitorImp& odeSolverMonitor,
                const double time )
    {
      const double nElements = odeSolverMonitor.numberOfElements_;
      const double ldt = tp.deltaT();

      std::vector<double> times(4);
      times[0] = maxNumDofs;
      times[1] = odeSolverMonitor.operatorTime_;
      times[2] = odeSolverMonitor.odeSolveTime_;
      times[3] = time;

      diagnostics_.write( tp.time() + ldt, ldt, nElements, times );
    }

    template< class TimeProviderImp, class OdeSolverMonitorImp, class AdaptationManagerImp >
    void write( TimeProviderImp& tp, const double maxNumDofs,
                const OdeSolverMonitorImp& odeSolverMonitor,
                const double time,
                const AdaptationManagerImp& adaptManager )
    {
      const double nElements = odeSolverMonitor.numberOfElements_;
      const double ldt = tp.deltaT();

      std::vector<double> times(6);
      times[0] = maxNumDofs;
      times[1] = odeSolverMonitor.operatorTime_;
      times[2] = odeSolverMonitor.odeSolveTime_;
      times[3] = adaptManager.adaptationTime();
      times[4] = adaptManager.loadBalanceTime();
      times[5] = time;

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
    template <typename ... Args>
    void write(Args ... a) const {};

    void finalize() const {};
  };

}
}

#endif
