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

    template< class TimeProviderImp, class DiscreteSolutionImp, class OdeSolverMonitorImp, class TimerImp, class AdaptationHandlerImp >
    void step( TimeProviderImp& tp, const DiscreteSolutionImp& solution,
               const OdeSolverMonitorImp& odeSolverMonitor,
               const TimerImp& timer,
               AdaptationHandlerImp& adaptHandler )
    {
      double maxNumDofs = solution.space().blockMapper().maxNumDofs() * solution.space().localBlockSize;
      const double nElements = odeSolverMonitor.numberOfElements_;
      const double ldt = tp.deltaT();

      std::vector< double > times;
      times.push_back( maxNumDofs );
      times.push_back( odeSolverMonitor.operatorTime_ );
      times.push_back( odeSolverMonitor.odeSolveTime_ );
      if( adaptHandler.adaptive() )
      {
        times.push_back( adaptHandler.adaptationTime() );
        times.push_back( adaptHandler.loadBalanceTime() );
      }
      times.push_back( timer.elapsed() );

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
    void step( Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };

}
}

#endif
