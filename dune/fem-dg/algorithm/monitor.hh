#ifndef DUNE_FEMDG_ALGORITHM_MONITOR_HH
#define DUNE_FEMDG_ALGORITHM_MONITOR_HH

#include <iostream>

#include <dune/fem/solver/timeprovider.hh>

namespace Dune
{
namespace Fem
{

  //! solver statistics and further info
  //! such as newton iterations and iterations of the linear solver
  class SolverMonitor
  {
  public:
    typedef std::pair< std::string, double > DoublePairType;
    typedef std::pair< std::string, int    > IntPairType;

    double gridWidth;
    double avgTimeStep;
    double minTimeStep;
    double maxTimeStep;
    size_t timeSteps;
    int newton_iterations;
    int ils_iterations;
    int total_newton_iterations;
    int total_ils_iterations;
    int max_newton_iterations;
    int max_ils_iterations;
    int operator_calls;
    int total_operator_calls;

    unsigned long elements ;

    SolverMonitor()
    : gridWidth( 0 ),
      avgTimeStep( 0 ),
      minTimeStep( std::numeric_limits<double>::max() ),
      maxTimeStep( 0 ),
      timeSteps( 0 ),
      newton_iterations( 0 ),
      ils_iterations( 0 ),
      total_newton_iterations( 0 ),
      total_ils_iterations( 0 ),
      max_newton_iterations( 0 ),
      max_ils_iterations( 0 ),
      operator_calls( 0 ),
      total_operator_calls( 0 ),
      elements( 0 )
    {}

    void setTimeStepInfo( const Dune::Fem::TimeProviderBase& tp )
    {
      if( tp.timeStepValid() )
      {
        const double ldt = tp.deltaT();
        // calculate time step info
        minTimeStep  = std::min( ldt, minTimeStep );
        maxTimeStep  = std::max( ldt, maxTimeStep );
        avgTimeStep += ldt;

        timeSteps = tp.timeStep() + 1;

        // accumulate iteration information
        total_newton_iterations += newton_iterations;
        total_ils_iterations += ils_iterations;
        total_operator_calls += operator_calls;
      }
    }

    void finalize( const double h = 0, const unsigned long el = 0 )
    {
      avgTimeStep /= double( timeSteps );
      gridWidth = h;
      elements = el;
    }

    std::vector< DoublePairType > doubleValues() const
    {
      std::vector< DoublePairType > values;
      values.reserve( 3 );
      values.push_back( DoublePairType("avg dt", avgTimeStep ) );
      values.push_back( DoublePairType("min dt", minTimeStep ) );
      values.push_back( DoublePairType("max dt", maxTimeStep ) );
      // causes warning: moving a local object in a return statement prevents copy elision [-Wpessimizing-move]
      // return std::move( values );
      return values;
    }

    std::vector< IntPairType > intValues() const
    {
      std::vector< IntPairType > values;
      values.reserve( 4 );
      values.push_back( IntPairType("Newton", total_newton_iterations ) );
      values.push_back( IntPairType("ILS", total_ils_iterations ) );
      values.push_back( IntPairType("max{Newton/linS}", max_newton_iterations ) );
      values.push_back( IntPairType("max{ILS/linS}", max_ils_iterations ) );
      values.push_back( IntPairType("OCs:", total_operator_calls ) );
      // causes warning: moving a local object in a return statement prevents copy elision [-Wpessimizing-move]
      // return std::move( values );
      return values;
    }

    void dump( std::ostream& out ) const
    {
      out << "SolverMonitorInfo (timesteps):" << std::endl;
      out << "min  dt: " << minTimeStep << std::endl;
      out << "max  dt: " << maxTimeStep << std::endl;
      out << "aver dt: " << avgTimeStep << std::endl;

      if( max_newton_iterations > 0 || max_ils_iterations > 0 )
      {
        out << "SolverMonitorInfo (iterations):" << std::endl;
        out << "newton   max: " << max_newton_iterations << std::endl;
        out << "newton   tot: " << total_newton_iterations << std::endl;
        out << "ils      max: " << max_ils_iterations << std::endl;
        out << "ils      tot: " << total_ils_iterations << std::endl;
        out << "op calls tot: " << total_operator_calls << std::endl;
      }
      out << std::endl;
    }
  };

  // ostream operator for SolverMonitor
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor& monitor )
  {
    monitor.dump( out );
    return out;
  }

} // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_MONITOR_HH
