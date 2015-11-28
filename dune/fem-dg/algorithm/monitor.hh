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
  template< int numSubMonitors >
  class SolverMonitor
  {
  public:
    typedef std::pair< std::string, double > DoublePairType;
    typedef std::pair< std::string, int    > IntPairType;

    double* gridWidth;
    double* avgTimeStep;
    double* minTimeStep;
    double* maxTimeStep;
    size_t* timeSteps;
    int* newton_iterations;
    int* ils_iterations;
    int* total_newton_iterations;
    int* total_ils_iterations;
    int* max_newton_iterations;
    int* max_ils_iterations;
    int* operator_calls;
    int* total_operator_calls;
    unsigned long* elements;

    SolverMonitor()
    {
      gridWidth_.fill( 0 );
      avgTimeStep_.fill( 0 );
      minTimeStep_.fill( std::numeric_limits<double>::max() );
      maxTimeStep_.fill( 0 );
      timeSteps_.fill( 0 );
      newton_iterations_.fill( 0 );
      ils_iterations_.fill( 0 );
      total_newton_iterations_.fill( 0 );
      total_ils_iterations_.fill( 0 );
      max_newton_iterations_.fill( 0 );
      max_ils_iterations_.fill( 0 );
      elements_.fill( 0 );
      operator_calls_.fill( 0 );
      total_operator_calls_.fill( 0 );

      setMonitorNumber( 0 );
    }

    void setTimeStepInfo( const Fem::TimeProviderBase& tp )
    {
      const double ldt = tp.deltaT();
      // calculate time step info
      minTimeStep_[ monitorNumber_ ]  = std::min( ldt, minTimeStep_[ monitorNumber_ ] );
      maxTimeStep_[ monitorNumber_ ]  = std::max( ldt, maxTimeStep_[ monitorNumber_ ] );
      avgTimeStep_[ monitorNumber_ ] += ldt;

      timeSteps_[ monitorNumber_ ] = tp.timeStep() + 1;

      // accumulate iteration information
      total_newton_iterations_[ monitorNumber_ ] += newton_iterations_[ monitorNumber_ ];
      total_ils_iterations_[ monitorNumber_ ] += ils_iterations_[ monitorNumber_ ];
      total_operator_calls_[ monitorNumber_ ] += operator_calls_[ monitorNumber_ ];
    }

    void finalize( const double h = 0, const unsigned long el = 0 )
    {
      avgTimeStep_[ monitorNumber_ ] /= double( timeSteps_[ monitorNumber_ ] );
      gridWidth_[ monitorNumber_ ] = h;
      elements_[ monitorNumber_ ] = el;
    }

    std::vector< DoublePairType > doubleValues() const
    {
      std::vector< DoublePairType > values;
      values.reserve( 3 );
      values.push_back( DoublePairType("avg dt", avgTimeStep_[ monitorNumber_ ] ) );
      values.push_back( DoublePairType("min dt", minTimeStep_[ monitorNumber_ ] ) );
      values.push_back( DoublePairType("max dt", maxTimeStep_[ monitorNumber_ ] ) );
      return std::move( values );
    }

    std::vector< IntPairType > intValues() const
    {
      std::vector< IntPairType > values;
      values.reserve( 4 );
      values.push_back( IntPairType("Newton", total_newton_iterations_[ monitorNumber_ ] ) );
      values.push_back( IntPairType("ILS", total_ils_iterations_[ monitorNumber_ ] ) );
      values.push_back( IntPairType("max{Newton/linS}", max_newton_iterations_[ monitorNumber_ ] ) );
      values.push_back( IntPairType("max{ILS/linS}", max_ils_iterations_[ monitorNumber_ ] ) );
      values.push_back( IntPairType("OCs:", total_operator_calls_[ monitorNumber_ ] ) );
      return std::move( values );
    }

    void dump( std::ostream& out, int monitorNumber = -1 ) const
    {
      if( monitorNumber == -1 )
        monitorNumber = monitorNumber_;

      out << "SolverMonitorInfo (timesteps):" << std::endl;
      out << "min  dt: " << minTimeStep[monitorNumber] << std::endl;
      out << "max  dt: " << maxTimeStep[monitorNumber] << std::endl;
      out << "aver dt: " << avgTimeStep[monitorNumber] << std::endl;

      if( max_newton_iterations > 0 || max_ils_iterations > 0 )
      {
        out << "SolverMonitorInfo (iterations):" << std::endl;
        out << "newton   max: " << max_newton_iterations_[monitorNumber] << std::endl;
        out << "newton   tot: " << total_newton_iterations_[monitorNumber] << std::endl;
        out << "ils      max: " << max_ils_iterations_[monitorNumber] << std::endl;
        out << "ils      tot: " << total_ils_iterations_[monitorNumber] << std::endl;
        out << "op calls tot: " << total_operator_calls_[monitorNumber] << std::endl;
      }
      out << std::endl;
    }

    void setMonitorNumber( int number )
    {
      monitorNumber_ = number;

      gridWidth = &gridWidth_[ monitorNumber_];
      avgTimeStep = &avgTimeStep_[ monitorNumber_];
      minTimeStep = &minTimeStep_[ monitorNumber_];
      maxTimeStep = &maxTimeStep_[ monitorNumber_];
      timeSteps = &timeSteps_[ monitorNumber_];
      newton_iterations = &newton_iterations_[ monitorNumber_];
      ils_iterations = &ils_iterations_[ monitorNumber_];
      total_newton_iterations = &total_newton_iterations_[ monitorNumber_];
      total_ils_iterations = &total_ils_iterations_[ monitorNumber_];
      max_newton_iterations = &max_newton_iterations_[ monitorNumber_];
      max_ils_iterations = &max_ils_iterations_[ monitorNumber_];
      elements = &elements_[ monitorNumber_];
      operator_calls = &operator_calls_[ monitorNumber_];
      total_operator_calls = &total_operator_calls_[ monitorNumber_];
   }

    void monitorNumber() const
    {
      return monitorNumber_;
    }


    private:

    int monitorNumber_;
    std::array<double, numSubMonitors> gridWidth_;
    std::array<double, numSubMonitors> avgTimeStep_;
    std::array<double, numSubMonitors> minTimeStep_;
    std::array<double, numSubMonitors> maxTimeStep_;
    std::array<size_t, numSubMonitors> timeSteps_;
    std::array<int, numSubMonitors> newton_iterations_;
    std::array<int, numSubMonitors> ils_iterations_;
    std::array<int, numSubMonitors> total_newton_iterations_;
    std::array<int, numSubMonitors> total_ils_iterations_;
    std::array<unsigned long, numSubMonitors> elements_;
    std::array<int, numSubMonitors> max_newton_iterations_;
    std::array<int, numSubMonitors> max_ils_iterations_;
    std::array<int, numSubMonitors> operator_calls_;
    std::array<int, numSubMonitors> total_operator_calls_;
  };

  // ostream operator for SolverMonitor
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<1>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<2>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<3>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<4>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<5>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<6>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<7>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<8>& monitor ) { monitor.dump( out ); return out; }
  inline std::ostream& operator << ( std::ostream& out, const SolverMonitor<9>& monitor ) { monitor.dump( out ); return out; }

} // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_MONITOR_HH
