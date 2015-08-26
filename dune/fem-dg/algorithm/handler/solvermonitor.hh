#ifndef FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH
#define FEMDG_ALGORITHM_SOLVERMONITORHANDLER_HH

#include <string>
#include <dune/fem-dg/algorithm/monitor.hh>

namespace Dune
{
namespace Fem
{

  template< class... StepperArg >
  class CombinedDefaultSolverMonitorHandler;


  template< class StepperHead, class ... StepperArg >
  class CombinedDefaultSolverMonitorHandler< StepperHead, StepperArg... >
  {
  public:
    typedef std::tuple< typename std::add_pointer< StepperHead >::type,
      typename std::add_pointer< StepperArg >::type... >                               StepperTupleType;
    typedef typename StepperHead::GridType                                             GridType;

    template< int i >
    struct SetTimeStepInformation
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
         std::get< i >( tuple )->monitor().monitor().setTimeStepInfo( args... );
      }
    };
    template< int i >
    struct Finalize
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
         std::get< i >( tuple )->monitor().finalize( args... );
      }
    };
    template< int i >
    struct StepPrintNewtonIterations
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::cout << *std::get< i >( tuple )->monitor().monitor().newton_iterations;
        if( i != std::tuple_size<Tuple>::value - 1 )
          std::cout << " | " << std::endl;
      }
    };
    template< int i >
    struct StepPrintILSIterations
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::cout << *std::get< i >( tuple )->monitor().monitor().ils_iterations;
        if( i != std::tuple_size<Tuple>::value )
          std::cout << " | " << std::endl;
      }
    };
    template< int i >
    struct StepPrintOperatorCalls
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::cout << *std::get< i >( tuple )->monitor().monitor().operator_calls;
        if( i != std::tuple_size<Tuple>::value )
          std::cout << " | " << std::endl;
      }
    };
    CombinedDefaultSolverMonitorHandler( StepperTupleType& tuple )
      : tuple_( tuple )
    {}


    void stepPrint()
    {
      std::cout << ",  Newton: ";
      ForLoop< StepPrintNewtonIterations, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
      std::cout << "  ILS: ";
      ForLoop< StepPrintILSIterations, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
      std::cout << "  OC: ";
      ForLoop< StepPrintOperatorCalls, 0, sizeof ... ( StepperArg ) >::apply( tuple_ );
      std::cout << std::endl;
    }

    template< class TimeProviderImp >
    void step( TimeProviderImp& tp )
    {
      ForLoop< SetTimeStepInformation, 0, sizeof ... ( StepperArg ) >::apply( tuple_, tp );
    }

    void finalize( const double gridWidth, const double gridSize )
    {
      ForLoop< Finalize, 0, sizeof ... ( StepperArg ) >::apply( tuple_, gridWidth, gridSize );
    }

  private:
    StepperTupleType&  tuple_;
  };


  template<>
  class CombinedDefaultSolverMonitorHandler<>
  {
  public:
    template< class ... Args >
    CombinedDefaultSolverMonitorHandler ( Args && ... ) {}

    template <class ... Args>
    void stepPrint(Args&& ... ) const {};

    template <class ... Args>
    void step(Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};
  };


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

    template< class OdeSolverMonitorImp >
    void step( OdeSolverMonitorImp& odeSolverMonitor )
    {
      *solverMonitor_.newton_iterations     = odeSolverMonitor.newtonIterations_;
      *solverMonitor_.ils_iterations        = odeSolverMonitor.linearSolverIterations_;
      *solverMonitor_.max_newton_iterations = odeSolverMonitor.maxNewtonIterations_ ;
      *solverMonitor_.max_ils_iterations    = odeSolverMonitor.maxLinearSolverIterations_;
      *solverMonitor_.operator_calls        = odeSolverMonitor.spaceOperatorCalls_;
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
      *solverMonitor_.elements = gridSize;

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
    typedef SolverMonitor<1>     SolverMonitorType;

    NoSolverMonitorHandler( const std::string keyPrefix = "" )
      : monitor_()
    {}

    template <class ... Args>
    void stepPrint(Args&& ... ) const {};

    template <class ... Args>
    void step(Args&& ... ) const {};

    template <class ... Args>
    void finalize(Args&& ... ) const {};

    template <class ... Args>
    SolverMonitorType& monitor(Args&& ... ) { return monitor_; };

  private:
    SolverMonitorType monitor_;

  };

}
}

#endif
