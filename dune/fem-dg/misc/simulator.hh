#ifndef DUNE_FEM_DG_SIMULATOR_HH
#define DUNE_FEM_DG_SIMULATOR_HH

// in dbug mode also enable FieldVector checking and dune devel mode
#ifndef NDEBUG
#define DUNE_ISTL_WITH_CHECKING
#define DUNE_DEVEL_MODE
#endif

// -1 means higher order FV
#if POLORDER == -1
#define HIGHER_ORDER_FV
#undef POLORDER
#define POLORDER 0
#endif

#ifdef ONLY_ONE_P
#define MIN_POLORD POLORDER
#define MAX_POLORD POLORDER
#else
#define MIN_POLORD 1
#define MAX_POLORD 4
#endif

#include <memory>
#include <dune/common/version.hh>

// streams for backup
#include <dune/fem-dg/misc/streams.hh>

#if defined USE_BASEFUNCTIONSET_CODEGEN || defined BASEFUNCTIONSET_CODEGEN_GENERATE
#define USE_FEMDG_BASISFUNCTIONSET
#endif

#ifdef USE_FEMDG_BASISFUNCTIONSET
#include <dune/fem-dg/main/default.hh>
#ifdef NEWBASEFCT_CACHING
#include <dune/fem-dg/main/codegen2.hh>
#else
#include <dune/fem-dg/main/codegen.hh>
#endif
#endif

#include <dune/fem-dg/pass/threadhandle.hh>
#if defined USE_SMP_PARALLEL
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem-dg/pass/threadpass.hh>

#if HAVE_DUNE_FEM_DG
#define NSMOD_USE_SMP_PARALLEL
#endif
#endif

#ifdef NEWBASEFCT_CACHING
#include <dune/fem-dg/main/codegen2.hh>
#else
#include <dune/fem-dg/main/codegen.hh>
#endif

#include <dune/fem-dg/stepper/base.hh>

#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/flops.hh>

namespace Dune
{
  struct FlopStartObject
  {
    FlopStartObject()
    {
      // initialize counters for master thread before all others
      runThread() ;
    }
    void runThread() const
    {
      Dune::Fem::FlopCounter::start();
    }
  };

  struct FlopStopObject
  {
    void runThread() const
    {
      Dune::Fem::FlopCounter::stop();
    }
  };



template <int polynomialOrder, class ProblemTraits>
inline void simulate(const ProblemTraits& problem)
{
  int polOrder = 1;
  polOrder = Dune::Fem::Parameter :: getValue("femdg.polynomialOrder", polOrder );

  // if polOrder and polynomialOrder differ don't do anything.
  if( polOrder != polynomialOrder ) return ;

  // get number of desired threads (default is 1)
  const int numThreads = Dune::Fem::Parameter::getValue< int >("fem.parallel.numberofthreads", 1);
  Dune :: Fem :: ThreadManager :: setMaxNumberThreads( numThreads );

  Dune::Fem::FemEoc::clear();

  const bool countFlops = Dune::Fem::Parameter::getValue< bool >("femdg.flopcounter", false );

  // if flop count is enabled count floating point operations (PAPI needed)
  // start flop counters for all threads
  if( countFlops )
  {
    FlopStartObject startObj ;
    Dune::Fem::ThreadHandle::run( startObj );
  }

  typedef Dune::GridSelector :: GridType GridType;

  // typedef ProblemCreator< GridType, polynomialOrder > ProblemTraits;

  // return type of initializeGrid is Dune::GridPtr, use release such that memory of GridPtr is released
  std::unique_ptr< GridType > gridptr( problem.initializeGrid().release() );

  typedef typename ProblemTraits :: template Stepper< polynomialOrder > :: Type StepperType;
  std::unique_ptr< StepperType > stepper( new StepperType( *gridptr, problem.moduleName() ) );

  // new method, the ProblemGenerator simply creates the stepper
  compute( *stepper );

  // stop flop counters for all threads
  if( countFlops )
  {
    FlopStopObject stopObj ;
    Dune::Fem::ThreadHandle::run( stopObj );
    // print results
    Dune::Fem::FlopCounter::print( std::cout );
  }
} // end simulate

  template <int polOrd>
  struct SimulatePolOrd
  {
    template <class ProblemTraits>
    static void apply( const ProblemTraits& problem )
    {
      simulate< polOrd > ( problem );
    }
  };

  struct Simulator
  {
    template <class ProblemTraits>
    static void run( const ProblemTraits& problem )
    {
      Dune::ForLoop< SimulatePolOrd, MIN_POLORD, MAX_POLORD > :: apply( problem );
    }
  };

} // namespace Dune
#endif // #ifndef DUNE_FEM_DG_SIMULATOR_HH
