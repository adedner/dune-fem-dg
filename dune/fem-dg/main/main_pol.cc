#if defined GRIDDIM
#ifndef CODEDIM
#define CODEDIM GRIDDIM
#endif
#endif

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

#include <config.h>

#include <memory>
#include <dune/common/version.hh>

#include "main.hh"

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

// problem dependent
#include <problemcreator.hh>

#ifndef NEW_STEPPER_SELECTOR_USED
#warning "Old method of stepper selection detected. Please consult dune/fem-dg/test and update your code"
#if defined EULER
#include <dune/fem-dg/stepper/advectionstepper.hh>
#define Stepper AdvectionStepper
#else
#include <dune/fem-dg/stepper/advectiondiffusionstepper.hh>
#define Stepper AdvectionDiffusionStepper
#endif
#endif

#include <dune/fem/misc/femeoc.hh>

#include <dune/fem/misc/flops.hh>

#if POLORDER == 0
#define LOOPSPACE DG_P0
#elif POLORDER == 1
#define LOOPSPACE DG_P1
#elif POLORDER == 2
#define LOOPSPACE DG_P2
#elif POLORDER == 3
#define LOOPSPACE DG_P3
#elif POLORDER == 4
#define LOOPSPACE DG_P4
#elif POLORDER == 5
#define LOOPSPACE DG_P5
#elif POLORDER == 6
#define LOOPSPACE DG_P6
#elif POLORDER == 7
#define LOOPSPACE DG_P7
#elif POLORDER == 8
#define LOOPSPACE DG_P8
#endif

namespace LOOPSPACE {

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

  void simulate()
  {
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

    // old method
#ifndef NEW_STEPPER_SELECTOR_USED
    typedef ProblemGenerator< GridType > ProblemTraits;

    // Note to me: problem description is for FemEOC
    const std::string advFlux  = ProblemTraits :: advectionFluxName();
    const std::string diffFlux = ProblemTraits :: diffusionFluxName();

    // return type of initializeGrid is Dune::GridPtr, use release such that memory of GridPtr is released
    std::unique_ptr< GridType > gridptr( ProblemTraits :: initializeGrid( advFlux + diffFlux ).release() );

    // the old legacy way
    typedef Stepper< GridType, ProblemTraits, POLORDER > StepperType ;
    StepperType stepper( *gridptr );

    // new method, the ProblemGenerator simply creates the stepper
    compute( stepper );
#else
    typedef ProblemCreator< GridType > ProblemTraits;

    // return type of initializeGrid is Dune::GridPtr, use release such that memory of GridPtr is released
    std::unique_ptr< GridType > gridptr( ProblemTraits :: initializeGrid().release() );

    typedef ProblemTraits :: StepperType StepperType;
    std::unique_ptr< StepperType > stepper( new StepperType( *gridptr, ProblemTraits::moduleName() ) );

    // new method, the ProblemGenerator simply creates the stepper
    compute( *stepper );
#endif

    // stop flop counters for all threads
    if( countFlops )
    {
      FlopStopObject stopObj ;
      Dune::Fem::ThreadHandle::run( stopObj );
      // print results
      Dune::Fem::FlopCounter::print( std::cout );
    }
  }

} // end namespace LOOPSPACE
