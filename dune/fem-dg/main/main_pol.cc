#if defined GRIDDIM
#ifndef CODEDIM
#define CODEDIM GRIDDIM
#endif
#endif

// in dbug mode also enable FieldVector checking and dune devel mode
#ifndef DNDEBUG 
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

#include <dune/common/version.hh>

#include <dune/fem/misc/threadmanager.hh>
#include <dune/fem-dg/pass/partitioner.hh>
#include <dune/fem-dg/pass/threadpass.hh>

#if defined USE_SMP_PARALLEL 
#if HAVE_DUNE_FEM_DG 
#define NSMOD_USE_SMP_PARALLEL
#endif
#endif

#include <dune/fem-dg/main/codegen.hh>

#include <dune/fem-dg/stepper/base.hh>

// problem dependent
#include <problemcreator.hh>

#if defined NS_ELLIPTIC_OPERATOR 
#include <dune/fem-dg/stepper/ellipt.hh>
#elif defined EULER 
#include <dune/fem-dg/stepper/advectionstepper.hh>
#else
#include <dune/fem-dg/stepper/stepper.hh>
#endif

#include <dune/fem/space/common/allgeomtypes.hh>

#include <dune/grid/io/visual/grapedatadisplay.hh>

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

  void simulate()
  {
    FemEoc::clear();

    typedef Dune::GridSelector :: GridType GridType;
    typedef ProblemGenerator< GridType > ProblemTraits;

    // ProblemType is a Dune::Function that evaluates to $u_0$ and also has a
    // method that gives you the exact solution.
    //typedef NSProblemType< GridType > ProblemType;
    //ProblemType problem;

    // Note to me: problem description is for FemEOC
    const std::string advFlux  = ProblemTraits :: advectionFluxName();
    const std::string diffFlux = ProblemTraits :: diffusionFluxName();

    // use problem specific initialize method since some problems do different things
    // there, e.g. poisson 
    Dune::GridPtr<GridType> gridptr = ProblemTraits :: initializeGrid( advFlux + diffFlux );

    // get grid reference 
    GridType & grid = *gridptr;

#if defined NS_ELLIPTIC_OPERATOR 
    EllipticAlgorithm<GridType, 
                      ProblemTraits, 
                      POLORDER> stepper( grid );
#else 
    Stepper<GridType, 
            ProblemTraits, 
            POLORDER> stepper( grid );
#endif

    compute( stepper );
  } 

} // end namespace LOOPSPACE
