#ifndef FEMDG_CHECKPOINTING_RPOBLEMCREATOR_HH
#define FEMDG_CHECKPOINTING_RPOBLEMCREATOR_HH
#include <config.h>

#ifndef POLORDER
#define POLORDER 2
#endif

// system includes
#include <string>

// dune includes
#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/stepper/base.hh>

// local includes
#include "problem.hh"
// this implements the Stepper
#include "checkpointing.hh"

template< class GridType >
struct ProblemCreator
{
  typedef Dune :: U0< GridType > ProblemType;
  template< class GridPart >
  struct Traits
  {
    typedef ProblemType InitialDataType;
  };

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    std::string description ("noflux");
    // use default implementation
    return initialize< GridType > ( description );
  }

  static ProblemType* problem()
  {
    // choice of explicit or implicit ode solver
    return new ProblemType ();
  }

  // this should be ok but could lead to a henn-egg problem
  typedef CheckPointingStepper< GridType, ProblemCreator<GridType>, POLORDER > StepperType;
};

#define NEW_STEPPER_SELECTOR_USED
#endif // FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
