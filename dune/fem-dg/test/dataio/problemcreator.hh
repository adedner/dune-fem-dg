#ifndef FEMDG_CHECKPOINTING_RPOBLEMCREATOR_HH
#define FEMDG_CHECKPOINTING_RPOBLEMCREATOR_HH
#include <config.h>

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
struct ProblemGenerator 
{
  template <class Traits> 
  struct Stepper
  {
    typedef CheckPointingStepper< GridType, Traits, POLORDER > Type;
  };

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
};

#endif // FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
