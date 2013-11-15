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

template< class GridType > 
struct ProblemGenerator 
{
  typedef Dune :: U0< GridType > ProblemType;
  template< class GridPart >
  struct Traits
  {
    typedef ProblemType InitialDataType;
  };

  static inline std::string advectionFluxName()
  {
    return "noflux";
  }

  static inline std::string diffusionFluxName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>               
  initializeGrid( const std::string description )
  {
    // use default implementation 
    return initialize< GridType > ( description );
  }

  static ProblemType* problem()
  {
    // choice of explicit or implicit ode solver
    return new ProblemType ();
  }
};

// this implements the Stepper
#include "checkpointing.hh"

#endif // FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
