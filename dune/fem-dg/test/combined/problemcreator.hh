#ifndef COMBINED_STEPPERCREATOR_HH
#define COMBINED_STEPPERCREATOR_HH
#include <config.h>

//fluxes
#include <dune/fem-dg/operator/fluxes/noflux.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

//discrete models
#include <dune/fem-dg/operator/dg/discretemodelcommon.hh>

// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

//stepper
#include "combinedstepper.hh"

//problem creator
#include <dune/fem-dg/test/advdiff/problemcreator.hh>
#include <dune/fem-dg/test/navierstokes/problemcreator.hh>


template< class GridType >
struct ProblemCreator
{
  typedef ProblemCreator< GridType > ProblemCreator1;
  typedef ProblemCreator< GridType > ProblemCreator2;

  //=======================================================================================================
  // inherited class from AlgorithmBase
  typedef CombinedStepper< ProblemCreator1, ProblemCreator2 > StepperType;
  //========================================================================================================


  static inline std::string moduleName()
  {
    return "atherosclerosis";
  }


  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    // use default implementation
    return initialize< GridType > ( moduleName() );
    // TODO do some crazy stuff, for example call gridgen
    // or read other important grid helper files
  }

};

#define NEW_STEPPER_SELECTOR_USED

#endif // FEMHOWTO_HEATSTEPPER_HH
