#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

// dune-fem includes
#include <dune/fem/io/parameter.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>

// local includes
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

// overload default stepper traits
#include <dune/fem-dg/stepper/advectiondiffusionstepper.hh>

#include "problem.hh"
#include "nsmodel.hh"
#include "stokesmodel.hh"

template< class GridType >
struct IncompressibleNavierStokesProblemCreator
{
  typedef Dune :: EvolutionProblemInterface<
                      Dune::Fem::FunctionSpace< double, double, GridType::dimension,
                      DIMRANGE>,
                      false > ProblemType;
  // define problem type here if interface should be avoided

  template< class GridPart >
  struct Traits
  {
    typedef ProblemType                                     InitialDataType;
    typedef NavierStokesModel< GridPart, InitialDataType >  ModelType;

    //typedef LLFFlux< ModelType >                      FluxType;
    typedef UpwindFlux< ModelType > FluxType;

    // choice of diffusion flux (see diffusionflux.hh for methods)
     static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId
       =  Dune :: method_general ;
  };


  static inline std::string advectionFluxName()
  {
    return "LLF";
  }


  static inline std::string diffusionFluxName()
  {
#ifdef EULER
    return "";
#elif (defined PRIMALDG)
    return Dune::Fem::Parameter::getValue< std::string >("dgdiffusionflux.method");
#else
    return "LDG";
#endif
  }

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    std::string description( advectionFluxName() + " " + diffusionFluxName() );
    // use default implementation
    return initialize< GridType > ( description );
  }

  static ProblemType* problem()
  {
    return new NavierStokesProblemDefault< GridType > ();
  }

  template <int polynomialOrder>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef AdvectionDiffusionStepper< GridType, IncompressibleNavierStokesProblemCreator<GridType>, polynomialOrder > Type;
  };

};

#ifndef COMBINED_PROBLEM_CREATOR
#define ProblemCreator IncompressibleNavierStokesProblemCreator
#endif

#define NEW_STEPPER_SELECTOR_USED
#endif
