#ifndef FEMHOWTO_HEATSTEPPER_HH
#define FEMHOWTO_HEATSTEPPER_HH
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

#include "problems/problem.hh"
#include "problems/problemQuasiHeatEqn.hh"
#include "problems/pulse.hh"
#include "problems/sin.hh"
#include "problems/deformationalflow.hh"
#include "models.hh"


template< class GridType >
struct ProblemCreator
{
  typedef Dune :: EvolutionProblemInterface<
                      Dune::Fem::FunctionSpace< double, double, GridType::dimension,
                      DIMRANGE>,
                      false > ProblemType;
  // define problem type here if interface should be avoided

  template< class GridPart >
  struct Traits
  {
    typedef ProblemType  InitialDataType;
    typedef HeatEqnModel< GridPart, InitialDataType > ModelType;

    // TODO: flux does not work yet
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
    // choice of explicit or implicit ode solver
    static const std::string probString[]  = { "heat" ,"quasi", "pulse", "sin" };
    const int probNr = Dune::Fem::Parameter::getEnum( "problem", probString, 0 );
    if( probNr == 0 )
      return new Dune :: U0< GridType, DIMRANGE > ();
    else if ( probNr == 1 )
      return new Dune :: QuasiHeatEqnSolution< GridType, DIMRANGE > ();
    else if ( probNr == 2 )
      return new Dune :: Pulse< GridType, DIMRANGE > ();
    else if ( probNr == 3 )
      return new Dune :: U0Sin< GridType, DIMRANGE > ();
    else
    {
      abort();
      return 0;
    }
  }

  // this should be ok but could lead to a henn-egg problem
  typedef AdvectionDiffusionStepper< GridType, ProblemCreator<GridType>, POLORDER > StepperType;
};

#define NEW_STEPPER_SELECTOR_USED
#endif // FEMHOWTO_HEATSTEPPER_HH
