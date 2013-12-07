#ifndef FEMHOWTO_HEATSTEPPER_HH
#define FEMHOWTO_HEATSTEPPER_HH
#include <config.h>

// dune-fem includes
#include <dune/fem/io/parameter.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>

// local includes
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>
#include <dune/fem-dg/stepper/base.hh>

#include "problem.hh"
#include "problemQuasiHeatEqn.hh"
#include "deformationalflow.hh"
#include "models.hh"


template< class GridType > 
struct ProblemGenerator 
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
    typedef LLFFlux< ModelType > FluxType;

    //typedef UpwindFlux< ModelType > FluxType;
    // choice of diffusion flux (see diffusionflux.hh for methods)
     static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId 
       =  Dune :: method_general ;
  };


  static inline std::string advectionFluxName()
  {
#if (FLUX==1)
    return "LLF";
#elif (FLUX==2)
    return "HLL(Dennis)";
#elif (FLUX==3)
    return "HLLC(Dennis)";
#elif (FLUX==4)
    return "HLL2C";
#elif (FLUX==5)
    return "HLL2";
#elif (FLUX==6)
    return "HLLEM(Mhd)";
#endif
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


  static inline Dune::GridPtr<GridType> 
  initializeGrid( const std::string description )
  {
    // use default implementation 
    return initialize< GridType > ( description );
  }

  static ProblemType* problem()
  {
    // choice of explicit or implicit ode solver
    //static const std::string probString[]  = { "heat" ,"quasi", "deformational" };
    //int probNr = Dune::Fem::Parameter::getEnum( "femhowto.problem", probString, 0 );
    //if( probNr == 0 ) 
    //  return new Dune :: U0< GridType > ();
    //else if ( probNr == 1 ) 
    //  return new Dune :: QuasiHeatEqnSolution< GridType > ();
    //else if ( probNr == 2 ) 
      return new Dune :: DeformationalFlow< GridType > ();
    //else 
    //{
    //  abort();
    //  return 0;
    //}
  }
};

#include "steppertraits.hh"

#if ADVECTION && DIFFUSION 
#include <dune/fem-dg/stepper/advectiondiffusionstepper.hh>
#else 
#include <dune/fem-dg/stepper/advectionstepper.hh>
#endif

#endif // FEMHOWTO_HEATSTEPPER_HH
