#ifndef FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
#define FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH

#ifdef SPGRID_COUNT_FLOPS
//#define COUNT_FLOPS
#warning "FLOP/s counting enabled!"
#endif

#include <config.h>

#ifndef POLORDER
#define POLORDER 1
#endif

// system includes
#include <string>

#include "passtraits.hh"

// dune includes
#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>

// local includes
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>
#include "nswaves.hh"
#include "ns_model.hh"

#include <dune/fem-dg/stepper/advectiondiffusionstepper.hh>


template< class GridType >
struct ProblemCreator
{
  static const int polynomialOrder =  POLORDER ;
  typedef NSWaves< GridType > ProblemType;

  template< class GridPart >
  struct Traits
  {
    typedef ProblemType      InitialDataType;
    typedef Dune::NSModel< GridPart, InitialDataType > ModelType;
    // choice of diffusion flux (see diffusionflux.hh for methods)
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId
           = Dune :: method_general ;

// for header check
#ifndef FLUX
#define FLUX 1
#endif

// ******************************** NUMERICAL FLUX *****************************
#if (FLUX==1)
#warning "FLUX: LLF"
    typedef LLFFlux< ModelType > FluxType;
#elif (FLUX==2)
#warning "FLUX: HLL (Dennis)"
    typedef HLLNumFlux< ModelType > FluxType;
#elif (FLUX==3)
#warning "FLUX: HLLC (Dennis)"
    typedef HLLCNumFlux< ModelType > FluxType;
#elif (FLUX==4)
#warning "FLUX: HLL2C"
    typedef HLL2CFlux< ModelType > FluxType;
#elif (FLUX==5)
#warning "FLUX: HLL2"
    typedef HLL2Flux< ModelType > FluxType;
#elif (FLUX==6)
#warning "FLUX: HLLEM (Mhd)"
    typedef HLLEMNumFlux< ModelType > FluxType;
#else
#error "Set the flag FLUX! See Makefile.am for details!"
#endif

// ****************************** END NUMERICAL FLUX ***************************
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
    return new ProblemType ();
  }

  // this should be ok but could lead to a henn-egg problem
  typedef AdvectionDiffusionStepper< GridType, ProblemCreator< GridType>, polynomialOrder > StepperType;
};

#define NEW_STEPPER_SELECTOR_USED
#endif // FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
