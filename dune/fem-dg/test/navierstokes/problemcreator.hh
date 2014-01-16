#ifndef FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
#define FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
#include <config.h>

// system includes
#include <string>

#include "passtraits.hh"

// dune includes
#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>

// local includes
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>
#include "problemtype.hh"
#include "ns_model.hh"

#include <dune/fem-dg/stepper/base.hh>


template< class GridType > 
struct ProblemGenerator 
{
  typedef NSProblemType ProblemType;

  template< class GridPart >
  struct Traits
  {
    typedef ProblemType InitialDataType;
    typedef Dune::NSModel< GridPart, InitialDataType > ModelType;
    // choice of diffusion flux (see diffusionflux.hh for methods)
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId 
           = Dune :: method_general ;

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

#endif // FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
