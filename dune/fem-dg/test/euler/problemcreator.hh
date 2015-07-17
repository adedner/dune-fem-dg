#ifndef FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
#define FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH

#include <config.h>

#ifndef POLORDER
#define POLORDER 1
#endif

#include <dune/fem/io/parameter.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// get a numerical flux
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

#include "eulermodel.hh"
#include "problems.hh"

#include <dune/fem-dg/stepper/advectionstepper.hh>

template <class GridType>
struct EulerProblemCreator
{
  typedef ProblemBase< GridType > ProblemType ;

  template <class GridPart>
  struct Traits
  {
    typedef ProblemType  InitialDataType;
    typedef Dune::Fem::EulerModel< GridPart, ProblemType > ModelType;

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
  typedef LLFFlux< ModelType > OriginalFluxType;
#elif (FLUX==2)
#warning "FLUX: HLL"
  typedef HLLNumFlux< ModelType > OriginalFluxType;
#elif (FLUX==3)
#warning "FLUX: HLLC"
  typedef HLLCNumFlux< ModelType > OriginalFluxType;
#else
#error "Use flag FLUX=1 for LLF or flag FLUX=2 for HLL or flag FLUX=3 for HLLC!"
#endif

  typedef OriginalFluxType FluxType;
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

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    std::string description( advectionFluxName() );
    // use default implementation
    return initialize< GridType >( description );
  }

  static ProblemType* problem()
  {
    const std::string problemNames []
        = { "sod" , "withman", "withmansmooth", "smooth1d" , "ffs" , "diffraction" , "shockbubble", "p123" };

    const int problemNumber = Dune :: Fem :: Parameter :: getEnum ( "problem" , problemNames );

    if( problemNames[ problemNumber ] == "sod" )
    {
      return new U0Sod< GridType > ( );
    }
    else if( problemNames[ problemNumber ] == "smooth1d" )
    {
      return new U0Smooth1D< GridType > ();
    }
    else if( problemNames[ problemNumber ] == "ffs" )
    {
      return new U0FFS< GridType > ();
    }
    /*
    else if( problemNames[ problemNumber ] == "sod" )
    {
      return new SodProblem< GridType > ();
    }
    else if( problemNames[ problemNumber ] == "sod" )
    {
      return new SodProblem< GridType > ();
    }
    */
    else if( problemNames[ problemNumber ] == "p123" )
    {
      return new U0P123< GridType >();
    }

    std::cerr << "Error: Problem " << problemNames[ problemNumber ]
              << " not implemented." << std::endl;

    // choice of explicit or implicit ode solver
    return new U0Smooth1D< GridType > ();
  }
  template <int polynomialOrder>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef AdvectionStepper< GridType, EulerProblemCreator< GridType> , polynomialOrder > Type;
  };
};

#ifndef COMBINED_PROBLEM_CREATOR
#define ProblemCreator EulerProblemCreator
#endif


#define NEW_STEPPER_SELECTOR_USED
#endif // FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH
