#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_SUB_STOKES_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_SUB_STOKES_HH
#include <config.h>

#ifndef GRIDDIM
#define GRIDDIM 2
#endif

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- SUB-ALGORITHMS -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/examples/stokes/stokesalgorithm.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "stokesmodel.hh"

//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/configurator.hh>

#include "subpoissonalgorithmcreator.hh"

namespace Dune
{
namespace Fem
{

  template< class AC >
  struct SubStokesAlgorithmCreator
  {
    //only for inner creators needed
    typedef SubPoissonAlgorithmCreator<AC>                               SubCreatorType;

    typedef typename AC::GridType                                        GridType;
    typedef typename AC::GridParts                                       HostGridPartType;
    typedef HostGridPartType                                             GridPartType;

    typedef NavierStokesProblemInterface< GridType >                     ProblemInterfaceType;

    typedef typename ProblemInterfaceType::StokesProblemType::FunctionSpaceType
                                                                         FunctionSpaceType;

    typedef StokesModel< GridPartType, ProblemInterfaceType, std::tuple<> > ModelType;

    static inline std::string moduleName() { return "";}

    static ProblemInterfaceType* problem()
    {
      return new NavierStokesProblem< GridType, NavierStokesProblemDefault > ();
    }

    template< int polOrd >
    struct DiscreteTraits
    {
      static const int redPolOrd = polOrd-0;
    public:
      typedef typename AC::template DiscreteFunctions< FunctionSpaceType, redPolOrd >   DiscreteFunctionType;

      typedef typename AC::template ExtraParameters< ModelType, std::tuple< AC,AC>, polOrd,polOrd >  ExtraParameters;

      class Operator
      {
        typedef typename AC::template DefaultAssembTraits< ModelType, FunctionSpaceType, redPolOrd, FunctionSpaceType, redPolOrd /*,ExtraParameters*/ >
                                                                                        OpTraits;
      public:
        typedef StokesAssembler< OpTraits, AC::template Containers,
                                   typename SubCreatorType::template DiscreteTraits< polOrd >::DiscreteFunctionType,
                                   DiscreteFunctionType >                               AssemblerType;
        //typedef StokesAssembler< AssTraits, OpTraits >                             AssemblerType;
        //the following typedef is not needed by stokes algorithm atm
        //typedef typename AssemblerType::MatrixType                               type;

        //typedef DGAdvectionDiffusionOperator< RhsOpTraits >                        RhsType;
      };

      struct Solver
      {
        typedef UzawaSolver< typename Operator::AssemblerType,typename SubCreatorType::template Algorithm<polOrd> >
                                                                                        type;
      };

      static_assert( (int)FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

    private:
      typedef typename SubCreatorType::template DiscreteTraits<polOrd>::ErrorEstimatorType  PoissonErrorEstimatorType;
      typedef typename SubCreatorType::template DiscreteTraits<polOrd>::PAdaptivityType     PoissonPAdaptivityType;

      typedef StokesSigmaEstimator< PoissonErrorEstimatorType, typename Operator::AssemblerType > StokesSigmaEstimatorType;
      typedef typename StokesSigmaEstimatorType::StokesErrorEstimatorType                   StokesErrorEstimatorType;

      typedef StokesPAdaptivity< PoissonPAdaptivityType,typename DiscreteFunctionType::DiscreteFunctionSpaceType,
                                 polOrd, StokesSigmaEstimatorType >                         StokesPAdaptivityType;
    public:

      typedef SubSolverMonitor< SolverMonitor >                                             SolverMonitorType;
      typedef SubDiagnostics< Diagnostics >                                                 DiagnosticsType;
      typedef StokesPAdaptIndicator< StokesPAdaptivityType, ModelType >                     AdaptIndicatorType;
      typedef SubDataWriter< SolutionOutput<DiscreteFunctionType>, ExactSolutionOutput<DiscreteFunctionType> >
                                                                                     DataWriterType;
    };

    template <int polOrd>
    using Algorithm = SubStokesAlgorithm< GridType, SubStokesAlgorithmCreator<AC>, SubCreatorType, polOrd >;
  };

}
}

#endif
