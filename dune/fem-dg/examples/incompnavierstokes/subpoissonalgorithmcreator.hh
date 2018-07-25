#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_SUB_POISSON_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_SUB_POISSON_HH
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
//--------- ADAPTATION -----------------------
#include <dune/fem-dg/operator/adaptation/poissonestimator.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- SUB-ALGORITHM -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
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

namespace Dune
{
namespace Fem
{

  template< class AC >
  struct SubPoissonAlgorithmCreator
  {

    typedef typename AC::GridType                         GridType;
    typedef typename AC::GridParts                        HostGridPartType;
    typedef HostGridPartType                              GridPartType;

    // define problem type here if interface should be avoided
    typedef typename NavierStokesProblemInterface< GridType >::PoissonProblemType  ProblemInterfaceType;

    typedef typename ProblemInterfaceType::FunctionSpaceType       FunctionSpaceType;

    typedef PoissonModel< GridPartType, ProblemInterfaceType, std::tuple<> > ModelType;

    static inline std::string moduleName() { return "";}

    static ProblemInterfaceType* problem()
    {
      return new typename NavierStokesProblem< GridType, NavierStokesProblemDefault >::PoissonProblemType();
    }

    template< int polOrd >
    struct DiscreteTraits
    {
      typedef typename AC::template DiscreteFunctions< FunctionSpaceType, polOrd > DiscreteFunctionType;

      typedef typename AC::template ExtraParameters< ModelType, std::tuple< AC,AC>, polOrd,polOrd >  ExtraParameters;

      class Operator
      {
        typedef typename AC::template DefaultAssembTraits< ModelType, FunctionSpaceType, polOrd, FunctionSpaceType, polOrd/*, ExtraParameters*/ >
                                                                                   OpTraits;

        typedef typename AC::template DefaultOpTraits< PoissonModel< GridPartType, typename ModelType::ProblemType, std::tuple</*__t*/> >,
                                                       FunctionSpaceType, polOrd, std::tuple<> /*ExtraParameters*/ >
                                                                                   RhsOpTraits;
      public:
        typedef typename AC::template Operators< OpTraits >                        AssemblerType;
        typedef typename AssemblerType::MatrixType                                 type;

        typedef DGAdvectionDiffusionOperator< RhsOpTraits >                        RhsType;
      };

      struct Solver
      {
        typedef typename AC::template LinearSolvers< DiscreteFunctionType, true >  type;
      };

      typedef typename AC::template WrappedDiscreteFunctions< SigmaDiscreteFunctionSelector, DiscreteFunctionType >
                                                                                   SigmaDiscreteFunctionType;

      typedef ErrorEstimator< DiscreteFunctionType, SigmaDiscreteFunctionType, typename Operator::AssemblerType >
                                                                                   ErrorEstimatorType;
      typedef PoissonSigmaEstimator< ErrorEstimatorType >                          SigmaEstimatorType;

      typedef PAdaptivity<typename DiscreteFunctionType::DiscreteFunctionSpaceType, polOrd, SigmaEstimatorType >
                                                                                   PAdaptivityType;

      typedef PAdaptIndicator< PAdaptivityType, ModelType >                        AdaptIndicatorType;
      // typedef NoPAdaptIndicator                                                   AdaptIndicatorType;

      typedef SubSolverMonitor< SolverMonitor >                                    SolverMonitorType;
      typedef SubDiagnostics< Diagnostics >                                        DiagnosticsType;
      typedef SubDataWriter< SolutionOutput<DiscreteFunctionType>, ExactSolutionOutput<DiscreteFunctionType> >
                                                                                   DataWriterType;
    };

    template <int polOrd>
    using Algorithm = SubEllipticAlgorithm< GridType, SubPoissonAlgorithmCreator<AC>, polOrd >;
  };

}
}



#endif
