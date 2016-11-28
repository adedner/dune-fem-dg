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
#include <dune/fem/misc/l2norm.hh>
//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>
//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/examples/stokes/stokesalgorithm.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
#include "incompnavierstokesalgorithm.hh"
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//#include "../stokes/problems.hh"
//--------- MODELS --------------------------
#include "stokesmodel.hh"
#include "models.hh"
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

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                 ProblemType;
      typedef ProblemInterfaceType                                 InitialDataType;
      typedef PoissonModel< GridPartType, InitialDataType, false > ModelType;

      template< class Solution, class Model, class ExactFunction, class SigmaFunction>
      static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
      {}
    };

    static inline std::string moduleName() { return "";}

    static ProblemInterfaceType* problem()
    {
      return new typename NavierStokesProblem< GridType, NavierStokesProblemDefault >::PoissonProblemType();
    }

    template< int polOrd >
    struct DiscreteTraits
    {
      typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                              DFSpaceType;
    public:
      typedef typename AC::template DiscreteFunctions< DFSpaceType >          DiscreteFunctionType;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >      IOTupleType;
      typedef std::tuple<>                                                    ExtraParameterTuple;

      class Operator
      {
        typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd, AnalyticalTraits >
                                                                              OpTraits;
        typedef typename AC::template RhsAnalyticalTraits< AnalyticalTraits, PoissonModel< GridPartType, typename AnalyticalTraits::InitialDataType, true > >
                                                                              RhsAnalyticalTraitsType;

        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, RhsAnalyticalTraitsType, ExtraParameterTuple >
                                                                              RhsOpTraits;
      public:
        typedef typename AC::template Operators< OpTraits >                   AssemblerType;
        typedef typename AssemblerType::MatrixType                            type;

        typedef DGAdvectionDiffusionOperator< RhsOpTraits >                   RhsType;
      };

      struct Solver
      {
        typedef typename AC::template LinearSolvers< DFSpaceType, true >      type;
      };

    private:
      //small helper class
      template< class SigmaDFSpaceType > struct SigmaFunctionChooser
      { typedef typename AC::template DiscreteFunctions< SigmaDFSpaceType > type; };
    public:
      typedef typename SigmaDiscreteFunctionSelector< DiscreteFunctionType, SigmaFunctionChooser >::type SigmaDiscreteFunctionType;

      typedef ErrorEstimator< DiscreteFunctionType, SigmaDiscreteFunctionType, typename Operator::AssemblerType >
                                                                              ErrorEstimatorType;
      typedef PoissonSigmaEstimator< ErrorEstimatorType >                     SigmaEstimatorType;

      typedef PAdaptivity<DFSpaceType, polOrd, SigmaEstimatorType >           PAdaptivityType;

      typedef PAdaptIndicator< PAdaptivityType, ProblemInterfaceType >        AdaptIndicatorType;
      // typedef NoPAdaptIndicator                                               AdaptIndicatorType;

      typedef SubSolverMonitor< SolverMonitor >                               SolverMonitorType;
      typedef SubDiagnostics< Diagnostics >                                   DiagnosticsType;
    };

    template <int polOrd>
    using Algorithm = SubEllipticAlgorithm< GridType, SubPoissonAlgorithmCreator<AC>, polOrd >;
  };

}
}



#endif
