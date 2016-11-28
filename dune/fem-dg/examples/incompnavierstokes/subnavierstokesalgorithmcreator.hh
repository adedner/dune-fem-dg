#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_SUB_NAVIER_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_SUB_NAVIER_HH
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
  struct SubNavierStokesAlgorithmCreator
  {

    typedef typename AC::GridType                         GridType;
    typedef typename AC::GridParts                        HostGridPartType;
    typedef HostGridPartType                              GridPartType;

    // define problem type here if interface should be avoided
    typedef typename NavierStokesProblemInterface< GridType >::NavierStokesProblemType ProblemInterfaceType;

    typedef typename ProblemInterfaceType::FunctionSpaceType             FunctionSpaceType;

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                       ProblemType;
      typedef ProblemInterfaceType                                       InitialDataType;
      typedef NavierStokesModel< GridPartType, InitialDataType, false >  ModelType;

      template< class Solution, class Model, class ExactFunction, class TimeProvider >
      static void addEOCErrors ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
      {
        static L2EOCError l2EocError( "$L^2$-Error");
        l2EocError.add( tp, u, model, f );
      }
    };

    static inline std::string moduleName() { return ""; }

    static ProblemInterfaceType* problem()
    {
      return new typename NavierStokesProblem< GridType, NavierStokesProblemDefault >::NavierStokesProblemType();
    }

    template< int polOrd >
    struct DiscreteTraits
    {
      typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                         DFSpaceType;
      typedef std::tuple<>                                                               ExtraParameterTuple;
    public:
      typedef typename AC::template DiscreteFunctions< DFSpaceType >                     DiscreteFunctionType;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                 IOTupleType;

      class Operator
      {
        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                         OpTraits;

        typedef typename AC::template RhsAnalyticalTraits< AnalyticalTraits, NavierStokesModel< GridPartType, typename AnalyticalTraits::InitialDataType, true > >
                                                                                         RhsAnalyticalTraitsType;
        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, RhsAnalyticalTraitsType, ExtraParameterTuple >
                                                                                         RhsOpTraits;

      public:
        typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::full >    type;
        typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::expl >    ExplicitType;
        typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::impl >    ImplicitType;

        typedef DGAdvectionDiffusionOperator< RhsOpTraits >                              RhsType;
      };

      struct Solver
      {
        typedef typename AC::template LinearSolvers< DFSpaceType >                       LinearSolverType;
        typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                      type;
      };

    private:
      typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                         OpTraits;
      typedef DGAdaptationIndicatorOperator< OpTraits >                                  IndicatorType;
      typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >  GradientIndicatorType ;
    public:

      typedef AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
      typedef SubSolverMonitor< SolverMonitor >                                          SolverMonitorType;
      typedef SubDiagnostics< Diagnostics >                                              DiagnosticsType;
      typedef ExactSolutionOutput< DiscreteFunctionType >                                AdditionalOutputType;
    };

    template <int polOrd>
    using Algorithm = SubAdvectionDiffusionAlgorithm< GridType, SubNavierStokesAlgorithmCreator<AC>, polOrd >;

  };

}
}



#endif
