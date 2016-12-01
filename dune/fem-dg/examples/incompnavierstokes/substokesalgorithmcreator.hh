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


#include "subpoissonalgorithmcreator.hh"

namespace Dune
{
namespace Fem
{

  template< class AC >
  struct SubStokesAlgorithmCreator
  {
    //only for inner creators needed
    typedef SubPoissonAlgorithmCreator<AC>                              SubCreatorType;

    typedef typename AC::GridType                                       GridType;
    typedef typename AC::GridParts                                      HostGridPartType;
    typedef HostGridPartType                                            GridPartType;

    typedef NavierStokesProblemInterface< GridType >                    ProblemInterfaceType;

    typedef typename ProblemInterfaceType::StokesProblemType::FunctionSpaceType
                                                                        FunctionSpaceType;

    struct AnalyticalTraits
    {
      typedef ProblemInterfaceType                                      ProblemType;
      typedef ProblemInterfaceType                                      InitialDataType;
      typedef StokesModel< GridPartType, InitialDataType, false >       ModelType;

      template< class Solution, class Model, class ExactFunction >
      static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f )
      {}
    };

    static inline std::string moduleName() { return "";}

    static ProblemInterfaceType* problem()
    {
      return new NavierStokesProblem< GridType, NavierStokesProblemDefault > ();
    }

    template< int polOrd >
    struct DiscreteTraits
    {
    private:
      typedef typename SubCreatorType::template DiscreteTraits< polOrd >           PoissonDiscreteTraits;
      typedef typename PoissonDiscreteTraits::DiscreteFunctionType                 VelDiscreteFunctionType;
      typedef typename SubCreatorType::FunctionSpaceType                           VelFunctionSpaceType;
      typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                   DFSpaceType;
      typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, VelFunctionSpaceType>
                                                                                   VelDFSpaceType;
    public:
      typedef typename AC::template DiscreteFunctions< DFSpaceType >               DiscreteFunctionType;

      typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >           IOTupleType;

      class Operator
      {
        typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd, AnalyticalTraits >
                                                                                   OpTraits;

        typedef AssemblerTraitsList< std::tuple< VelDiscreteFunctionType, DiscreteFunctionType >, AC::template Containers > AssTraits;

        typedef typename AC::template RhsAnalyticalTraits< AnalyticalTraits, StokesModel< GridPartType, typename AnalyticalTraits::InitialDataType, true > >
                                                                              RhsAnalyticalTraitsType;
        typedef typename AC::template DefaultOpTraits< VelDFSpaceType, polOrd, RhsAnalyticalTraitsType >
                                                                              RhsOpTraits;
      public:
        typedef StokesAssembler< OpTraits, AC::template Containers,VelDiscreteFunctionType,DiscreteFunctionType >
                                                                                   AssemblerType;
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

      static_assert( (int)DFSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

    private:
      typedef typename SubCreatorType::template DiscreteTraits<polOrd>::ErrorEstimatorType  PoissonErrorEstimatorType;
      typedef typename SubCreatorType::template DiscreteTraits<polOrd>::SigmaEstimatorType  PoissonSigmaEstimatorType;
      typedef typename SubCreatorType::template DiscreteTraits<polOrd>::PAdaptivityType     PoissonPAdaptivityType;

      typedef StokesSigmaEstimator< PoissonErrorEstimatorType, typename Operator::AssemblerType > StokesSigmaEstimatorType;
      typedef typename StokesSigmaEstimatorType::StokesErrorEstimatorType             StokesErrorEstimatorType;

      typedef StokesPAdaptivity< PoissonPAdaptivityType,
                                 DFSpaceType, polOrd, StokesSigmaEstimatorType >       StokesPAdaptivityType;
    public:

      typedef SubSolverMonitor< SolverMonitor >                                     SolverMonitorType;
      typedef SubDiagnostics< Diagnostics >                                         DiagnosticsType;
      typedef StokesPAdaptIndicator< StokesPAdaptivityType, ProblemInterfaceType >    AdaptIndicatorType;
    };

    template <int polOrd>
    using Algorithm = SubStokesAlgorithm< GridType, SubStokesAlgorithmCreator<AC>, SubCreatorType, polOrd >;
  };

}
}

#endif
