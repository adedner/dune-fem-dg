#ifndef FEMDG_NAVIERSTOKESSTEPPER_HH
#define FEMDG_NAVIERSTOKESSTEPPER_HH
#include <config.h>

#ifndef GRIDDIM
#define GRIDDIM 2
#endif

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM + 2
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/l1eocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/configurator.hh>


namespace Dune
{
namespace Fem
{

  template< class GridSelectorGridType >
  struct NavierStokesAlgorithmCreator
  {
    typedef AlgorithmConfigurator< GridSelectorGridType,
                                   Galerkin::Enum::dg,
                                   Adaptivity::Enum::yes,
                                   DiscreteFunctionSpaces::Enum::legendre,
                                   //Solver::Enum::istl,
                                   Solver::Enum::fem,
                                   AdvectionLimiter::Enum::unlimited,
                                   Matrix::Enum::matrixfree,
                                   AdvectionFlux::Enum::llf,
                                   DiffusionFlux::Enum::primal > ACNavier;

    template< class AC >
    struct SubNavierStokesAlgorithmCreator
    {
      typedef typename AC::GridType                                 GridType;
      typedef typename AC::GridParts                                HostGridPartType;
      typedef HostGridPartType                                      GridPartType;

      // define problem type here if interface should be avoided
      typedef NSWaves< GridType >                                   ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType      FunctionSpaceType;

      typedef NSModel< GridType, ProblemInterfaceType >             ModelType;

      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem() { return new ProblemInterfaceType(); }


      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctions< FunctionSpaceType, polOrd >       DiscreteFunctionType;

        class Operator
        {
          typedef typename AC::template DefaultOpTraits< ModelType, FunctionSpaceType, polOrd >
                                                                                           OpTraits;
        public:
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::full >    type;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::expl >    ExplicitType;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::impl >    ImplicitType;
        };

        struct Solver
        {
          typedef typename AC::template LinearSolvers< DiscreteFunctionType >              LinearSolverType;
          typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                      type;
        };

      private:
        typedef typename AC::template DefaultOpTraits< ModelType, FunctionSpaceType, polOrd >
                                                                                           OpTraits;
        typedef DGAdaptationIndicatorOperator< OpTraits >                                  IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename ModelType::ProblemType >         GradientIndicatorType ;
      public:

        typedef AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
        typedef SubSolverMonitor< SolverMonitor >                                          SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                              DiagnosticsType;
        typedef SubDataWriter< SolutionOutput<DiscreteFunctionType>, ExactSolutionOutput<DiscreteFunctionType> >
                                                                                           DataWriterType;
      };

      template <int polOrd>
      using Algorithm = SubAdvectionDiffusionAlgorithm< GridType, SubNavierStokesAlgorithmCreator<AC>, polOrd >;

    };

    template <int polOrd>
    using Algorithm = EvolutionAlgorithm< polOrd, SubNavierStokesAlgorithmCreator<ACNavier> >;


    static inline std::string moduleName() { return ""; }

    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef typename SubNavierStokesAlgorithmCreator<ACNavier>::GridType GridType;
      //Discrete Functions
      typedef typename SubNavierStokesAlgorithmCreator<ACNavier>::template DiscreteTraits<polOrd>::DiscreteFunctionType
                                                                           DFType;

      //Item1
      typedef _t< SubEvolutionContainerItem >                              Steady;
      typedef std::tuple< Steady >                                         Item1TupleType;

      //Item2
      typedef _t< EmptyContainerItem >                                     Empty;
      typedef std::tuple< std::tuple< Empty > >                            Item2TupleType;


      //Sub (discrete function argument ordering)
      typedef std::tuple<__0 >                                             AdvDiffOrder;

      typedef std::tuple< AdvDiffOrder >                                   SubOrderRowType;
      typedef SubOrderRowType                                              SubOrderColType;

      //external params lists
      typedef ExtraArg<>                                                   ExtraType;

      //Global container
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, ExtraType, DFType >
                                                                           GlobalContainerType;

      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( moduleName(), gridptr );
    }
  };

}
}

#endif // FEMHOWTO_HEATSTEPPER_HH
