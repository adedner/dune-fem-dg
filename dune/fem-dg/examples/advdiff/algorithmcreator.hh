#ifndef FEMHOWTO_HEATSTEPPER_HH
#define FEMHOWTO_HEATSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>
#include <dune/fem-dg/algorithm/monitor.hh>

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

//undef the following to prescribe an external discrete velocity.
//This will not produce any physical results, yet. Just an example
//#define VELO

  /**
   *  \brief problem creator for an advection diffusion problem
   */
  template< class GridSelectorGridType >
  struct AdvectionDiffusionAlgorithmCreator
  {
    typedef AlgorithmConfigurator< GridSelectorGridType,
                                   Galerkin::Enum::default_,
                                   Adaptivity::Enum::default_,
                                   DiscreteFunctionSpaces::Enum::default_, //legendre
                                   Solver::Enum::default_,
                                   AdvectionLimiter::Enum::default_,
                                   Matrix::Enum::default_,
                                   AdvectionFlux::Enum::upwind,
                                   DiffusionFlux::Enum::primal > ACAdvDiff;
                                   // DiffusionFlux::Enum::local > ACAdvDiff;

    template< class AC >
    struct SubAdvectionDiffusionAlgorithmCreator
    {
      typedef typename AC::GridType                         GridType;
      typedef typename AC::GridParts                        HostGridPartType;
      typedef HostGridPartType                              GridPartType;

      // define problem type here if interface should be avoided
      typedef EvolutionProblemInterface< typename AC::template FunctionSpaces<DIMRANGE> >
                                                                        ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType          FunctionSpaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                    ProblemType;
        typedef ProblemInterfaceType                                    InitialDataType;
        typedef HeatEqnModel< GridType, InitialDataType, std::tuple<
#ifdef VELO
          __t
#endif
          > >                                                           ModelType;

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
        return AnalyticalAdvDiffProblemCreator<FunctionSpaceType,GridType>::apply();
      }


      template< int polOrd >
      struct DiscreteTraits
      {
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, FunctionSpaceType>
                                                                                           DFSpaceType;

#ifdef VELO
        typedef typename AC::template FunctionSpaces<GRIDDIM>                              VelFunctionSpaceType;
        typedef typename AC::template DiscreteFunctionSpaces< GridPartType, polOrd, VelFunctionSpaceType>
                                                                                           VelDFSpaceType;
#endif

      public:

        typedef typename AC::template DiscreteFunctions< DFSpaceType >                     DiscreteFunctionType;
#ifdef VELO
        typedef typename AC::template DiscreteFunctions< VelDFSpaceType >                  VelDiscreteFunctionType;
        typedef std::tuple< VelDiscreteFunctionType >                                      ExtraParameters;
#else
        typedef std::tuple<>                                                               ExtraParameters;
#endif

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                 IOTupleType;

        class Operator
        {
          typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameters >
                                                                                           OpTraits;
        public:
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::full >    type;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::expl >    ExplicitType;
          typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::impl >    ImplicitType;
        };

        struct Solver
        {
          typedef typename AC::template LinearSolvers< DFSpaceType >                       LinearSolverType;
          typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                      type;
        };

      private:
        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameters >
                                                                                           OpTraits;
        typedef DGAdaptationIndicatorOperator< OpTraits >                                  IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >  GradientIndicatorType ;
      public:

        typedef AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
        typedef SubSolverMonitor< SolverMonitor >                                          SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                              DiagnosticsType;
        typedef ExactSolutionOutput< DiscreteFunctionType >                                AdditionalOutputType;
      };

      template< int polOrd >
      using Algorithm = SubAdvectionDiffusionAlgorithm< GridType, SubAdvectionDiffusionAlgorithmCreator<AC>, polOrd >;

    };

    template< int polOrd >
    using Algorithm = EvolutionAlgorithm< polOrd, SubAdvectionDiffusionAlgorithmCreator<ACAdvDiff> >;


    static inline std::string moduleName() { return ""; }

    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef typename SubAdvectionDiffusionAlgorithmCreator<ACAdvDiff>::GridType GridType;
      //Discrete Functions
      typedef typename SubAdvectionDiffusionAlgorithmCreator<ACAdvDiff>::template DiscreteTraits<polOrd>::DiscreteFunctionType
                                                                                  DFType;
#ifdef VELO
      typedef typename SubAdvectionDiffusionAlgorithmCreator<ACAdvDiff>::template DiscreteTraits<polOrd>::VelDiscreteFunctionType
                                                                       VelDFType;
#endif

      //Item1
      typedef _t< SubEvolutionContainerItem >                                     Steady;
#ifdef VELO
      typedef _t< SolutionContainerItem >                                         Solut;
      typedef std::tuple< Steady, Solut >                                         Item1TupleType;
#else
      typedef std::tuple< Steady >                                                Item1TupleType;
#endif

      //Item2
      typedef _t< EmptyContainerItem >                                            Empty;
#ifdef VELO
      typedef std::tuple< std::tuple< Empty,Empty >,
                          std::tuple< Empty,Empty > >                             Item2TupleType;
#else
      typedef std::tuple< std::tuple< Empty > >                                   Item2TupleType;

#endif


      //Sub (discrete function argument ordering)
#ifdef VELO
      typedef std::tuple<__0, __1 >                                               AdvDiffOrder;
#else
      typedef std::tuple<__0 >                                                    AdvDiffOrder;
#endif

      typedef std::tuple< AdvDiffOrder >                                          SubOrderRowType;
      typedef SubOrderRowType                                                     SubOrderColType;

      //external params lists
#ifdef VELO
      typedef ExtraArg<_e< SolutionSelect, __0, __1 > >                           ExtraType;
#else
      typedef ExtraArg< >                                                          ExtraType;
#endif

      //Global container
#ifdef VELO
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, ExtraType, DFType, VelDFType >
#else
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, ExtraType, DFType >
#endif
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
