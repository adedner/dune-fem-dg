#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
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

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/test/stokes/stokesalgorithm.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
#include "../stokes/problems.hh"
//--------- MODELS --------------------------
#include "stokesmodel.hh"
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>
//--------- HANDLER -------------------------
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>

namespace Dune
{
namespace Fem
{

  template< class GridImp >
  struct IncompressibleNavierStokesProblemCreator
  {

  //=======================================
  //     STOKES PROBLEM CREATOR
  //=======================================

    struct SubStokesProblemCreator
    {

      struct SubPoissonProblemCreator
      {
        typedef GridImp                                         GridType;
        typedef Fem::DGAdaptiveLeafGridPart< GridType >         HostGridPartType;
        typedef HostGridPartType                                GridPartType;

        // define problem type here if interface should be avoided
        typedef NavierStokesProblemDefault< GridType >                 ProblemInterfaceType;

        typedef typename ProblemInterfaceType::FunctionSpaceType       FunctionSpaceType;

        struct AnalyticalTraits
        {
          typedef ProblemInterfaceType                                   ProblemType;
          typedef ProblemInterfaceType                                   InitialDataType;
          typedef StokesModel< GridPartType, InitialDataType, false >    ModelType;

          template< class Solution, class Model, class ExactFunction, class SigmaFunction>
          static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
          {}
        };

        static inline std::string moduleName() { return "";}

        static ProblemInterfaceType* problem()
        {
          return new NavierStokesProblemDefault< GridType > ();
        }


        //Stepper Traits
        template< int polOrd >
        struct DiscreteTraits
        {
          static const SolverType solverType = istl;
        private:
          static const DiscreteFunctionSpaceIdentifier::id spaceId = DiscreteFunctionSpaceIdentifier::legendre;
          static const GalerkinIdentifier::id dgId = GalerkinIdentifier::dg;
          static const bool symmetricSolver = true;
        public:
          typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, spaceId, dgId >::type    DiscreteFunctionSpaceType;
          typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                          DiscreteFunctionType;
          typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                      JacobianOperatorType;

          typedef std::tuple<>                                                                        ExtraParameterTuple;

          typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

          class Operator
          {
            friend DiscreteTraits;
            friend SolverType;
            typedef DGAdvectionFlux< typename AnalyticalTraits::ModelType, AdvectionFluxIdentifier::none > AdvectionFluxType;
            typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceType, typename AnalyticalTraits::ModelType, DGDiffusionFluxIdentifier::general >    DiffusionFluxType;

            typedef DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                           DiscreteFunctionType, AdvectionFluxType, DiffusionFluxType, ExtraParameterTuple >  OperatorTraitsType;

            typedef DGAdvectionDiffusionOperator< OperatorTraitsType >                            AssemblyOperatorType;
            typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>               SolversType;



            struct RhsAnalyticalTraits
            {
              typedef ProblemInterfaceType                                   ProblemType;
              typedef ProblemInterfaceType                                   InitialDataType;
              typedef StokesModel< GridPartType, InitialDataType, true >     ModelType;
            };

            typedef DGAdvectionFlux< typename RhsAnalyticalTraits::ModelType, AdvectionFluxIdentifier::none > RhsAdvectionFluxType;
            typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceType, typename RhsAnalyticalTraits::ModelType, DGDiffusionFluxIdentifier::general >    RhsDiffusionFluxType;
            typedef DefaultOperatorTraits< GridPartType, polOrd, RhsAnalyticalTraits,
                                           DiscreteFunctionType, RhsAdvectionFluxType, RhsDiffusionFluxType, ExtraParameterTuple,
                                           FunctionSpaceType >                                       RhsOperatorTraitsType;


          public:
            typedef DGPrimalMatrixAssembly< AssemblyOperatorType >                                AssemblerType;
            typedef typename SolversType::LinearOperatorType                                      type;
            typedef DGAdvectionDiffusionOperator< RhsOperatorTraitsType >                         RhsType;
          };

          struct Solver
          {
            typedef typename Operator::SolversType::LinearInverseOperatorType                           type;
          };

          typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                                SolverMonitorHandlerType;
          typedef Fem::SubDiagnosticsHandler< Diagnostics >                                         DiagnosticsHandlerType;
        };

        template <int polOrd>
        struct Stepper
        {
          // this should be ok but could lead to a henn-egg problem
          typedef Fem::EllipticAlgorithm< GridType, SubPoissonProblemCreator, polOrd > Type;
        };
      };

      typedef typename SubPoissonProblemCreator::GridType             GridType;
      typedef typename SubPoissonProblemCreator::HostGridPartType     HostGridPartType;
      typedef typename SubPoissonProblemCreator::GridPartType         GridPartType;

      typedef NavierStokesProblemDefault< GridType >                  ProblemInterfaceType;

      typedef typename ProblemInterfaceType::PressureFunctionSpaceType FunctionSpaceType;

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

      static ProblemInterfaceType* problem() { return new NavierStokesProblemDefault< GridType > (); }

      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        typedef typename SubPoissonProblemCreator::template DiscreteTraits< polOrd >          PoissonDiscreteTraits;
        static const SolverType solverType = PoissonDiscreteTraits::solverType;
        static const DiscreteFunctionSpaceIdentifier::id spaceId = DiscreteFunctionSpaceIdentifier::legendre;
        static const GalerkinIdentifier::id dgId = GalerkinIdentifier::dg;
        static const bool symmetricSolver = true;
        typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, spaceId, dgId >::type             DiscreteFunctionSpaceType;
      public:
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

        typedef std::tuple<>                                                                        ExtraParameterTuple;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

        class Operator
        {
          friend DiscreteTraits;
          friend SolverType;
          typedef DGAdvectionFlux< typename AnalyticalTraits::ModelType, AdvectionFluxIdentifier::none > AdvectionFluxType;
          typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceType, typename AnalyticalTraits::ModelType, DGDiffusionFluxIdentifier::general >    DiffusionFluxType;

          typedef DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                         DiscreteFunctionType, AdvectionFluxType, DiffusionFluxType, ExtraParameterTuple,
                                         typename SubPoissonProblemCreator::FunctionSpaceType >       OperatorTraitsType;

          struct RhsAnalyticalTraits
          {
            typedef ProblemInterfaceType                                   ProblemType;
            typedef ProblemInterfaceType                                   InitialDataType;
            typedef StokesModel< GridPartType, InitialDataType, true >     ModelType;
          };

          typedef DGAdvectionFlux< typename RhsAnalyticalTraits::ModelType, AdvectionFluxIdentifier::none > RhsAdvectionFluxType;
          typedef DefaultOperatorTraits< GridPartType, polOrd, RhsAnalyticalTraits,
                                         DiscreteFunctionType, RhsAdvectionFluxType, ExtraParameterTuple,
                                         typename SubPoissonProblemCreator::FunctionSpaceType >        RhsOperatorTraitsType;

        public:

          typedef StokesAssembler< typename PoissonDiscreteTraits::DiscreteFunctionType,
                                   DiscreteFunctionType,
                                   OperatorTraitsType>                                                AssemblerType;
          typedef typename PoissonDiscreteTraits::Operator::type                                      type;
          //typedef DGAdvectionDiffusionOperator< RhsOperatorTraitsType >                             RhsType;
          typedef void                         RhsType;
        };

        struct Solver
        {
          typedef UzawaSolver< typename PoissonDiscreteTraits::DiscreteFunctionType, DiscreteFunctionType, typename Operator::AssemblerType,
                               typename PoissonDiscreteTraits::Solver::type >                   type;
        };

        static_assert( (int)DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

        typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                    SolverMonitorHandlerType;
        typedef Fem::SubDiagnosticsHandler< Diagnostics >                             DiagnosticsHandlerType;
      };

      template <int polOrd>
      struct Stepper
      {
        typedef Fem::StokesAlgorithm< GridType, SubStokesProblemCreator, SubPoissonProblemCreator, polOrd > Type;
      };
    };


  //=======================================
  //     NAVIER-STOKES PROBLEM CREATOR
  //=======================================

    struct SubNavierStokesProblemCreator
    {

      typedef GridImp                                         GridType;
      typedef Fem::DGAdaptiveLeafGridPart< GridType >         HostGridPartType;
      typedef HostGridPartType                                GridPartType;

      // define problem type here if interface should be avoided
      typedef NavierStokesProblemDefault< GridType >                       ProblemInterfaceType;

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

      static ProblemInterfaceType* problem() { return new ProblemInterfaceType(); }

      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        static const SolverType solverType = fem ;
        static const DiscreteFunctionSpaceIdentifier::id spaceId = DiscreteFunctionSpaceIdentifier::legendre;
        static const GalerkinIdentifier::id dgId = GalerkinIdentifier::dg;
        static const AdvectionLimiterIdentifier::id advLimitId = AdvectionLimiterIdentifier::unlimited;
        typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, spaceId, dgId >::type             DiscreteFunctionSpaceType;
      public:
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
        typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

        typedef std::tuple<> ExtraParameterTuple;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

        class Operator
        {
          friend DiscreteTraits;
          typedef DGAdvectionFlux< typename AnalyticalTraits::ModelType, AdvectionFluxIdentifier::llf > AdvectionFluxType;
          typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceType, typename AnalyticalTraits::ModelType, DGDiffusionFluxIdentifier::general >    DiffusionFluxType;

          typedef DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits, DiscreteFunctionType, AdvectionFluxType, DiffusionFluxType, ExtraParameterTuple >
                                                                                                      OperatorTraitsType;

          static const int hasAdvection = AnalyticalTraits::ModelType::hasAdvection;
          static const int hasDiffusion = AnalyticalTraits::ModelType::hasDiffusion;
          typedef AdvectionDiffusionOperators< OperatorTraitsType, hasAdvection, hasDiffusion, advLimitId > AdvectionDiffusionOperatorType;

          struct RhsAnalyticalTraits
          {
            typedef ProblemInterfaceType                                       ProblemType;
            typedef ProblemInterfaceType                                       InitialDataType;
            typedef NavierStokesModel< GridPartType, InitialDataType, true >   ModelType;
          };

          typedef DGAdvectionFlux< typename AnalyticalTraits::ModelType, AdvectionFluxIdentifier::llf > RhsAdvectionFluxType;
          typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceType, typename AnalyticalTraits::ModelType, DGDiffusionFluxIdentifier::general >    RhsDiffusionFluxType;

          typedef DefaultOperatorTraits< GridPartType, polOrd, RhsAnalyticalTraits,
                                         DiscreteFunctionType, RhsAdvectionFluxType, RhsDiffusionFluxType, ExtraParameterTuple >     RhsOperatorTraitsType;


        public:
          typedef typename AdvectionDiffusionOperatorType::FullOperatorType                           type;
          typedef typename AdvectionDiffusionOperatorType::ExplicitOperatorType                       ExplicitType;
          typedef typename AdvectionDiffusionOperatorType::ImplicitOperatorType                       ImplicitType;
          typedef DGDiffusionOperator< RhsOperatorTraitsType >                                        RhsType;
        };

        struct Solver
        {
          // type of linear solver for implicit ode
          typedef Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >                  BasicLinearSolverType;

          typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                                 type;
        };

      private:
        typedef DGAdaptationIndicatorOperator< typename Operator::OperatorTraitsType, Operator::hasAdvection, Operator::hasDiffusion >
                                                                                                      IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >             GradientIndicatorType ;
      public:

        typedef Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
        typedef Fem::SubSolverMonitorHandler< Fem::SolverMonitor >                              SolverMonitorHandlerType;
        typedef Fem::SubDiagnosticsHandler< Diagnostics >                                       DiagnosticsHandlerType;
       // typedef Fem::ExactSolutionOutputHandler< DiscreteFunctionType >                         AdditionalOutputHandlerType;
      };

      template <int polOrd>
      struct Stepper
      {
        typedef Fem::AdvectionDiffusionAlgorithm< GridType, SubNavierStokesProblemCreator, polOrd > Type;
      };

    };


    template <int polOrd>
    struct Stepper
    {
      typedef Fem::EvolutionAlgorithm< polOrd, SubStokesProblemCreator, SubNavierStokesProblemCreator, SubStokesProblemCreator > Type;
    };

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return Fem::DefaultGridInitializer< GridType >::initialize(); }

  };

}
}
#endif
