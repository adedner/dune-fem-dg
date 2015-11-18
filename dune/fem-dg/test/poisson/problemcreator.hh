#ifndef FEMDG_POISSONSTEPPER_HH
#define FEMDG_POISSONSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

//--------- HANDLER --------------------------------
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
#include "gridinitializer.hh"
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/elliptic.hh>
#include <dune/fem-dg/algorithm/steadystate.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/problemcreatorselector.hh>

namespace Dune
{
namespace Fem
{

  template< class GridImp >
  struct PoissonProblemCreator
  {

    struct SubPoissonProblemCreator
    {
      typedef GridImp                                         GridType;
      //typedef typename GridPartChooser<GridType,galerkinType>::Type     HostGridPartType;
      typedef Fem::DGAdaptiveLeafGridPart< GridType >       HostGridPartType;
      typedef HostGridPartType                                GridPartType;

      // define problem type here if interface should be avoided
      typedef ProblemInterface< Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE > >
                                                                    ProblemInterfaceType;

      typedef typename ProblemInterfaceType::FunctionSpaceType      FunctionSpaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                ProblemType;
        typedef ProblemInterfaceType                                InitialDataType;
        typedef PoissonModel< GridPartType, InitialDataType >       ModelType;

        template< class Solution, class Model, class ExactFunction, class SigmaFunction>
        static void addEOCErrors ( Solution &u, Model &model, ExactFunction &f, SigmaFunction& sigma )
        {
          static L2EOCError l2EocError( "$L^2$-Error" );
          l2EocError.add( u, f );
          static DGEOCError dgEocError( "DG-Error" );
          dgEocError.add( u, f );
          static H1EOCError sigmaEocError( "sigma-norm" );
          sigmaEocError.add( sigma, f );
        }

      };

      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem()
      {
        int probNr = Fem::Parameter::getValue< int > ( "problem" );
        return new PoissonProblem< GridType, DIMRANGE > ( probNr );
      }


      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
        static const SolverType solverType = istl;
        static const DiscreteFunctionSpaceIdentifier::id spaceId = DiscreteFunctionSpaceIdentifier::legendre;
        static const GalerkinIdentifier::id dgId = GalerkinIdentifier::dg;
        static const bool symmetricSolver = true;
        typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polOrd, spaceId, dgId >::type             DiscreteFunctionSpaceType;
      public:
        typedef typename DiscreteFunctionSelector< DiscreteFunctionSpaceType, solverType >::type                          DiscreteFunctionType;
        typedef typename DiscreteFunctionSelector< DiscreteFunctionSpaceType, solverType >::jacobian                      JacobianOperatorType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                          IOTupleType;

        typedef std::tuple<>                                                                        ExtraParameterTuple;

        class Operator
        {
          friend DiscreteTraits;
          friend SolverType;
          typedef DGAdvectionFlux< typename AnalyticalTraits::ModelType, AdvectionFluxIdentifier::upwind > AdvectionFluxType;
          typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceType, typename AnalyticalTraits::ModelType, DGDiffusionFluxIdentifier::general >    DiffusionFluxType;

          typedef DefaultOperatorTraits< GridPartType, polOrd, AnalyticalTraits,
                                         DiscreteFunctionType, AdvectionFluxType, DiffusionFluxType, ExtraParameterTuple > OperatorTraitsType;

          typedef DGAdvectionDiffusionOperator< OperatorTraitsType >                                  AssemblyOperatorType;
          typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver>                     SolversType;
        public:
          typedef DGPrimalMatrixAssembly< AssemblyOperatorType >                                      AssemblerType;
          typedef typename SolversType::LinearOperatorType                                            type;
        };

        struct Solver
        {
          typedef typename Operator::SolversType::LinearInverseOperatorType                           type;
        };

      private:
        //typedef DGAdaptationIndicatorOperator< OperatorTraitsType, true, true >                       IndicatorType;
        //typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >             GradientIndicatorType ;
      public:

        //typedef Fem::AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
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

    template <int polOrd>
    struct Stepper
    {
      typedef Fem::SteadyStateAlgorithm< polOrd, SubPoissonProblemCreator > Type;
    };

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return Fem::DefaultGridInitializer< GridType >::initialize(); }



  };

}
}

#endif // FEMHOWTO_HEATSTEPPER_HH
