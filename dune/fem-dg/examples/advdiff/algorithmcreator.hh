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
  template< class DomainFunction, class RangeFunction, template< class, class > class MatrixObjectImp >
  struct LinearOperatorInterface
  : public MatrixObjectImp< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
    public Fem::AssembledOperator< DomainFunction, RangeFunction >
  {
    typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
    typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
    typedef MatrixObjectImp< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType > BaseType;

    static constexpr bool assembled = true;

    using BaseType::apply;

    LinearOperatorInterface( const std::string & ,
                             const DomainSpaceType &domainSpace,
                             const RangeSpaceType &rangeSpace ) :
      BaseType( domainSpace, rangeSpace )
    {}

    virtual void operator()( const DomainFunction &arg, RangeFunction &dest ) const
    {
      BaseType()( arg, dest );
    }

    const BaseType &systemMatrix() const
    {
      return *this;
    }

    BaseType &systemMatrix()
    {
      return *this;
    }

    void communicate()
    {}
  };



  /**
   *  \brief problem creator for an advection diffusion problem
   */
  template< class GridImp >
  struct AdvectionDiffusionAlgorithmCreator
  {

    struct SubAdvectionDiffusionAlgorithmCreator
    {
      typedef AlgorithmConfigurator< GridImp,
                                     Galerkin::Enum::default_,
                                     Adaptivity::Enum::default_,
                                     DiscreteFunctionSpaces::Enum::default_, //legendre
                                     Solver::Enum::fem,
                                     AdvectionLimiter::Enum::default_,
                                     Matrix::Enum::default_,
                                     AdvectionFlux::Enum::upwind,
                                     DiffusionFlux::Enum::primal > AC;
                                     // DiffusionFlux::Enum::local > AC;

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
        typedef HeatEqnModel< GridPartType, InitialDataType >           ModelType;

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
        typedef std::tuple<>                                                               ExtraParameterTuple;
      public:
        typedef typename AC::template DiscreteFunctions< DFSpaceType >                     DiscreteFunctionType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                 IOTupleType;

        //class Operator
        //{
        //  typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
        //                                                                                   OpTraits;
        //public:
        //  typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::full >    type;
        //  typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::expl >    ExplicitType;
        //  typedef typename AC::template Operators< OpTraits,OperatorSplit::Enum::impl >    ImplicitType;
        //};

        class Operator
        {
          public:
          typedef typename AC::template DefaultAssembTraits< DFSpaceType, DFSpaceType, polOrd, AnalyticalTraits >
                                                                                   OpTraits;
          typedef typename AC::template Operators< OpTraits >                      AssemblerType;

          typedef AssembledPoissonSpaceOperator< typename Operator::OpTraits /*OP*/ >  type;
          typedef type                                                                 ExplicitType;
          typedef type                                                                 ImplicitType;
        };

        struct Solver
        {

          typedef typename AC::template Containers< DFSpaceType, DFSpaceType >             LinearContainerType;
          typedef typename AC::template LinearSolvers< DFSpaceType >                       LinearSolverType;

          //typedef DGHelmholtzJacobianOperator< LinearContainerType >                       RealLinearContainerType;
          typedef typename Operator::type::JacobianOperatorType                            RealLinearContainerType;

          typedef NewtonInverseOperator< RealLinearContainerType, LinearSolverType >           NonLinearSolverType;
          //NewtonInverseOperator< ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType >,
          //                       ISTLGMResOp< DiscreteFunctionType, JacobianOperatorType > >

          typedef RungeKuttaSolver< typename Operator::type,
                                    typename Operator::type,
                                    typename Operator::type,
                                    NonLinearSolverType,
                                    DuneODE::RungeKuttaOperatorAlgorithm<typename Operator::type> >   type;


#if 0


          //typedef typename AC::template LinearSolvers< DFSpaceType >                       LinearSolverType;
          //typedef DuneODE::OdeSolverInterface< DiscreteFunctionType >                      type;

          typedef typename AC::template Containers< DFSpaceType, DFSpaceType >               LinearContainerType;

          template< class LinOpType >
          using LinearSolverType = Dune::Fem::OEMBICGSTABOp< DiscreteFunctionType, LinOpType >;



          //typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinOpType;
          //typedef Dune::Fem::Operator< DiscreteFunctionType >                          LinOpType;
          //typedef LinearOperatorInterface< DiscreteFunctionType, DiscreteFunctionType, Dune::Fem::SparseRowLinearOperator >
          //                                                                                   LinOpType;

          typedef RungeKuttaSolver< typename Operator::type,
                                    typename Operator::type,
                                    typename Operator::type,
                                    LinearSolverType<LinearContainerType > >          type;
                                    //LinearSolverType<typename Operator::type::JacobianOperatorType > >          type;
          //typedef RungeKuttaSolver< typename Operator::type,
          //                          typename Operator::ExplicitType,
          //                          typename Operator::ImplicitType,
          //                          LinearSolverType<LinOpType> >                      type;
//                                    LinearSolverType<LinearContainerType> >            type;
#endif
        };

      private:
        typedef typename AC::template DefaultOpTraits< DFSpaceType, polOrd, AnalyticalTraits, ExtraParameterTuple >
                                                                                           OpTraits;
        typedef DGAdaptationIndicatorOperator< OpTraits >                                  IndicatorType;
        typedef Estimator< DiscreteFunctionType, typename AnalyticalTraits::ProblemType >  GradientIndicatorType ;
      public:

        //typedef AdaptIndicator< IndicatorType, GradientIndicatorType >                     AdaptIndicatorType;
        typedef AdaptIndicator<>                                                           AdaptIndicatorType;
        typedef SubSolverMonitor< SolverMonitor >                                          SolverMonitorType;
        typedef SubDiagnostics< Diagnostics >                                              DiagnosticsType;
        typedef ExactSolutionOutput< DiscreteFunctionType >                                AdditionalOutputType;
      };

      template< int polOrd >
      using Algorithm = SubAdvectionDiffusionAlgorithm< GridType, SubAdvectionDiffusionAlgorithmCreator, polOrd >;

    };

    template< int polOrd >
    using Algorithm = EvolutionAlgorithm< polOrd, UncoupledSubAlgorithms, SubAdvectionDiffusionAlgorithmCreator >;

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return DefaultGridInitializer< GridType >::initialize(); }


  };

}
}

#endif // FEMHOWTO_HEATSTEPPER_HH
