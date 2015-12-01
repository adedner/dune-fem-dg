#ifndef FEMDG_PROBLEMCREATOR_CONFIGURATOR_HH
#define FEMDG_PROBLEMCREATOR_CONFIGURATOR_HH

#include "problemcreatorselector.hh"

namespace Dune
{
namespace Fem
{


  ////DiffusionFlux 	BR, CDG, CDG2...
////AdvectionFlux 	LLF, HLLC...
////GalerkinType 	CG, DG
////Solver 	istl, fem...
////polOrder 	0, 1, 2...
////GridDim 	1,2,3...
////Adaptation 	no, yes (Estimator...)
////DimRange 	1,2,3..
////Assembly 	matrix-free/matrix-based
////Codegen 	yes/no
////Parallel 	yes/no


 //     template< class AnalyticalTraitsImp, class NewModelImp >
 //     struct RhsAnalyticalTraits
 //     {
 //       typedef typename AnalyticalTraitsImp::ProblemInterfaceType    ProblemType;
 //       typedef typename AnalyticalTraitsImp::InitialDataType         InitialDataType;
 //       typedef typename AnalyticalTraitsImp::NewModelImp             ModelType;
 //     };
 //   public:

 //     template< class AnalyticalTraitsImp, class AdvectionFluxIdentifierImp, class DiffusionFluxIdentifierImp >
 //     using Traits = DefaultOperatorTraits< polOrd, AnalyticalTraits, DiscreteFunctionImp,
 //                                           AdvectionFluxes< AnalyticalTraitsImp, AdvectionFluxIdentifierImp >,
 //                                           DiffusionFluxes< AnalyticaltraitsImp, typename DiscreteFunctionImp::DiscreteFunctionSpaceType, DiffusionFluxIdentifierImp >,
 //                                           typename Config::template DiffusionFluxes< DFSpaceType >,
 //                                           ExtraParameterTuple >;


 //     typedef DefaultOperatorTraits< GridPartType, polOrd, RhsAnalyticalTraits,

  template< class GridImp,
            Galerkin::Enum dgId,
            DiscreteFunctionSpaces::Enum spaceId,
            Solver::Enum solverId,
            AdvectionLimiter::Enum advLimitId,
            Matrix::Enum matrixId,
            class AdvectionFluxIdentifierImp,
            class DiffusionFluxIdentifierImp >
  class AlgorithmConfigurator
  {
    template< class AnalyticalTraitsImp, class NewModelImp >
    struct RhsAnalyticalTraits
    {
      typedef typename AnalyticalTraitsImp::ProblemInterfaceType    ProblemType;
      typedef typename AnalyticalTraitsImp::InitialDataType         InitialDataType;
      typedef typename AnalyticalTraitsImp::NewModelImp             ModelType;
    };

  public:
    typedef GridImp                                           GridType;

    using GridParts = typename GridPartSelector< GridType, dgId >::type;

    template< int dimRange >
    using FunctionSpaces = Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, dimRange >;

    template< class GridPartImp, int polOrd, class FunctionSpaceImp >
    using DiscreteFunctionSpaces = typename DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrd, spaceId, dgId >::type;

    template< class DFSpace >
    using DiscreteFunctions = typename SolverSelector< solverId, false, DFSpace >::DiscreteFunctionType;

    template< class ModelImp, class AdvectionFluxIdentifier = AdvectionFluxIdentifierImp >
    using AdvectionFluxes = typename AdvectionFluxSelector< AdvectionFluxIdentifier >::template type<ModelImp>;

    template< class ModelImp, class DFSpace, class DiffusionFluxIdentifier = DiffusionFluxIdentifierImp >
    using DiffusionFluxes = typename DiffusionFluxSelector< DiffusionFluxIdentifier >::template type< DFSpace, ModelImp>;

    template< class DomainDFSpace,
              class RangeDFSpace,
              int polOrd,
              class AnalyticalTraitsImp,
              class AdvectionFluxIdentifier = AdvectionFluxIdentifierImp,
              class DiffusionFluxIdentifier = DiffusionFluxIdentifierImp >
    using DefaultAssembTraits = DefaultAssemblerTraits< polOrd,
                                                        AnalyticalTraitsImp,
                                                        typename SolverSelector< solverId, false, DomainDFSpace, RangeDFSpace >::LinearOperatorType,
                                                        AdvectionFluxes< typename AnalyticalTraitsImp::ModelType, AdvectionFluxIdentifier >,
                                                        DiffusionFluxes< typename AnalyticalTraitsImp::ModelType, DomainDFSpace, DiffusionFluxIdentifier >,
                                                        DiscreteFunctions< DomainDFSpace >,
                                                        DiscreteFunctions< RangeDFSpace > >;
    template< class DomainDFSpace,
              int polOrd,
              class AnalyticalTraitsImp,
              class ExtraParameterTupleImp = std::tuple<>,
              class AdaptationIndicatorFunctionSpaceImp = typename DomainDFSpace::FunctionSpaceType,
              class AdvectionFluxIdentifier = AdvectionFluxIdentifierImp,
              class DiffusionFluxIdentifier = DiffusionFluxIdentifierImp >
    using DefaultOpTraits = DefaultOperatorTraits< polOrd,
                                                   AnalyticalTraitsImp,
                                                   DiscreteFunctions< DomainDFSpace >,
                                                   AdvectionFluxes< typename AnalyticalTraitsImp::ModelType, AdvectionFluxIdentifier >,
                                                   DiffusionFluxes< typename AnalyticalTraitsImp::ModelType, DomainDFSpace, DiffusionFluxIdentifier >,
                                                   ExtraParameterTupleImp,
                                                   AdaptationIndicatorFunctionSpaceImp >;

    //Operator/Assembler
    template< class OpTraits, OperatorSplit::Enum opSplit = OperatorSplit::Enum::full, Matrix::Enum matrix = matrixId >
    using Operators = typename OperatorSelector< OpTraits, advLimitId, opSplit, matrix >::type;

    //Matrix Containers
    template< class DomainDFSpace, class RangeDFSpace  >
    using Containers = typename SolverSelector< solverId, false, DomainDFSpace, RangeDFSpace >::LinearOperatorType;

    //Solver
    template< class DFSpace, bool symmetric = false >
    using LinearSolvers = typename SolverSelector< solverId, symmetric, DFSpace >::LinearInverseOperatorType;


    ////Real Jacobian for Runge-Kutta
    //template< class OperatorImp >
    //using JacobianOperators = DGHelmholtzOperator< OperatorImp >; //i.e. JacobianImplicitOperators

    //template< class OperatorImp >
    //using JacobianOperators = DGHelmholtzOperator< DGPrimalMatrixAssembly >; //i.e. JacobianImplicitOperators

    //general, rungekutta, rungekutta_pardg, rungekutta_assembled
    //template< class DiscreteFunctionImp >
    //using OdeSolvers = DuneODE::OdeSolverInterface< DiscreteFunctionImp >;

    //using OdeSolvers = RungeKuttaSolver< OperatorImp, ExplicitOperatorImp, ImplicitOperatorImp, BasicLinearSolverImp >;
    //using OdeSolvers = RungeKuttaSolver< OperatorImp, ExplicitOperatorImp, ImplicitOperatorImp, JacobianImplicitOperatorImp, BasicLinearSolverImp >;

  };


}
}
#endif
