#ifndef FEMDG_PROBLEMCREATOR_CONFIGURATOR_HH
#define FEMDG_PROBLEMCREATOR_CONFIGURATOR_HH

#include "algorithmcreatorselector.hh"

namespace Dune
{
namespace Fem
{


  /**
   * \brief Convenience class for AlgorithmCreator classes.
   *
   * The module dune-fem-dg needs a lot of templates and typedefs which might confuse some users.
   * Furthermore, redundancy of templates expressions which are not updated correctly
   * may lead to nasty compiler errors.
   *
   * This class protects the user from this situation and enables all available implementations.
   *
   * \tparam GridImp type of the grid
   * \tparam Galerkin::Enum enum of the galerkin type
   * \tparam Adaptivity::Enum enum indicating whether we have adaptivity or not
   * \tparam DiscreteFunctionSpaces::Enum enum defining a discrete function space
   * \tparam Solver::Enum enum defining a solver
   * \tparam AdvectionLimiter::Enum enum defining the limiting of the advection operator
   * \tparam Matrix::Enum enum describing whether to assemble or not
   * \tparam AdvectionFlux::Enum enum describing the chosen numerical advection flux,
     \tparam PrimalDiffusionFlux::Enum enum describing the chosen primal numerical diffusion flux
   */
  template< class GridImp,
            Galerkin::Enum dgId,
            Adaptivity::Enum adap,
            DiscreteFunctionSpaces::Enum spaceId,
            Solver::Enum solverId,
            AdvectionLimiter::Enum advLimitId,
            Matrix::Enum matrixId,
            AdvectionFlux::Enum advFluxId,
            PrimalDiffusionFlux::Enum diffFluxId,
            LocalDiffusionFlux::Enum localDiffFluxId = LocalDiffusionFlux::Enum::general,
            Formulation::Enum formId = Formulation::Enum::primal >
  class AlgorithmConfigurator
  {
    template< class AnalyticalTraitsImp, class NewModelImp >
    struct NewRhsAnalyticalTraits
    {
      typedef typename AnalyticalTraitsImp::ProblemType     ProblemType;
      typedef typename AnalyticalTraitsImp::InitialDataType InitialDataType;
      typedef NewModelImp                                   ModelType;
    };

  public:
    typedef GridImp                                           GridType;

    using GridParts = typename GridPartSelector< GridType, dgId, adap >::type;

    template< int dimRange >
    using FunctionSpaces = Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, dimRange >;

    template< class GridPartImp, int polOrd, class FunctionSpaceImp >
    using DiscreteFunctionSpaces = typename DiscreteFunctionSpaceSelector< FunctionSpaceImp, GridPartImp, polOrd, spaceId, dgId >::type;

    template< class DFSpace >
    using DiscreteFunctions = typename SolverSelector< solverId, false, DFSpace >::DiscreteFunctionType;

    template< class ModelImp, AdvectionFlux::Enum id = advFluxId >
    using AdvectionFluxes = DGAdvectionFlux< ModelImp, id >;

    template< class ModelImp, class DFSpace, PrimalDiffusionFlux::Enum id = diffFluxId, LocalDiffusionFlux::Enum did = localDiffFluxId >
    using DiffusionFluxes = typename DiffusionFluxSelector< ModelImp, DFSpace, id, did, formId >::type;

    template< class DomainDFSpace,
              class RangeDFSpace,
              int polOrd,
              class AnalyticalTraitsImp,
              AdvectionFlux::Enum advId = advFluxId,
              PrimalDiffusionFlux::Enum diffId = diffFluxId,
              LocalDiffusionFlux::Enum localDiffId = localDiffFluxId >
    using DefaultAssembTraits = DefaultAssemblerTraits< polOrd,
                                                        AnalyticalTraitsImp,
                                                        typename SolverSelector< solverId, false, DomainDFSpace, RangeDFSpace >::LinearOperatorType,
                                                        AdvectionFluxes< typename AnalyticalTraitsImp::ModelType, advId >,
                                                        DiffusionFluxes< typename AnalyticalTraitsImp::ModelType, DomainDFSpace, diffId, localDiffId >,
                                                        DiscreteFunctions< DomainDFSpace >,
                                                        DiscreteFunctions< RangeDFSpace > >;
    template< class DomainDFSpace,
              int polOrd,
              class AnalyticalTraitsImp,
              class ExtraParameterTupleImp = std::tuple<>,
              class AdaptationIndicatorFunctionSpaceImp = typename DomainDFSpace::FunctionSpaceType,
              AdvectionFlux::Enum advId = advFluxId,
              PrimalDiffusionFlux::Enum diffId = diffFluxId,
              LocalDiffusionFlux::Enum localDiffId = localDiffFluxId >
    using DefaultOpTraits = DefaultOperatorTraits< polOrd,
                                                   AnalyticalTraitsImp,
                                                   DiscreteFunctions< DomainDFSpace >,
                                                   AdvectionFluxes< typename AnalyticalTraitsImp::ModelType, advId >,
                                                   DiffusionFluxes< typename AnalyticalTraitsImp::ModelType, DomainDFSpace, diffId, localDiffId >,
                                                   ExtraParameterTupleImp,
                                                   AdaptationIndicatorFunctionSpaceImp >;

    template< class AnalyticalTraitsImp, class NewModelImp >
    using RhsAnalyticalTraits = NewRhsAnalyticalTraits< AnalyticalTraitsImp, NewModelImp >;

    //Operator/Assembler
    template< class OpTraits, OperatorSplit::Enum opSplit = OperatorSplit::Enum::full, Matrix::Enum matrix = matrixId >
    using Operators = typename OperatorSelector< OpTraits, formId, advLimitId, opSplit, matrix >::type;

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
