#ifndef DUNE_STEPPERTRAITS_HH
#define DUNE_STEPPERTRAITS_HH

#warning "Overloaded StepperTraits"

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/solver/odesolver.hh>
#include <dune/fem/solver/pardginverseoperators.hh>

#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>

#include "passtraits.hh"
#include "../stokes/stokesalgorithm.hh"

template <class GridImp,
          class ProblemTraits,
          int polOrd,
          class ExtraParameterTuple = std::tuple< > >
struct StepperTraits
{
  static const int polynomialOrder = polOrd ;

  // type of Grid
  typedef GridImp                                   GridType;

  typedef typename ProblemTraits :: StokesProblemTraits StokesProblemTraits;
  typedef StokesAlgorithm< GridImp, StokesProblemTraits, polynomialOrder >  StokesAlgorithmType;
  typedef typename StokesAlgorithmType :: DiscreteFunctionType    VelocityFunctionType;

  // Choose a suitable GridView
  typedef Dune :: Fem :: DGAdaptiveLeafGridPart< GridType >       HostGridPartType;
  //typedef Dune :: Fem :: LeafGridPart< GridType >       HostGridPartType;
  typedef HostGridPartType  GridPartType ;
  //typedef Dune :: Fem :: IdGridPart< HostGridPartType >       GridPartType;

private:
  typedef typename StokesProblemTraits :: template Traits< GridPartType, true > RhsLaplaceModelTraits;

  // traits for the operator class
  struct RhsLaplaceOperatorTraits :
    public Dune::PassTraits< RhsLaplaceModelTraits, polynomialOrder, RhsLaplaceModelTraits::ModelType::dimRange, istl >
  {
    static const int limiterPolynomialOrder = polynomialOrder;
    typedef std::tuple<> ExtraParameterTupleType;
  };
public:

  typedef typename ProblemTraits :: template Traits< GridPartType >   ModelTraits;

  // traits for the operator class
  struct OperatorTraits :
    public Dune::PassTraits< ModelTraits, polynomialOrder == -1 ? 0 : polynomialOrder, ModelTraits::ModelType::dimRange, fem >
  {
    static const int limiterPolynomialOrder = polynomialOrder == -1 ? 1 : polynomialOrder;
    typedef std::tuple< VelocityFunctionType*, VelocityFunctionType* > ExtraParameterTupleType;
  };

  // problem dependent types
  typedef typename OperatorTraits :: InitialDataType   InitialDataType;
  typedef typename OperatorTraits :: ModelType         ModelType;
  typedef typename OperatorTraits :: FluxType          FluxType;
  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId
    = OperatorTraits :: PrimalDiffusionFluxId ;

private:
  typedef Dune :: DGAdvectionDiffusionOperator< OperatorTraits >  DgType;

  typedef typename ProblemTraits :: template Traits< GridPartType, false >  RhsModelTraits;

  // traits for the operator class
  struct RhsOperatorTraits :
    public Dune::PassTraits< RhsModelTraits, polynomialOrder, RhsModelTraits::ModelType::dimRange, istl >
  {
    static const int limiterPolynomialOrder = polynomialOrder;
    typedef std::tuple< VelocityFunctionType* > ExtraParameterTupleType;
  };

public:

  typedef DgType    FullOperatorType;
  typedef DgType    ImplicitOperatorType;
  typedef DgType    ExplicitOperatorType;

  typedef Dune :: DGAdvectionDiffusionOperator< RhsOperatorTraits >  RhsOperatorType;
  typedef Dune :: DGDiffusionOperator< RhsLaplaceOperatorTraits >     RhsLaplaceOperatorType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename DgType :: DestinationType                         DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename DgType :: SpaceType                               DiscreteSpaceType;

  // The ODE Solvers
  typedef DuneODE :: OdeSolverInterface< DiscreteFunctionType > OdeSolverType ;

  typedef DiscreteFunctionType  IndicatorType;

  // type of restriction/prolongation projection for adaptive simulations
  typedef Dune :: Fem :: RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;

  // type of IOTuple
  typedef Dune::tuple< DiscreteFunctionType*, DiscreteFunctionType*, IndicatorType* >  IOTupleType;

  // type of linear solver for implicit ode
  typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >  LinearInverseOperatorType;
  typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >  LinearInverseImplicitOperatorType;
};

#endif
