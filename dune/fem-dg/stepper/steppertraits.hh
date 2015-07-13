#ifndef DUNE_STEPPERTRAITS_HH
#define DUNE_STEPPERTRAITS_HH

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/solver/odesolver.hh>
#include <dune/fem/solver/pardginverseoperators.hh>

#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>

template <class GridImp,
          class ProblemTraits,
          int polOrd,
          class ExtraParameterTuple = std::tuple<> >
struct StepperTraits
{
  static const int polynomialOrder = polOrd ;

  // type of Grid
  typedef GridImp                                   GridType;

  // Choose a suitable GridView
  typedef Dune :: Fem :: DGAdaptiveLeafGridPart< GridType >       HostGridPartType;
  //typedef Dune :: Fem :: LeafGridPart< GridType >       HostGridPartType;
  typedef HostGridPartType  GridPartType ;
  //typedef Dune :: Fem :: IdGridPart< HostGridPartType >       GridPartType;

  typedef typename ProblemTraits :: template Traits< GridPartType >   ModelTraits;

  // traits for the operator class
  struct OperatorTraits :
    public Dune::PassTraits< ModelTraits, polynomialOrder == -1 ? 0 : polynomialOrder, ModelTraits::ModelType::dimRange >
  {
    static const int limiterPolynomialOrder = polynomialOrder == -1 ? 1 : polynomialOrder;
    typedef Dune::Fem::FiniteVolumeSpace< typename ModelTraits::ModelType::Traits,
            GridPartType, 0 > SpaceType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< SpaceType > VeloType;
    //typedef std::tuple< VeloType* > ExtraParameterTupleType;

    typedef ExtraParameterTuple ExtraParameterTupleType;
  };

  // problem dependent types
  typedef typename OperatorTraits :: InitialDataType   InitialDataType;
  typedef typename OperatorTraits :: ModelType         ModelType;
  typedef typename OperatorTraits :: FluxType          FluxType;
  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId
    = OperatorTraits :: PrimalDiffusionFluxId ;

private:
  // The DG Operator (using 2 Passes)
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
#ifdef LIMITER
  #if (not defined EULER) and (defined FLUXDG)
  #warning "DGAdvectionDiffusionOperator: using LIMITER."
    typedef Dune :: DGLimitedAdvectionDiffusionOperator< OperatorTraits >  DgType;
  #else
  #warning "DGAdvectionDiffusionOperator: LIMITER can NOT be used. Not supported -> LIMITER, no EULER, no FLUXDG."
    typedef Dune :: DGAdvectionDiffusionOperator< OperatorTraits >  DgType;
  #endif
  #ifndef HIGHER_ORDER_FV
  #warning "DGAdvectionOperator: using LIMITER."
    typedef Dune :: DGLimitedAdvectionOperator< OperatorTraits >    DgAdvectionType;
  #else
  #warning "DGAdvectionOperator: using HIGHER ORDER FV."
    typedef typename ProblemTraits :: template Traits< GridPartType, -1 > FVOperatorTraits;
    typedef Dune :: DGLimitedAdvectionOperator< FVOperatorTraits >  DgAdvectionType;
  #endif
#else // no LIMITER
#warning "No limiter is applied to the numerical solution !!"
  typedef Dune :: DGAdvectionDiffusionOperator< OperatorTraits >  DgType;
  typedef Dune :: DGAdvectionOperator< OperatorTraits >           DgAdvectionType;
#endif
  typedef Dune :: DGDiffusionOperator< OperatorTraits >           DgDiffusionType;

  // default is that both are enabled
  template < bool advection, bool diffusion >
  struct OperatorChooser
  {
    typedef DgType          FullOperatorType;
    typedef DgDiffusionType ImplicitOperatorType;
    typedef DgAdvectionType ExplicitOperatorType;
  };

  // advection operator only, i.e. linear advection equation
  template < bool advection >
  struct OperatorChooser< advection, false >
  {
    typedef DgAdvectionType  FullOperatorType;
    typedef FullOperatorType ImplicitOperatorType;
    typedef FullOperatorType ExplicitOperatorType;
  };

  // diffusion operator only, i.e. for the heat equation
  template < bool diffusion >
  struct OperatorChooser< false, diffusion >
  {
    typedef DgDiffusionType  FullOperatorType;
    typedef FullOperatorType ImplicitOperatorType;
    typedef FullOperatorType ExplicitOperatorType;
  };
public:

  static const bool advection = ModelType :: hasAdvection ;
  static const bool diffusion = ModelType :: hasDiffusion ;

  typedef OperatorChooser< advection, diffusion > OperatorChooserType ;

  typedef typename OperatorChooserType :: FullOperatorType      FullOperatorType;
  typedef typename OperatorChooserType :: ImplicitOperatorType  ImplicitOperatorType;
  typedef typename OperatorChooserType :: ExplicitOperatorType  ExplicitOperatorType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename DgType :: DestinationType                         DiscreteFunctionType;

  // The indicator function in case of limiting
  typedef typename DgAdvectionType :: IndicatorType                  IndicatorType;

  // ... as well as the Space type
  typedef typename DgType :: SpaceType                               DiscreteSpaceType;

  // The ODE Solvers                                                         /*@LST1S@*/
  typedef DuneODE :: OdeSolverInterface< DiscreteFunctionType > OdeSolverType ;

  // type of restriction/prolongation projection for adaptive simulations
  typedef Dune :: Fem :: RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;

  // type of IOTuple
  typedef Dune::tuple< DiscreteFunctionType*, DiscreteFunctionType*, IndicatorType* >  IOTupleType;

  // type of linear solver for implicit ode
  typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >  LinearInverseOperatorType;
  typedef Dune::Fem::ParDGGeneralizedMinResInverseOperator< DiscreteFunctionType >  LinearInverseImplicitOperatorType;

};

#endif
