#ifndef DUNE_STEPPERTRAITS_HH
#define DUNE_STEPPERTRAITS_HH

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/solver/odesolver.hh>
#include <dune/fem/solver/pardginverseoperators.hh>

#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>

template <class GridImp,
          class ProblemTraits, 
          int polOrd>
struct StepperTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                   GridType;

  // Choose a suitable GridView
  typedef Dune :: Fem :: DGAdaptiveLeafGridPart< GridType >       HostGridPartType;
  typedef HostGridPartType  GridPartType ;
  //typedef Dune :: Fem :: IdGridPart< HostGridPartType >       GridPartType;

  // problem dependent types 
  typedef typename ProblemTraits :: template Traits< GridPartType > :: InitialDataType   InitialDataType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: ModelType         ModelType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: FluxType          FluxType;
  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId 
    = ProblemTraits :: template Traits< GridPartType > ::PrimalDiffusionFluxId ;

private:
  // The DG Operator (using 2 Passes)
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
#ifdef LIMITER
  #if (not defined EULER) and (defined FLUXDG)
  #warning "DGAdvectionDiffusionOperator: using LIMITER."
    typedef Dune :: DGLimitedAdvectionDiffusionOperator< ModelType, FluxType,
                          DiffusionFluxId,  polynomialOrder >            DgType;
  #else
  #warning "DGAdvectionDiffusionOperator: LIMITER can NOT be used. Not supported -> LIMITER, no EULER, no FLUXDG."
    typedef Dune :: DGAdvectionDiffusionOperator< ModelType, FluxType,
                          DiffusionFluxId,  polynomialOrder >            DgType;
  #endif
  #ifndef HIGHER_ORDER_FV 
  #warning "DGAdvectionOperator: using LIMITER."
    typedef Dune :: DGLimitedAdvectionOperator< ModelType, FluxType,
                                 DiffusionFluxId, polynomialOrder >      DgAdvectionType;
  #else 
  #warning "DGAdvectionOperator: using HIGHER ORDER FV."
    typedef Dune :: DGLimitedAdvectionOperator< ModelType, FluxType,
                                 DiffusionFluxId, -1 >      DgAdvectionType;
  #endif
#else // no LIMITER 
#warning "No limiter is applied to the numerical solution !!"
  typedef Dune :: DGAdvectionDiffusionOperator< ModelType, FluxType,
                        DiffusionFluxId,  polynomialOrder >            DgType;
  typedef Dune :: DGAdvectionOperator< ModelType, FluxType,
                               DiffusionFluxId, polynomialOrder >      DgAdvectionType;
#endif                                       
  typedef Dune :: DGDiffusionOperator< ModelType, FluxType,
                               DiffusionFluxId, polynomialOrder >      DgDiffusionType;

public:
  typedef DgType          FullOperatorType;
  typedef DgDiffusionType ImplicitOperatorType;
  typedef DgAdvectionType ExplicitOperatorType;

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
