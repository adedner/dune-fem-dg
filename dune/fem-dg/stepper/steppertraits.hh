#ifndef DUNE_STEPPERTRAITS_HH
#define DUNE_STEPPERTRAITS_HH

#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/solver/odesolver.hh>

#include <dune/fem-dg/pass/threadpass.hh>

// DefaultL2Projection
// -------------------

template< class DType, class RType >
class DefaultL2Projection
: public Dune::Fem::L2Projection< DType, RType >
{
  typedef Dune::Fem::L2Projection< DType, RType > BaseType;

public:
  DefaultL2Projection ()
  : BaseType( -1, !Dune::NonBlockingCommParameter::nonBlockingCommunication() )
  {}
};


// StepperTraits
// -------------

template< class GridImp, class ProblemTraits, int polOrd >
struct StepperTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                   GridType;

  typedef typename ProblemTraits::GridPartType GridPartType;

  // problem dependent types 
  typedef typename ProblemTraits :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: FluxType         FluxType;
  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId 
    = ProblemTraits :: template Traits< GridPartType > ::PrimalDiffusionFluxId ;

  // The DG Operator (using 2 Passes)
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
#ifdef LIMITER
  #if (not defined EULER) and (defined FLUXDG)
  #warning "DGAdvectionDiffusionOperator: using LIMITER."
    typedef Dune :: DGLimitedAdvectionDiffusionOperator< ModelType, FluxType,
                          DiffusionFluxId,  polynomialOrder >            DgType; /*@LST1E@*/
  #else
  #warning "DGAdvectionDiffusionOperator: LIMITER can NOT be used. Not supported -> LIMITER, no EULER, no FLUXDG."
    typedef Dune :: DGAdvectionDiffusionOperator< ModelType, FluxType,
                          DiffusionFluxId,  polynomialOrder >            DgType; /*@LST1E@*/
  #endif
  #ifndef HIGHER_ORDER_FV 
  #warning "DGAdvectionOperator: using LIMITER."
    typedef Dune :: DGLimitedAdvectionOperator< ModelType, FluxType,
                                 DiffusionFluxId, polynomialOrder >      DgAdvectionType; /*@LST1E@*/
  #else 
  #warning "DGAdvectionOperator: using HIGHER ORDER FV."
    typedef Dune :: DGLimitedAdvectionOperator< ModelType, FluxType,
                                 DiffusionFluxId, -1 >      DgAdvectionType; /*@LST1E@*/
  #endif
#else // no LIMITER 
#warning "No limiter is applied to the numerical solution !!"
  typedef Dune :: DGAdvectionDiffusionOperator< ModelType, FluxType,
                        DiffusionFluxId,  polynomialOrder >            DgType; /*@LST1E@*/
  typedef Dune :: DGAdvectionOperator< ModelType, FluxType,
                               DiffusionFluxId, polynomialOrder >      DgAdvectionType; /*@LST1E@*/
#endif                                       
  typedef Dune :: DGDiffusionOperator< ModelType, FluxType,
                               DiffusionFluxId, polynomialOrder >      DgDiffusionType; /*@LST1E@*/

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

  // initial data projection
  typedef typename InitialDataType::TimeDependentFunctionType TimeDependentFunctionType;
  typedef DefaultL2Projection< TimeDependentFunctionType, DiscreteFunctionType > InitialProjectionType;
};

#endif // #ifndef DUNE_STEPPERTRAITS_HH
