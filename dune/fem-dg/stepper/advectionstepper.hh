#ifndef NAVIER_STOKES_STEPPER_HH
#define NAVIER_STOKES_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>

// local includes
#include "stepperbase.hh"

template <class GridImp,
          class ProblemTraits, 
          int polynomialOrder>             
struct Stepper 
  : public StepperBase< GridImp, ProblemTraits, polynomialOrder >
{
  typedef StepperBase< GridImp, ProblemTraits, polynomialOrder > BaseType ;

  // type of Grid
  typedef typename BaseType :: GridType                 GridType;

  // Choose a suitable GridView
  typedef typename BaseType :: GridPartType             GridPartType;

  // initial data type 
  typedef typename BaseType :: InitialDataType          InitialDataType;

  // An analytical version of our model
  typedef typename BaseType :: ModelType                 ModelType;

  // The flux for the discretization of advection terms
  typedef typename BaseType :: FluxType                  FluxType;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef typename BaseType :: DgAdvectionType           DgType;
  typedef DgType  DgAdvectionType;
  typedef DgType  DgDiffusionType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename BaseType :: DiscreteFunctionType      DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename BaseType :: DiscreteSpaceType         DiscreteSpaceType;

  // The ODE Solvers
  typedef typename BaseType :: OdeSolverType     OdeSolverType;

  typedef typename BaseType :: TimeProviderType        TimeProviderType;
  typedef typename BaseType :: AdaptationManagerType   AdaptationManagerType;
  typedef typename BaseType :: AdaptationHandlerType   AdaptationHandlerType;

  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId =
    BaseType::Traits::DiffusionFluxId ;

  // advection = true , diffusion = false 
  typedef Dune :: DGAdaptationIndicatorOperator< ModelType, FluxType,
            DiffusionFluxId, polynomialOrder, true, false >  DGIndicatorType;

  typedef Estimator< DiscreteFunctionType, InitialDataType > GradientIndicatorType ;

  using BaseType :: grid_;
  using BaseType :: gridPart_;
  using BaseType :: convectionFlux_ ;
  using BaseType :: problem;
  using BaseType :: adaptationHandler_ ;
  using BaseType :: solution_ ;
  using BaseType :: adaptive_ ;
  using BaseType :: adaptationParameters_;

  Stepper( GridType& grid ) :
    BaseType( grid ),
    dgAdvectionOperator_(gridPart_, convectionFlux_),
    dgIndicator_( gridPart_, convectionFlux_ ),
    gradientIndicator_( solution_, problem() )
  {
  }

  virtual OdeSolverType* createOdeSolver(TimeProviderType& tp) 
  {
    if( adaptive_ )
    {
      if( ! adaptationHandler_ && adaptationParameters_.aposterioriIndicator() )
      {
        adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
        dgIndicator_.setAdaptationHandler( *adaptationHandler_ );
      }
    }

    typedef SmartOdeSolver< DgAdvectionType, DgAdvectionType, DgAdvectionType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgAdvectionOperator_, 
                              dgAdvectionOperator_,
                              dgAdvectionOperator_ );
  }

  //! return overal number of grid elements 
  virtual size_t gridSize() const
  {
    // is adaptation handler exists use the information to avoid global comm
    if( adaptationHandler_ )
      return adaptationHandler_->globalNumberOfElements() ;

    size_t grSize  = dgAdvectionOperator_.numberOfElements();
    return grid_.comm().sum( grSize );
  }

  //! call limiter (only if dgAdvectionOperator_ is DGLimitedAdvectionOperator)
  void limitSolution() 
  { 
    dgAdvectionOperator_.limit( solution_ );
  }

  //! estimate and mark solution 
  virtual void initialEstimateMarkAdapt( AdaptationManagerType& am )
  {
    doEstimateMarkAdapt( dgIndicator_, gradientIndicator_, am, true );
  }

  //! estimate and mark solution 
  virtual void estimateMarkAdapt( AdaptationManagerType& am )
  {
    doEstimateMarkAdapt( dgIndicator_, gradientIndicator_, am, false );
  }

protected:
  DgAdvectionType         dgAdvectionOperator_;
  DGIndicatorType         dgIndicator_;
  GradientIndicatorType   gradientIndicator_;
};
#endif // FEMHOWTO_STEPPER_HH
