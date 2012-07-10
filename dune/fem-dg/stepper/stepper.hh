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
  typedef typename BaseType :: DgType                    DgType;
  typedef typename BaseType :: DgAdvectionType           DgAdvectionType;
  typedef typename BaseType :: DgDiffusionType           DgDiffusionType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename BaseType :: DiscreteFunctionType      DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename BaseType :: DiscreteSpaceType         DiscreteSpaceType;

  // The ODE Solvers
  typedef typename BaseType :: OdeSolverType     OdeSolverType;

  typedef typename BaseType :: TimeProviderType       TimeProviderType;
  typedef typename BaseType :: AdaptationManagerType  AdaptationManagerType;
  typedef typename BaseType :: AdaptationHandlerType  AdaptationHandlerType; 

  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId =
    BaseType::Traits::DiffusionFluxId ;

  // advection = true , diffusion = true
  typedef Dune :: DGAdaptationIndicatorOperator< ModelType, FluxType,
            DiffusionFluxId, polynomialOrder, true, true >  DGIndicatorType;

  typedef Estimator< DiscreteFunctionType, InitialDataType > GradientIndicatorType ;

  using BaseType :: grid_;
  using BaseType :: gridPart_;
  using BaseType :: space;
  using BaseType :: convectionFlux_ ;
  using BaseType :: problem;
  using BaseType :: solution_;
  using BaseType :: adaptationHandler_ ;
  using BaseType :: adaptationParameters_;
  using BaseType :: adaptive_ ;

  Stepper( GridType& grid ) :
    BaseType( grid ),
    dgOperator_( gridPart_, convectionFlux_ ),
    dgAdvectionOperator_( gridPart_, convectionFlux_ ),
    dgDiffusionOperator_( gridPart_, convectionFlux_ ),
    dgIndicator_( gridPart_, convectionFlux_ ),
    gradientIndicator_( solution_, problem() )
  {
  }

  //! return overal number of grid elements 
  virtual size_t gridSize() const 
  {
    // is adaptation handler exists use the information to avoid global comm
    if( adaptationHandler_ ) 
      return adaptationHandler_->globalNumberOfElements() ;

    // one of them is not zero, 
    // use int because the unintialized size_t is the largest 
    size_t advSize  = dgAdvectionOperator_.numberOfElements();
    size_t diffSize = dgDiffusionOperator_.numberOfElements();
    size_t grSize   = std::max( advSize, diffSize );
    grSize = std::max( dgOperator_.numberOfElements(), grSize );
    return grid_.comm().sum( grSize );
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

    typedef SmartOdeSolver< DgType, DgAdvectionType, DgDiffusionType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgOperator_, 
                              dgAdvectionOperator_,
                              dgDiffusionOperator_ );
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
  DgType                  dgOperator_;
  DgAdvectionType         dgAdvectionOperator_;
  DgDiffusionType         dgDiffusionOperator_;
  DGIndicatorType         dgIndicator_;
  GradientIndicatorType   gradientIndicator_;
};
#endif // FEMHOWTO_STEPPER_HH
