#ifndef DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH
#define DUNE_FEMDG_ADVECTIONDIFFUSION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include "stepperbase.hh"

template <class GridImp,
          class ProblemTraits,
          int polynomialOrder>
struct AdvectionDiffusionStepper
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
  typedef typename BaseType :: FullOperatorType               FullOperatorType;
  typedef typename BaseType :: ExplicitOperatorType           ExplicitOperatorType;
  typedef typename BaseType :: ImplicitOperatorType           ImplicitOperatorType;

  typedef typename BaseType :: LinearInverseOperatorType LinearInverseOperatorType;

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

  // gradient estimator
  typedef Estimator< DiscreteFunctionType, InitialDataType > GradientIndicatorType ;

  // type of 64bit unsigned integer
  typedef typename BaseType :: UInt64Type  UInt64Type;

  using BaseType :: grid_;
  using BaseType :: gridPart_;
  using BaseType :: space;
  using BaseType :: problem;
  using BaseType :: adaptationHandler_ ;
  using BaseType :: adaptationParameters_;
  using BaseType :: adaptive ;
  using BaseType :: doEstimateMarkAdapt ;
  using BaseType :: name ;

  AdvectionDiffusionStepper( GridType& grid, const std::string name = "" ) :
    BaseType( grid, name ),
    dgOperator_( gridPart_, problem(), name  ),
    dgAdvectionOperator_( gridPart_, problem(), name ),
    dgDiffusionOperator_( gridPart_, problem(), name ),
    dgIndicator_( gridPart_, problem() ),
    gradientIndicator_( space(), problem() )
  {
  }

  //! return overal number of grid elements
  virtual UInt64Type gridSize() const
  {
    // is adaptation handler exists use the information to avoid global comm
    if( adaptationHandler_ )
    {
      UInt64Type globalElements = adaptationHandler_->globalNumberOfElements() ;
      if( Dune::Fem::Parameter::verbose () )
      {
        std::cout << "grid size (sum,min,max) = ( "
          << globalElements << " , "
          << adaptationHandler_->minNumberOfElements() << " , "
          << adaptationHandler_->maxNumberOfElements() << ")" << std::endl;
      }
      return globalElements;
    }

    // one of them is not zero,
    size_t advSize     = dgAdvectionOperator_.numberOfElements();
    size_t diffSize    = dgDiffusionOperator_.numberOfElements();
    size_t dgIndSize   = gradientIndicator_.numberOfElements();
    size_t dgSize      = dgOperator_.numberOfElements();
    UInt64Type grSize  = std::max( std::max(advSize, dgSize ), std::max( diffSize, dgIndSize ) );
    double minMax[ 2 ] = { double(grSize), 1.0/double(grSize) } ;
    grid_.comm().max( &minMax[ 0 ], 2 );
    if( Dune::Fem::Parameter :: verbose () )
    {
      std::cout << "grid size (min,max) = ( " << size_t(1.0/minMax[ 1 ]) << " , " << size_t(minMax[ 0 ]) << ")" << std::endl;
    }
    return grid_.comm().sum( grSize );
  }

  virtual OdeSolverType* createOdeSolver(TimeProviderType& tp)
  {
    // create adaptation handler in case of apost indicator
    if( adaptive() )
    {
      if( ! adaptationHandler_ && adaptationParameters_.aposterioriIndicator() )
      {
        adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
        dgIndicator_.setAdaptation( *adaptationHandler_ );
      }
    }

    // create ODE solver
    typedef RungeKuttaSolver< FullOperatorType, ExplicitOperatorType, ImplicitOperatorType,
                              LinearInverseOperatorType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgOperator_,
                              dgAdvectionOperator_,
                              dgDiffusionOperator_,
                              name() );
  }

  //! estimate and mark solution
  virtual void initialEstimateMarkAdapt( )
  {
    doEstimateMarkAdapt( dgIndicator_, gradientIndicator_, true );
  }

  //! estimate and mark solution
  virtual void estimateMarkAdapt( )
  {
    doEstimateMarkAdapt( dgIndicator_, gradientIndicator_, false );
  }

  const ModelType& model() const { return dgOperator_.model(); }

protected:
  FullOperatorType        dgOperator_;
  ExplicitOperatorType    dgAdvectionOperator_;
  ImplicitOperatorType    dgDiffusionOperator_;
  DGIndicatorType         dgIndicator_;
  GradientIndicatorType   gradientIndicator_;
};
#endif // FEMHOWTO_STEPPER_HH
