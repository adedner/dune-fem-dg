#ifndef DUNE_FEMDG_ALGORITHM_ADVECTION_STEPPER_HH
#define DUNE_FEMDG_ALGORITHM_ADVECTION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include <dune/fem-dg/algorithm/evolution.hh>


namespace Dune
{
namespace Fem
{


  template <class GridImp,
            class ProblemTraits,
            int polynomialOrder,
            class ExtraParameterTuple = std::tuple<> >
  struct AdvectionStepper
    : public EvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder, ExtraParameterTuple >
  {
    typedef EvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder, ExtraParameterTuple > BaseType ;

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
    typedef typename BaseType :: ExplicitOperatorType           FullOperatorType;
    typedef FullOperatorType  ExplicitOperatorType;
    typedef FullOperatorType  ImplicitOperatorType;

    typedef typename BaseType :: LinearInverseOperatorType LinearInverseOperatorType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType :: DiscreteFunctionType      DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType :: DiscreteSpaceType         DiscreteSpaceType;

    // The ODE Solvers
    typedef typename BaseType :: OdeSolverType             OdeSolverType;

    typedef typename BaseType :: TimeProviderType          TimeProviderType;
    typedef typename BaseType :: AdaptationManagerType     AdaptationManagerType;
    typedef typename BaseType :: AdaptationHandlerType     AdaptationHandlerType;
    typedef typename BaseType :: IndicatorType             IndicatorType;


    typedef typename FullOperatorType :: ExtraParameterTupleType  ExtraParameterTupleType;

    static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId =
      BaseType::Traits::DiffusionFluxId ;

    // advection = true , diffusion = false
    typedef Dune :: DGAdaptationIndicatorOperator< typename BaseType::OperatorTraits, true, false >  DGIndicatorType;

    // gradient estimator
    typedef Estimator< DiscreteFunctionType, InitialDataType > GradientIndicatorType ;

    // type of 64bit unsigned integer
    typedef typename BaseType :: UInt64Type  UInt64Type;

    using BaseType :: grid_;
    using BaseType :: gridPart_;
    using BaseType :: problem;
    using BaseType :: adaptationHandler_ ;
    using BaseType :: space ;
    using BaseType :: solution ;
    using BaseType :: adaptive ;
    using BaseType :: adaptParam_;
    using BaseType :: doEstimateMarkAdapt ;
    using BaseType :: name ;

    // constructor
    AdvectionStepper( GridType& grid, const std::string name = "",
                      ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
      BaseType( grid, name ),
      dgAdvectionOperator_(gridPart_, problem(), tuple, name ),
      dgIndicator_( gridPart_, problem(), tuple, name ),
      gradientIndicator_( space(), problem(), adaptParam_ )
    {
    }

    virtual OdeSolverType* createOdeSolver(TimeProviderType& tp)
    {
      if( adaptive() )
      {
        if( ! adaptationHandler_ && adaptParam_.aposterioriIndicator() )
        {
          adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
          dgIndicator_.setAdaptation( *adaptationHandler_ );
        }
      }

      //typedef SmartOdeSolver< ExplicitOperatorType, ExplicitOperatorType, ExplicitOperatorType > OdeSolverImpl;

      typedef RungeKuttaSolver< ExplicitOperatorType, ExplicitOperatorType, ExplicitOperatorType,
                                LinearInverseOperatorType > OdeSolverImpl;
      return new OdeSolverImpl( tp, dgAdvectionOperator_,
                                dgAdvectionOperator_,
                                dgAdvectionOperator_,
                                name() );
    }

    //! return overal number of grid elements
    virtual uint64_t gridSize() const
    {
      // is adaptation handler exists use the information to avoid global comm
      if( adaptationHandler_ )
        return adaptationHandler_->globalNumberOfElements() ;

      size_t  advSize   = dgAdvectionOperator_.numberOfElements();
      size_t  dgIndSize = gradientIndicator_.numberOfElements();
      uint64_t grSize   = std::max( advSize, dgIndSize );
      return grid_.comm().sum( grSize );
    }

    //! call limiter (only if dgAdvectionOperator_ is DGLimitedAdvectionOperator)
    void limitSolution()
    {
      dgAdvectionOperator_.limit( solution() );
    }

    // return indicator pointer for data output, see above
    IndicatorType* indicator ()
    {
      return dgAdvectionOperator_.indicator();
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

    const ModelType& model() const { return dgAdvectionOperator_.model(); }

  protected:
    ExplicitOperatorType    dgAdvectionOperator_;
    DGIndicatorType         dgIndicator_;
    GradientIndicatorType   gradientIndicator_;
  };
}
}
#endif // FEMHOWTO_STEPPER_HH