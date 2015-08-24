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
            int polynomialOrder >
  struct AdvectionStepper
    : public SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder >
  {
    typedef SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder > BaseType ;

    // type of Grid
    typedef typename BaseType::GridType                        GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                    GridPartType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                       ModelType;

    // The DG space operator
    // The first operator is sum of the other two
    // The other two are needed for semi-implicit time discretization
    typedef typename BaseType::OperatorType::FullType          FullOperatorType;
    typedef FullOperatorType                                   ExplicitOperatorType;
    typedef FullOperatorType                                   ImplicitOperatorType;

    typedef typename BaseType::BasicLinearSolverType           BasicLinearSolverType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType            DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType       DiscreteFunctionSpaceType;

    // The ODE Solvers
    typedef typename BaseType::OdeSolverType                   OdeSolverType;

    typedef typename BaseType::TimeProviderType                TimeProviderType;

    typedef typename FullOperatorType::ExtraParameterTupleType ExtraParameterTupleType;

    // type of 64bit unsigned integer
    typedef typename BaseType::UInt64Type                      UInt64Type;

    using BaseType::grid_;
    using BaseType::gridPart_;
    using BaseType::problem;
    using BaseType::space ;
    using BaseType::solution ;
    using BaseType::name ;
    using BaseType::adaptHandler_;
    using BaseType::solutionLimiterHandler_;

    // constructor
    AdvectionStepper( GridType& grid, const std::string name = "",
                      ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
      BaseType( grid, name ),
      dgAdvectionOperator_(gridPart_, problem(), tuple, name )
    {
       adaptHandler_.setIndicator( problem(), tuple );
       solutionLimiterHandler_.setLimiter( &dgAdvectionOperator_ );
    }

    virtual OdeSolverType* createOdeSolver(TimeProviderType& tp)
    {
      adaptHandler_.setAdaptation( tp );

      typedef RungeKuttaSolver< ExplicitOperatorType, ExplicitOperatorType, ExplicitOperatorType,
                                BasicLinearSolverType > OdeSolverImpl;
      return new OdeSolverImpl( tp, dgAdvectionOperator_,
                                dgAdvectionOperator_,
                                dgAdvectionOperator_,
                                name() );
    }

    //! return overal number of grid elements
    virtual uint64_t gridSize() const
    {
      adaptHandler_.globalNumberOfElements();

      size_t  advSize   = dgAdvectionOperator_.numberOfElements();
      size_t  dgIndSize = adaptHandler_.numberOfElements();
      uint64_t grSize   = std::max( advSize, dgIndSize );
      return grid_.comm().sum( grSize );
    }

    const ModelType& model() const { return dgAdvectionOperator_.model(); }

  protected:
    ExplicitOperatorType    dgAdvectionOperator_;
  };
}
}
#endif // FEMHOWTO_STEPPER_HH
