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

    // initial data type
    typedef typename BaseType::ProblemType                    ProblemType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                       ModelType;

    // The DG space operator
    // The first operator is sum of the other two
    // The other two are needed for semi-implicit time discretization
    typedef typename BaseType::OperatorType::type             FullOperatorType;
    typedef typename BaseType::OperatorType::ExplicitType     ExplicitOperatorType;
    typedef typename BaseType::OperatorType::ImplicitType     ImplicitOperatorType;

    typedef typename BaseType::SolverType::BasicLinearSolverType BasicLinearSolverType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType            DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType       DiscreteFunctionSpaceType;

    typedef typename BaseType::IOTupleType                     IOTupleType;

    // The ODE Solvers
    typedef typename BaseType::SolverType::type                SolverType;

    typedef typename BaseType::TimeProviderType                TimeProviderType;

     // type of 64bit unsigned integer
    typedef typename BaseType::UInt64Type                      UInt64Type;

    typedef typename BaseType::ExtraParameterTupleType    ExtraParameterTupleType;

    typedef typename BaseType::AdaptIndicatorType             AdaptIndicatorType;

    using BaseType::grid_;
    using BaseType::gridPart_;
    using BaseType::space ;
    using BaseType::solution ;
    using BaseType::problem;
    using BaseType::name ;
    using BaseType::limitSolution;

    // constructor
    AdvectionStepper( GridType& grid, const std::string name = "",
                      ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
      BaseType( grid, name ),
      tuple_( ),
      advectionOperator_(gridPart_, problem(), tuple, name ),
      adaptIndicator_( solution(), problem(), tuple_, name )
    {}

    virtual AdaptIndicatorType* adaptIndicator()
    {
      return &adaptIndicator_;
    }

    virtual void limit()
    {
      if( limitSolution() )
        advectionOperator_.limit( *limitSolution() );
    }

    //! return overal number of grid elements
    virtual UInt64Type gridSize() const
    {
      int globalElements = adaptIndicator_.globalNumberOfElements();
      if( globalElements > 0 )
        return globalElements;

      // one of them is not zero,
      size_t  advSize   = advectionOperator_.numberOfElements();
      size_t  dgIndSize = adaptIndicator_.numberOfElements();
      UInt64Type grSize   = std::max( advSize, dgIndSize );
      return grid_.comm().sum( grSize );
    }

    virtual SolverType* createSolver(TimeProviderType& tp)
    {
      adaptIndicator_.setAdaptation( tp );

      typedef RungeKuttaSolver< ExplicitOperatorType, ExplicitOperatorType, ExplicitOperatorType,
                                BasicLinearSolverType > SolverImpl;
      return new SolverImpl( tp, advectionOperator_,
                             advectionOperator_,
                             advectionOperator_,
                             name() );
    }

   const ModelType& model() const { return advectionOperator_.model(); }

  protected:
    ExtraParameterTupleType tuple_;


    ExplicitOperatorType    advectionOperator_;
    mutable AdaptIndicatorType       adaptIndicator_;
  };
}
}
#endif // FEMHOWTO_STEPPER_HH
