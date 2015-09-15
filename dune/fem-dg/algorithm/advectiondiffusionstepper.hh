#ifndef DUNE_FEMDG_ALGORITHM_ADVECTIONDIFFUSION_STEPPER_HH
#define DUNE_FEMDG_ALGORITHM_ADVECTIONDIFFUSION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include <dune/fem-dg/algorithm/evolution.hh>

#include <dune/fem-dg/operator/dg/primaloperator.hh>

#include <dune/fem-dg/algorithm/handler/adapt.hh>

namespace Dune
{
namespace Fem
{

  template <class GridImp,
            class ProblemTraits,
            int polynomialOrder >
  struct AdvectionDiffusionStepper
    : public SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder >
  {

    typedef SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder > BaseType ;

    // type of Grid
    typedef typename BaseType::GridType                       GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                   GridPartType;

    // initial data type
    typedef typename BaseType::ProblemType                    ProblemType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                      ModelType;

    // The DG space operator
    // The first operator is sum of the other two
    // The other two are needed for semi-implicit time discretization
    typedef typename BaseType::OperatorType::type             OperatorType;
    typedef typename BaseType::OperatorType::ExplicitType     ExplicitOperatorType;
    typedef typename BaseType::OperatorType::ImplicitType     ImplicitOperatorType;

    typedef typename BaseType::SolverType::BasicLinearSolverType BasicLinearSolverType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType           DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;

    typedef typename BaseType::IOTupleType                    IOTupleType;

    // The ODE Solvers
    typedef typename BaseType::SolverType::type               SolverType;

    typedef typename BaseType::TimeProviderType               TimeProviderType;

    // type of 64bit unsigned integer
    typedef typename BaseType::UInt64Type                     UInt64Type;

    typedef typename BaseType::ExtraParameterTupleType        ExtraParameterTupleType;

    typedef typename BaseType::AdaptIndicatorType             AdaptIndicatorType;

    using BaseType::grid_;
    using BaseType::gridPart_;
    using BaseType::space;
    using BaseType::solution;
    using BaseType::problem;
    using BaseType::name;

    AdvectionDiffusionStepper( GridType& grid,
                               const std::string name = "",
                               ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
      BaseType( grid, name ),
      //vSpace_( gridPart_ ),
      //velo_( "velocity", vSpace_ ),
      //tuple_( &velo_ ),
      tuple_( ),
      operator_( gridPart_, problem(), tuple_, name ),
      advectionOperator_( gridPart_, problem(), tuple_, name ),
      diffusionOperator_( gridPart_, problem(), tuple_, name ),
      adaptIndicator_( solution(), problem(), tuple_, name )
    {}

    virtual AdaptIndicatorType* adaptIndicator()
    {
      return adaptIndicator_.value();
    }

    //! return overal number of grid elements
    virtual UInt64Type gridSize() const
    {
      int globalElements = adaptIndicator_ ? adaptIndicator_.globalNumberOfElements() : 0;
      if( globalElements > 0 )
        return globalElements;

      // one of them is not zero,
      size_t advSize     = advectionOperator_.numberOfElements();
      size_t diffSize    = diffusionOperator_.numberOfElements();
      size_t dgIndSize   = adaptIndicator_ ? adaptIndicator_.numberOfElements() : diffSize;
      size_t dgSize      = operator_.numberOfElements();
      UInt64Type grSize  = std::max( std::max(advSize, dgSize ), std::max( diffSize, dgIndSize ) );
      double minMax[ 2 ] = { double(grSize), 1.0/double(grSize) } ;
      grid_.comm().max( &minMax[ 0 ], 2 );
      if( Dune::Fem::Parameter::verbose () )
      {
        std::cout << "grid size (min,max) = ( " << size_t(1.0/minMax[ 1 ]) << " , " << size_t(minMax[ 0 ]) << ")" << std::endl;
      }
      return grid_.comm().sum( grSize );
    }

    virtual SolverType* createSolver(TimeProviderType& tp)
    {
      if( adaptIndicator_ )
        adaptIndicator_.setAdaptation( tp );

      // create ODE solver
      typedef RungeKuttaSolver< OperatorType, ExplicitOperatorType, ImplicitOperatorType,
                                BasicLinearSolverType > SolverImpl;
      return new SolverImpl( tp, operator_,
                             advectionOperator_,
                             diffusionOperator_,
                             name() );
    }

    const ModelType& model() const { return operator_.model(); }

  protected:
    //typename OperatorTraits::SpaceType vSpace_;
    //typename OperatorTraits::VeloType  velo_;
    ExtraParameterTupleType tuple_;

    OperatorType                         operator_;
    ExplicitOperatorType                 advectionOperator_;
    ImplicitOperatorType                 diffusionOperator_;
    mutable AdaptIndicatorOptional<AdaptIndicatorType> adaptIndicator_;
  };
}
}
#endif // FEMHOWTO_STEPPER_HH
