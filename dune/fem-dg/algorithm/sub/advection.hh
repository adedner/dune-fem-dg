#ifndef DUNE_FEMDG_ALGORITHM_ADVECTION_STEPPER_HH
#define DUNE_FEMDG_ALGORITHM_ADVECTION_STEPPER_HH

#include <dune/fem-dg/misc/memory.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include <dune/fem-dg/algorithm/sub/evolution.hh>

#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Sub-Algorithm modeling an advection equation.
   *
   * \ingroup SubAlgorithms
   */
  template <class GridImp,
            class ProblemTraits,
            int polynomialOrder >
  struct SubAdvectionAlgorithm
    : public SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder >
  {
    typedef SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder > BaseType ;

    // type of Grid
    typedef typename BaseType::GridType                        GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                    GridPartType;

    // initial data type
    typedef typename BaseType::ProblemType                     ProblemType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                       ModelType;

    // The DG space operator
    // The first operator is sum of the other two
    // The other two are needed for semi-implicit time discretization
    typedef typename BaseType::OperatorType::type              FullOperatorType;

    typedef typename BaseType::SolverType::LinearSolverType    LinearSolverType;

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

    typedef typename FullOperatorType::ExtraParameterTupleType ExtraParameterTupleType;

    typedef typename BaseType::AdaptIndicatorType              AdaptIndicatorType;

    using BaseType::grid_;
    using BaseType::gridPart_;
    using BaseType::solution ;
    using BaseType::problem;
    using BaseType::name ;
    using BaseType::limitSolution;

    typedef typename BaseType::ContainerType                   ContainerType;

    // constructor
    SubAdvectionAlgorithm( GridType& grid, ContainerType& container,
                           ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
      BaseType( grid, container ),
      tuple_( ),
      advectionOperator_( Std::make_unique< FullOperatorType >( gridPart_, problem(), tuple_, name() ) ),
      adaptIndicator_( Std::make_unique< AdaptIndicatorOptional<AdaptIndicatorType> >( solution(), problem(), tuple_, name() ) )
    {}

    virtual AdaptIndicatorType* adaptIndicator () override
    {
      assert( adaptIndicator_ );
      return adaptIndicator_->value();
    }

    virtual void limit () override
    {
      if( limitSolution() )
        advectionOperator_->limit( *limitSolution() );
    }

    //! return overal number of grid elements
    virtual UInt64Type gridSize () const override
    {
      assert( advectionOperator_ );
      assert( adaptIndicator_ );

      int globalElements = adaptIndicator_ ? adaptIndicator_->globalNumberOfElements() : 0;
      if( globalElements > 0 )
        return globalElements;

      // one of them is not zero,
      size_t  advSize   = advectionOperator_->numberOfElements();
      size_t  dgIndSize = *adaptIndicator_ ? adaptIndicator_->numberOfElements() : advSize;
      UInt64Type grSize   = std::max( advSize, dgIndSize );
      return grid_.comm().sum( grSize );
    }

    virtual std::shared_ptr< SolverType > doCreateSolver ( TimeProviderType& tp ) override
    {
      assert( advectionOperator_ );
      assert( adaptIndicator_ );

      if( adaptIndicator_ )
        adaptIndicator_->setAdaptation( tp );

      typedef RungeKuttaSolver< FullOperatorType, FullOperatorType, FullOperatorType,
                                LinearSolverType > SolverImpl;
      return std::make_shared< SolverImpl >( tp, *advectionOperator_,
                                             *advectionOperator_,
                                             *advectionOperator_,
                                             name() );
    }

   const ModelType& model () const override { assert( advectionOperator_ ); return advectionOperator_->model(); }

  protected:
    ExtraParameterTupleType tuple_;


    std::unique_ptr< FullOperatorType >    advectionOperator_;
    mutable std::unique_ptr< AdaptIndicatorOptional<AdaptIndicatorType> > adaptIndicator_;
  };
}
}
#endif // FEMHOWTO_STEPPER_HH
