#ifndef DUNE_FEMDG_ALGORITHM_ADVECTION_STEPPER_HH
#define DUNE_FEMDG_ALGORITHM_ADVECTION_STEPPER_HH

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

    // The DG space operator
    // The first operator is sum of the other two
    // The other two are needed for semi-implicit time discretization
    typedef typename BaseType::OperatorType::type              FullOperatorType;

    // type of advective flux implementation
    typedef typename FullOperatorType :: AdvectionFluxType     AdvectionFluxType;

    typedef typename BaseType::SolverType::LinearSolverType    LinearSolverType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType            DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType       DiscreteFunctionSpaceType;

    // The ODE Solvers
    typedef typename BaseType::SolverType::type                SolverType;

    typedef typename BaseType::TimeProviderType                TimeProviderType;

     // type of 64bit unsigned integer
    typedef typename BaseType::UInt64Type                      UInt64Type;

    typedef typename BaseType::AdaptIndicatorType              AdaptIndicatorType;

  private:
    using BaseType::model_;
  public:
    using BaseType::solution ;
    using BaseType::name ;
    using BaseType::limitSolution;

    typedef typename BaseType::ContainerType                   ContainerType;

    template< class ContainerImp, class ExtraArgsImp >
    SubAdvectionAlgorithm( const std::shared_ptr< ContainerImp >& cont,
                           const std::shared_ptr< ExtraArgsImp >& extra )
    : BaseType( cont, extra ),
      numFlux_( model_ ),
      advectionOperator_( std::make_unique< FullOperatorType >( solution().space().gridPart(), model_, numFlux_, extra, name() ) ),
      adaptIndicator_( std::make_unique< AdaptIndicatorOptional<AdaptIndicatorType> >( solution(), model_, numFlux_, extra, name() ) )
    {
    }


    virtual AdaptIndicatorType* adaptIndicator () override
    {
      assert( adaptIndicator_ );
      return adaptIndicator_->value();
    }

    virtual void limit () override
    {
      if( limitSolution() )
        advectionOperator_->applyLimiter( *limitSolution() );
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
      return solution().space().gridPart().grid().comm().sum( grSize );
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

  protected:
    AdvectionFluxType                      numFlux_;
    std::unique_ptr< FullOperatorType >    advectionOperator_;
    mutable std::unique_ptr< AdaptIndicatorOptional<AdaptIndicatorType> > adaptIndicator_;
  };
}
}
#endif // FEMHOWTO_STEPPER_HH
