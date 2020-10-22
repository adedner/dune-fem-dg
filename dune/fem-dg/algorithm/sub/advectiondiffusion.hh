#ifndef DUNE_FEMDG_ALGORITHM_ADVECTIONDIFFUSION_STEPPER_HH
#define DUNE_FEMDG_ALGORITHM_ADVECTIONDIFFUSION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include <dune/fem-dg/algorithm/sub/evolution.hh>

#include <dune/fem-dg/algorithm/caller/adapt.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Sub-Algorithm modeling an advection-diffusion equation
   *
   * \ingroup SubAlgorithms
   */
  template <class GridImp,
            class ProblemTraits,
            int polynomialOrder >
  struct SubAdvectionDiffusionAlgorithm
    : public SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder >
  {

    typedef SubEvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder > BaseType ;

    // type of Grid
    typedef typename BaseType::GridType                       GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                   GridPartType;

    // The DG space operator
    // The first operator is sum of the other two
    // The other two are needed for semi-implicit time discretization
    typedef typename BaseType::OperatorType::type             OperatorType;
    typedef typename BaseType::OperatorType::ExplicitType     ExplicitOperatorType;
    typedef typename BaseType::OperatorType::ImplicitType     ImplicitOperatorType;

    typedef typename OperatorType :: AdvectionFluxType        AdvectionFluxType;

    typedef typename BaseType::SolverType::LinearSolverType   LinearSolverType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType           DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;

    // The ODE Solvers
    typedef typename BaseType::SolverType::type               SolverType;

    typedef typename BaseType::TimeProviderType               TimeProviderType;

    // type of 64bit unsigned integer
    typedef typename BaseType::UInt64Type                     UInt64Type;

    typedef typename BaseType::AdaptIndicatorType             AdaptIndicatorType;

    using BaseType::grid;
    using BaseType::solution;
    using BaseType::name;

  protected:
    using BaseType::model_;
  public:
    typedef typename BaseType::ContainerType                  ContainerType;

    template< class ContainerImp, class ExtraArgsImp >
    SubAdvectionDiffusionAlgorithm( const std::shared_ptr< ContainerImp >& cont,
                                    const std::shared_ptr< ExtraArgsImp >& extra )
    : BaseType( cont, extra ),
      numFlux_( model_ ),
      operator_( std::make_unique< OperatorType >( solution().space().gridPart(), model_, numFlux_, extra, name() ) ),
      advectionOperator_( std::make_unique< ExplicitOperatorType >( solution().space().gridPart(), model_, numFlux_, extra, name() ) ),
      diffusionOperator_( std::make_unique< ImplicitOperatorType >( solution().space().gridPart(), model_, numFlux_, extra, name() ) ),
      adaptIndicator_( std::make_unique< AdaptIndicatorOptional<AdaptIndicatorType> >( solution(), model_, numFlux_, extra, name() ) )
    {}

    virtual AdaptIndicatorType* adaptIndicator () override
    {
      assert( adaptIndicator_ );
      return adaptIndicator_->value();
    }

    //! return overal number of grid elements
    virtual UInt64Type gridSize () const override
    {
      assert( operator_ );
      assert( advectionOperator_ );
      assert( diffusionOperator_ );
      assert( adaptIndicator_ );

      int globalElements = *adaptIndicator_ ? adaptIndicator_->globalNumberOfElements() : 0;
      if( globalElements > 0 )
        return globalElements;

      // one of them is not zero,
      size_t advSize     = advectionOperator_->numberOfElements();
      size_t diffSize    = diffusionOperator_->numberOfElements();
      size_t dgIndSize   = *adaptIndicator_ ? adaptIndicator_->numberOfElements() : diffSize;
      size_t dgSize      = operator_->numberOfElements();
      UInt64Type grSize  = std::max( std::max(advSize, dgSize ), std::max( diffSize, dgIndSize ) );
      double minMax[ 2 ] = { double(grSize), 1.0/double(grSize) } ;
      grid().comm().max( &minMax[ 0 ], 2 );
      if( Dune::Fem::Parameter::verbose () )
      {
        std::cout << "grid size per core (min,max) = ( " << size_t(1.0/minMax[ 1 ]) << " , " << size_t(minMax[ 0 ]) << ")" << std::endl;
      }
      return grid().comm().sum( grSize );
    }

    virtual std::shared_ptr< SolverType > doCreateSolver ( TimeProviderType& tp ) override
    {
      assert( operator_ );
      assert( advectionOperator_ );
      assert( diffusionOperator_ );
      assert( adaptIndicator_ );

      if( adaptIndicator_ )
        adaptIndicator_->setAdaptation( tp );

      // create ODE solver
      typedef RungeKuttaSolver< OperatorType, ExplicitOperatorType, ImplicitOperatorType,
                                LinearSolverType > SolverImpl;
      return std::make_shared< SolverImpl >( tp, *operator_, *advectionOperator_, *diffusionOperator_, name() );
    }

  protected:
    AdvectionFluxType                         numFlux_;
    std::unique_ptr< OperatorType >           operator_;
    std::unique_ptr< ExplicitOperatorType >   advectionOperator_;
    std::unique_ptr< ImplicitOperatorType >   diffusionOperator_;
    mutable std::unique_ptr< AdaptIndicatorOptional<AdaptIndicatorType> > adaptIndicator_;
  };
}
}
#endif // FEMHOWTO_STEPPER_HH
