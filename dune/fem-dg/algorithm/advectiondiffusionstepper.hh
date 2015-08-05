#ifndef DUNE_FEMDG_ALGORITHM_ADVECTIONDIFFUSION_STEPPER_HH
#define DUNE_FEMDG_ALGORITHM_ADVECTIONDIFFUSION_STEPPER_HH

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// local includes
#include <dune/fem-dg/algorithm/evolution.hh>

#include <dune/fem-dg/operator/dg/primaloperator.hh>

namespace Dune
{
namespace Fem
{

  template <class GridImp,
            class ProblemTraits,
            int polynomialOrder >
  struct AdvectionDiffusionStepper
    : public EvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder >
  {
    typedef EvolutionAlgorithm< GridImp, ProblemTraits, polynomialOrder > BaseType ;

    // type of Grid
    typedef typename BaseType :: GridType                       GridType;

    // Choose a suitable GridView
    typedef typename BaseType :: GridPartType                   GridPartType;

    // initial data type
    typedef typename BaseType :: ProblemType                    ProblemType;

    // An analytical version of our model
    typedef typename BaseType :: ModelType                      ModelType;

    // The DG space operator
    // The first operator is sum of the other two
    // The other two are needed for semi-implicit time discretization
    typedef typename BaseType :: FullOperatorType               FullOperatorType;
    typedef typename BaseType :: ExplicitOperatorType           ExplicitOperatorType;
    typedef typename BaseType :: ImplicitOperatorType           ImplicitOperatorType;

    typedef typename BaseType :: BasicLinearSolverType          BasicLinearSolverType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType :: DiscreteFunctionType           DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType :: DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;

    // The ODE Solvers
    typedef typename BaseType :: OdeSolverType                  OdeSolverType;

    typedef typename BaseType :: TimeProviderType               TimeProviderType;

    typedef typename BaseType::OperatorTraits                   OperatorTraits;

    typedef typename std::remove_pointer<typename std::tuple_element<0,typename BaseType::IndicatorTupleType>::type>::type  DGIndicatorType;
    typedef typename std::remove_pointer<typename std::tuple_element<1,typename BaseType::IndicatorTupleType>::type>::type GradientIndicatorType ;

    // type of 64bit unsigned integer
    typedef typename BaseType :: UInt64Type                     UInt64Type;

    typedef typename OperatorTraits :: ExtraParameterTupleType  ExtraParameterTupleType;

    using BaseType :: grid_;
    using BaseType :: gridPart_;
    using BaseType :: space;
    using BaseType :: problem;
    using BaseType :: name ;
    using BaseType :: adaptHandler_ ;

    AdvectionDiffusionStepper( GridType& grid,
                               const std::string name = "",
                               ExtraParameterTupleType tuple = ExtraParameterTupleType() ) :
      BaseType( grid, name ),
      //vSpace_( gridPart_ ),
      //velo_( "velocity", vSpace_ ),
      //tuple_( &velo_ ),
      tuple_( ),
      dgOperator_( gridPart_, problem(), tuple_, name ),
      dgAdvectionOperator_( gridPart_, problem(), tuple_, name ),
      dgDiffusionOperator_( gridPart_, problem(), tuple_, name ),
      dgIndicator_( gridPart_, problem(), tuple_, name ),
      gradientIndicator_( space(), problem(), adaptHandler_.params() )
    {
      adaptHandler_.setIndicator( &dgIndicator_, &gradientIndicator_ );
    }

    //! return overal number of grid elements
    virtual UInt64Type gridSize() const
    {
      int globalElements = adaptHandler_.globalNumberOfElements();
      if( globalElements > 0 )
        return globalElements;

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
      adaptHandler_.setAdaptation( tp, dgIndicator_ );

      // create ODE solver
      typedef RungeKuttaSolver< FullOperatorType, ExplicitOperatorType, ImplicitOperatorType,
                                BasicLinearSolverType > OdeSolverImpl;
      return new OdeSolverImpl( tp, dgOperator_,
                                dgAdvectionOperator_,
                                dgDiffusionOperator_,
                                name() );
    }

    const ModelType& model() const { return dgOperator_.model(); }

  protected:
    //typename OperatorTraits::SpaceType vSpace_;
    //typename OperatorTraits::VeloType  velo_;
    ExtraParameterTupleType tuple_;

    FullOperatorType        dgOperator_;
    ExplicitOperatorType    dgAdvectionOperator_;
    ImplicitOperatorType    dgDiffusionOperator_;
    DGIndicatorType         dgIndicator_;
    GradientIndicatorType   gradientIndicator_;
  };

}
}
#endif // FEMHOWTO_STEPPER_HH
