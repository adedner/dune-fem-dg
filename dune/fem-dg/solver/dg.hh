#ifndef DUNE_FEM_DG_OPERATOR_HH
#define DUNE_FEM_DG_OPERATOR_HH

// system includes
#include <string>

// dune-fem includes
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/dg/operatortraits.hh>
#include <dune/fem-dg/operator/dg/primaloperator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>
#include <dune/fem-dg/models/modelwrapper.hh>
#include <dune/fem-dg/misc/algorithmcreatorselector.hh>

namespace Dune
{
namespace Fem
{

  // DG solver operator
  //---------------------

  template < class DestinationImp,
             class AdvectionModel,
             class DiffusionModel >
  class DGOperator : public Fem::SpaceOperatorInterface< DestinationImp >
  {
  public:
    typedef DestinationImp   DestinationType;
    typedef typename DestinationType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    static const int polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder ;

    typedef typename DiscreteFunctionSpaceType :: GridPartType   GridPartType;
    typedef typename GridPartType::GridType                      GridType;

    typedef typename GridType :: CollectiveCommunicationType     CollectiveCommunicationType;
    typedef TimeProvider< CollectiveCommunicationType >          TimeProviderType;

    typedef ModelImplementationWrapper< AdvectionModel >         ProblemType;
    typedef ModelWrapper< GridType, ProblemType >                ModelType;

    //enum { limiterId   =  AdvectionLimiter::Enum::limited };
    //enum { formId      =  Formulation::Enum::primal };
    //enum { advFluxId   =  AdvectionFlux::Enum::llf };
    //enum { diffFluxId  =  DiffusionFlux::Enum::primal };
    //enum { solverId    =  Solver::Enum::fem };
    static constexpr bool symmetric  =  false ;
    static constexpr bool matrixfree =  true  ;

    typedef DGAdvectionFlux< ModelType, AdvectionFlux::Enum::llf >       AdvectionFluxType;
    typedef typename DiffusionFluxSelector< ModelType, DiscreteFunctionSpaceType, DiffusionFlux::Enum::primal, Formulation::Enum::primal >::type  DiffusionFluxType;

    typedef DefaultOperatorTraits< ModelType, DestinationType, AdvectionFluxType, DiffusionFluxType >  OpTraits;

    typedef typename AdvectionDiffusionOperatorSelector< OpTraits, Formulation::Enum::primal, AdvectionLimiter::Enum::limited > :: FullOperatorType
      DGOperatorType ;

    // solver selection, available fem, istl, petsc, ...
    typedef typename SolverSelector< Solver::Enum::fem, symmetric, matrixfree > :: type  LinearSolverType ;

    typedef DuneODE::OdeSolverInterface< DestinationType >      OdeSolverInterfaceType;
    // type of runge kutta solver
    typedef RungeKuttaSolver< DGOperatorType, DGOperatorType, DGOperatorType,
                              LinearSolverType > RKSolverType;

    typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

    static std::string name() { return std::string("DGOperator"); }

    DGOperator( const DiscreteFunctionSpaceType& space,
                const AdvectionModel &advectionModel,
                const DiffusionModel &diffusionModel )
      : space_( space ),
        extra_(),
        tp_( space_.gridPart().comm() ),
        model_(advectionModel),
        dgOperator_( space.gridPart(), model_, extra_, name() ),
        rkSolver_( tp_, dgOperator_, dgOperator_, dgOperator_, name() )
    {}

    //! evaluate the operator
    void operator()( const DestinationType& arg, DestinationType& dest ) const
    {
      dest.assign( arg );
      solve( dest );
      // dgOperator_( arg, dest );
    }

    const DiscreteFunctionSpaceType& space () const { return space_; }
    const DiscreteFunctionSpaceType& domainSpace () const { return space_; }
    const DiscreteFunctionSpaceType& rangeSpace () const { return space_; }

    void solve( DestinationType& dest ) const
    {
      rkSolver_.solve( dest, monitor_ );
    }

    void setTimeStepSize( const double dt )
    {
      tp_.provideTimeStepEstimate( dt );
    }

  protected:
    const DiscreteFunctionSpaceType&      space_;

    std::tuple<>                          extra_;
    TimeProviderType                      tp_;
    mutable MonitorType                   monitor_;
    ModelType                             model_;
    DGOperatorType                        dgOperator_;
    mutable RKSolverType                  rkSolver_;
  };

}
}
#endif
