#ifndef DUNE_FEM_DG_OPERATOR_HH
#define DUNE_FEM_DG_OPERATOR_HH

// system includes
#include <string>

// dune-fem includes
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/dg/primaloperator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>
#include <dune/fem-dg/models/modelwrapper.hh>
#include <dune/fem-dg/misc/configurator.hh>
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
  class DGOperator : public Fem::SpaceOperatorInterface< typename Traits::DestinationType >
  {
  public:
    typedef DestinationImp   DestinationType;
    typedef typename DestinationType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    static const int polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder ;

    typedef typename DiscreteFunctionSpaceType :: GridPartType   GridPartType;
    typedef typename GridPartType                                GridType;

    typedef ModelImplementationWrapper< AdvectionModel >         ProblemType;
    typedef ModelWrapper< GridType, ProblemType >                ModelType;

    typedef DefaultOpTraits< ModelType, FunctionSpaceType, polynomialOrder >  OpTraits;

    enum { limiter = AdvectionLimiter::Enum::limited };
    enum { form    =  Formulation::Enum::primal };

    typedef typename AdvectionDiffusionOperatorSelector< OpTraits, form, limiter > :: FullOperatorType DGOperatorType ;

    // type of runge kutta solver
    typedef RungeKuttaSolver< FullOperatorType, FullOperatorType, FullOperatorType,
                              LinearSolverType > RKSolverType;

    typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

    static std::string name() const { return std::string("DGOperator"); }

    DGOperator( const DiscreteFunctionSpaceType& space )
      : extra_(),
        tp_(),
        model_(),
        dgOperator( space.gridPart(), model_, extra_, name() ),
        rkSolver_( tp_, dgOperator_, dgOperator_, dgOperator_, name() )
    {}

    //! evaluate the operator
    void operator()( const DestinationType& arg, DestinationType& dest ) const
    {
      dest.assign( arg );
      solve( dest );
      // dgOperator_( arg, dest );
    }

    void solve( DestinationType& dest ) const
    {
      rkSolver_.solve( dest, monitor_ );
    }

    void setTimeStepSize( const double dt )
    {
      tp_.provideTimeStepEstimate( dt );
    }

  protected:
    std::tuple<>          extra_;
    TimeProviderType      tp_;
    mutable MonitorType   monitor_;
    ModelType             model_;
    DGOperatorType        dgOperator_;
    mutable RKSolverType  rkSolver_;
  };

}
}
#endif
