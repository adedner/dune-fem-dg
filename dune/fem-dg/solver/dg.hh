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
             class DiffusionModel,
             class Additional>
  class DGOperator : public Fem::SpaceOperatorInterface< DestinationImp >
  {
  public:
    static const Solver::Enum solverId             = Additional::solverId;
    static const Formulation::Enum formId          = Additional::formId;
    static const AdvectionLimiter::Enum limiterId  = Additional::limiterId;
    // for valid advection fluxes see dune/fem-dg/operator/fluxes/advection/parameters.hh
    static const AdvectionFlux::Enum advFluxId     = Additional::advFluxId;
    // for valid diffusion fluxes see dune/fem-dg/operator/fluxes/diffusion/parameters.hh
    static const DiffusionFlux::Enum diffFluxId    = Additional::diffFluxId;
    typedef DestinationImp   DestinationType;
    typedef typename DestinationType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    static const int polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder ;

    typedef typename DiscreteFunctionSpaceType :: GridPartType    GridPartType;
    typedef typename GridPartType::GridType                       GridType;

    typedef typename GridType :: CollectiveCommunication          CollectiveCommunicationType;
    typedef TimeProvider< CollectiveCommunicationType >           TimeProviderType;

    typedef ModelWrapper< GridType, AdvectionModel, DiffusionModel, Additional >  ModelType;

    static constexpr bool symmetric  =  false ;
    static constexpr bool matrixfree =  true  ;

    typedef DGAdvectionFlux< ModelType, advFluxId >       AdvectionFluxType;
    typedef typename DiffusionFluxSelector< ModelType, DiscreteFunctionSpaceType, diffFluxId, formId >::type  DiffusionFluxType;

    typedef DefaultOperatorTraits< ModelType, DestinationType, AdvectionFluxType, DiffusionFluxType >  OpTraits;

    typedef DGLimitedAdvectionOperator< OpTraits, Additional::hasAdvection, Additional::hasDiffusion > DGOperatorType;
    //typedef typename AdvectionDiffusionOperatorSelector< OpTraits, formId, limiterId > :: FullOperatorType
    //  DGOperatorType ;

    // typedef typename AdvectionDiffusionOperatorSelector< OpTraits, formId, limiterId > :: ExplicitOperatorType  DGOperatorType ;

    // solver selection, available fem, istl, petsc, ...
    typedef typename MatrixFreeSolverSelector< solverId, symmetric > :: template LinearInverseOperatorType< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType >  LinearSolverType ;

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
        model_( advectionModel, diffusionModel ),
        dgOperator_( space.gridPart(), model_, extra_, name() ),
        rkSolver_( tp_, dgOperator_, dgOperator_, dgOperator_, name() ),
        initialized_( false )
    {
      std::string keyPrefix("femdg.stepper.");
      const double maxTimeStep =
        Dune::Fem::Parameter::getValue< double >( keyPrefix +  "maxtimestep", std::numeric_limits<double>::max());
      fixedTimeStep_ = Dune::Fem::Parameter::getValue< double >( keyPrefix + "fixedtimestep" , 0.0 );
      // start first time step with prescribed fixed time step
      // if it is not 0 otherwise use the internal estimate
      tp_.provideTimeStepEstimate(maxTimeStep);

      // adjust fixed time step with timeprovider.factor()
      fixedTimeStep_ /= tp_.factor() ;
      if ( fixedTimeStep_ > 1e-20 )
        tp_.init( fixedTimeStep_ );
      else
        tp_.init();

      std::cout << "cfl = " << double(tp_.factor()) << std::endl;
    }

    const DiscreteFunctionSpaceType& space () const { return space_; }
    const DiscreteFunctionSpaceType& domainSpace () const { return space_; }
    const DiscreteFunctionSpaceType& rangeSpace () const { return space_; }

    //! evaluate the operator
    void operator()( const DestinationType& arg, DestinationType& dest ) const
    {
      dest.assign( arg );
      solve( dest );
    }

    void limit( DestinationType &u) { dgOperator_.limit(u); }

    void solve( DestinationType& dest ) const
    {
      if( !initialized_ )
      {
        const double dt = tp_.timeStepEstimate();
        rkSolver_.initialize( dest );
        initialized_ = true;
        tp_.provideTimeStepEstimate( dt );
      }

      // next time step is prescribed by fixedTimeStep
      if ( fixedTimeStep_ > 1e-20 )
        tp_.next( fixedTimeStep_ );
      else
        tp_.next();

      rkSolver_.solve( dest, monitor_ );
    }

    void setTimeStepSize( const double dt )
    {
      fixedTimeStep_ = dt ;
      tp_.provideTimeStepEstimate( dt );
    }

    double deltaT() const { return tp_.deltaT(); }

  protected:
    const DiscreteFunctionSpaceType&      space_;

    std::tuple<>                          extra_;
    mutable TimeProviderType              tp_;
    mutable MonitorType                   monitor_;
    ModelType                             model_;
    DGOperatorType                        dgOperator_;
    mutable RKSolverType                  rkSolver_;
    mutable double                        fixedTimeStep_ ;
    mutable bool                          initialized_;
  };

} // end namespace Fem
} // end namespace Dune
#endif
