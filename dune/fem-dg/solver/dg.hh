#ifndef DUNE_FEM_DG_OPERATOR_HH
#define DUNE_FEM_DG_OPERATOR_HH

// system includes
#include <string>

// dune-fem includes
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#include <dune/fem-dg/algorithm/evolution.hh>

// dune-fem-dg includes
#include <dune/fem-dg/algorithm/evolution.hh>
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
#include <dune/fem-dg/operator/dg/primaloperator.hh>
#include <dune/fem-dg/solver/rungekuttasolver.hh>
#include <dune/fem-dg/models/modelwrapper.hh>
#include <dune/fem-dg/misc/algorithmcreatorselector.hh>

#ifdef EULER_WRAPPER_TEST
#include <dune/fem-dg/models/additional.hh>
#endif

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
#ifdef EULER_WRAPPER_TEST
  class DGOperator : public DuneODE :: OdeSolverInterface< DestinationImp >
#else
  class DGOperator : public Fem::SpaceOperatorInterface< DestinationImp >
#endif
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

#ifdef EULER_WRAPPER_TEST
    typedef U0Sod< GridType > Problem;
#endif

    typedef ModelWrapper< GridType, AdvectionModel, DiffusionModel, Additional
#ifdef EULER_WRAPPER_TEST
        , Problem
#endif
      >  ModelType;

    static constexpr bool symmetric  =  false ;
    static constexpr bool matrixfree =  true  ;

    typedef DGAdvectionFlux< ModelType, advFluxId >       AdvectionFluxType;
    typedef typename DiffusionFluxSelector< ModelType, DiscreteFunctionSpaceType, diffFluxId, formId >::type  DiffusionFluxType;

    typedef DefaultOperatorTraits< ModelType, DestinationType, AdvectionFluxType, DiffusionFluxType >  OpTraits;

    typedef AdvectionDiffusionOperatorSelector< OpTraits, formId, limiterId > OperatorSelectorType ;

    typedef typename OperatorSelectorType :: FullOperatorType      FullOperatorType;
    typedef typename OperatorSelectorType :: ExplicitOperatorType  ExplicitOperatorType;
    typedef typename OperatorSelectorType :: ImplicitOperatorType  ImplicitOperatorType;

    // solver selection, available fem, istl, petsc, ...
    typedef typename MatrixFreeSolverSelector< solverId, symmetric > :: template LinearInverseOperatorType< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType >  LinearSolverType ;

    typedef DuneODE::OdeSolverInterface< DestinationType >      OdeSolverInterfaceType;
    // type of runge kutta solver
    typedef RungeKuttaSolver< FullOperatorType, ExplicitOperatorType, ImplicitOperatorType,
                              LinearSolverType > RKSolverType;

    typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

    static std::string name() { return std::string(""); }

    DGOperator( const DiscreteFunctionSpaceType& space,
                const AdvectionModel &advectionModel,
                const DiffusionModel &diffusionModel,
                const TimeSteppingParameters& param = TimeSteppingParameters() )
      : space_( space ),
        extra_(),
        tpPtr_( new TimeProviderType(space_.gridPart().comm()) ),
        tp_( *tpPtr_ ),
        model_( advectionModel, diffusionModel ),
        fullOperator_( space.gridPart(), model_, extra_, name() ),
        explOperator_( space.gridPart(), model_, extra_, name() ),
        implOperator_( space.gridPart(), model_, extra_, name() ),
        rkSolver_( tp_, fullOperator_, explOperator_, implOperator_, name() ),
        initialized_( false )
    {
      const double maxTimeStep = param.maxTimeStep();
      fixedTimeStep_ = param.fixedTimeStep();
      // start first time step with prescribed fixed time step
      // if it is not 0 otherwise use the internal estimate
      tp_.provideTimeStepEstimate(maxTimeStep);

      // adjust fixed time step with timeprovider.factor()
      fixedTimeStep_ /= tp_.factor() ;
      if ( fixedTimeStep_ > 1e-20 )
        tp_.init( fixedTimeStep_ );
      else
        tp_.init();

      std::cout << "cfl = " << double(tp_.factor()) << " " << tp_.time() << std::endl;
    }

    DGOperator( TimeProviderType& tp,
                const DiscreteFunctionSpaceType& space,
                const AdvectionModel &advectionModel,
                const DiffusionModel &diffusionModel,
                const TimeSteppingParameters& param = TimeSteppingParameters())
      : space_( space ),
        extra_(),
        tp_( tp ),
        model_( advectionModel, diffusionModel ),
        fullOperator_( space.gridPart(), model_, extra_, name() ),
        explOperator_( space.gridPart(), model_, extra_, name() ),
        implOperator_( space.gridPart(), model_, extra_, name() ),
        rkSolver_( tp_, fullOperator_, explOperator_, implOperator_, name() ),
        initialized_( false )
    {
      std::cout << "cfl = " << double(tp_.factor()) << " " << tp_.time() << std::endl;
    }

    virtual void initialize( const DestinationType& dest )
    {
      checkInitialize( dest );
    }

    virtual void description( std::ostream&) const {}

    const DiscreteFunctionSpaceType& space () const { return space_; }
    const DiscreteFunctionSpaceType& domainSpace () const { return space_; }
    const DiscreteFunctionSpaceType& rangeSpace () const { return space_; }

#ifndef EULER_WRAPPER_TEST
    //! evaluate the operator
    void operator()( const DestinationType& arg, DestinationType& dest ) const
    {
      dest.assign( arg );
      solve( dest );
    }
#endif

    void limit( DestinationType &u) const { explOperator_.limit(u); }

#ifdef EULER_WRAPPER_TEST
    void solve( DestinationType& dest, MonitorType& )
    {
      solve( dest );
    }
#endif

    void checkInitialize( const DestinationType& dest ) const
    {
      if( !initialized_ )
      {
        //const double dt = tp_.timeStepEstimate();
        rkSolver_.initialize( dest );
        initialized_ = true;
        // tp_.provideTimeStepEstimate( dt );
      }
    }

    void solve( DestinationType& dest ) const
    {
      checkInitialize( dest );

      //std::cout << "Start operator solve!" << std::endl;
      Dune :: Timer timer;
      // limit( dest );
      rkSolver_.solve( dest, monitor_ );

      //std::cout << "Operator solve called: " << timer.elapsed() << std::endl;

      if( tpPtr_ )
      {
        // next time step is prescribed by fixedTimeStep
        if ( fixedTimeStep_ > 1e-20 )
          tp_.next( fixedTimeStep_ );
        else
          tp_.next();
      }

       //       // reset time step estimate
       // tp.provideTimeStepEstimate( maxTimeStep );

#ifndef EULER_WRAPPER_TEST
       // return limited solution, to be discussed
       limit( dest );
#endif
    }

    void setTimeStepSize( const double dt )
    {
      fixedTimeStep_ = dt ;
      fixedTimeStep_ /= tp_.factor() ;
      tp_.provideTimeStepEstimate( dt );
    }

    double deltaT() const { return tp_.deltaT(); }

  protected:
    const DiscreteFunctionSpaceType&      space_;

    std::tuple<>                          extra_;
    std::unique_ptr< TimeProviderType >   tpPtr_;
    TimeProviderType&                     tp_;
    mutable MonitorType                   monitor_;
    ModelType                             model_;
    mutable FullOperatorType              fullOperator_;
    mutable ExplicitOperatorType          explOperator_;
    mutable ImplicitOperatorType          implOperator_;
    mutable RKSolverType                  rkSolver_;
    mutable double                        fixedTimeStep_ ;
    mutable bool                          initialized_;
  };

} // end namespace Fem
} // end namespace Dune
#endif
