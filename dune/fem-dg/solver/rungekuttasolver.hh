#ifndef DUNE_FEM_DG_RUNGEKUTTA_HH
#define DUNE_FEM_DG_RUNGEKUTTA_HH

#include <limits>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/solver/rungekutta/explicit.hh>
#include <dune/fem/solver/rungekutta/implicit.hh>
#include <dune/fem/solver/rungekutta/semiimplicit.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/solver/pardginverseoperators.hh>

#include <dune/fem/operator/dghelmholtz.hh>
#include <dune/fem-dg/solver/smartodesolver.hh>
#include <dune/fem-dg/misc/parameterkey.hh>

namespace Dune
{
namespace Fem
{


  /**
   * \brief Runge-Kutta solver.
   *
   * \ingroup Solvers
   */
  template <class Operator,
            class ExplicitOperator,
            class ImplicitOperator,
            class LinearInverseOperator>
  class RungeKuttaSolver :
    public DuneODE :: OdeSolverInterface< typename Operator :: DestinationType >
  {
    typedef Operator          OperatorType;
    typedef ExplicitOperator  ExplicitOperatorType;
    typedef ImplicitOperator  ImplicitOperatorType;

    typedef DuneODE :: OdeSolverInterface< typename Operator :: DestinationType > BaseType;

    template <class AdvOp, class DiffOp>
    struct IsSame
    {
      static bool check(const AdvOp&, const DiffOp& )
      {
        return false;
      }
    };

    template <class AdvOp>
    struct IsSame< AdvOp, AdvOp>
    {
      static bool check(const AdvOp& a, const AdvOp& d)
      {
        return true;
      }
    };

  public:
    typedef typename OperatorType :: DestinationType DestinationType ;
    typedef DestinationType  DiscreteFunctionType;

    typedef DuneODE :: OdeSolverInterface< DestinationType >        OdeSolverInterfaceType;
    typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;


    typedef Dune::Fem::Operator< DestinationType, DestinationType > HelmHoltzOperatorType ;

    typedef std::pair< OdeSolverInterfaceType* ,  HelmHoltzOperatorType* > solverpair_t ;

    typedef SmartOdeSolverParameters ParameterType;

    /////////////////////////////////////////////////////////////////////////
    //  ODE solvers from dune-fem/dune/fem/solver/rungekutta
    /////////////////////////////////////////////////////////////////////////
    template < class Op, class DF, bool pardgOdeSolver >
    struct OdeSolverSelection
    {
      template < class OdeParameter >
      static solverpair_t
      createExplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps, const OdeParameter& param, const std::string& name = ""  )
      {
        typedef DuneODE :: ExplicitRungeKuttaSolver< DiscreteFunctionType >          ExplicitOdeSolverType;
        return solverpair_t( new ExplicitOdeSolverType( op, tp, rkSteps ), nullptr );
      }

      template < class OdeParameter >
      static solverpair_t
      createImplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps, const OdeParameter& param, const std::string& name = "" )
      {
#ifdef COUNT_FLOPS
        return solverpair_t();
#else
        typedef Dune::Fem::DGHelmholtzOperator< Op >  HelmholtzOperatorType;
        HelmholtzOperatorType* helmOp = new HelmholtzOperatorType( op );

        typedef Dune::Fem::NewtonInverseOperator<
                      typename HelmholtzOperatorType::JacobianOperatorType,
                      LinearInverseOperator > NonlinearInverseOperatorType;

        typedef DuneODE::ImplicitRungeKuttaSolver< HelmholtzOperatorType,
                      NonlinearInverseOperatorType > ImplicitOdeSolverType;

        typedef typename NonlinearInverseOperatorType::ParametersType NonlinParametersType;

        return solverpair_t(new ImplicitOdeSolverType( *helmOp, tp, rkSteps, param, NonlinParametersType( ParameterKey::generate( name, "fem.solver.newton." ) ) ), helmOp );
#endif
      }

      template < class ExplOp, class ImplOp, class OdeParameter >
      static solverpair_t
      createSemiImplicitSolver( ExplOp& explOp, ImplOp& implOp, Fem::TimeProviderBase& tp,
                                const int rkSteps, const OdeParameter& param, const std::string& name = "" )
      {
#ifdef COUNT_FLOPS
        return solverpair_t();
#else
        typedef Dune::Fem::DGHelmholtzOperator< ImplOp >  HelmholtzOperatorType;
        HelmholtzOperatorType* helmOp = new HelmholtzOperatorType( implOp );

        typedef Dune::Fem::NewtonInverseOperator<
                      typename HelmholtzOperatorType::JacobianOperatorType,
                      LinearInverseOperator > NonlinearInverseOperatorType;

        typedef DuneODE::SemiImplicitRungeKuttaSolver< ExplicitOperatorType,
                      HelmholtzOperatorType, NonlinearInverseOperatorType > SemiImplicitOdeSolverType ;


        typedef typename NonlinearInverseOperatorType::ParametersType NonlinParametersType;


        return solverpair_t(new SemiImplicitOdeSolverType( explOp, *helmOp, tp, rkSteps, param, NonlinParametersType( ParameterKey::generate( name, "fem.solver.newton." ) ) ), helmOp );
#endif
      }

    };

    /////////////////////////////////////////////////////////////////////////
    //  parDG ODE solvers based on double*
    /////////////////////////////////////////////////////////////////////////
    template < class Op >
    struct OdeSolverSelection< Op, Fem :: AdaptiveDiscreteFunction< typename Op::SpaceType >, true >
    {
      typedef Fem :: AdaptiveDiscreteFunction< typename Operator::SpaceType >  DiscreteFunctionType ;
      // old ode solver based on double* from pardg
      typedef DuneODE :: ImplicitOdeSolver< DiscreteFunctionType >       ImplicitOdeSolverType;
      typedef DuneODE :: ExplicitOdeSolver< DiscreteFunctionType >       ExplicitOdeSolverType;
      typedef DuneODE :: SemiImplicitOdeSolver< DiscreteFunctionType >   SemiImplicitOdeSolverType;

      template < class OdeParameter >
      static solverpair_t
      createExplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps, const OdeParameter& param, const std::string& name = ""  )
      {
        return solverpair_t(new ExplicitOdeSolverType( op, tp, rkSteps ), nullptr );
      }

      template < class OdeParameter >
      static solverpair_t
      createImplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps, const OdeParameter& param, const std::string& name = "" )
      {
        return solverpair_t(new ImplicitOdeSolverType( op, tp, rkSteps, param), nullptr );
      }

      template <class ExplOp, class ImplOp, class OdeParameter >
      static solverpair_t
      createSemiImplicitSolver( ExplOp& explOp, ImplOp& implOp, Fem::TimeProviderBase& tp, const int rkSteps, const OdeParameter& param, const std::string& name = "" )
      {
        return solverpair_t(new SemiImplicitOdeSolverType( explOp, implOp, tp, rkSteps, param), nullptr );
      }

    };

    static const bool useParDGSolvers = false ;
    typedef OdeSolverSelection< OperatorType, DestinationType, useParDGSolvers >  OdeSolversType ;

    using BaseType :: solve;
  protected:
    OperatorType&           operator_;
    Fem::TimeProviderBase&  timeProvider_;

    const std::string      name_;
    ExplicitOperatorType&  explicitOperator_;
    ImplicitOperatorType&  implicitOperator_;
    const SmartOdeSolverParameters* param_;

    OdeSolverInterfaceType* explicitSolver_;
    OdeSolverInterfaceType* odeSolver_;
    HelmHoltzOperatorType * helmholtzOperator_;

    const double explFactor_ ;
    const int verbose_ ;
    const int rkSteps_ ;
    const int odeSolverType_ ;
    int imexCounter_ , exCounter_;
    int minIterationSteps_, maxIterationSteps_ ;
    bool useImex_ ;
  public:
    RungeKuttaSolver( Fem::TimeProviderBase& tp,
                      OperatorType& op,
                      ExplicitOperatorType& advOp,
                      ImplicitOperatorType& diffOp,
                      const std::string name = "" )
     : operator_( op ),
       timeProvider_( tp ),
       name_( name ),
       explicitOperator_( advOp ),
       implicitOperator_( diffOp ),
       param_( new SmartOdeSolverParameters( ParameterKey::generate( name_, "fem.ode." ) ) ),
       explicitSolver_( 0 ),
       odeSolver_( 0 ),
       helmholtzOperator_( 0 ),
       explFactor_( param_->explicitFactor() ),
       verbose_( param_->verbose() ),
       rkSteps_( param_->obtainRungeKuttaSteps( operator_.space().order() + 1 ) ),
       odeSolverType_( param_->obtainOdeSolverType() ),
       imexCounter_( 0 ), exCounter_ ( 0 ),
       minIterationSteps_( std::numeric_limits< int > :: max() ),
       maxIterationSteps_( 0 ),
       useImex_( odeSolverType_ > 1 )
    {
      solverpair_t solver( nullptr, nullptr ) ;
      // create implicit or explicit ode solver
      if( odeSolverType_ == 0 )
      {
        solver = OdeSolversType :: createExplicitSolver( operator_, tp, rkSteps_, *param_, name_ );
      }
      else if (odeSolverType_ == 1)
      {
        solver = OdeSolversType :: createImplicitSolver( operator_, tp, rkSteps_, *param_, name_ );
      }
      else if( odeSolverType_ > 1 )
      {
        // make sure that advection and diffusion operator are different
        if( IsSame< ExplicitOperatorType, ImplicitOperatorType >::check( explicitOperator_, implicitOperator_ ) )
        {
          DUNE_THROW(Dune::InvalidStateException,"Advection and Diffusion operator are the same, therefore IMEX cannot work!");
        }

        solver = OdeSolversType :: createSemiImplicitSolver( explicitOperator_, implicitOperator_, tp, rkSteps_, *param_, name_ );

        // IMEX+
        if( odeSolverType_ == 3 )
          explicitSolver_ = OdeSolversType :: createExplicitSolver( operator_, tp, rkSteps_, *param_, name_ ).first;
      }
      else
      {
        DUNE_THROW(NotImplemented,"Wrong ODE solver selected");
      }

      odeSolver_ = solver.first;
      helmholtzOperator_ = solver.second;
    }

    //! destructor
    ~RungeKuttaSolver()
    {
      delete param_;           param_ = 0;
      delete odeSolver_;       odeSolver_ = 0;
      delete explicitSolver_ ; explicitSolver_ = 0;
      delete helmholtzOperator_; helmholtzOperator_ = 0;
    }

    //! initialize method
    void initialize( const DestinationType& U )
    {
      if( explicitSolver_ )
      {
        explicitSolver_->initialize( U );
      }
      assert( odeSolver_ );
      odeSolver_->initialize( U );
    }

    void getAdvectionDiffsionTimeSteps( double& advStep, double& diffStep ) const
    {
      if( odeSolverType_ > 1 )
      {
        double steps[ 2 ] = { explicitOperator_.timeStepEstimate(),
                              implicitOperator_.timeStepEstimate() };
        // get global min
        Dune :: Fem :: MPIManager :: comm().min( &steps[ 0 ], 2 );

        advStep  = steps[ 0 ];
        diffStep = steps[ 1 ];
      }
    }

    //! solver the ODE
    void solve( DestinationType& U ,
                MonitorType& monitor )
    {
      // take CPU time of solution process
      Dune::Timer timer ;

      // switch upwind direction
      operator_.switchupwind();
      explicitOperator_.switchupwind();
      implicitOperator_.switchupwind();

      // reset compute time counter
      resetComputeTime();

      double maxAdvStep  = std::numeric_limits< double > :: max();
      double maxDiffStep = std::numeric_limits< double > :: max();

      if( explicitSolver_ && ! useImex_ )
      {
        explicitSolver_->solve( U, monitor );
        ++exCounter_ ;
      }
      else
      {
        assert( odeSolver_ );
        odeSolver_->solve( U, monitor );

        ++imexCounter_ ;

        const int iterationSteps = monitor.newtonIterations_ * monitor.linearSolverIterations_ ;
        minIterationSteps_ = std::min( minIterationSteps_, iterationSteps );
        maxIterationSteps_ = std::max( maxIterationSteps_, iterationSteps );

        if( verbose_ == 3 && MPIManager::rank()<=0 )
        {
          // get advection and diffusion time step
          getAdvectionDiffsionTimeSteps( maxAdvStep, maxDiffStep );

          double factor = explFactor_ ;
          //if( averageIterationSteps > 0 )
          //  factor *= averageIterationSteps / (rkSteps_ + 1 ) ;
          std::cout << maxAdvStep << " a | d " << maxDiffStep << "  factor: " << factor
            << "  " << minIterationSteps_ << " min | max " << maxIterationSteps_ << "  use imex = " << useImex_ << std::endl;
        }
      }

      if( explicitSolver_ )
      {
        // get advection and diffusion time step
        getAdvectionDiffsionTimeSteps( maxAdvStep, maxDiffStep );

        const int averageIterationSteps = (minIterationSteps_ + maxIterationSteps_) / 2;
        double factor = explFactor_ ;
        if( averageIterationSteps > 0 )
          factor *= double(rkSteps_ + 1) / double(averageIterationSteps) ;

        // if true solve next time step with semi implicit solver
        useImex_ = ( maxDiffStep < (factor * maxAdvStep) ) ;

        if( verbose_ == 3 && MPIManager::rank()<=0 )
        {
          std::cout << maxAdvStep << " a | d " << maxDiffStep << "  factor: " << factor
            << "  " << minIterationSteps_ << " min | max " << maxIterationSteps_
            << "  use imex = " << useImex_ << "  ex steps: " << exCounter_ << std::endl;
        }

        // make sure the correct time step is used for the explicit solver
        //if( ! useImex_ )
        //  timeProvider_.provideTimeStepEstimate( operator_.timeStepEstimate() ) ;
      }

      // store needed time
      monitor.odeSolveTime_     = timer.elapsed();
      monitor.operatorTime_     = operatorTime();
      monitor.numberOfElements_ = numberOfElements();
    }

    //! return CPU time needed for the operator evaluation
    double operatorTime() const
    {
      if( useImex_ )
      {
        return explicitOperator_.computeTime() +
               implicitOperator_.computeTime() ;
      }
      else
        return operator_.computeTime();
    }

    //! return number of elements meat during operator evaluation
    size_t numberOfElements() const
    {
      if( useImex_ )
        return explicitOperator_.numberOfElements();
      else
        return operator_.numberOfElements();
    }

    void description(std::ostream&) const {}

    // gather information from the space operator, the time integratior
    // and the problem to output before each table in tex file
    std::string description() const
    {
      std::string latexInfo;

      if ((odeSolverType_==0) || (odeSolverType_==1))
        latexInfo = operator_.description();
      else
        latexInfo = explicitOperator_.description()
                    + implicitOperator_.description();

      std::stringstream odeInfo;
      odeSolver_->description( odeInfo );

      latexInfo += odeInfo.str() + "\n";
      std::stringstream info;
      info << "Regular  Solver used: " << imexCounter_ << std::endl;
      info << "Explicit Solver used: " << exCounter_ << std::endl;
      latexInfo += info.str();

      return latexInfo;
    }

  protected:
    void resetComputeTime() const
    {
      // this will reset the internal time counters
      operator_.computeTime() ;
      explicitOperator_.computeTime();
      implicitOperator_.computeTime();
    }
  }; // end RungeKuttaSolver

} // end namespace
} // end namespace Dune
#endif
