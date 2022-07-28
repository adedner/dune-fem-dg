#ifndef DUNE_FEM_DG_RUNGEKUTTA_HH
#define DUNE_FEM_DG_RUNGEKUTTA_HH

#include <limits>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/solver/rungekutta/explicit.hh>
#include <dune/fem/solver/rungekutta/implicit.hh>
#include <dune/fem/solver/rungekutta/semiimplicit.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/solver/timeprovider.hh>

#include <dune/fem/operator/dghelmholtz.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/timesteppingparameter.hh>

namespace Dune
{
namespace Fem
{
  /**
   * \brief Parameter class for selection of ode solver.
   *
   * \ingroup ParameterClass
   */
  struct RungeKuttaSolverParameters
    : public DuneODE::ImplicitRungeKuttaSolverParameters
  {
    /**
     * \brief Constructor
     *
     * \param[in] keyPrefix key prefix for parameter file.
     */
    RungeKuttaSolverParameters( const std::string keyPrefix = "fem.ode.",
                                const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : DuneODE::ImplicitRungeKuttaSolverParameters( keyPrefix, parameter )
    {}

    using DuneODE::ImplicitRungeKuttaSolverParameters::keyPrefix_;
    using DuneODE::ImplicitRungeKuttaSolverParameters::parameter;

    /**
     * \brief returns a factor for the decision whether to use IMEX Runge-Kutta
     * oder Explicit Runge-Kutta scheme.
     *
     * This parameter is needed for the "IMEX+" scheme where for each time step
     * a solver is chosen regarding the following rule
     * - Implicit/Explicit Runge-Kutta, if \f$ c_\text{diff} < c_\text{factor} c_\text{adv} \f$
     * - Explicit Runge-Kutts, otherwise
     *
     * where \f$ c_\text{diff} \f$ is the maximal time step estimate of the diffusion term,
     * \f$ c_\text{adv} \f$ is the maximal time step estimate of the advection term
     * and \f$ c_\text{factor} \f$ is the factor returned by this function.
     */
    virtual double explicitFactor() const
    {
      return parameter().getValue< double >( keyPrefix_ + "explicitfactor" , 1.0 );
    }

    /**
     * \brief Clones the object.
     */
    RungeKuttaSolverParameters* clone() const
    {
      return new RungeKuttaSolverParameters( *this );
    }

    /**
     * \brief Chooses the type of the Runge-Kutta scheme.
     *
     * Possible values are
     * | Scheme                           | Description       | Integer value |
     * | -------------------------------- | ----------------- | ------------- |
     * | Explicit RK                      | "EX"              | 0             |
     * | Implicit RK                      | "IM"              | 1             |
     * | Implicit/Explicit RK             | "EX"              | 2             |
     * | IMEX/EX RK                       | "IMEX+"           | 3             |
     *
     *
     * \returns An integer describing the selected Runge-Kutta scheme.
     */
    virtual int obtainOdeSolverType () const
    {
      // we need this choice of explicit or implicit ode solver
      // defined here, so that it can be used later in two different
      // methods
      static const std::string odeSolver[]  = { "EX", "IM" ,"IMEX", "IMEX+"  };
      std::string key( keyPrefix_ + "odesolver" );
      return parameter().getEnum( key, odeSolver, 0 );
    }

    /**
     * \brief returns the number of Runge-Kutta steps
     *
     * \param[in] defaultRKOrder The default value if not specified in parameter file.
     * \returns The number of Runge-Kutta steps.
     */
    virtual int obtainRungeKuttaSteps( const int defaultRKOrder ) const
    {
      return parameter().getValue< int > ( keyPrefix_+"order", defaultRKOrder );
    }
  };



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
    typedef RungeKuttaSolver< Operator, ExplicitOperator, ImplicitOperator, LinearInverseOperator > ThisType;

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
        return &a == &d;
      }
    };

  public:
    typedef typename OperatorType :: DestinationType DestinationType ;
    typedef DestinationType  DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType :: GridType GridType;
    typedef typename GridType :: Communication            CommunicationType;
    typedef Dune::Fem::TimeProvider< CommunicationType >  TimeProviderType;

    typedef DuneODE :: OdeSolverInterface< DestinationType >        OdeSolverInterfaceType;
    typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

    typedef DiscreteFunctionType DomainFunctionType;
    typedef DiscreteFunctionType RangeFunctionType;

    typedef Dune::Fem::Operator< DestinationType, DestinationType > HelmHoltzOperatorType ;

    typedef std::pair< OdeSolverInterfaceType* ,  HelmHoltzOperatorType* > solverpair_t ;

    typedef RungeKuttaSolverParameters ParameterType;

    /////////////////////////////////////////////////////////////////////////
    //  ODE solvers from dune-fem/dune/fem/solver/rungekutta
    /////////////////////////////////////////////////////////////////////////
    template < class Op, class DF, bool pardgOdeSolver >
    struct OdeSolverSelection
    {
      template < class OdeParameter >
      static solverpair_t
      createExplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps, const OdeParameter& param,
                            const Dune::Fem::ParameterReader &parameter,
                            const std::string& name = ""  )
      {
        typedef DuneODE :: ExplicitRungeKuttaSolver< DiscreteFunctionType >          ExplicitOdeSolverType;
        return solverpair_t( new ExplicitOdeSolverType( op, tp, rkSteps ), nullptr );
      }

      template < class OdeParameter >
      static solverpair_t
      createImplicitSolver( Op& op, Fem::TimeProviderBase& tp, const int rkSteps, const OdeParameter& param,
                            const Dune::Fem::ParameterReader &parameter,
                            const std::string& name = "" )
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

        typedef typename NonlinearInverseOperatorType::ParameterType NonlinParametersType;

        return solverpair_t(new ImplicitOdeSolverType( *helmOp, tp, rkSteps, param, NonlinParametersType( parameter ) ),
                            helmOp );
#endif
      }

      template < class ExplOp, class ImplOp, class OdeParameter >
      static solverpair_t
      createSemiImplicitSolver( ExplOp& explOp, ImplOp& implOp, Fem::TimeProviderBase& tp,
                                const int rkSteps, const OdeParameter& param,
                                const Dune::Fem::ParameterReader &parameter,
                                const std::string& name = "" )
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


        typedef typename NonlinearInverseOperatorType::ParameterType NonlinParametersType;


        return solverpair_t(new SemiImplicitOdeSolverType( explOp, *helmOp, tp, rkSteps, param, NonlinParametersType( parameter ) ),
                            helmOp );
#endif
      }

    };

    static const bool useParDGSolvers = false ;
    typedef OdeSolverSelection< OperatorType, DestinationType, useParDGSolvers >  OdeSolversType ;

    using BaseType :: solve;
  protected:
    OperatorType&           operator_;
    Fem::TimeProviderBase&  timeProvider_;
    std::unique_ptr< TimeProviderType > tpPtr_;

    const std::string      name_;
    ExplicitOperatorType&  explicitOperator_;
    ImplicitOperatorType&  implicitOperator_;
    std::unique_ptr< const RungeKuttaSolverParameters > param_;

    std::unique_ptr< OdeSolverInterfaceType > explicitSolver_;
    std::unique_ptr< OdeSolverInterfaceType > odeSolver_;
    std::unique_ptr< HelmHoltzOperatorType  > helmholtzOperator_;

    const   double explFactor_ ;
    mutable double fixedTimeStep_ ;
    const int verbose_ ;
    const int rkSteps_ ;
    const int odeSolverType_ ;
    int imexCounter_ , exCounter_;
    int minIterationSteps_, maxIterationSteps_ ;
    bool useImex_ ;
    mutable bool initialized_;
    mutable MonitorType monitor_;

  public:
    RungeKuttaSolver( Fem::TimeProviderBase& tp,
                      OperatorType& op,
                      ExplicitOperatorType& advOp,
                      ImplicitOperatorType& diffOp,
                      const std::string name = "",
                      const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
     : operator_( op ),
       timeProvider_( tp ),
       name_( name ),
       explicitOperator_( advOp ),
       implicitOperator_( diffOp ),
       param_( new RungeKuttaSolverParameters( ParameterKey::generate( name_, "fem.ode." ), parameter ) ),
       explicitSolver_(),
       odeSolver_(),
       helmholtzOperator_(),
       explFactor_( param_->explicitFactor() ),
       verbose_( param_->verbose() ),
       rkSteps_( param_->obtainRungeKuttaSteps( operator_.space().order() + 1 ) ),
       odeSolverType_( param_->obtainOdeSolverType() ),
       imexCounter_( 0 ), exCounter_ ( 0 ),
       minIterationSteps_( std::numeric_limits< int > :: max() ),
       maxIterationSteps_( 0 ),
       useImex_( odeSolverType_ > 1 ),
       initialized_( false )
    {
      solverpair_t solver( nullptr, nullptr ) ;
      // create implicit or explicit ode solver
      if( odeSolverType_ == 0 )
      {
        solver = OdeSolversType :: createExplicitSolver( operator_, tp, rkSteps_, *param_, parameter, name_ );
      }
      else if (odeSolverType_ == 1)
      {
        solver = OdeSolversType :: createImplicitSolver( operator_, tp, rkSteps_, *param_, parameter, name_ );
      }
      else if( odeSolverType_ > 1 )
      {
        // make sure that advection and diffusion operator are different
        if( IsSame< ExplicitOperatorType, ImplicitOperatorType >::check( explicitOperator_, implicitOperator_ ) )
        {
          DUNE_THROW(Dune::InvalidStateException,"Advection and Diffusion operator are the same, therefore IMEX cannot work!");
        }

        solver = OdeSolversType :: createSemiImplicitSolver( explicitOperator_, implicitOperator_, tp, rkSteps_, *param_, parameter, name_ );

        // IMEX+
        if( odeSolverType_ == 3 )
        {
          explicitSolver_.reset( OdeSolversType :: createExplicitSolver( operator_, tp, rkSteps_, *param_, parameter, name_ ).first );
        }
      }
      else
      {
        DUNE_THROW(NotImplemented,"Wrong ODE solver selected");
      }

      odeSolver_.reset( solver.first );
      helmholtzOperator_.reset( solver.second );
    }

    RungeKuttaSolver( OperatorType& op,
                      ExplicitOperatorType& advOp,
                      ImplicitOperatorType& diffOp,
                      const int  odeSolverType,
                      const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
     : RungeKuttaSolver( *(new TimeProviderType( op.space().gridPart().comm(), parameter )),
                         op, advOp, diffOp, "", parameter )
    {
      tpPtr_.reset( static_cast< TimeProviderType* > (&timeProvider_) );
      const TimeSteppingParameters param("femdg.stepper.",parameter);
      const double maxTimeStep = param.maxTimeStep();
      // start first time step with prescribed fixed time step
      // if it is not 0 otherwise use the internal estimate
      tpPtr_->provideTimeStepEstimate(maxTimeStep);
      // adjust fixed time step with timeprovider.factor()
      // SAD: how is fixedTimeStep_ initialized here?
      fixedTimeStep_ = 0;
      fixedTimeStep_ /= tpPtr_->factor() ;
      if ( fixedTimeStep_ > 1e-20 )
        tpPtr_->init( fixedTimeStep_ );
      else
        tpPtr_->init();
      tpPtr_->init();

      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout << "cfl = " << double(tpPtr_->factor()) << ", T_0 = " << tpPtr_->time() << std::endl;
      }
    }

    void initialize( const DestinationType& U )
    {
      if ( ! initialized_)
      {
        if( explicitSolver_ )
        {
          std::abort();
          explicitSolver_->initialize( U );
        }
        assert( odeSolver_ );
        odeSolver_->initialize( U );

        if( tpPtr_ )
        {
          if ( fixedTimeStep_ > 1e-20 )
            tpPtr_->init( fixedTimeStep_ );
          else
            tpPtr_->init();
        }
        initialized_ = true;
      }
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

    void operator () (const DestinationType& U, DestinationType& dest ) const
    {
      dest.assign( U );
      const_cast< ThisType& > (*this).solve( dest );
    }

    void setTimeStepSize( const double dt )
    {
      const double factor_1 = tpPtr_ ? 1.0 / tpPtr_->factor() : 1.0;
      fixedTimeStep_  = dt * factor_1 ;
      timeProvider_.provideTimeStepEstimate( fixedTimeStep_ );
    }

    double deltaT() const { return timeProvider_.deltaT(); }
    double time()   const { return timeProvider_.time(); }

    //! solver the ODE
    void step( DestinationType& U ) const
    {
      const_cast< ThisType& > (*this).solve( U, monitor_ );
    }

    //! solver the ODE
    void solve( DestinationType& U )
    {
      solve( U, monitor_ );
    }

    //! solver the ODE
    void solve( DestinationType& U ,
                MonitorType& monitor )
    {
      initialize( U );

      // make sure the current time step is valid
      assert( timeProvider_.timeStepValid() );

      // take CPU time of solution process
      Dune::Timer timer ;

      // switch upwind direction
      //operator_.switchupwind();
      //explicitOperator_.switchupwind();
      //implicitOperator_.switchupwind();

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

      // if TimeProvider was created locally, it needs to be updated here
      if( tpPtr_ )
      {
        // next time step is prescribed by fixedTimeStep
        if ( fixedTimeStep_ > 1e-20 )
          tpPtr_->next( fixedTimeStep_ );
        else
          tpPtr_->next();

        // apply Limiter to make solution physical (only when tpPtr_ was set)
        // first check if explicitSolver is set (IMEX case)
        if( explicitSolver_ && explicitOperator_.hasLimiter() )
          explicitOperator_.applyLimiter( U );
        else if ( operator_.hasLimiter() )
          operator_.applyLimiter( U );
      }
    }

    //! return CPU time needed for the operator evaluation
    double operatorTime() const
    {
      if( useImex_ )
      {
        return 0.0;
          //explicitOperator_.computeTime() +
          //     implicitOperator_.computeTime() ;
      }
      else
        return 0.0;//operator_.computeTime();
    }

    //! return number of elements meat during operator evaluation
    size_t numberOfElements() const
    {
      if( useImex_ )
        return 0; //explicitOperator_.numberOfElements();
      else
        return 0; //operator_.numberOfElements();
    }

    void description(std::ostream&) const {}

    // gather information from the space operator, the time integratior
    // and the problem to output before each table in tex file
    std::string description() const
    {
      std::string latexInfo;

      /*
      if ((odeSolverType_==0) || (odeSolverType_==1))
        latexInfo = operator_.description();
      else
        latexInfo = explicitOperator_.description()
                    + implicitOperator_.description();
                    */

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
      //operator_.computeTime() ;
      //explicitOperator_.computeTime();
      //implicitOperator_.computeTime();
    }
  }; // end RungeKuttaSolver

  template <class DestinationType>
  using SimpleRungeKuttaSolver = RungeKuttaSolver<
          Dune::Fem::SpaceOperatorInterface< DestinationType >,
          Dune::Fem::SpaceOperatorInterface< DestinationType >,
          Dune::Fem::SpaceOperatorInterface< DestinationType >,
          Dune::Fem::KrylovInverseOperator< DestinationType > >;

} // end namespace
} // end namespace Dune
#endif
