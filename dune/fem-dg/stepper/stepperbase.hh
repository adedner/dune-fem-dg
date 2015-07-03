#ifndef FEMHOWTO_DGSTEPPERBASE_HH
#define FEMHOWTO_DGSTEPPERBASE_HH
#include <config.h>

#ifdef LOCALDEBUG
static double sum_ = 0.;
static double sum2_ = 0.;
static double localMaxRatio_ = 0.;
static double localMinRatio_ = 1e+100;
static double maxRatioOfSums = 0.;
static double minRatioOfSums = 1e+100;
#endif


#ifndef NDEBUG
// enable fvector and fmatrix checking
#define DUNE_ISTL_WITH_CHECKING
#endif

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/projection/l2projection.hh>

#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>

// include local header files
#include <dune/fem-dg/stepper/baseevolution.hh>
#include <dune/fem-dg/stepper/steppertraits.hh>

#include <dune/fem-dg/misc/diagnostics.hh>
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>

using namespace Dune;


template <class GridImp,
          class ProblemTraits,
          int polynomialOrder>
struct StepperBase
  : public AlgorithmBase< StepperTraits< GridImp, ProblemTraits, polynomialOrder> >
{
  // my traits class
  typedef StepperTraits< GridImp, ProblemTraits, polynomialOrder> Traits ;
  typedef AlgorithmBase< Traits > BaseType;


  // type of Grid
  typedef typename Traits :: GridType                 GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType             GridPartType;

  // initial data type
  typedef typename Traits :: InitialDataType          InitialDataType;

  // model type
  typedef typename Traits :: ModelType                ModelType ;

  // type of convective numerical flux
  typedef typename Traits :: FluxType                 FluxType ;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef typename Traits :: FullOperatorType              FullOperatorType;
  typedef typename Traits :: ExplicitOperatorType          ExplicitOperatorType;
  typedef typename Traits :: ImplicitOperatorType          ImplicitOperatorType;

  // type of linear solver for implicit ode
  typedef typename Traits :: LinearInverseOperatorType    LinearInverseOperatorType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename Traits :: DiscreteFunctionType     DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename Traits :: DiscreteSpaceType        DiscreteSpaceType;

  // the indicator function type (for limiting only)
  typedef typename Traits :: IndicatorType            IndicatorType;

  // The ODE Solvers
  typedef typename Traits :: OdeSolverType            OdeSolverType;
  typedef typename OdeSolverType :: MonitorType       OdeSolverMonitorType ;

  typedef typename BaseType :: TimeProviderType       TimeProviderType;

  typedef AdaptationHandler< GridType,
           typename DiscreteSpaceType::FunctionSpaceType >  AdaptationHandlerType;

  // type of run time diagnostics
  typedef Diagnostics                                 DiagnosticsType;

  // type of most simple check pointer
  typedef Dune::Fem::CheckPointer< GridType >         CheckPointerType;

  // type of solver monitor
  typedef typename BaseType :: SolverMonitorType      SolverMonitorType;

  // type of solver monitor
  typedef typename BaseType :: IOTupleType            IOTupleType;

  // type of data writer
  typedef Dune::Fem::DataWriter< GridType, IOTupleType >    DataWriterType;

  // type of parameter class
  typedef Dune::Fem::Parameter                        ParameterType ;

  // the restriction and prolongation projection
  typedef typename Traits :: RestrictionProlongationType RestrictionProlongationType;

  // type of adaptation manager
  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType ;

  using BaseType :: grid_ ;
  using BaseType :: limitSolution ;

  // constructor taking grid
  StepperBase(GridType& grid, const std::string name = "" ) :
    BaseType( grid, name  ),
    gridPart_( grid_ ),
    space_( gridPart_ ),
    solution_( "U", space() ),
    additionalVariables_( ParameterType :: getValue< bool >("femhowto.additionalvariables", false) ?
        new DiscreteFunctionType("A", space() ) : 0 ),
    problem_( ProblemTraits::problem() ),
    adaptationHandler_( 0 ),
    diagnostics_( true ),
    checkPointer_( 0 ),
    overallTimer_(),
    eocId_( Fem::FemEoc::addEntry(std::string("$L^2$-error")) ),
    odeSolver_( 0 ),
    rp_( solution_ ),
    adaptationManager_( grid_, rp_ ),
    adaptationParameters_( ),
    dataWriter_( 0 )
  {
    // set refine weight
    rp_.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );

  }

  //! destructor
  virtual ~StepperBase()
  {
    delete checkPointer_;
    checkPointer_ = 0;
    delete odeSolver_;
    odeSolver_ = 0;
    delete dataWriter_;
    dataWriter_ = 0 ;
    delete problem_ ;
    problem_ = 0;
    delete adaptationHandler_ ;
    adaptationHandler_ = 0;
    delete additionalVariables_;
    additionalVariables_ = 0;
  }

  virtual const DiscreteSpaceType& space() const { return space_ ; }

  double gridWidth () const { return Dune::Fem::GridWidth::calcGridWidth( gridPart_ ); }

  // return reference to discrete function holding solution
  DiscreteFunctionType& solution() { return solution_; }

  IOTupleType dataTuple()
  {
    // tuple with additionalVariables
    return IOTupleType( &solution_, additionalVariables_, indicator() );
  }

  bool checkDofsValid ( TimeProviderType& tp, const int loop ) const
  {
    return solution_.dofsValid();
  }

  void writeData( TimeProviderType& tp, const bool writeAnyway = false )
  {
    if( dataWriter_ )
    {
      const bool reallyWrite = writeAnyway ? true : dataWriter_->willWrite( tp );
      if( reallyWrite && additionalVariables_ )
      {
        setupAdditionalVariables( tp, solution(), model(), *additionalVariables_ );
      }

      //todo: check it
      dataWriter_->write( tp );
    }

    // possibly write a grid solution
    problem().postProcessTimeStep( grid_, solution(), tp.time() );
  }

  bool adaptive () const { return adaptationManager_.adaptive(); }

  // function creating the ode solvers
  virtual OdeSolverType* createOdeSolver( TimeProviderType& ) = 0;

  CheckPointerType& checkPointer( TimeProviderType& tp ) const
  {
    // create check point if not exsistent
    if( ! checkPointer_ )
      checkPointer_ = new CheckPointerType( grid_, tp );

    return *checkPointer_;
  }

  // restore data if checkpoint file is given
  virtual bool restoreFromCheckPoint( TimeProviderType& tp )
  {
    // add solution to persistence manager for check pointing
    Dune::Fem::persistenceManager << solution_ ;

    std::string checkPointRestartFile = checkPointRestartFileName();

    // if check file is non-zero a restart is performed
    if( checkPointRestartFile.size() > 0 )
    {
      // restore data
      checkPointer( tp ).restoreData( grid_, checkPointRestartFile );
      return false;
    }

    return true ;
  }

  void writeDiagnostics( TimeProviderType& tp ) const
  {
    const double ldt = tp.deltaT();
    const int maxNumDofs = space().blockMapper().maxNumDofs() * space().localBlockSize;

    // write times to run file
    diagnostics_.write( tp.time() + ldt, ldt,
                    odeSolverMonitor_.numberOfElements_,    // number of elements
                    maxNumDofs,                             // number of dofs per element (max)
                    odeSolverMonitor_.operatorTime_,        // time for operator evaluation
                    odeSolverMonitor_.odeSolveTime_,        // ode solver
                    adaptationManager_.adaptationTime(),    // time for adaptation
                    adaptationManager_.loadBalanceTime(),   // time for load balance
                    overallTimer_.elapsed());               // time step overall time
  }

  // write checkpoint data and also run time diagnostics
  virtual void writeCheckPoint( TimeProviderType& tp ) const
  {
    // write data checkpoint (see datawriter.hh)
    checkPointer( tp ).write( tp );
  }

  //! returns data prefix for EOC loops ( default is loop )
  virtual std::string dataPrefix() const
  {
    return problem_->dataPrefix();
  }

  // gather information from the space operator, the time integratior
  // and the problem to output before each table in tex file
  std::string description() const
  {
    std::string latexInfo = odeSolver_->description();
    latexInfo +=  problem_->description() + "\n\n";
    return latexInfo;
  }

  // before first step, do data initialization
  void initializeStep( TimeProviderType& tp, const int loop )
  {
    DiscreteFunctionType& U = solution_;

    if( odeSolver_ == 0 ) odeSolver_ = this->createOdeSolver( tp );
    assert( odeSolver_ );

    if( dataWriter_ == 0 )
    {
      // copy data tuple
      dataTuple_ = dataTuple () ;
      dataWriter_ = new DataWriterType( grid_, dataTuple_, tp,
        EocDataOutputParameters( loop, problem_->dataPrefix() ) );
    }
    assert( dataWriter_ );

    typedef typename InitialDataType :: TimeDependentFunctionType
      TimeDependentFunctionType;

    // communication is needed when blocking communication is used
    // but has to be avoided otherwise (because of implicit solver)
    const bool doCommunicate = ! NonBlockingCommParameter :: nonBlockingCommunication ();

    // create L2 projection
    Fem :: L2Projection< TimeDependentFunctionType,
        DiscreteFunctionType > l2pro( 2 * U.space().order(), doCommunicate );

    // L2 project initial data
    l2pro( problem().fixedTimeFunction( tp.time() ), U );

    // ode.initialize applies the DG Operator once to get an initial
    // estimate on the time step. This does not change the initial data u.
    odeSolver_->initialize( U );
  }


  void step(TimeProviderType& tp,
            SolverMonitorType& monitor )
  {
    DiscreteFunctionType& U = solution();

    // reset overall timer
    overallTimer_.reset();

    // solve ODE
    assert(odeSolver_);
    odeSolver_->solve( U, odeSolverMonitor_ );

    // limit solution if necessary
    limitSolution ();

    // copy information to solver monitor
    monitor.newton_iterations     = odeSolverMonitor_.newtonIterations_;
    monitor.ils_iterations        = odeSolverMonitor_.linearSolverIterations_;
    monitor.max_newton_iterations = odeSolverMonitor_.maxNewtonIterations_ ;
    monitor.max_ils_iterations    = odeSolverMonitor_.maxLinearSolverIterations_;
    monitor.operator_calls        = odeSolverMonitor_.spaceOperatorCalls_;

    // set time step size to monitor
    monitor.setTimeStepInfo( tp );

#ifdef LOCALDEBUG
    maxRatioOfSums = std::max( maxRatioOfSums, std::abs(sum_/sum2_) );
    minRatioOfSums = std::min( minRatioOfSums, std::abs(sum_/sum2_) );

    std::cout <<"localMaxRatioPerTimeStep: " <<localMaxRatio_ <<std::endl;
    std::cout <<"localMinRatioPerTimeStep: " <<localMinRatio_ <<std::endl;
    std::cout <<"maxRatioOfSumsPerTimeStep: " <<maxRatioOfSums <<std::endl;
    std::cout <<"minRatioOfSumsPerTimeStep: " <<minRatioOfSums <<std::endl;

    sum_ = 0.;
    sum2_ = 0.;
#endif

    writeDiagnostics( tp );
  }

  inline double error(TimeProviderType& tp, DiscreteFunctionType& u)
  {
    Fem :: L2Norm< GridPartType > l2norm( u.space().gridPart(), 2*u.space().order()+4 );
    return l2norm.distance( problem().fixedTimeFunction( tp.time() ), u );
  }

  void finalizeStep(TimeProviderType& tp)
  {
    DiscreteFunctionType& u = solution();
    // write run file (in writeatonce mode)
    diagnostics_.flush();

    bool doFemEoc = problem().calculateEOC( tp, u, eocId_ );

    // ... and print the statistics out to a file
    if( doFemEoc )
      Fem::FemEoc::setErrors(eocId_, error(tp, u ));

#ifdef LOCALDEBUG
    std::cout <<"maxRatioOfSums: " <<maxRatioOfSums <<std::endl;
    std::cout <<"minRatioOfSums: " <<minRatioOfSums <<std::endl;
#endif

    // delete ode solver
    delete odeSolver_;
    odeSolver_ = 0;

    delete dataWriter_;
    dataWriter_ = 0 ;

    delete adaptationHandler_;
    adaptationHandler_ = 0;
  }                                                       /*@LST1E@*/

  const InitialDataType& problem() const
  {
    assert( problem_ );
    return *problem_;
  }

  InitialDataType& problem()
  {
    assert( problem_ );
    return *problem_;
  }

  // return reference to additional variables
  virtual IndicatorType* indicator()  { return 0; }

  virtual const ModelType& model() const = 0 ;

protected:
  template <class IndicatorOperator, class GradientIndicator>
  void doEstimateMarkAdapt( const IndicatorOperator& dgIndicator,
                            GradientIndicator& gradientIndicator,
                            const bool initialAdaptation = false )
  {
    if( adaptationManager_.adaptive() )
    {
      // get grid sequence before adaptation
      const int sequence = space().sequence();

      if( adaptationHandler_ )
      {
        // call operator once to calculate indicator
        dgIndicator.evaluateOnly( solution_ );

        // do marking and adaptation
        adaptationHandler_->adapt( adaptationManager_, initialAdaptation );
      }
      else if( adaptationParameters_.gradientBasedIndicator() )
      {
        gradientIndicator.estimateAndMark( solution_ );
        adaptationManager_.adapt();
      }
      else if( adaptationParameters_.shockIndicator() )
      {
        // marking has been done by limiter
        adaptationManager_.adapt();
      }

      // if grid has changed then limit solution again
      if( sequence != space().sequence() )
      {
        limitSolution();
      }
    }
  }

  ///////////////////////////////////////////
  //  instances
  ///////////////////////////////////////////

  GridPartType         gridPart_;
  DiscreteSpaceType    space_;

  // the solution
  DiscreteFunctionType   solution_;
  DiscreteFunctionType*  additionalVariables_;

  // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  InitialDataType*  problem_;
  // Initial flux for advection discretization (UpwindFlux)
  AdaptationHandlerType*  adaptationHandler_;

  // diagnostics file
  mutable DiagnosticsType diagnostics_;
  // check point writer
  mutable CheckPointerType* checkPointer_;

  Timer                   overallTimer_;
  double                  odeSolve_;
  const unsigned int      eocId_;
  OdeSolverType*          odeSolver_;
  OdeSolverMonitorType    odeSolverMonitor_;
  int                     odeSolverType_;

  RestrictionProlongationType rp_;

  AdaptationManagerType   adaptationManager_;
  AdaptationParameters    adaptationParameters_;

  IOTupleType             dataTuple_ ;
  DataWriterType*         dataWriter_ ;
};
#endif // FEMHOWTO_STEPPER_HH
