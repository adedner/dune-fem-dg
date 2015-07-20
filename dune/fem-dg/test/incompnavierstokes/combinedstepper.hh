#ifndef COMBINED_DGSTEPPERBASE_HH
#define COMBINED_DGSTEPPERBASE_HH
#include <config.h>

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
#include <dune/fem-dg/solver/rungekuttasolver.hh>

// include local header files
#include <dune/fem-dg/operator/adaptation/estimator.hh>
#include <dune/fem-dg/stepper/baseevolution.hh>
#include "combinedsteppertraits.hh"

#include <dune/fem-dg/misc/diagnostics.hh>
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>


using namespace Dune;

template< class Stepper1, class Stepper2 >
class CombinedSolverMonitor
{
};

template< class StokesProblemCreator,
          class AdvectionDiffusionProblemCreator >
struct CombinedStepper
  : public AlgorithmBase< typename CombinedStepperTraits< typename StokesProblemCreator::StepperType, typename AdvectionDiffusionProblemCreator::StepperType >::GridType,
                          CombinedSolverMonitor< typename StokesProblemCreator::StepperType::SolverMonitorType, typename AdvectionDiffusionProblemCreator::StepperType::SolverMonitorType > >
{
  // my traits class
  typedef typename CombinedStepperTraits<
          typename StokesProblemCreator::StepperType,
          typename AdvectionDiffusionProblemCreator::StepperType > Traits;

  typedef typename StokesProblemCreator::StepperType StokesStepper;

  typedef typename StokesStepper :: PressureDiscreteFunctionType PressureDiscreteFunctionType;
  typedef typename StokesStepper :: VelocityDiscreteFunctionType VelocityDiscreteFunctionType;

  typedef std::tuple< VelocityDiscreteFunctionType*, PressureDiscreteFunctionType* > AdvectionDiffusionParameterType;

  // this should be ok but could lead to a henn-egg problem
  typedef AdvectionDiffusionStepper< GridType, AdvectionDiffusionProblemCreator, polynomialOrder, AdvectionDiffusionParameterType > Type;

  typedef StokesStepper              Stepper1;
  typedef AdvectionDiffusionStepper  Stepper2;

  typedef CombinedStepperTraits< Stepper1, Stepper2 > Traits ;
  typedef CombinedSolverMonitor< typename Stepper1::SolverMonitorType, typename Stepper2::SolverMonitorType > CombinedSolverMonitorType;
  typedef AlgorithmBase< typename Traits::GridType, CombinedSolverMonitorType > BaseType;


  // type of Grid
  typedef typename Traits :: GridType                    GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType                GridPartType;

  typedef typename BaseType :: TimeProviderType          TimeProviderType;

  // type of run time diagnostics
  typedef Diagnostics                                    DiagnosticsType;

  // type of most simple check pointer
  typedef Dune::Fem::CheckPointer< GridType >            CheckPointerType;

  // type of solver monitor
  typedef typename BaseType :: SolverMonitorType         SolverMonitorType;

  // type of solver monitor
  typedef typename Traits :: IOTupleType                 IOTupleType;

  // type of solver monitor
  typedef typename Traits :: RestrictionProlongationType RestrictionProlongationType;

  // type of adaptation manager
  typedef typename Traits :: AdaptationManagerType       AdaptationManagerType;

  // type of data writer
  typedef Dune::Fem::DataWriter< GridType, IOTupleType > DataWriterType;

  // type of parameter class
  typedef Dune::Fem::Parameter                           ParameterType ;

  using BaseType::timeStepTimer_;
  using BaseType::fixedTimeStep_;
  using BaseType::fixedTimeStepEocLoopFactor_;
  using BaseType :: grid_ ;
  using BaseType::gridWidth;
  using BaseType::gridSize;
  using BaseType::solve;

  // constructor taking grid
  CombinedStepper(GridType& grid, const std::string name = "" ) :
    BaseType( grid, name ),
    checkPointer_( 0 ),
    dataWriter_( 0 ),
    stepper1_( grid, ProblemCreator1::moduleName() ),
    stepper2_( grid, ProblemCreator2::moduleName() ),
    rp_( std::make_tuple( stepper1_.restrictionProlongation(), stepper2_.restrictionProlongation() ) ) ,
    adaptationManager_( 0 )
  {
    // set refine weight
  }


  //! destructor
  virtual ~CombinedStepper()
  {
    delete checkPointer_;
    checkPointer_ = 0;
    delete dataWriter_;
    dataWriter_ = 0 ;
  }

  IOTupleType dataTuple()
  {
    return std::tuple_cat( stepper1_.dataTuple(), stepper2_.dataTuple() );
  }

  bool checkDofsValid ( TimeProviderType& tp, const int loop ) const
  {
    return stepper1_.checkDofsValid( tp, loop ) && stepper2_.checkDofsValid( tp, loop );
  }

  void writeData( TimeProviderType& tp, const bool writeAnyway = false )
  {
    if( dataWriter_ )
    {
      const bool reallyWrite = writeAnyway ? true : dataWriter_->willWrite( tp );
      dataWriter_->write( tp );
    }
  }

  bool adaptive () const { return stepper1_.adaptive() && stepper2_.adaptive(); }

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
    Dune::Fem::persistenceManager << stepper1_.solution() ;
    Dune::Fem::persistenceManager << stepper2_.solution() ;

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
    //stepper1_.writeDiagnostics( tp );
    //stepper2_.writeDiagnostics( tp );
  }

  // write checkpoint data and also run time diagnostics
  virtual void writeCheckPoint( TimeProviderType& tp ) const
  {
    checkPointer( tp ).write( tp );
  }

  // gather information from the space operator, the time integratior
  // and the problem to output before each table in tex file
  std::string description() const
  {
    std::string latexInfo = stepper1_.description() + stepper2_.description();
    return latexInfo;
  }

  // before first step, do data initialization
  void initializeStep( TimeProviderType& tp, const int loop )
  {
    if( dataWriter_ == 0 )
    {
      // copy data tuple
      dataTuple_ = dataTuple ();
      EocParameters eocParam( ParameterKey::generate( "", "fem.eoc." ) );
      dataWriter_ = new DataWriterType( grid_, dataTuple_, tp,
          eocParam.dataOutputParameters( loop, stepper1_.problem().dataPrefix() ) );
    }
    assert( dataWriter_ );

    //draft:
    //stepper1_.setFunctionTuple( get<2>(tuple), get<4>(tuple) )
    //stepper1_.setFunctionTuple( get<1>(tuple), get<6>(tuple) )

    stepper1_.initializeStep( tp, loop );
    stepper2_.initializeStep( tp, loop );
  }

  void step(TimeProviderType& tp,
            SolverMonitorType& monitor )
  {
    // draft:
    // MoreFancyTimeProviderType tp1( tp )
    // FancyTimeProviderType     tp2( tp )
    // SubMonitorStep1      m1( monitor1 )
    // SubMonitorStep2      m2( monitor2 )
    //
    // while( error(stepper1_.solution - stepper1_.oldsolution()) < eps )
    // {
    //   stepper1_.step( tp1, m1 );
    //   stepper2_.inputTuple( stepper1_.outputTuple() );
    //   stepper2_.step( tp2, m2 );
    // }

    // solve stokes problem
    stokes_.solve();

    // solve advection-diffusion problem
    advectionDiffusion_.step( tp, monitor );

    // solve stokes problem
    stokes_.solve();
  }

  void finalizeStep(TimeProviderType& tp)
  {
    //stepper1_.finalizeStep( tp );
    //stepper2_.finalizeStep( tp );
  }                                                       /*@LST1E@*/

  //! estimate and mark solution
  virtual void initialEstimateMarkAdapt( )
  {
    //stepper1_.initialEstimateMarkAdapt();
    //stepper2_.initialEstimateMarkAdapt();
  }

  //! estimate and mark solution
  virtual void estimateMarkAdapt( )
  {
    //stepper1_.estimateMarkAdapt();
    //stepper2_.estimateMarkAdapt();
  }

  virtual AdaptationManagerType& adaptationManager()
  {
    if( !adaptationManager_ )
      adaptationManager_ = new AdaptationManagerType( grid_, rp_ );
    return *adaptationManager_;
  }


protected:
  template <class IndicatorOperator, class GradientIndicator>
  void doEstimateMarkAdapt( const IndicatorOperator& dgIndicator,
                            GradientIndicator& gradientIndicator,
                            const bool initialAdaptation = false )
  {
    //stepper1_.doEstimateMarkAdapt( dgIndicator, gradientIndicator, initialAdaptation );
    //stepper2_.doEstimateMarkAdapt( dgIndicator, gradientIndicator, initialAdaptation );
  }

  ///////////////////////////////////////////
  //  instances
  ///////////////////////////////////////////
  mutable CheckPointerType*    checkPointer_;

  IOTupleType                  dataTuple_ ;
  DataWriterType*              dataWriter_ ;

  StokesStepper                stokes_;
  AdvectionDiffusionStepper    advectionDiffusion_;
  RestrictionProlongationType  rp_;
  AdaptationManagerType*       adaptationManager_;
};

#endif // FEMHOWTO_STEPPER_HH
