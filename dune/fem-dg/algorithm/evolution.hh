#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH


#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

// helm holz operator wrapper
//#include <dune/fem/operator/dg/helmholtz.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/operator/dg/passtraits.hh>
#include <dune/fem/misc/femtimer.hh>

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem-dg/pass/threadpass.hh>
#include <dune/common/timer.hh>
#include "monitor.hh"
#include <dune/fem-dg/misc/diagnostics.hh>

namespace Dune
{

namespace Fem
{

  // EvolutionAlgorithmTraits
  // -------------------------

  template< class Grid,
            class ProblemTraits,
            int polOrder,
            class ExtraParameterTuple = std::tuple< > >
  struct EvolutionAlgorithmTraits
  {
    static const int polynomialOrder = polOrder;

    // type of Grid
    typedef Grid                                          GridType;

    typedef DGAdaptiveLeafGridPart< GridType >            HostGridPartType;

    typedef HostGridPartType                              GridPartType;

    typedef typename ProblemTraits::template AnalyticalTraits< HostGridPartType >  AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< HostGridPartType, polynomialOrder >  DiscreteTraits;

    // obtain the problem dependent types, analytical context
    typedef typename AnalyticalTraits::ModelType                   ModelType;
    typedef typename AnalyticalTraits::ProblemType                 ProblemType;
    typedef typename AnalyticalTraits::InitialDataType             InitialDataType;

    // type of discrete function space and discrete function

    typedef typename DiscreteTraits::InitialProjectorType          InitialProjectorType;

    // type of dg operator
    typedef typename DiscreteTraits::FullOperatorType              FullOperatorType;
    typedef typename DiscreteTraits::ImplicitOperatorType          ImplicitOperatorType;
    typedef typename DiscreteTraits::ExplicitOperatorType          ExplicitOperatorType;

    typedef typename ProblemTraits::FunctionSpaceType              FunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;
    typedef typename DiscreteTraits::IndicatorType                 IndicatorType;

    typedef typename DiscreteTraits::OperatorTraitsType            OperatorTraits;

    // wrap operator
    //typedef typename DgHelmHoltzOperatorType::JacobianOperatorType JacobianOperatorType;
    //typedef typename ProblemTraits::template Solver< DgHelmHoltzOperatorType >::Type OdeSolverType;

    typedef typename DiscreteTraits::OdeSolverType                 OdeSolverType;

    typedef typename DiscreteTraits::BasicLinearSolverType         BasicLinearSolverType;
    typedef typename DiscreteTraits::RestrictionProlongationType   RestrictionProlongationType;

    // type of IOTuple
    typedef Dune::tuple< DiscreteFunctionType * >                  IOTupleType;

    typedef SolverMonitor<1>                                       SolverMonitorType;

    // type of Data io
    typedef DataWriter< GridType, IOTupleType >                    DataWriterType;

    // error handling
    typedef typename AnalyticalTraits::EOCErrorIDs                 EOCErrorIDs;


    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                 AdaptationHandlerType;
    typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType >  AdaptationManagerType;
    typedef typename OdeSolverType::MonitorType                                    OdeSolverMonitorType;
    typedef Diagnostics                                                            DiagnosticsType;
  };


  // EvolutionAlgorithm
  // ------------------

  template< class Grid, class ProblemTraits, int polynomialOrder, class ExtraParameterTuple = std::tuple< >>
  class EvolutionAlgorithm
    : public EvolutionAlgorithmBase< EvolutionAlgorithmTraits< Grid, ProblemTraits, polynomialOrder, ExtraParameterTuple > >
  {
    typedef EvolutionAlgorithmTraits< Grid, ProblemTraits, polynomialOrder, ExtraParameterTuple > Traits;
    typedef EvolutionAlgorithmBase< Traits >                                                      BaseType;

  public:
    typedef typename BaseType::GridType GridType;
    typedef typename BaseType::IOTupleType IOTupleType;
    typedef typename BaseType::SolverMonitorType SolverMonitorType;
    typedef typename BaseType::TimeProviderType TimeProviderType;

    typedef typename Traits::HostGridPartType HostGridPartType;

    // An analytical version of our model
    typedef typename Traits::ModelType ModelType;
    typedef typename Traits::ProblemType ProblemType;

    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;

    typedef typename Traits::BasicLinearSolverType   BasicLinearSolverType;

    // The DG space operator
    typedef typename Traits::FullOperatorType     FullOperatorType;
    typedef typename Traits::ExplicitOperatorType ExplicitOperatorType;
    typedef typename Traits::ImplicitOperatorType ImplicitOperatorType;

    typedef typename Traits::OperatorTraits       OperatorTraits;
    // The ODE Solvers
    typedef typename Traits::OdeSolverType OdeSolverType;

    // type of initial interpolation
    typedef typename Traits::InitialProjectorType InitialProjectorType;

    // type of DataWriter
    typedef typename Traits::DataWriterType DataWriterType;

    // analytical Tratis
    typedef typename Traits::AnalyticalTraits AnalyticalTraits;

    // discrete Traits
    typedef typename Traits::DiscreteTraits DiscreteTraits;

    // error handling
    typedef typename Traits::EOCErrorIDs EOCErrorIDs;



    typedef typename Traits::AdaptationHandlerType        AdaptationHandlerType;
    typedef typename Traits::DiagnosticsType              DiagnosticsType;
    typedef typename Traits::OdeSolverMonitorType         OdeSolverMonitorType;
    typedef typename Traits::RestrictionProlongationType  RestrictionProlongationType;
    typedef typename Traits::AdaptationManagerType        AdaptationManagerType;

    using BaseType::grid_;
    using BaseType::adaptParam_ ;
    using BaseType::eocParam_;
    using BaseType::param_;
    using BaseType::limitSolution;

    EvolutionAlgorithm ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      solution_( "U_"+name, space() ),
      problem_( ProblemTraits::problem() ),
      model_( *problem_ ),
      adaptationHandler_( 0 ),
      diagnostics_( true ),
      overallTimer_(),
      eocIds_( AnalyticalTraits::initEoc() ),
      odeSolver_( 0 ),
      rp_( solution_ ),
      adaptationManager_( 0 ),
      dataWriter_( 0 )
    {
      // set refine weight
      rp_.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
    }


    //! destructor
    ~EvolutionAlgorithm ()
    {
      delete odeSolver_;
      odeSolver_ = 0;
      delete dataWriter_;
      dataWriter_ = 0 ;
      delete problem_ ;
      problem_ = 0;
      delete adaptationHandler_ ;
      adaptationHandler_ = 0;
    }

    // function creating the ode solvers
    virtual OdeSolverType* createOdeSolver( TimeProviderType& ) = 0;

    // return reference to the discrete function space
    const DiscreteFunctionSpaceType &space () const
    {
      return space_;
    }

    // return reference to discrete function holding solution
    DiscreteFunctionType &solution ()
    {
      return solution_;
    }

    IOTupleType dataTuple ()
    {
      return IOTupleType( &solution_ );
    }


    void writeData ( TimeProviderType &tp, const bool writeAnyway = false )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        dataWriter_->write( tp );
      }
    }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const
    {
      return problem().dataPrefix();
    }

    // before first step, do data initialization
    void initializeStep ( TimeProviderType &tp, int loop, SolverMonitorType &monitor )
    {
      DiscreteFunctionType& U = solution_;

      if( odeSolver_ == 0 ) odeSolver_ = this->createOdeSolver( tp );
      assert( odeSolver_ );

      if( dataWriter_ == 0 )
      {
        // copy data tuple
        dataTuple_ = dataTuple();
        dataWriter_ = new DataWriterType( grid_, dataTuple_, tp,
          eocParam_.dataOutputParameters( loop, problem_->dataPrefix() ) );
      }
      assert( dataWriter_ );

      typedef typename ProblemType::TimeDependentFunctionType
        TimeDependentFunctionType;

      // communication is needed when blocking communication is used
      // but has to be avoided otherwise (because of implicit solver)
      const bool doCommunicate = ! NonBlockingCommParameter :: nonBlockingCommunication ();

      // create projection
      InitialProjectorType projection( 2 * U.space().order(), doCommunicate );

      // L2 project initial data
      projection( problem().fixedTimeFunction( tp.time() ), U );

      // ode.initialize applies the DG Operator once to get an initial
      // estimate on the time step. This does not change the initial data u.
      odeSolver_->initialize( U );
    }


    void step ( TimeProviderType &tp,
                SolverMonitorType &monitor )
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
      *monitor.newton_iterations     = odeSolverMonitor_.newtonIterations_;
      *monitor.ils_iterations        = odeSolverMonitor_.linearSolverIterations_;
      *monitor.max_newton_iterations = odeSolverMonitor_.maxNewtonIterations_ ;
      *monitor.max_ils_iterations    = odeSolverMonitor_.maxLinearSolverIterations_;
      *monitor.operator_calls        = odeSolverMonitor_.spaceOperatorCalls_;

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
    }

    void finalizeStep ( TimeProviderType &tp )
    {
      DiscreteFunctionType& u = solution();
      // write run file (in writeatonce mode)
      diagnostics_.flush();

      AnalyticalTraits::addEOCErrors( eocIds_, tp, u, model(), problem() );

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

    }

    std::string description () const { return problem().description(); }

    ProblemType &problem ()
    {
      assert( problem_ );
      return *problem_;
    }

    const ProblemType &problem () const
    {
      assert( problem_ );
      return *problem_;
    }

    const ModelType &model () const { return model_; }
    ModelType &model () { return model_; }

    GridPartType &gridPart () { return gridPart_; }

    virtual AdaptationManagerType& adaptationManager()
    {
      if( !adaptationManager_ )
        adaptationManager_ = new AdaptationManagerType( grid_, rp_ );
      return *adaptationManager_;
    }

  protected:
    OdeSolverType &odeSolver ()
    {
      assert( odeSolver_ );
      return *odeSolver_;
    }

    template <class IndicatorOperator, class GradientIndicator>
    void doEstimateMarkAdapt( const IndicatorOperator& dgIndicator,
                              GradientIndicator& gradientIndicator,
                              const bool initialAdaptation = false )
    {
      if( adaptParam_.adaptive() )
      {
        // get grid sequence before adaptation
        const int sequence = space().sequence();

        if( adaptationHandler_ )
        {
          // call operator once to calculate indicator
          dgIndicator.evaluateOnly( solution_ );

          // do marking and adaptation
          adaptationHandler_->adapt( adaptationManager(), initialAdaptation );
        }
        else if( adaptParam_.gradientBasedIndicator() )
        {
          gradientIndicator.estimateAndMark( solution_ );
          adaptationManager().adapt();
        }
        else if( adaptParam_.shockIndicator() )
        {
          // marking has been done by limiter
          adaptationManager().adapt();
        }

        // if grid has changed then limit solution again
        if( sequence != space().sequence() )
        {
          limitSolution();
        }
      }
    }

    GridPartType              gridPart_;
    DiscreteFunctionSpaceType space_;

    // the solution
    DiscreteFunctionType     solution_;

    // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
    // method that gives you the exact solution.
    ProblemType*         problem_;
    ModelType            model_;
    // Initial flux for advection discretization (UpwindFlux)
    AdaptationHandlerType*  adaptationHandler_;

    // diagnostics file
    mutable DiagnosticsType diagnostics_;

    Dune::Timer             overallTimer_;
    double                  odeSolve_;
    EOCErrorIDs             eocIds_;
    OdeSolverType*          odeSolver_;
    OdeSolverMonitorType    odeSolverMonitor_;
    int                     odeSolverType_;

    RestrictionProlongationType rp_;

    AdaptationManagerType*   adaptationManager_;

    IOTupleType             dataTuple_ ;
    DataWriterType*         dataWriter_ ;
  };

} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
