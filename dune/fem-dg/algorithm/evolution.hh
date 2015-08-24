#ifndef DUNE_FEMDG_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_ALGORITHM_EVOLUTION_HH


#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/gridwidth.hh>

// helm holz operator wrapper
//#include <dune/fem/operator/dg/helmholtz.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/operator/dg/passtraits.hh>
#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>
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


#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/algorithm/handler/datawriter.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/solutionlimiter.hh>
#include <dune/fem-dg/algorithm/handler/adapt.hh>

namespace Dune
{

namespace Fem
{

  // EvolutionAlgorithmTraits
  // -------------------------

  template< class Grid,
            class ProblemTraits,
            int polOrder >
  struct EvolutionAlgorithmTraits
  {
    static const int polynomialOrder = polOrder;

    typedef ProblemTraits                                         ProblemTraitsType;

    // type of Grid
    typedef Grid                                          GridType;
    typedef typename ProblemTraits :: HostGridPartType    HostGridPartType;
    typedef typename ProblemTraits :: GridPartType        GridPartType;

    typedef typename ProblemTraits::AnalyticalTraits               AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polynomialOrder >  DiscreteTraits;

    // obtain the problem dependent types, analytical context
    typedef typename AnalyticalTraits::ModelType                   ModelType;
    typedef typename AnalyticalTraits::ProblemType                 ProblemType;
    typedef typename AnalyticalTraits::InitialDataType             InitialDataType;

    // type of discrete function space and discrete function
    typedef typename DiscreteTraits::InitialProjectorType          InitialProjectorType;

    // type of dg operator
    typedef typename DiscreteTraits::OperatorType                  OperatorType;

    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;

    typedef typename DiscreteTraits::ExtraParameterTuple           ExtraParameterTupleType;
    typedef typename DiscreteTraits::IOTupleType                   IOTupleType;

    // wrap operator
    typedef GridTimeProvider< GridType >                           TimeProviderType;

    typedef typename DiscreteTraits::OdeSolverType                 OdeSolverType;

    typedef typename DiscreteTraits::BasicLinearSolverType         BasicLinearSolverType;

    typedef typename DiscreteTraits::AdaptIndicatorType            AdaptIndicatorType;
    typedef typename DiscreteTraits::SolverMonitorHandlerType      SolverMonitorHandlerType;
    typedef typename DiscreteTraits::DiagnosticsHandlerType        DiagnosticsHandlerType;
  };



  template< class EvolutionAlgorithmTraits >
  class EvolutionAlgorithmBase;

  template< class Grid, class ProblemTraits, int polOrder >
  class EvolutionAlgorithm
    : public EvolutionAlgorithmBase< EvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder > >
  {
    typedef EvolutionAlgorithmTraits< Grid, ProblemTraits, polOrder >                    Traits;
    typedef EvolutionAlgorithmBase< Traits >                                             BaseType;
  public:
    EvolutionAlgorithm( Grid &grid, const std::string name = "" )
    : BaseType( grid, name  )
    {}
  };


  // EvolutionAlgorithm
  // ------------------

  template< class EvolutionAlgorithmTraits >
  class EvolutionAlgorithmBase
  {
    typedef EvolutionAlgorithmTraits                                                    Traits;

  public:
    typedef typename Traits::GridType                            GridType;
    typedef typename Traits::IOTupleType                         IOTupleType;
    typedef typename Traits::TimeProviderType                    TimeProviderType;

    typedef typename Traits::HostGridPartType                    HostGridPartType;

    // An analytical version of our model
    typedef typename Traits::ModelType                           ModelType;
    typedef typename Traits::ProblemType                         ProblemType;

    typedef typename Traits::GridPartType                        GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType           DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType                DiscreteFunctionType;

    typedef typename Traits::BasicLinearSolverType               BasicLinearSolverType;

    // The DG space operator
    typedef typename Traits::OperatorType                        OperatorType;

    // The ODE Solvers
    typedef typename Traits::OdeSolverType                       OdeSolverType;

    // type of initial interpolation
    typedef typename Traits::InitialProjectorType                InitialProjectorType;

    // analytical Tratis
    typedef typename Traits::AnalyticalTraits                    AnalyticalTraits;

    // discrete Traits
    typedef typename Traits::DiscreteTraits                      DiscreteTraits;

    typedef typename Traits::DiscreteTraits::GridExactSolutionType
                                                                 GridExactSolutionType;
    typedef typename Traits::ExtraParameterTupleType  ExtraParameterTupleType;

    typedef uint64_t                                              UInt64Type ;

    typedef DiscreteFunctionType                                  CheckPointDiscreteFunctionType;
    typedef DiscreteFunctionType                                  LimitDiscreteFunctionType;
    typedef DiscreteFunctionType                                  AdaptationDiscreteFunctionType;

    typedef typename Traits::AdaptIndicatorType                   AdaptIndicatorType;
    typedef typename Traits::DiagnosticsHandlerType               DiagnosticsHandlerType;
    typedef typename Traits::SolverMonitorHandlerType            SolverMonitorHandlerType;


    EvolutionAlgorithmBase ( GridType &grid, const std::string name = "" )
    : grid_( grid ),
      algorithmName_( name ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      solution_( "U_"+name, space() ),
      problem_( Traits::ProblemTraitsType::problem() ),
      model_( *problem_ ),
      exact_( "exact solution", space() ),
      dataTuple_( std::make_tuple( &solution(), &exact_ ) ),
      solverMonitorHandler_( "" ),
      diagnosticsHandler_(),
      overallTimer_(),
      odeSolver_()
    {}

    virtual const std::string name () { return algorithmName_; }

    GridType& grid () const { return grid_; }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    // return size of grid
    virtual UInt64Type gridSize () const { UInt64Type grSize = grid().size(0); return grid().comm().sum( grSize); }

    virtual bool checkDofsValid ( TimeProviderType& tp, const int loop  ) const { return solution_.dofsValid(); }

    // function creating the ode solvers
    virtual OdeSolverType* createOdeSolver ( TimeProviderType& ) = 0;

    // return reference to the discrete function space
    const DiscreteFunctionSpaceType& space () const { return space_; }

    // return reference to discrete function holding solution
    DiscreteFunctionType& solution () { return solution_; }

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType& monitor() { return solverMonitorHandler_; }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType& diagnostics() { return diagnosticsHandler_; }

    //LIMITING
    virtual void limit( const LimitDiscreteFunctionType& limitSolution ){}
    virtual LimitDiscreteFunctionType* limitSolution () { return &solution_; }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() { return (AdaptIndicatorType*)0; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return &solution_; }

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () { return &solution_; }

    //DATAWRITING
    virtual IOTupleType* dataTuple () { return &dataTuple_; }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const {  return problem().dataPrefix(); }

    virtual void initializeStep ( TimeProviderType &tp, int loop )
    {
      // project initial data
      const bool doCommunicate = ! NonBlockingCommParameter :: nonBlockingCommunication ();
      InitialProjectorType projection( 2 * solution().space().order(), doCommunicate );
      projection( problem().fixedTimeFunction( tp.time() ), solution() );

      // setup ode solver
      odeSolver_.reset( this->createOdeSolver( tp ) );
      assert( odeSolver_ );

      // initialize ode solver
      odeSolver_->initialize( solution() );
    }


    virtual void step ( TimeProviderType &tp )
    {
      overallTimer_.reset();

      // solve ODE
      assert(odeSolver_);
      typename OdeSolverType::MonitorType odeSolverMonitor;
      odeSolver_->solve( solution(), odeSolverMonitor );

      // diagnostic information
      solverMonitorHandler_.step( odeSolverMonitor );
      diagnosticsHandler_.addTimings( odeSolverMonitor.operatorTime_, odeSolverMonitor.odeSolveTime_ );
      //TODO: fix the following line
      //if( adaptIndicator() )
      //  diagnosticsHandler_.addTimings( adaptIndicator()->timings() );
      //diagnosticsHandler_.addTimings( timer_.timings() );
    }

    void finalizeStep ( TimeProviderType &tp )
    {
      // setup exact solution
      // TODO: Add this projection to addtional output writer, later
      Dune::Fem::DGL2ProjectionImpl::project( problem().exactSolution( tp.time() ), exact_ );

      // add eoc errors
      AnalyticalTraits::addEOCErrors( tp, solution(), model(), problem() );

      // delete ode solver
      odeSolver_.reset();
    }

    //! finalize problem, i.e. calculated EOC ...
    virtual void finalize ( const int eocloop ) {}

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


  protected:
    OdeSolverType &odeSolver ()
    {
      assert( odeSolver_ );
      return *odeSolver_;
    }

    GridType&                      grid_;
    std::string                    algorithmName_;
    GridPartType                   gridPart_;
    DiscreteFunctionSpaceType      space_;

    // the solution
    DiscreteFunctionType           solution_;

    // InitialDataType evaluates to $u_0$
    std::unique_ptr< ProblemType > problem_;
    ModelType                      model_;
    DiscreteFunctionType           exact_;
    IOTupleType                    dataTuple_;

    DiagnosticsHandlerType         diagnosticsHandler_;
    SolverMonitorHandlerType       solverMonitorHandler_;

    Dune::Timer                    overallTimer_;
    std::unique_ptr< OdeSolverType > odeSolver_;
  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
