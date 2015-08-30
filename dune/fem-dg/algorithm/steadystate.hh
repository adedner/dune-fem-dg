#ifndef DUNE_FEMDG_ALGORITHM_STEADYSTATEALGORITHM_HH
#define DUNE_FEMDG_ALGORITHM_STEADYSTATEALGORITHM_HH

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem-dg/algorithm/base.hh>

#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>


namespace Dune
{

namespace Fem
{

  template<class Grid,
            class ProblemTraits,
            int polOrd >
  struct SubSteadyStateTraits
  {
    enum { polynomialOrder = polOrd };

    // type of Grid
    typedef Grid                                                   GridType;
    typedef typename ProblemTraits::HostGridPartType               HostGridPartType;
    typedef typename ProblemTraits::GridPartType                   GridPartType;

    typedef typename ProblemTraits::AnalyticalTraits               AnalyticalTraits;
    typedef typename ProblemTraits::template DiscreteTraits< polynomialOrder >  DiscreteTraits;

    // obtain the problem dependent types, analytical context
    typedef typename AnalyticalTraits::ModelType                   ModelType;
    typedef typename AnalyticalTraits::ProblemType                 ProblemType;

    // type of discrete function space and discrete function
    typedef typename ProblemTraits::FunctionSpaceType              FunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteTraits::DiscreteFunctionType          DiscreteFunctionType;

    typedef typename DiscreteTraits::ExtraParameterTuple           ExtraParameterTuple;

    //typedef typename DiscreteTraits::OdeSolverType                 OdeSolverType;
    typedef typename DiscreteTraits::BasicLinearSolverType         BasicLinearSolverType;

    // tpye of jacobian operator used in the nested newton loops
    typedef typename DiscreteTraits::OperatorType                  OperatorType;

    typedef typename DiscreteTraits::AssemblerType                 AssemblerType;

    // type of IOTuple
    typedef typename DiscreteTraits::IOTupleType                   IOTupleType;

    //typedef typename DiscreteTraits::AdaptIndicatorType            AdaptIndicatorType;
    typedef typename DiscreteTraits::SolverMonitorHandlerType      SolverMonitorHandlerType;
    typedef typename DiscreteTraits::DiagnosticsHandlerType        DiagnosticsHandlerType;
  };


  template< class Grid, class ProblemTraits, int polOrder >
  class SubSteadyStateAlgorithm
  {
    typedef SubSteadyStateTraits< Grid, ProblemTraits, polOrder >    Traits;

  public:
    typedef typename Traits::GridType                           GridType;
    typedef typename Traits::IOTupleType                        IOTupleType;

    typedef typename Traits::HostGridPartType                     HostGridPartType;

    // initial data type
    typedef typename Traits::ProblemType                          ProblemType;

    // An analytical version of our model
    typedef typename Traits::ModelType                            ModelType;

    typedef typename Traits::GridPartType                         GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType            DiscreteFunctionSpaceType;
    typedef typename Traits::DiscreteFunctionType                 DiscreteFunctionType;

    // The DG space operator
    typedef typename Traits::OperatorType                         OperatorType;

    // type of steady state solver
    typedef typename Traits::BasicLinearSolverType                BasicLinearSolverType;

    // type of analytical traits
    typedef typename Traits::AnalyticalTraits                     AnalyticalTraits;

    // type of discrete traits
    typedef typename Traits::DiscreteTraits                       DiscreteTraits;

    typedef typename Traits::AssemblerType                        AssemblerType;

    typedef uint64_t                                              UInt64Type ;

    //typedef typename Traits::AdaptIndicatorType                   AdaptIndicatorType;
    typedef typename Traits::SolverMonitorHandlerType             SolverMonitorHandlerType;
    typedef typename Traits::DiagnosticsHandlerType               DiagnosticsHandlerType;

    SubSteadyStateAlgorithm ( GridType &grid, const std::string name = "", std::string solName = ""  )
      : grid_( grid ),
        algorithmName_( name ),
        problem_( ProblemTraits::problem() ),
        model_( problem() ),
        solverMonitorHandler_( "" ),
        gridPart_( grid ),
        space_( gridPart_ ),
        solution_( "solution-" + solName + name, space_ ),
        exactSolution_( "exact-solution-" + solName + name, space_ ),
        solverIterations_()
    {
      solution().clear();
    }

    virtual const std::string name () { return algorithmName_; }

    GridType& grid () const { return grid_; }

    //! return reference to discrete space
    DiscreteFunctionSpaceType& space () { return space_; }
    const DiscreteFunctionSpaceType& space () const { return space_; }

    virtual DiscreteFunctionType& solution () { return solution_; }

    virtual DiscreteFunctionType& rhs() = 0;

    std::string description () const { return problem().description(); }

    //! returns data prefix for EOC loops ( default is loop )
    virtual std::string dataPrefix () const
    {
      return problem().dataPrefix();
    }

    virtual SolverMonitorHandlerType& monitor()
    {
      return solverMonitorHandler_;
    }

    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const { return GridWidth::calcGridWidth( gridPart_ ); }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type grSize = grid().size( 0 );
      return grid().comm().sum( grSize );
    }

    virtual BasicLinearSolverType* createSolver( DiscreteFunctionType* rhs ) = 0;

    virtual void initialize ( const int loop )
    {
      //initialize solverMonitor
      solverMonitorHandler_.registerData( "GridWidth", solverMonitorHandler_.monitor().gridWidth, nullptr, true );
      solverMonitorHandler_.registerData( "Elements", solverMonitorHandler_.monitor().elements, nullptr, true );
      solverMonitorHandler_.registerData( "TimeSteps", solverMonitorHandler_.monitor().timeSteps, nullptr, true );
      solverMonitorHandler_.registerData( "ILS", solverMonitorHandler_.monitor().ils_iterations, &solverIterations_ );
      solverMonitorHandler_.registerData( "MaxILS", solverMonitorHandler_.monitor().max_ils_iterations );
    }

    virtual void preSolve( const int loop )
    {
      solution().clear();
      solver_.reset( this->createSolver( &rhs() ) );
    }

    virtual void solve ( const int loop )
    {
      (*solver_)( rhs(), solution() );
      solverIterations_ = solver_->iterations();
    }

    virtual void postSolve( const int loop )
    {
      solverMonitorHandler_.finalize( gridWidth(), gridSize() );
    }

    void finalize ( const int loop )
    {
      // add eoc errors
      AnalyticalTraits::addEOCErrors( solution(), model(), problem() );

      solver_.reset( nullptr );
    }

    //DATAWRITING
    IOTupleType dataTuple () { return std::make_tuple( &solution(), &exactSolution_ ); }

    const ProblemType &problem () const
    {
      assert( problem_ );
      return *problem_;
    }

    ProblemType &problem ()
    {
      assert( problem_ );
      return *problem_;
    }

    ModelType &model () { return model_; }
    const ModelType &model () const { return model_; }

  protected:
    GridType&   grid_;
    std::string algorithmName_;
    ProblemType *problem_;
    ModelType model_;

    SolverMonitorHandlerType       solverMonitorHandler_;

    GridPartType gridPart_;      // reference to grid part, i.e. the leaf grid
    DiscreteFunctionSpaceType space_;    // the discrete function space
    DiscreteFunctionType solution_;
    DiscreteFunctionType exactSolution_;

    std::unique_ptr< BasicLinearSolverType > solver_;
    int solverIterations_;
  };

}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
