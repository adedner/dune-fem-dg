#ifndef DUNE_FEMDG_ALGORITHM_INTERFACE_HH
#define DUNE_FEMDG_ALGORITHM_INTERFACE_HH

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/subevolution.hh>

#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>


namespace Dune
{

namespace Fem
{

  template< class Grid, class ProblemTraits, int polOrder >
  struct SubAlgorithmInterfaceTraits
  {

  private:
    CHECK_TYPEDEF_EXISTS( AdaptIndicatorType )
    CHECK_TYPEDEF_EXISTS( AdditionalOutputHandlerType )
    CHECK_TYPEDEF_EXISTS( SolverMonitorHandlerType )
    CHECK_TYPEDEF_EXISTS( DiagnosticsHandlerType )

    typedef ProblemTraits                                         Traits;

    typedef typename Traits::AnalyticalTraits                     AnalyticalTraits;
    typedef typename Traits::template DiscreteTraits< polOrder >  DiscreteTraits;

  public:

    typedef typename Traits::GridType                             GridType;
    typedef typename Traits::GridPartType                         GridPartType;
    typedef typename Traits::HostGridPartType                     HostGridPartType;

    typedef GridTimeProvider< GridType >                          TimeProviderType;

    typedef typename AnalyticalTraits::ModelType                  ModelType;
    typedef typename AnalyticalTraits::ProblemType                ProblemType;

    typedef typename DiscreteTraits::DiscreteFunctionType         DiscreteFunctionType;

    typedef uint64_t                                              UInt64Type ;

    typedef DiscreteFunctionType                                  CheckPointDiscreteFunctionType;
    typedef DiscreteFunctionType                                  LimitDiscreteFunctionType;
    typedef DiscreteFunctionType                                  AdaptationDiscreteFunctionType;

    typedef CovariantTuple< typename DiscreteTraits::IOTupleType >          IOTupleType;
    typedef typename AdaptIndicatorTypes< DiscreteTraits >::type            AdaptIndicatorType;
    typedef typename AdditionalOutputHandlerTypes< DiscreteTraits >::type   AdditionalOutputHandlerType;
    typedef typename SolverMonitorHandlerTypes< DiscreteTraits >::type      SolverMonitorHandlerType;
    typedef typename DiagnosticsHandlerTypes< DiscreteTraits >::type        DiagnosticsHandlerType;
  };


  template< class Grid, class ProblemTraits, int polOrder >
  class SubAlgorithmInterface
  {

  public:
    typedef SubAlgorithmInterfaceTraits< Grid, ProblemTraits, polOrder > Traits;

    typedef typename Traits::GridType                             GridType;
    typedef typename Traits::GridPartType                         GridPartType;
    typedef typename Traits::HostGridPartType                     HostGridPartType;

    typedef typename Traits::TimeProviderType                     TimeProviderType;


    typedef typename Traits::ModelType                            ModelType;
    typedef typename Traits::ProblemType                          ProblemType;

    typedef typename Traits::DiscreteFunctionType                 DiscreteFunctionType;

    typedef typename Traits::UInt64Type                           UInt64Type;

    typedef typename Traits::CheckPointDiscreteFunctionType       CheckPointDiscreteFunctionType;
    typedef typename Traits::LimitDiscreteFunctionType            LimitDiscreteFunctionType;
    typedef typename Traits::AdaptationDiscreteFunctionType       AdaptationDiscreteFunctionType;

    typedef typename Traits::IOTupleType                          IOTupleType;
    typedef typename Traits::AdaptIndicatorType                   AdaptIndicatorType;
    typedef typename Traits::DiagnosticsHandlerType               DiagnosticsHandlerType;
    typedef typename Traits::SolverMonitorHandlerType             SolverMonitorHandlerType;
    typedef typename Traits::AdditionalOutputHandlerType          AdditionalOutputHandlerType;


    SubAlgorithmInterface ( GridType &grid, const std::string name = "" )
    : grid_( grid ),
      algorithmName_( name ),
      problem_( ProblemTraits::problem() ),
      model_( problem() )
    {}

    virtual const std::string name () { return algorithmName_; }

    GridType& grid () const { return grid_; }

    // return grid width of grid
    virtual double gridWidth () const
    {
      GridWidthProvider< GridType > gwp( &grid_ );
      return gwp.gridWidth();
    }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type grSize = grid().size(0);
      return grid().comm().sum( grSize );
    }


    // return reference to discrete function holding solution
    virtual DiscreteFunctionType& solution () = 0;

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType* monitor() { return nullptr; }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType* diagnostics() { return nullptr; }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType* additionalOutput() { return nullptr; }

    //LIMITING
    virtual void limit(){}
    virtual LimitDiscreteFunctionType* limitSolution () { return nullptr; }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() { return nullptr; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return nullptr; }

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () { return nullptr; }

    //DATAWRITING
    virtual IOTupleType& dataTuple () = 0;

    virtual bool checkSolutionValid ( const int loop, TimeProviderType* tp ) const { return checkSolutionValid( loop, tp ); }
    virtual bool checkSolutionValid ( const int loop ) const { return true; }


    virtual void initialize ( int loop, TimeProviderType* tp ){ initialize( loop ); }
    virtual void initialize ( int loop ){}

    virtual void preSolve( int loop, TimeProviderType* tp ){ preSolve( loop ); }
    virtual void preSolve( int loop ){}

    virtual void solve ( int loop, TimeProviderType* tp ){ solve( loop ); }
    virtual void solve ( int loop ){}

    virtual void postSolve( int loop, TimeProviderType* tp ){ postSolve( loop ); }
    virtual void postSolve( int loop ){}

    virtual void finalize ( int loop, TimeProviderType* tp ){ finalize( loop ); }
    virtual void finalize ( int loop ){}


    virtual ProblemType& problem ()
    {
      assert( problem_ );
      return *problem_;
    }

    virtual const ProblemType& problem () const
    {
      assert( problem_ );
      return *problem_;
    }

    virtual const ModelType& model () const { return model_; }
    virtual ModelType& model () { return model_; }

  protected:
    GridType&    grid_;
    std::string  algorithmName_;
    ProblemType* problem_;
    ModelType    model_;

  };



}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_INTERFACE_HH
