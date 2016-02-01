#ifndef DUNE_FEMDG_CHECKPOINTING_STEPPER_HH
#define DUNE_FEMDG_CHECKPOINTING_STEPPER_HH


// local includes
#include <dune/fem-dg/algorithm/sub/evolution.hh>
#include <dune/fem/function/common/instationary.hh>
#include <dune/fem/function/common/rangegenerators.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem-dg/algorithm/handler/checkpoint.hh>
#include <dune/fem-dg/test/dataio/checkedcheckpointhandler.hh>

namespace Dune
{
namespace Fem
{

  // calculates || u-u_h ||_L2 including the ghost cells
  template <class DiscreteFunctionType>
  class L2ErrorNoComm
  {
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

  public:
    template<class FunctionType >
    double norm (FunctionType &f, DiscreteFunctionType &discFunc, int polOrd = -1 )
    {
      const DiscreteFunctionSpaceType &space = discFunc.space();

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

      typedef typename GridPartType :: template Codim< 0 > ::
        template Partition< Dune::All_Partition > :: IteratorType IteratorType ;

      typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

      if( polOrd < 0 ) polOrd = 2*space.order() + 4 ;

      RangeType ret (0.0);
      RangeType phi (0.0);

      double sum = 0.0;

      IteratorType it    = space.gridPart().template begin< 0, Dune::All_Partition > ();
      IteratorType endit = space.gridPart().template end< 0, Dune::All_Partition > ();

      for(; it != endit ; ++it)
      {
        Dune::Fem::CachingQuadrature<GridPartType,0> quad(*it, polOrd);
        LocalFuncType lf = discFunc.localFunction(*it);
        for( size_t qP = 0; qP < quad.nop(); ++qP )
        {
          double det = (*it).geometry().integrationElement(quad.point(qP));
          f.evaluate((*it).geometry().global(quad.point(qP)), ret);
          lf.evaluate(quad[qP],phi);
          RangeType diff = ret - phi ;
          sum += det * quad.weight(qP) * ( diff * diff );
        }
      }
      return std::sqrt( Dune::Fem::MPIManager::comm().sum( sum ) );
    }
  };





 /**
   *  \brief Algorithm for solving an instationary PDE.
   *
   *  \ingroup SubAlgorithms
   */
  template <class GridImp,
            class ProblemTraits,
            int polynomialOrder >
  class SubCheckPointingAlgorithm
    : public SubAlgorithmInterface< GridImp, ProblemTraits, polynomialOrder >
  {

    typedef SubAlgorithmInterface< GridImp, ProblemTraits, polynomialOrder > BaseType ;

  public:

    // type of Grid
    typedef typename BaseType::GridType                       GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                   GridPartType;

    // initial data type
    typedef typename BaseType::ProblemType                    ProblemType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                      ModelType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType           DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename BaseType::IOTupleType                    IOTupleType;

    typedef typename BaseType::TimeProviderType               TimeProviderType;

    // type of 64bit unsigned integer
    typedef typename BaseType::UInt64Type                     UInt64Type;

    typedef typename BaseType::AdaptIndicatorType             AdaptIndicatorType;
    typedef typename BaseType::DiagnosticsType                DiagnosticsType;
    typedef typename BaseType::SolverMonitorType              SolverMonitorType;
    typedef typename BaseType::AdditionalOutputType           AdditionalOutputType;
    typedef typename BaseType::LimitDiscreteFunctionType      LimitDiscreteFunctionType;
    typedef typename BaseType::CheckPointDiscreteFunctionType CheckPointDiscreteFunctionType;
    typedef typename BaseType::AdaptationDiscreteFunctionType AdaptationDiscreteFunctionType;

    typedef SubEvolutionContainer< DiscreteFunctionType >     ContainerType;

    typedef typename ProblemTraits::AnalyticalTraits          AnalyticalTraits;

    using BaseType::grid;
    using BaseType::name;
    using BaseType::problem;
    using BaseType::model;
    using BaseType::solution;
    using BaseType::gridSize;

    SubCheckPointingAlgorithm ( GridType &grid, const ContainerType& container )
    : BaseType( grid ),
      container_( container ),
      solution_( container_.solution() ),
      ioTuple_( std::make_tuple( &solution(), nullptr ) ),
      error_( 0.0 )
    {}

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () { return &solution(); }

    //DATAWRITING
    virtual IOTupleType& dataTuple () { return ioTuple_; }

    virtual DiscreteFunctionType& solution ()
    {
      assert( solution_ );
      return *solution_;
    }

    void checkCheckPointSolutionValid( TimeProviderType& tp )
    {
      // reset ghost cells to make sure we rely on the communication
      resetNonInterior( solution() );

      // communicate data first to check communication
      solution().communicate();

      // Compute L2 error of discretized solution ...
      double error = computeError( tp, solution() );

      std::cout << "Algorithm::consistencyCheck: L2-error after restore: " << error
                << "  stored value: " << error_ << std::endl;
      if( std::abs( error - error_ ) > 1e-14 )
      {
        std::cerr << "ERROR: backup/restore not consistent" << std::endl;
        //DUNE_THROW(Dune::InvalidStateException, "Error in backup/restore" );
      }
    }

  private:
    virtual void doInitialize ( const int loop, TimeProviderType& tp ) override
    {
      auto ftp = problem().fixedTimeFunction( tp.time() );
      GridFunctionAdapter< typename ProblemType::InstationaryFunctionType, GridPartType >
        adapter( "-exact", ftp, solution().gridPart(), solution().space().order()+2 );
      interpolate( adapter, solution());
    }

    virtual void doPreSolve ( const int loop, TimeProviderType& tp ) override
    {
      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout << "Try to write checkpoint: error = " << error_ << std::endl;
      }
    }

    virtual void doSolve ( const int loop, TimeProviderType& tp ) override
    {
      auto ftp = problem().fixedTimeFunction( tp.time() );
      GridFunctionAdapter< typename ProblemType::InstationaryFunctionType, GridPartType >
        adapter( "-exact", ftp, solution().gridPart(), solution().space().order()+2 );
      interpolate( adapter, solution() );

      // exchange data to ghost cells
      solution().communicate();

      // compute error for backup and restore (including ghost cells)
      error_ = computeError(tp, solution() );
    }

    virtual void doFinalize ( const int loop, TimeProviderType& tp ) override
    {
      // add eoc errors
      AnalyticalTraits::addEOCErrors( tp, solution(), model(), problem() );

    }

    double computeError(TimeProviderType& tp, DiscreteFunctionType& u)
    {
      L2ErrorNoComm< DiscreteFunctionType > l2norm;
      // Compute L2 error of discretized solution ...
      typedef Dune::Fem::InstationaryFunction< ProblemType, Dune::Fem::__InstationaryFunction::HoldReference > FunctionType;
      FunctionType function( problem(), tp.time() );
      return l2norm.norm( function, u );
    }

    // reset solution on ghost cells
    void resetNonInterior( DiscreteFunctionType& solution )
    {
      for( auto&& entity : entities(solution) )
        if( entity.partitionType() != Dune::InteriorEntity )
          solution.localFunction( entity ).clear();
    }

  protected:
    const ContainerType&                    container_;
    std::shared_ptr< DiscreteFunctionType > solution_;

    IOTupleType                             ioTuple_;
    double                                  error_;

  };


}
}

#endif // #ifndef DUNE_FEMDG_CHECKPOINTING_STEPPER_HH
