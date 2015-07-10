#ifndef DUNE_FEMDG_CHECKPOINTING_STEPPER_HH
#define DUNE_FEMDG_CHECKPOINTING_STEPPER_HH

#include <dune/fem-dg/misc/streams.hh>

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/solver/odesolver.hh>

// space and function
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

// include local header files
#include <dune/fem-dg/stepper/baseevolution.hh>
#include "problem.hh"

template <class GridImp, class ProblemTraits, int order>               /*@LST1S@*/
struct StepperTraits {
  // type of Grid
  typedef GridImp                                                    GridType;
  // Choose a suitable GridView
  //typedef Dune::Fem::LeafGridPart< GridType >              GridPartType;
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >              GridPartType;
  //typedef AdaptiveLeafGridPart< GridType >                         GridPartType;
  //typedef IdBasedLeafGridPart< GridType >                         GridPartType;

  // type of initial data
  typedef typename ProblemTraits :: template Traits< GridPartType > :: InitialDataType  InitialDataType;

  // type of function space
  typedef typename InitialDataType :: FunctionSpaceType  FunctionSpaceType;

  // ... as well as the Space type
  typedef Dune::Fem::DiscontinuousGalerkinSpace < FunctionSpaceType, GridPartType, order> DiscreteSpaceType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >           DiscreteFunctionType;

  // type of restriction/prolongation projection for adaptive simulations
  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >  RestrictionProlongationType;

  // fake indicator type
  typedef DiscreteFunctionType  IndicatorType ;

  // type of IOTuple
  typedef Dune::tuple< DiscreteFunctionType*, DiscreteFunctionType* >  IOTupleType;
};


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


template <class GridImp, class ProblemTraits, int order>
struct CheckPointingStepper : public AlgorithmBase< GridImp>
{
  // my traits class
  typedef StepperTraits< GridImp, ProblemTraits, order> Traits ;

  // my base class
  typedef AlgorithmBase < GridImp > BaseType;

  // type of Grid
  typedef typename Traits :: GridType                  GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType              GridPartType;

  // type of problem data
  typedef typename Traits :: InitialDataType           InitialDataType ;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename Traits :: DiscreteFunctionType      DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename Traits :: DiscreteSpaceType         DiscreteSpaceType;

  typedef typename BaseType :: TimeProviderType   TimeProviderType;

  typedef typename Traits :: IOTupleType     IOTupleType;

  typedef typename BaseType :: SolverMonitorType  SolverMonitorType;

  // type of most simple check pointer
  typedef Dune::Fem::CheckPointer< GridType >   CheckPointerType;

  using BaseType :: grid_;

  CheckPointingStepper( GridType& grid, const std::string name = "" )
  : BaseType ( grid, name ),
    gridPart_( grid ),
    space_( gridPart_ ),
    solution_( "solution-"+name, space() ),
    problem_( ProblemTraits::problem() ),
    eocId_( Dune::Fem::FemEoc::addEntry(std::string("$L^2$-error")) ),
    checkPointer_(),
    checkFile_( 0 )
  {
    Dune::Fem::persistenceManager << error_;
  }

  const DiscreteSpaceType& space() const { return space_ ; }

  const InitialDataType& problem () const { assert( problem_ ); return *problem_; }

  // return reference to discrete function holding solution
  DiscreteFunctionType& solution() { return solution_; }

  IOTupleType dataTuple()
  {
    // tuple with additionalVariables
    return IOTupleType( &solution_, (DiscreteFunctionType*) 0 );
  }

  CheckPointerType& checkPointer( TimeProviderType& tp ) const
  {
    // create check point if not exsistent
    if( ! checkPointer_ )
      checkPointer_.reset( new CheckPointerType( grid_, tp ) );

    return *checkPointer_;
  }

  // restore data, return true for new start
  bool restoreFromCheckPoint(TimeProviderType& tp )
  {
    // add solution to persistence manager for check pointing
    bool writeData = Dune::Fem::Parameter::getValue<bool>("fem.io.writedata", true );
    if( writeData )
    {
      Dune::Fem::persistenceManager << solution_ ;
    }

    std::string checkPointRestartFile = checkPointRestartFileName();

    // if check file is non-zero a restart is performed
    if( checkPointRestartFile.size() > 0 )
    {
      // restore data
      checkPointer( tp ).restoreData( grid_, checkPointRestartFile );

      // check consistency of check point
      consistencyCheck( tp, solution_ );

      // return false for no new start
      return false;
    }
    // do new start
    return true ;
  }

  // backup data
  void writeCheckPoint(TimeProviderType& tp) const
  {
    if( Dune::Fem::Parameter::verbose() )
    {
      std::cout << "Try to write checkpoint: error = " << error_ << std::endl;
    }
    checkPointer( tp ).write( tp );
  }

  // before first step, do data initialization
  void initializeStep(TimeProviderType& tp, const int loop )
  {
    Dune::Fem::L2Projection< InitialDataType, DiscreteFunctionType > l2pro;
    l2pro( problem(), solution_);
  }

  // solve ODE for one time step
  void step(TimeProviderType& tp, SolverMonitorType& monitor )
  {
    // do new projection
    typedef Dune::Fem::InstationaryFunction< InitialDataType, Dune::Fem::__InstationaryFunction::HoldReference > FunctionType;
    FunctionType function( problem(), tp.time() );
    Dune::Fem::L2Projection< FunctionType, DiscreteFunctionType > l2pro;
    l2pro(function, solution_);

    // exchange data to ghost cells
    solution_.communicate();

    // compute error for backup and restore (including ghost cells)
    error_ = computeError(tp, solution_ );
  }

  double computeError(TimeProviderType& tp, DiscreteFunctionType& u)
  {
    L2ErrorNoComm< DiscreteFunctionType > l2norm;
    // Compute L2 error of discretized solution ...
    typedef Dune::Fem::InstationaryFunction< InitialDataType, Dune::Fem::__InstationaryFunction::HoldReference > FunctionType;
    FunctionType function( problem(), tp.time() );
    return l2norm.norm( function, u );
  }

  // after last step, do EOC calculation
  void finalizeStep(TimeProviderType& tp)
  {
    // ... and print the statistics out to the eocOutputPath file
    Dune::Fem::FemEoc::setErrors(eocId_, computeError(tp, solution_ ) );
  }


  // reset solution on ghost cells
  void resetNonInterior( DiscreteFunctionType& solution )
  {
    typedef typename GridPartType :: template Codim< 0 > :: template
        Partition< Dune::All_Partition > :: IteratorType  IteratorType;

    typedef typename IteratorType :: Entity  EntityType ;

    IteratorType it    = solution.space().gridPart().template begin< 0, Dune::All_Partition > ();
    IteratorType endit = solution.space().gridPart().template end  < 0, Dune::All_Partition > ();

    for( ; it != endit; ++ it )
    {
      const EntityType& entity = * it ;
      if( entity.partitionType() != Dune::InteriorEntity )
      {
        solution.localFunction( entity ).clear();
      }
    }
  }

  // after last step, do EOC calculation
  void consistencyCheck(TimeProviderType& tp, DiscreteFunctionType& u)
  {
    // reset ghost cells to make sure we rely on the communication
    resetNonInterior( u );

    // communicate data first to check communication
    u.communicate();

    // Compute L2 error of discretized solution ...
    double error = computeError( tp, u );

    std::cout << "Stepper::consistencyCheck: L2-error after restore: " << error
              << "  stored value: " << error_ << std::endl;
    if( std::abs( error - error_ ) > 1e-14 )
    {
      std::cerr << "ERROR: backup/restore not consistent" << std::endl;
      //DUNE_THROW(Dune::InvalidStateException, "Error in backup/restore" );
    }
  }

protected:
  GridPartType          gridPart_;
  DiscreteSpaceType     space_;
  DiscreteFunctionType  solution_;

  // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  std::unique_ptr< const InitialDataType > problem_;
  // Initial flux for advection discretization (UpwindFlux)
  const unsigned int      eocId_;

  // check point writer
  mutable std::unique_ptr< CheckPointerType > checkPointer_;
  // name of checkpoint file
  const char* checkFile_;

  double error_;
};

#endif // #ifndef DUNE_FEMDG_CHECKPOINTING_STEPPER_HH
