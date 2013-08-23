#ifndef DUNE_FEM_HOWTO_CHECKPOINTING_HH
#define DUNE_FEM_HOWTO_CHECKPOINTING_HH


#if HAVE_SIONLIB && WANT_SIONLIB
#include <dune/fem/io/streams/sionlibstreams.hh>
namespace Dune {

struct PersistenceManagerTraits 
{
  typedef Fem :: SIONlibOutStream  BackupStreamType ;
  typedef Fem :: SIONlibInStream   RestoreStreamType ;
  static const bool singleBackupRestoreFile = true ;
};
#define FEM_PERSISTENCEMANAGERSTREAMTRAITS  PersistenceManagerTraits
}
#endif

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
#include <dune/fem-howto/baseevolution.hh>
#include "problem.hh" 

// approximation order
const int order   = POLORDER;

template <class GridImp, class InitialDataType>               /*@LST1S@*/
struct StepperTraits {
  // type of Grid
  typedef GridImp                                                    GridType;
  // Choose a suitable GridView
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >              GridPartType;
  //typedef AdaptiveLeafGridPart< GridType >                         GridPartType;
  //typedef IdBasedLeafGridPart< GridType >                         GridPartType;

  // type of function space 
  typedef Dune::Fem::FunctionSpace< typename  GridType :: ctype, double, 
                               GridType :: dimensionworld, 1 >  FunctionSpaceType;

  // ... as well as the Space type
  typedef Dune::Fem::DiscontinuousGalerkinSpace < FunctionSpaceType, GridPartType, order> DiscreteSpaceType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType >           DiscreteFunctionType;

  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >  RestrictionProlongationType;
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



template <class GridImp, class InitialDataType>               /*@LST1S@*/
struct Stepper : public AlgorithmBase< StepperTraits< GridImp, InitialDataType> >
{
  // my traits class 
  typedef StepperTraits< GridImp, InitialDataType> Traits ;

  // my base class 
  typedef AlgorithmBase < Traits > BaseType;

  // type of Grid
  typedef typename Traits :: GridType                  GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType              GridPartType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename Traits :: DiscreteFunctionType      DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename Traits :: DiscreteSpaceType         DiscreteSpaceType;

  typedef typename BaseType :: TimeProviderType   TimeProviderType;

  // type of most simple check pointer 
  typedef Dune::Fem::CheckPointer< GridType >   CheckPointerType;

  using BaseType :: grid_;
  using BaseType :: gridPart_;

  Stepper(GridType& grid, const InitialDataType& problem, const char* checkFile) :
    BaseType ( grid ),
    problem_(problem),
    eocId_( Dune::Fem::FemEoc::addEntry(std::string("$L^2$-error")) ),
    checkPointer_( 0 ),
    checkFile_( checkFile )
  {
    Dune::Fem::persistenceManager << error_;
  }                                                                        /*@LST1E@*/

  ~Stepper() 
  {
    removeObj();
  }

  void removeObj() 
  {
    delete checkPointer_;
    checkPointer_ = 0;
  }

  void finalize(DiscreteFunctionType&) 
  {
    removeObj();
  }

  CheckPointerType& checkPointer( TimeProviderType& tp ) const
  {
    // create check point if not exsistent 
    if( ! checkPointer_ ) 
      checkPointer_ = new CheckPointerType( grid_, tp );

    return *checkPointer_;
  }

  // restore data 
  void restoreFromCheckPoint(TimeProviderType& tp,
                             DiscreteFunctionType& solution) 
  {
    // add solution to persistence manager for check pointing 
    Dune::Fem::persistenceManager << solution ;

    // if check file is non-zero a restart is performed 
    if( checkFile_ ) 
    {
      // restore data 
      checkPointer( tp ).restoreData( grid_, checkFile_ );
  
      // check consistency of check point 
      consistencyCheck( tp, solution );
    }
  }

  // backup data  
  void writeCheckPoint(TimeProviderType& tp, DiscreteFunctionType & solution ) const 
  {
    if( Dune::Fem::Parameter::verbose() )
    {
      std::cout << "Try to write checkpoint: error = " << error_ << std::endl;
    }
    checkPointer( tp ).write( tp );
  }

  // before first step, do data initialization 
  void initializeStep(TimeProviderType& tp, DiscreteFunctionType& u)
  {
    Dune::Fem::L2Projection< InitialDataType, DiscreteFunctionType > l2pro;
    l2pro(problem_, u);
  }

  // solve ODE for one time step 
  void step(TimeProviderType& tp, DiscreteFunctionType& u) 
  {
    // do new projection 
    typedef Dune::TimeDependentFunction< InitialDataType > FunctionType;
    FunctionType function( problem_, tp.time() );
    Dune::Fem::L2Projection< FunctionType, DiscreteFunctionType > l2pro;
    l2pro(function, u); 

    // exchange data to ghost cells
    u.communicate();

    // compute error for backup and restore (including ghost cells)
    error_ = computeError(tp, u );
  }

  double computeError(TimeProviderType& tp, DiscreteFunctionType& u)
  {
    L2ErrorNoComm< DiscreteFunctionType > l2norm;
    // Compute L2 error of discretized solution ...
    typedef Dune::TimeDependentFunction< InitialDataType > FunctionType;
    FunctionType function( problem_, tp.time() );
    return l2norm.norm( function, u );
  }

  // after last step, do EOC calculation 
  void finalizeStep(TimeProviderType& tp, DiscreteFunctionType& u)
  {
    // ... and print the statistics out to the eocOutputPath file
    Dune::Fem::FemEoc::setErrors(eocId_, computeError(tp, u ) );
  }

  // reset solution on ghost cells 
  void resetNonInterior( DiscreteFunctionType& solution ) 
  {
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
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
      //std::cerr << "ERROR: backup/restore not consistent" << std::endl;
      DUNE_THROW(Dune::InvalidStateException, "Error in backup/restore" );
    }
  }
private:
  // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  const InitialDataType  &problem_;
  // Initial flux for advection discretization (UpwindFlux)
  const unsigned int      eocId_;

  // check point writer 
  mutable CheckPointerType* checkPointer_;
  // name of checkpoint file 
  const char* checkFile_;

  double error_;
};

#endif // #ifndef DUNE_FEM_HOWTO_CHECKPOINTING_HH
