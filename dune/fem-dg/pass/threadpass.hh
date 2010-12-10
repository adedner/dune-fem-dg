#ifndef DUNE_THREADPASS_HH
#define DUNE_THREADPASS_HH

//- system includes 
#include <dune/fem/misc/utility.hh>

#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/selection.hh>
// #include "discretemodel.hh"
#include <dune/fem/pass/dgmodelcaller.hh>

#include <dune/fem/solver/timeprovider.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

#include <dune/fem/space/common/allgeomtypes.hh> 
#include <dune/fem/space/common/arrays.hh> 
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/selection.hh>

#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/misc/threaditerator.hh>

namespace Dune {

  template < class InnerPass > 
  class ThreadPass :
    public LocalPass< typename InnerPass :: DiscreteModelType,
                      typename InnerPass :: PreviousPassType, 
                      InnerPass :: passId > 
  {
    typedef ThreadPass< InnerPass > ThisType;
  public:
    typedef InnerPass InnerPassType;
    typedef typename InnerPass :: DiscreteModelType  DiscreteModelType;
    typedef typename InnerPass :: PreviousPassType   PreviousPassType;

    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelType, PreviousPassType, InnerPass :: passId>  BaseType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef typename EntityType :: EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Iterator over the space
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType:: BaseFunctionSetType
      BaseFunctionSetType; 

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef CombinedSelector< ThisType , SelectorType > CombinedSelectorType;
    typedef DGDiscreteModelCaller< DiscreteModelType 
                                   , ArgumentType 
                                   , CombinedSelectorType
                                 > DiscreteModelCallerType;

    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };

    // type of local id set 
    typedef typename GridPartType::IndexSetType IndexSetType; 

    typedef Fem::ThreadIterator< DiscreteFunctionSpaceType > ThreadIteratorType;

  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order 
    ThreadPass(DiscreteModelType& problem, 
               PreviousPassType& pass, 
               const DiscreteFunctionSpaceType& spc,
               const int volumeQuadOrd = -1,
               const int faceQuadOrd = -1) :
      BaseType(pass, spc),
      iterators_( spc ),
      problems_( Fem::ThreadManager::maxThreads() ),
      passes_( Fem::ThreadManager::maxThreads() ),
      firstCall_( true )
    {
      for(int i=0; i<Fem::ThreadManager::maxThreads(); ++i)
      {
        // use serparate discrete problem for each thread 
        problems_[ i ] = new DiscreteModelType( problem );
        // create dg passes, the last bool disables communication in the pass itself
        passes_[ i ]   = new InnerPassType( *problems_[ i ], pass, spc, volumeQuadOrd, faceQuadOrd, false );
      }

#ifndef NDEBUG
      if( Parameter :: verbose() )
        std::cout << "Thread Pass initialized\n";
#endif
    }

    virtual ~ThreadPass () 
    {
      for(int i=0; i<Fem::ThreadManager::maxThreads(); ++i)
      {
        delete passes_[ i ];
        delete problems_[ i ];
      }
    }
   
    //! print tex info
    void printTexInfo(std::ostream& out) const 
    {
      BaseType::printTexInfo(out);
      pass( 0 ).printTexInfo(out);
    }

    //! Estimate for the timestep size 
    double timeStepEstimateImpl() const 
    {
      double dtMin = pass( 0 ).timeStepEstimateImpl();
      const int maxThreads = Fem::ThreadManager::maxThreads();
      for( int i = 1; i < maxThreads ; ++i)
      {
        dtMin = std::min( dtMin, pass( i ).timeStepEstimateImpl() );
      }
      return dtMin;
    }

  protected:  
    //! returns true for flux evaluation if neighbor 
    //! is on same thread as entity  
    struct NBChecker 
    {
      const ThreadIteratorType& storage_;
      const int myThread_; 
      NBChecker( const ThreadIteratorType& st ) 
        : storage_( st ),
          myThread_( Fem::ThreadManager::thread() )
      {}

      // returns true if niehhbor can be updated 
      bool operator () ( const EntityType& en, const EntityType& nb ) const 
      {
        return myThread_ == storage_.thread( nb );
      }
    };

    InnerPass& pass( const int thread ) const 
    {
      assert( (int) passes_.size() > thread );
      return *( passes_[ thread ] );
    }

  public:  
    //! switch upwind direction
    void switchUpwind() 
    {
      const int maxThreads = Fem::ThreadManager::maxThreads();
      for(int i=0; i<maxThreads; ++i ) 
        problems_[ i ]->switchUpwind(); 
    }

    //! overload compute method to use thread iterators 
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      // get stopwatch 
      Timer timer;

      // clear destination 
      dest.clear();

      // for the first call we only run on one thread to avoid 
      // clashes with the singleton storages for quadratures 
      // and base function caches etc.
      // after one grid traversal everything should be set up 
      if( firstCall_ ) 
      {
        //! get pass for my thread  
        InnerPassType& pass = *(passes_[ 0 ]);
        //std::cout << "Thread Pass :: First Call !!! " << this << std::endl;

        // pepare 
        pass.prepare( arg, dest );

        // Iterator is of same type as the space iterator 
        typedef typename DiscreteFunctionSpaceType :: IteratorType Iterator;
        const Iterator endit = iterators_.space().end();
        for (Iterator it = iterators_.space().begin(); it != endit; ++it)
        {
          pass.applyLocal( *it );
        }

        // finalize pass 
        pass.finalize(arg, dest);

        // set tot false since first call has been done
        firstCall_ = false ;
      }
      else 
      {
        // update thread iterators in case grid changed 
        iterators_.update();

        // call prepare before parallel area 
        const int maxThreads = Fem::ThreadManager::maxThreads();
        for(int i=0; i<maxThreads; ++i ) 
        {
          // prepare pass (make sure pass doesn't clear dest, this will conflict)
          pass( i ).prepare( arg, dest );
        }

        /////////////////////////////////////////////////
        // BEGIN PARALLEL REGION 
        /////////////////////////////////////////////////
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          //! get pass for my thread  
          InnerPassType& myPass = pass( Fem::ThreadManager::thread() );

          // create NB checker 
          NBChecker nbChecker( iterators_ );

          // set current time to my pass
          myPass.setTime( this->time() );

          // Iterator is of same type as the space iterator 
          typedef typename ThreadIteratorType :: IteratorType Iterator;
          const Iterator endit = iterators_.end();
          for (Iterator it = iterators_.begin(); it != endit; ++it)
          {
            assert( iterators_.thread( *it ) == Fem::ThreadManager::thread() );
            myPass.applyLocal( *it, nbChecker );
          }

          // finalize pass (make sure communication is done in case of thread parallel
          // program, this would give conflicts)
          myPass.finalize(arg, dest);
        } 
        /////////////////////////////////////////////////
        // END PARALLEL REGION 
        /////////////////////////////////////////////////

      } // end if first call 

      // accumulate time 
      this->computeTime_ += timer.elapsed()/((double) Fem::ThreadManager::maxThreads());

      // communicate calculated function 
      dest.communicate();
    }

    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // we use prepare of internal operator
      abort();
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      // we use finalize of internal operator
      abort();
    }

    void applyLocal( const EntityType& en) const 
    {
      // we use applyLocal of internal operator
      abort();
    }
  private:
    ThreadPass();
    ThreadPass(const ThreadPass&);
    ThreadPass& operator=(const ThreadPass&);

  protected:  
    mutable ThreadIteratorType iterators_;
    std::vector< DiscreteModelType* > problems_; 
    std::vector< InnerPassType* > passes_;
    mutable bool firstCall_;
  };
//! @}  
} // end namespace Dune

#endif
