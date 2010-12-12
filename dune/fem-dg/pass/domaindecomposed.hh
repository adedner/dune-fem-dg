#ifndef DUNE_FEM_THREADITERATOR2_HH
#define DUNE_FEM_THREADITERATOR2_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/misc/threadmanager.hh>

#include "threadfilter.hh"
#include "communicator.hh"

namespace Dune {

  namespace Fem {

    /** \brief Thread iterator */
    template <class DiscreteFunctionSpace>  
    class DomainDecomposedIterator
    {
      typedef DiscreteFunctionSpace SpaceType;
    public:  
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename SpaceType :: GridType      GridType;
      typedef typename SpaceType :: GridPartType  GridPartType;
      typedef typename SpaceType :: IndexSetType  IndexSetType;
      typedef ThreadFilter<GridPartType> FilterType;
      typedef FilteredGridPart< GridPartType, FilterType > FilteredGridPartType;

#ifdef _OPENMP
      typedef typename FilteredGridPartType :: template Codim< 0 > :: IteratorType
        IteratorType;
#else 
      typedef typename SpaceType :: IteratorType IteratorType ;
#endif

      typedef typename IteratorType :: Entity EntityType ;
      typedef typename SpaceType :: GridType :: template Codim<0> ::
        EntityPointer EntityPointer ;

    protected:  
      const SpaceType& space_ ;
      const IndexSetType& indexSet_;

#ifdef _OPENMP
      int sequence_;
      std::vector< FilteredGridPartType* > filteredGridParts_;
      MutableArray< int > threadNum_;
#endif
    public:  
      //! contructor creating thread iterators 
      explicit DomainDecomposedIterator( const SpaceType& spc )
        : space_( spc ),
          indexSet_( space_.indexSet() )
#ifdef _OPENMP
        , sequence_( -1 )  
        , filteredGridParts_( Fem :: ThreadManager :: maxThreads() )
#endif
      {
#ifdef _OPENMP
        FilterType filter( space_.gridPart(), threadNum_ );
        GridType& grid = const_cast< GridType& > (space_.grid());
        for(int i=0; i< Fem :: ThreadManager :: maxThreads(); ++i )
        {
          filteredGridParts_[ i ] = new FilteredGridPartType( grid, filter );
          // set thread number of grid part 
          filteredGridParts_[ i ]->filter().setThread( i );
        }

        threadNum_.setMemoryFactor( 1.1 ); 
#endif
        update();
      }

      //! return reference to space 
      const SpaceType& space() const { return space_; }

      //! update internal list of iterators 
      void update() 
      {
#ifdef _OPENMP
        // if grid got updated also update iterators 
        if( sequence_ != space_.sequence() )
        {
          if( ! ThreadManager :: singleThreadMode() ) 
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          // get maximal number of threads possible 
          const size_t maxThreads = ThreadManager :: maxThreads() ;

          // create partitioner 
          Partitioner< GridType > db( space_.gridPart() , maxThreads );
          // do partitioning 
          db.serialPartition( false );

          // get end iterator
          typedef typename SpaceType :: IteratorType SpaceIteratorType;
          const SpaceIteratorType endit = space_.end();

          // get size for index set 
          const size_t size = indexSet_.size( 0 );

          // resize threads storage 
          threadNum_.resize( size );
          // set all values to default value 
          for(size_t i = 0; i<size; ++i) threadNum_[ i ] = -1;

          // just for diagnostics 
          std::vector< int > counter( maxThreads , 0 );

          for(SpaceIteratorType it = space_.begin(); it != endit; ++it ) 
          {
            const EntityType& entity  = * it;
            const int rank = db.getRank( entity );
            //std::cout << "Got rank = " << rank << "\n";
            threadNum_[ indexSet_.index( entity ) ] = rank ; 
            ++counter[ rank ];
          }

          // update sequence number 
          sequence_ = space_.sequence();

          if( Parameter :: verbose() )
          {
            for(size_t i = 0; i<maxThreads; ++i ) 
              std::cout << "T[" << i << "] = " << counter[ i ] << std::endl;
          }
            
          //for(size_t i = 0; i<size; ++i ) 
          //  std::cout << threadNum_[ i ] << std::endl;
        }
#endif
      }

      //! return begin iterator for current thread 
      IteratorType begin() const 
      {
#ifdef _OPENMP
        return filteredGridParts_[ ThreadManager :: thread() ]->template begin< 0 > ();
#else 
        return space_.begin();
#endif
      }

      //! return end iterator for current thread 
      IteratorType end() const 
      {
#ifdef _OPENMP
        return filteredGridParts_[ ThreadManager :: thread() ]->template end< 0 > ();
#else 
        return space_.end();
#endif
      }

      //! return thread number this entity belongs to 
      int index(const EntityType& entity ) const 
      {
        return indexSet_.index( entity );
      }

      //! return thread number this entity belongs to 
      int thread(const EntityType& entity ) const 
      {
#ifdef _OPENMP
        assert( (size_t) threadNum_.size() > indexSet_.index( entity ) );
        // NOTE: this number can also be negative for ghost elements or elements
        // that do not belong to the set covered by the space iterators 
        return threadNum_[ indexSet_.index( entity ) ];
#else 
        return 0;
#endif
      }
    };
  }
}

#endif
