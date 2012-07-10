#ifndef DUNE_FEM_DG_THREADITERATOR2_HH
#define DUNE_FEM_DG_THREADITERATOR2_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/threadmanager.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>

#include <dune/fem-dg/pass/threadfilter.hh>
#include <dune/fem-dg/pass/partitioner.hh>

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
      typedef FilteredGridPart< GridPartType, FilterType, false > FilteredGridPartType ;

#ifdef USE_SMP_PARALLEL
      typedef typename FilteredGridPartType :: template Codim< 0 > :: IteratorType
        IteratorType;
#else 
      typedef typename SpaceType :: IteratorType IteratorType ;
#endif

      typedef typename IteratorType :: Entity EntityType ;
      typedef typename GridPartType :: template Codim<0> ::
        EntityPointerType EntityPointer ;

    protected:  
      const SpaceType& space_ ;
      const IndexSetType& indexSet_;

#ifdef USE_SMP_PARALLEL
      int sequence_;
      std::vector< FilteredGridPartType* > filteredGridParts_;

      typedef typename FilterType :: ThreadArrayType ThreadArrayType;
      ThreadArrayType threadNum_;
#endif
      const bool verbose_ ;

    public:  
      //! contructor creating thread iterators 
      explicit DomainDecomposedIterator( const SpaceType& spc )
        : space_( spc ),
          indexSet_( space_.indexSet() )
#ifdef USE_SMP_PARALLEL
        , sequence_( -1 )  
        , filteredGridParts_( Fem :: ThreadManager :: maxThreads() )
#endif
        , verbose_( Parameter::verbose() && 
                    Parameter::getValue<bool>("femdg.threads.verbose", false ) )
      {
#ifdef USE_SMP_PARALLEL
        for(int thread=0; thread < Fem :: ThreadManager :: maxThreads(); ++thread )
        {
          // i is the thread number of this filter 
          filteredGridParts_[ thread ] 
            = new FilteredGridPartType( space_.gridPart(), 
                                        FilterType( space_.gridPart(), threadNum_, thread ) );
        }

        threadNum_.setMemoryFactor( 1.1 ); 
#endif
        update();
      }

      //! return reference to space 
      const SpaceType& space() const { return space_; }

      //! return filter for given thread 
      const FilterType& filter( const int thread ) const 
      {
        return filteredGridParts_[ thread ]->filter();
      }

      //! update internal list of iterators 
      void update() 
      {
#ifdef USE_SMP_PARALLEL
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
          ThreadPartitioner< GridPartType > db( space_.gridPart() , maxThreads );
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

          {
            // just for diagnostics 
            std::vector< int > counter( maxThreads , 0 );

            int numInteriorElems = 0;
            for(SpaceIteratorType it = space_.begin(); it != endit; ++it, ++numInteriorElems ) 
            {
              const EntityType& entity  = * it;
              const int rank = db.getRank( entity );
              assert( rank >= 0 );
              //std::cout << "Got rank = " << rank << "\n";
              threadNum_[ indexSet_.index( entity ) ] = rank ; 
              ++counter[ rank ];
            }

            // update sequence number 
            sequence_ = space_.sequence();

            if( verbose_ )
            {
              std::cout << "DomainDecomposedIterator: sequence = " << sequence_ << " size = " << numInteriorElems << std::endl;
              for(size_t i = 0; i<maxThreads; ++i ) 
                std::cout << "DomainDecomposedIterator: T[" << i << "] = " << counter[ i ] << std::endl;
            }
              
#ifndef NDEBUG
            // check that all interior elements have got a valid thread number 
            for(SpaceIteratorType it = space_.begin(); it != endit; ++it ) 
            {
              assert( threadNum_[ indexSet_.index( *it ) ] >= 0 );
            }
#endif
            //for(size_t i = 0; i<size; ++i ) 
            //{
            //  //std::cout << threadNum_[ i ] << std::endl;
            //}
          }
        }
#endif
      }

      //! return begin iterator for current thread 
      IteratorType begin() const 
      {
#ifdef USE_SMP_PARALLEL
        return filteredGridParts_[ ThreadManager :: thread() ]->template begin< 0 > ();
#else 
        return space_.begin();
#endif
      }

      //! return end iterator for current thread 
      IteratorType end() const 
      {
#ifdef USE_SMP_PARALLEL
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
#ifdef USE_SMP_PARALLEL
        assert( (int)threadNum_.size() > (int) indexSet_.index( entity ) );
        // NOTE: this number can also be negative for ghost elements or elements
        // that do not belong to the set covered by the space iterators 
        return threadNum_[ indexSet_.index( entity ) ];
#else 
        return 0;
#endif
      }
    };

    /** \brief Thread iterator */
    template <class DiscreteFunctionSpace>  
    class DomainDecomposedIteratorStorage
    {
      typedef DiscreteFunctionSpace SpaceType;
    public:  
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

      typedef typename SpaceType :: GridPartType  GridPartType;
      typedef typename SpaceType :: IndexSetType  IndexSetType;

      typedef DomainDecomposedIterator< DiscreteFunctionSpace >  DomainIterator;
      typedef typename DomainIterator :: FilterType    FilterType ;
      typedef typename DomainIterator :: IteratorType  IteratorType;

      typedef typename IteratorType :: Entity EntityType ;

    private:
      struct IteratorFactory
      {
        struct Key
        {
          const DiscreteFunctionSpaceType& space_;
          const IndexSetType& indexSet_;
          Key(const DiscreteFunctionSpaceType& space)
           : space_( space ), indexSet_( space_.indexSet() )
          {}

          bool operator ==( const Key& other ) const
          {
            // compare grid pointers 
            return (&indexSet_) == (& other.indexSet_ );
          }
          const DiscreteFunctionSpaceType& space() const { return space_; }
        };

        typedef DomainIterator ObjectType;
        typedef Key KeyType;

        inline static ObjectType *createObject ( const KeyType &key )
        {
          return new ObjectType( key.space() );
        }

        inline static void deleteObject ( ObjectType *object )
        {
          delete object;
        }
      };


     typedef typename IteratorFactory :: KeyType KeyType;
     typedef SingletonList
      < KeyType, DomainIterator, IteratorFactory > IndexSetProviderType;

    protected:  
      DomainIterator& iterators_;

    public:  
      //! contructor creating thread iterators 
      explicit DomainDecomposedIteratorStorage( const DiscreteFunctionSpaceType& spc )
        : iterators_( IndexSetProviderType::getObject( KeyType( spc ) ) )
      {
        update();
      }

      ~DomainDecomposedIteratorStorage() 
      {
        IndexSetProviderType::removeObject( iterators_ );
      }

      //! return reference to space 
      const SpaceType& space() const { return iterators_.space(); }

      //! return filter for given thread 
      const FilterType& filter( const int thread ) const 
      {
        return iterators_.filter( thread );
      }

      //! update internal list of iterators 
      void update() 
      {
        iterators_.update();
      }

      //! return begin iterator for current thread 
      IteratorType begin() const 
      {
        return iterators_.begin();
      }

      //! return end iterator for current thread 
      IteratorType end() const 
      {
        return iterators_.end();
      }

      //! return thread number this entity belongs to 
      int index(const EntityType& entity ) const 
      {
        return iterators_.index( entity );
      }

      //! return thread number this entity belongs to 
      int thread(const EntityType& entity ) const 
      {
        return iterators_.thread( entity );
      }
    };
  }
}

#endif
