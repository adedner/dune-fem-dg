#ifndef DUNE_THREADPARTITIONER_HH
#define DUNE_THREADPARTITIONER_HH

//- system includes 
#include <string>
#include <list>
#include <map>
#include <vector>

#include <dune/common/mpihelper.hh>
#include <dune/grid/alugrid.hh>

#include <dune/grid/io/visual/grapegriddisplay.hh>

#if ALU3DGRID_PARALLEL

#warning "Using the Partitioner"

#define ALU2DGRID_READALL_GRIDS 

namespace Dune {

template < class GridPartImp >
class ThreadPartitioner 
{

protected:  
  typedef GridPartImp GridPartType;
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType;
  typedef typename LocalIdSetType :: IdType IdType;

  typedef typename GridPartType :: IndexSetType  IndexSetType;
  typedef typename IndexSetType :: IndexType IndexType;

protected:  
  typedef ALU3DSPACE LoadBalancer LoadBalancerType;
  typedef typename LoadBalancerType :: DataBase DataBaseType;

  // type of communicator 
  typedef ALU3DSPACE MpAccessLocal MPAccessInterfaceType;
  typedef ALU3DSPACE MpAccessMPI   MPAccessImplType;
  mutable MPAccessImplType  mpAccess_;

  DataBaseType db_;

  const GridPartType& gridPart_;
  const IndexSetType& indexSet_;

  typedef typename GridPartType :: template Codim<0> :: EntityType         EntityType;
  typedef typename GridPartType :: template Codim<0> :: EntityPointerType  EntityPointerType;

  // load balancer bounds 
  const double ldbOver_ ;
  const double ldbUnder_;
  const int pSize_ ;

  std::vector< int > index_;
  int indexCounter_ ;

  std::vector< int > partition_;

public:
  ThreadPartitioner( const GridPartType& gridPart, const int pSize )
    : mpAccess_( MPIHelper::getCommunicator() ),
      db_ (),
      gridPart_( gridPart )
    , indexSet_( gridPart_.indexSet() )
    , ldbOver_(1.2)
    , ldbUnder_(0.0)
    , pSize_( pSize )
    , index_( indexSet_.size( 0 ), -1 )
    , indexCounter_( 0 )
  {
    calculateGraph( gridPart_ );
  }

  int getIndex( const size_t idx ) 
  {
    assert( idx < index_.size() );
    if( index_[ idx ] < 0 ) index_[ idx ] = indexCounter_ ++ ;
    return index_[ idx ] ;
  }

  int getIndex( const EntityType& entity ) 
  {
    return getIndex( indexSet_.index( entity ) );
  }

  int index( const EntityType& entity ) const 
  {
    const size_t idx = indexSet_.index( entity );
    assert( idx < index_.size() );
    assert( index_[ idx ] >= 0 );
    return index_[ idx ];
  }

  void calculateGraph( const GridPartType& gridPart )
  {
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType Iterator;
    const Iterator end = gridPart.template end<0> ();
    for(Iterator it = gridPart.template begin<0> (); it != end; ++it )
    {
      const EntityType& entity = *it;
      //assert( entity.partitionType() == InteriorEntity );
      ldbUpdateVertex ( entity,
                        gridPart.ibegin( entity ),
                        gridPart.iend( entity ),
                        db_ );
    }
  }

  template <class IntersectionIteratorType>
  void ldbUpdateVertex ( const EntityType & entity, 
                         const IntersectionIteratorType& ibegin,
                         const IntersectionIteratorType& iend,
                         DataBaseType & db )
  {
    int weight = 1; // a least weight 1 for macro element 

    // calculate weight, which is number of children 
    {
      const int mxl = gridPart_.grid().maxLevel();
      if( mxl > entity.level() ) 
      {
        typedef typename EntityType :: HierarchicIterator HierIt; 
        HierIt endit = entity.hend( mxl ); 
        for(HierIt it = entity.hbegin( mxl ); it != endit; ++it)
        {
          ++weight;
        }
      }

      enum { dim = GridType :: dimension };

#ifdef GRAPHVERTEX_WITH_CENTER
      double center[ 3 ] = { 0.0, 0.0, 0.0 };
      {
        typedef typename GridType :: ctype ctype;
        const GenericReferenceElement< ctype, dim > & refElem =
          GenericReferenceElements< ctype, dim>::general(entity.type());

        FieldVector<ctype, dim> c = 
          entity.geometry().global( refElem.position(0,0) );

        for( int i = 0; i < dim ; ++ i ) 
        {
          center[ i ] = c[ i ];
        }
      }
#endif

      db.vertexUpdate ( typename LoadBalancerType :: 
           GraphVertex ( getIndex( entity ), weight
#ifdef GRAPHVERTEX_WITH_CENTER
             , center
#endif
             ) ) ;
    }
    
    // set weight for faces (to be revised)
    updateFaces( entity, ibegin, iend, weight, db );   
  }

  template <class IntersectionIteratorType>
  void updateFaces(const EntityType& en,
                   IntersectionIteratorType nit,
                   const IntersectionIteratorType endit,
                   const int weight,
                   DataBaseType & db) 
  {
    for( ; nit != endit; ++nit )
    {
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;
      const IntersectionType& intersection = *nit;
      if(intersection.neighbor())
      {
        EntityPointerType ep = intersection.outside();
        const EntityType& nb = *ep;

        if( nb.partitionType() == InteriorEntity )
        {
          const int eid = getIndex( en );
          const int nid = getIndex( nb );
          if( eid < nid ) 
          {
            db.edgeUpdate ( 
                typename LoadBalancerType :: GraphEdge ( eid, nid, weight )
                );
          }
        }
      }
    }
  }

public:
  bool repartition() 
  {
    return db_.repartition( mpAccess_,  DataBaseType :: METIS_PartGraphKway );
  }

  bool serialPartition(const bool useKway = true ) 
  {
#ifdef ITERATORS_WITHOUT_MYALLOC
    if( pSize_ > 1 ) 
    {
      if( useKway ) 
        partition_ = db_.repartition( mpAccess_, DataBaseType :: METIS_PartGraphKway, pSize_ );
      else 
        partition_ = db_.repartition( mpAccess_, DataBaseType :: METIS_PartGraphRecursive, pSize_ );
      assert( partition_.size() > 0 );
      /*
      std::vector< int > counter( pSize_ , 0 );
      for( size_t i =0; i<partition_.size(); ++i) 
      {
        std::cout << "part[" << i << "] = " << partition_[ i ]  << endl;
        ++counter[  partition_[ i ]  ];
      }

      for( int i=0; i<pSize_; ++ i) 
        std::cout << counter[ i ] << " counter \n";
      */
      return partition_.size() > 0;
    }
    else 
#endif 
    {
      partition_.resize( indexSet_.size( 0 ) );
      for( size_t i =0; i<partition_.size(); ++i )
        partition_[ i ] = 0;
      return false ;
    }

  }

  bool checkPartition() const 
  {
    bool neu = false ;
    {
      const int np = mpAccess_.psize ();
      // Kriterium, wann eine Lastneuverteilung vorzunehmen ist:
      // 
      // load  - eigene ElementLast
      // mean  - mittlere ElementLast
      // nload - Lastverh"altnis

      double load = db_.accVertexLoad () ;
      std :: vector < double > v (mpAccess_.gcollect (load)) ;
      double mean = std::accumulate (v.begin (), v.end (), 0.0) / double (np) ;

      for (std :: vector < double > :: iterator i = v.begin () ; i != v.end () ; ++i)
        neu |= (*i > mean ? (*i > (ldbOver_ * mean) ? true : false) : 
            (*i < (ldbUnder_ * mean) ? true : false)) ;
    }

    int val = (neu) ? 1 : 0; 
    // get global maximum 
    int v2  = mpAccess_.gmax( val ); 

    return (v2 == 1) ? true : false;
  }

  int getDestination( const int idx ) const 
  {
    return db_.getDestination ( idx ) ;
  }

  std::set < int, std::less < int > > scan() const 
  {
    return db_.scan();
  }

  int getRank( const EntityType& entity ) const 
  {
    assert( (int) partition_.size() > index( entity ) );
    return partition_[ index( entity ) ];
  }

  bool validEntity( const EntityType& entity, const int rank ) const 
  {
    return getRank( entity ) == rank;
  }

};

} // end namespace Dune 
#else 
#warning "Parallel ALUGrid not available"
#endif // ALU3DGRID_PARALLEL
#endif
