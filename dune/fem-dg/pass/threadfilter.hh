#ifndef DUNE_THREADFILTER_HH
#define DUNE_THREADFILTER_HH

//- System includes
#include <vector>
#include <cassert>

//- Dune includes
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/filteredgrid.hh>

#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {

//***************************************************************************
// 
// Thread Filter 
//
//***************************************************************************
  template <class GridPartType>
  class ThreadFilter : 
    public FilterDefaultImplementation<
      DefaultFilterTraits<ThreadFilter<GridPartType>,GridPartType> >
  {
    typedef typename GridPartType :: GridType GridType;
  public:
    typedef DefaultFilterTraits<ThreadFilter<GridPartType>, GridPartType> Traits;
    typedef FilterInterface<DefaultFilterTraits<ThreadFilter<GridPartType>,GridPartType> > BaseType;

    using BaseType :: has0Entity;
    
    typedef typename BaseType::FilterType FilterType;
    typedef typename BaseType::EntityCodim0Type EntityCodim0Type;
    typedef typename BaseType::EntityPointerCodim0Type EntityPointerCodim0Type;

  private:
    typedef typename GridType::Traits::template Codim<0>::Geometry GeometryType;
    enum{dim = GeometryType::dimension};
    typedef typename Dune::FieldVector<typename GridType::ctype, dim> FieldVectorType;
    typedef typename GridType :: ctype coordType;
    typedef typename GridPartType :: IndexSetType IndexSetType;

    const GridPartType& gridPart_;
    const IndexSetType& indexSet_;
    const MutableArray< int >& threadNum_ ;
    int thread_;
  public:
    ThreadFilter(const GridPartType& gridPart, const MutableArray< int > & threadNum )
     : gridPart_(gridPart),
       indexSet_( gridPart.indexSet() ),
       threadNum_( threadNum ),
       thread_( 0 )
    {
    }
    
    //! copy constructor 
    ThreadFilter(const ThreadFilter& org) : 
      gridPart_(org.gridPart_) ,
      indexSet_(org.indexSet_) , 
      threadNum_( org.threadNum_ ),
      thread_( org.thread_ )
    {
    }

    void setThread ( const int thread )
    {
      thread_ = thread ;
    }

    //! return true if entity is in FFS domain 
    inline bool has0Entity(const EntityPointerCodim0Type & e) const 
    {
      return has0Entity( *e );
    }

    //! return true if entity is in FFS domain 
    template <class EntityType>
    inline bool has0Entity(const EntityType & e) const 
    {
      return threadNum_[ indexSet_.index(e) ] == thread_ ;
    }

    // true if all intersections are interior intersections 
    template <class IntersectionType>
    inline bool interiorIntersection(const IntersectionType & it) const 
    {
      return false;
    }
    
    template <class IntersectionType>
    inline bool intersectionBoundary(const IntersectionType & it) const 
    {
      return it.boundary();
    }
    
    template <class IntersectionType>
    inline int intersectionBoundaryId(const IntersectionType & it) const 
    {
      // default value for lower region 
      return ( it.boundary () ) ? it.boundaryId() : 111;
    }
    
    template <class IntersectionType>
    inline bool intersectionNeighbor(const IntersectionType & it) const 
    {
      return it.neighbor();
    }
  }; // end ThreadFilter

}  // end namespace Dune
#endif
