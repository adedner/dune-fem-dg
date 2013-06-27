#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH

// C++ includes
#include <cstddef>
#include <vector>

// dune-common includes
#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

// dune-fem includes
#include <dune/fem/misc/functor.hh>
#include <dune/fem/quadrature/caching/registry.hh>
#include <dune/fem/quadrature/cachingpointlist.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

namespace Dune
{

  namespace Fem
  {

    // CachingShapeFunctionSet
    // -----------------------

    template< class ShapeFunctionSet >
    class CachingShapeFunctionSet
    : private QuadratureStorageRegistry::StorageInterface
    {
      typedef CachingShapeFunctionSet< ShapeFunctionSet > ThisType;

    public:
      typedef typename ShapeFunctionSet::FunctionSpaceType FunctionSpaceType;
      
      typedef typename ShapeFunctionSet::DomainType DomainType;
      typedef typename ShapeFunctionSet::RangeType RangeType;
      typedef typename ShapeFunctionSet::JacobianRangeType JacobianRangeType;
      typedef typename ShapeFunctionSet::HessianRangeType HessianRangeType;

      typedef MutableArray< MutableArray<RangeType> >         RangeVectorType;
      typedef MutableArray< MutableArray<JacobianRangeType> > JacobianRangeVectorType;

      typedef std::vector< RangeVectorType >          ValueCacheVectorType ;
      typedef std::vector< JacobianRangeVectorType >  JacobianCacheVectorType;

      explicit CachingShapeFunctionSet ( const GeometryType &type,
                                         const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : type_( type ),
        shapeFunctionSet_( shapeFunctionSet ),
        valueCaches_(),
        jacobianCaches_(),
        localRangeCache_( ThreadManager::maxThreads() ),
        localJacobianCache_( ThreadManager::maxThreads() )
      {
        QuadratureStorageRegistry::registerStorage( *this );
      }

      ~CachingShapeFunctionSet ();

      int order () const { return shapeFunctionSet_.order(); }

      std::size_t size () const
      {
        return shapeFunctionSet_.size();
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.evaluateEach( x, functor );
      }

      template< class Quadrature, class Functor >
      void evaluateEach ( const QuadraturePointWrapper< Quadrature > &x, Functor functor ) const
      {
        const bool cacheable = Conversion< Quadrature, CachingInterface >::exists;
        evaluateEach( x.quadrature(), x.point(), functor, integral_constant< bool, cacheable >() );
      }
     
      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.jacobianEach( x, functor );
      }

      template< class Quadrature, class Functor >
      void jacobianEach ( const QuadraturePointWrapper< Quadrature > &x, Functor functor ) const
      {
        const bool cacheable = Conversion< Quadrature, CachingInterface >::exists;
        jacobianEach( x.quadrature(), x.point(), functor, integral_constant< bool, cacheable >() );
      }

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.hessianEach( x, functor );
      }

      GeometryType type () const DUNE_DEPRECATED {  return type_; }
      
      GeometryType geometryType () const DUNE_DEPRECATED {  return type(); }

      template< class Quad, bool cacheable > 
      struct ReturnCache 
      {
        static unsigned int maxCachingPoint( const Quad& quad ) 
        {
          unsigned int maxCp = 0; 
          const unsigned int nop = quad.nop(); 
          for( unsigned int i=0; i<nop; ++i ) 
          {
            const unsigned int cp = quad.cachingPoint( i );
            maxCp = std::max( cp, maxCp );
          }
          return maxCp;
        }

        static const RangeVectorType& 
        ranges( const ThisType& shapeFunctionSet, 
                const Quad& quad,     
                const ValueCacheVectorType&,
                RangeVectorType& storage ) 
        {
          // evaluate all basis functions and multiply with dof value 
          const unsigned int nop  = quad.nop(); 
          const unsigned int size = shapeFunctionSet.size();

          // check caching points 
          assert( maxCachingPoint( quad ) < nop );

          // make sure cache has the appropriate size 
          storage.resize( nop );

          typedef MutableArray<RangeType>         RangeVector;
          for( unsigned int qp = 0 ; qp < nop; ++ qp ) 
          {
            const int cacheQp = quad.cachingPoint( qp );
            storage[ cacheQp ].resize( size );
            AssignFunctor< RangeVector > funztor( storage[ cacheQp ] );
            shapeFunctionSet.evaluateEach( quad[ qp ], funztor );
          }
          return storage;
        }

        static const JacobianRangeVectorType&
        jacobians( const ThisType& shapeFunctionSet, 
                   const Quad& quad,
                   const JacobianCacheVectorType&,
                   JacobianRangeVectorType& storage ) 
        {
          // evaluate all basis functions and multiply with dof value 
          const unsigned int nop  = quad.nop(); 
          const unsigned int size = shapeFunctionSet.size();

          // check caching points 
          assert( maxCachingPoint( quad ) < nop );

          // make sure cache has the appropriate size 
          storage.resize( nop );

          typedef MutableArray<JacobianRangeType> JacobianRangeVector;
          for( unsigned int qp = 0 ; qp < nop; ++ qp ) 
          {
            const int cacheQp = quad.cachingPoint( qp );
            storage[ cacheQp ].resize( size );
            AssignFunctor< JacobianRangeVector > funztor( storage[ cacheQp ] );
            shapeFunctionSet.jacobianEach( quad[ qp ], funztor );
          }
          return storage;
        }
      };

      template< class Quad > 
      struct ReturnCache< Quad, true >  
      {
        static const RangeVectorType& 
        ranges( const ThisType& shapeFunctionSet, 
                const Quad& quad,
                const ValueCacheVectorType& cache,
                const RangeVectorType& ) 
        {
          return cache[ quad.id() ];
        }

        static const JacobianRangeVectorType&  
        jacobians( const ThisType& shapeFunctionSet, 
                   const Quad& quad,
                   const JacobianCacheVectorType& cache,
                   const JacobianRangeVectorType& ) 
        {
          return cache[ quad.id() ];
        }
      };

      template < class QuadratureType > 
      const RangeVectorType& rangeCache( const QuadratureType& quadrature ) const 
      {
        return ReturnCache< QuadratureType, Conversion< QuadratureType, CachingInterface >::exists > :: 
          ranges( *this, quadrature, valueCaches_, localRangeCache_[ ThreadManager::thread() ] );
      }

      template < class QuadratureType > 
      const JacobianRangeVectorType& jacobianCache( const QuadratureType& quadrature ) const 
      {
        return ReturnCache< QuadratureType, Conversion< QuadratureType, CachingInterface >::exists > :: 
          jacobians( *this, quadrature, jacobianCaches_, localJacobianCache_[ ThreadManager::thread() ] );
      }

    private:
      template< class Quadrature, class Functor >
      void evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, false > ) const
      {
        evaluateEach( quadrature.point( pt ), functor );
      }

      template< class Quadrature, class Functor >
      void evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, true > ) const;

      template< class Quadrature, class Functor >
      void jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, false > ) const
      {
        jacobianEach( quadrature.point( pt ), functor );
      }

      template< class Quadrature, class Functor >
      void jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                          integral_constant< bool, true > ) const;


      void cacheQuadrature( std::size_t id, std::size_t codim, std::size_t size );

      template< class PointVector >
      void cachePoints ( std::size_t id, const PointVector &points );

      // prohibit copying and assignment
      CachingShapeFunctionSet ( const ThisType & );
      const ThisType &operator= ( const ThisType & );

      GeometryType type_;
      ShapeFunctionSet shapeFunctionSet_;
      ValueCacheVectorType valueCaches_;
      JacobianCacheVectorType jacobianCaches_;

      mutable ValueCacheVectorType     localRangeCache_ ;
      mutable JacobianCacheVectorType  localJacobianCache_;
    };



    // Implementation of CachingShapeFunctionSet
    // -----------------------------------------

    template< class ShapeFunctionSet >
    inline CachingShapeFunctionSet< ShapeFunctionSet >::~CachingShapeFunctionSet ()
    {
      QuadratureStorageRegistry::unregisterStorage( *this );
    }


    template< class ShapeFunctionSet >
    template< class Quadrature, class Functor >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::evaluateEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                       integral_constant< bool, true > ) const
    {
      typedef MutableArray<RangeType>         RangeVector;
      const RangeVector& cache = valueCaches_[ quadrature.id() ][ quadrature.cachingPoint( pt ) ];
      const std::size_t numShapeFunctions = size();
      for( std::size_t i = 0; i < numShapeFunctions; ++i )
        functor( i, cache[ i ] );
    }


    template< class ShapeFunctionSet >
    template< class Quadrature, class Functor >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::jacobianEach ( const Quadrature &quadrature, std::size_t pt, Functor functor,
                       integral_constant< bool, true > ) const
    {
      typedef MutableArray<JacobianRangeType>       JacobianRangeVector;
      const JacobianRangeVector& cache = jacobianCaches_[ quadrature.id() ][ quadrature.cachingPoint( pt ) ];

      const std::size_t numShapeFunctions = size();
      for( std::size_t i = 0; i < numShapeFunctions; ++i )
        functor( i, cache[ i ] );
    }


    template< class ShapeFunctionSet >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::cacheQuadrature( std::size_t id, std::size_t codim, std::size_t size )
    {
      if( id >= valueCaches_.size() ) 
      {
        valueCaches_.resize( id+10 );
        jacobianCaches_.resize( id+10 );
      }

      if( valueCaches_[ id ].size() == 0 )
      {
        typedef typename FunctionSpaceType::DomainFieldType ctype;
        const int dim = FunctionSpaceType::dimDomain;
        switch( codim )
        {
        case 0:
          cachePoints( id, PointProvider< ctype, dim, 0 >::getPoints( id, type_ ) );
          break;

        case 1:
          cachePoints( id, PointProvider< ctype, dim, 1 >::getPoints( id, type_ ) );
          break;

        default:
          DUNE_THROW( NotImplemented, "Caching for codim > 1 not implemented." );
        }
      }
    }


    template< class ShapeFunctionSet >
    template< class PointVector >
    inline void CachingShapeFunctionSet< ShapeFunctionSet >
      ::cachePoints ( std::size_t id, const PointVector &points )
    {
      const std::size_t numShapeFunctions = size();
      const std::size_t numPoints = points.size();

      valueCaches_[ id ].resize( numPoints );
      jacobianCaches_[ id ].resize( numPoints );

      for( size_t pt=0; pt<numPoints; ++ pt )
      {
        valueCaches_[ id ][ pt ].resize( numShapeFunctions );
        jacobianCaches_[ id ][ pt ].resize( numShapeFunctions );
      }

      typedef MutableArray<RangeType>         RangeVector;
      typedef MutableArray<JacobianRangeType> JacobianRangeVector;
      for( std::size_t pt = 0; pt < numPoints; ++pt )
      {
        evaluateEach( points[ pt ], AssignFunctor< RangeVector >( valueCaches_[ id ][ pt ] ) );
        jacobianEach( points[ pt ], AssignFunctor< JacobianRangeVector >( jacobianCaches_[ id ][ pt ] ) );
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CACHING_HH
