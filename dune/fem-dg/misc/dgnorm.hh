#ifndef DUNE_FEMDG_DGNORM_HH
#define DUNE_FEMDG_DGNORM_HH

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/compatibility.hh>

namespace Dune
{

namespace Fem
{

  template< class GridPart >
  class DGNorm
  : public LPNormBase< GridPart, DGNorm< GridPart> >
  {
    typedef DGNorm< GridPart > ThisType;
    typedef LPNormBase< GridPart, DGNorm< GridPart> > BaseType;

  public:
    typedef GridPart GridPartType;

    template< class Function >
    struct FunctionJacobianSquare;

  protected:
    typedef typename BaseType::GridIteratorType GridIteratorType;

    typedef typename GridIteratorType::Entity EntityType;
    typedef typename EntityType::Geometry Geometry;

    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType :: Intersection IntersectionType ;
    typedef typename IntersectionType::Geometry IntersectionGeometryType;

    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
    typedef ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  public:
    typedef Integrator< QuadratureType > IntegratorType;


    using BaseType::gridPart;
    using BaseType::comm;

  public:
    explicit DGNorm ( const GridPartType &gridPart );
    DGNorm ( const ThisType &other );

    template< class DiscreteFunctionType >
    typename DiscreteFunctionType::RangeFieldType
    norm ( const DiscreteFunctionType &u ) const;

    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    typename UDiscreteFunctionType::RangeFieldType
    distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const;

    template< class UDiscreteFunctionType,
              class VDiscreteFunctionType,
              class ReturnType >
    inline void
    distanceLocal ( const EntityType& entity, const unsigned int order,
                    const UDiscreteFunctionType &u,
                    const VDiscreteFunctionType &v,
                    ReturnType& sum ) const ;

    template< class UDiscreteFunctionType,
              class ReturnType >
    inline void
    normLocal ( const EntityType& entity, const unsigned int order,
                    const UDiscreteFunctionType &u,
                    ReturnType& sum ) const ;

  private:
    // prohibit assignment
    ThisType operator= ( const ThisType &other );
  };



  // DGNorm::FunctionJacobianSquare
  // ------------------------------

  template< class GridPart >
  template< class Function >
  struct DGNorm< GridPart >::FunctionJacobianSquare
  {
    typedef Function FunctionType;

    typedef typename FunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > RangeType;

  public:
    explicit FunctionJacobianSquare ( const FunctionType &function )
    : function_( function )
    {}

    template< class Point >
    void evaluate ( const Point &x, RangeType &ret ) const
    {
      const int dimRange = FunctionType::RangeType::dimension;

      typename FunctionType::RangeType phi;
      function_.evaluate( x, phi );
      ret[ 0 ] = phi * phi;

      typename FunctionType::JacobianRangeType grad;
      function_.jacobian( x, grad );
      for( int i = 0; i < dimRange; ++i )
        ret[ 0 ] += (grad[ i ] * grad[ i ]);
    }

  private:
    const FunctionType &function_;
  };



  // Implementation of DG Norm
  // -------------------------

  template< class GridPart >
  inline DGNorm< GridPart >::DGNorm ( const GridPartType &gridPart )
  : BaseType( gridPart )
  {}



  template< class GridPart >
  inline DGNorm< GridPart >::DGNorm ( const ThisType &other )
  : BaseType( other )
  {}


  template< class GridPart >
  template< class DiscreteFunctionType >
  inline typename DiscreteFunctionType::RangeFieldType
  DGNorm< GridPart >::norm ( const DiscreteFunctionType &u ) const
  {
    typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > ReturnType ;

    ReturnType sum = BaseType :: forEach( u, ReturnType( 0 ) );

    // return result, e.g. sqrt of calculated sum
    return sqrt( comm().sum( sum[ 0 ] ) );
  }

  template< class GridPart >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType >
  inline typename UDiscreteFunctionType::RangeFieldType
  DGNorm< GridPart >::distance ( const UDiscreteFunctionType &u,
                                 const VDiscreteFunctionType &v ) const
  {
    typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > ReturnType ;

    ReturnType sum = BaseType :: forEach( u, v, ReturnType( 0 ) );

    // return result, e.g. sqrt of calculated sum
    return sqrt( comm().sum( sum[ 0 ] ) );
  }

  template< class GridPart >
  template< class UDiscreteFunctionType,
            class VDiscreteFunctionType,
            class ReturnType >
  inline void
  DGNorm< GridPart >::distanceLocal ( const EntityType& entity, const unsigned int order,
                                      const UDiscreteFunctionType &u,
                                      const VDiscreteFunctionType &v,
                                      ReturnType& sum ) const
  {
    typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
    typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;
    ULocalFunctionType ulocal = u.localFunction( entity );
    VLocalFunctionType vlocal = v.localFunction( entity );

    typedef typename L2Norm< GridPart >::template FunctionDistance< ULocalFunctionType, VLocalFunctionType >
      LocalDistanceType;

    IntegratorType integrator( order );

    LocalDistanceType dist( ulocal, vlocal );
    FunctionJacobianSquare< LocalDistanceType > dist2( dist );

    integrator.integrateAdd( entity, dist2, sum );

    unsigned int enIdx = gridPart().indexSet().index(entity);
    const Geometry& geometry = entity.geometry();

    double jumpTerm = 0;
    {
      const IntersectionIteratorType endiit = gridPart().iend( entity );
      for ( IntersectionIteratorType iit = gridPart().ibegin( entity ); iit != endiit ; ++ iit )
      {
        const IntersectionType& intersection = *iit ;
        if( intersection.neighbor() )
        {
          const EntityType& neighbor = make_entity( intersection.outside() );
          const Geometry& geometryNb = neighbor.geometry();


          unsigned int nbIdx = gridPart().indexSet().index(neighbor);
          unsigned int nbOrder = std::max( uint(2 * u.space().order( neighbor )) , order );
          if( (enIdx < nbIdx) || (neighbor.partitionType() != Dune::InteriorEntity) )
          {
            typedef typename IntersectionType :: Geometry IntersectionGeometry;
            const IntersectionGeometry intersectionGeometry = intersection.geometry();

            const double intersectionArea = intersectionGeometry.volume();
            const double heInverse = intersectionArea / std::min( geometry.volume(), geometryNb.volume() );
            ULocalFunctionType ulocalNb = u.localFunction( neighbor ); // local u on neighbor element
            VLocalFunctionType vlocalNb = v.localFunction( neighbor ); // local u on neighbor element
            LocalDistanceType distNb( ulocalNb, vlocalNb );

            FaceQuadratureType quadInside ( gridPart(), intersection, nbOrder, FaceQuadratureType::INSIDE  );
            FaceQuadratureType quadOutside( gridPart(), intersection, nbOrder, FaceQuadratureType::OUTSIDE );
            const size_t numQuadraturePoints = quadInside.nop();
            for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
            {
              const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
              typename LocalDistanceType::RangeType  distIn(0),distOut(0),jump(0);
              dist.evaluate( quadInside[ pt ], distIn );
              distNb.evaluate( quadOutside[ pt ], distOut );
              jump = distIn - distOut;
              double weight = quadInside.weight( pt )*heInverse * intersectionGeometry.integrationElement( x );
              jumpTerm += (jump*jump) * weight;
            }
          }
        }
      }
    }
    sum[0] += jumpTerm;
  }


  template< class GridPart >
  template< class DiscreteFunctionType, class ReturnType >
  inline void
  DGNorm< GridPart >::normLocal ( const EntityType& entity, const unsigned int order,
                                  const DiscreteFunctionType &u,
                                  ReturnType& sum ) const
  {
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
    // evaluate norm locally

    IntegratorType integrator( order );

    LocalFunctionType ulocal = u.localFunction( entity );
    FunctionJacobianSquare< LocalFunctionType > ulocal2( ulocal );

    integrator.integrateAdd( entity, ulocal2, sum );

    unsigned int enIdx = gridPart().indexSet().index(entity);
    const Geometry& geometry = entity.geometry();

    double jumpTerm = 0;
    {
      const IntersectionIteratorType endiit = gridPart().iend( entity );
      for ( IntersectionIteratorType iit = gridPart().ibegin( entity ); iit != endiit ; ++ iit )
      {
        const IntersectionType& intersection = *iit ;
        if( intersection.neighbor() )
        {
          const EntityType& neighbor = make_entity( intersection.outside() );
          const Geometry& geometryNb = neighbor.geometry();
          unsigned int nbIdx = gridPart().indexSet().index(neighbor);
          unsigned int nbOrder = std::max( uint(2 * u.space().order( neighbor )) , order );
          if( (enIdx < nbIdx) || (neighbor.partitionType() != Dune::InteriorEntity) )
          {
            const double intersectionArea = intersection.geometry().volume();
            const double heInverse = intersectionArea / std::min( geometry.volume(), geometryNb.volume() );
            LocalFunctionType ulocalNb = u.localFunction( neighbor ); // local u on neighbor element

            FaceQuadratureType quadInside( gridPart(), intersection, nbOrder, FaceQuadratureType::INSIDE );
            FaceQuadratureType quadOutside( gridPart(), intersection, nbOrder, FaceQuadratureType::OUTSIDE );
            const size_t numQuadraturePoints = quadInside.nop();
            for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
            {
              const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
              typename LocalFunctionType::RangeType  distIn(0),distOut(0),jump(0);
              ulocal.evaluate( quadInside[ pt ], distIn );
              ulocalNb.evaluate( quadOutside[ pt ], distOut );
              jump = distIn - distOut;
              double weight = quadInside.weight( pt )*heInverse * intersection.geometry().integrationElement( x );
              jumpTerm += (jump*jump) * weight;
            }
          }
        }
      }
    }
    sum[0] += jumpTerm;
  }

}

using Fem :: DGNorm ;

} // end namespace Dune

#endif // #ifndef DUNE_FEM_DGNORM_HH
