#ifndef DUNE_FEMDG_DGNORM_HH
#define DUNE_FEMDG_DGNORM_HH

#include <dune/fem/misc/l2norm.hh>

namespace Dune
{

namespace Fem
{

  template< class GridPart >
  class DGNorm
  : public IntegralBase< GridPart, DGNorm< GridPart> >
  {
    typedef DGNorm< GridPart > ThisType;
    typedef IntegralBase< GridPart, DGNorm< GridPart> > BaseType;

  public:
    typedef GridPart GridPartType;

    template< class Function >
    struct FunctionJacobianSquare;

  protected:

    typedef typename BaseType::EntityType                   EntityType;
    typedef typename EntityType::Geometry                   Geometry;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType         IntersectionType;
    typedef typename IntersectionType::Geometry             IntersectionGeometryType;

    typedef CachingQuadrature< GridPartType, 0 >            QuadratureType;
    typedef ElementQuadrature< GridPartType, 1 >            FaceQuadratureType;
  public:
    typedef Integrator< QuadratureType >                    IntegratorType;


    using BaseType::gridPart;
    using BaseType::comm;

  public:

    /** \brief constructor
     *    \param gridPart     specific gridPart for selection of entities
     *    \param order        order of integration quadrature (default = 2*space.order())
     *    \param communicate  if true global (over all ranks) norm is computed (default = true)
     */
    explicit DGNorm ( const GridPartType &gridPart, const unsigned int order = 0, const bool communicate = true );
    DGNorm ( const ThisType &other );

    //! || u ||_H1 + || [ u ] ||_L2(\Gamma) on given set of entities (partition set)
    template< class DiscreteFunctionType, class PartitionSet >
    typename DiscreteFunctionType::RangeFieldType
    norm ( const DiscreteFunctionType &u, const PartitionSet &partition ) const;

    //! || u ||_H1 + || [ u ] ||_L2(\Gamma_I) on interior partition entities (partition set)
    template< class DiscreteFunctionType >
    typename DiscreteFunctionType::RangeFieldType
    norm ( const DiscreteFunctionType &u ) const
    {
      return norm( u, Partitions::interior );
    }

    //! || u - v ||_H1 + || [ u - v ] ||_L2(\Gamma) on given set of entities (partition set)
    template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
    typename UDiscreteFunctionType::RangeFieldType
    distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v, const PartitionSet &partition ) const;

    //! || u - v ||_H1 + || [ u - v ] ||_L2(\Gamma_I) on interior partition entities
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    typename UDiscreteFunctionType::RangeFieldType
    distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
    {
      return distance( u, v, Partitions::interior );
    }

    template< class ULocalFunction,
              class VLocalFunction,
              class ReturnType >
    inline void
    distanceLocal ( const EntityType& entity, const unsigned int order,
                    const ULocalFunction &u,
                    const VLocalFunction &v,
                    ReturnType& sum ) const ;

    template< class ULocalFunction,
              class ReturnType >
    inline void
    normLocal ( const EntityType& entity, const unsigned int order,
                const ULocalFunction &u, ReturnType& sum ) const ;

  private:
    // prohibit assignment
    ThisType operator= ( const ThisType &other );
    const unsigned int order_;
    const bool communicate_;
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
  inline DGNorm< GridPart >::DGNorm ( const GridPartType &gridPart, unsigned int order, const bool communicate )
  : BaseType( gridPart ),
    order_( order ),
    communicate_( BaseType::checkCommunicateFlag( communicate ) )
  {}



  template< class GridPart >
  inline DGNorm< GridPart >::DGNorm ( const ThisType &other )
  : BaseType( other ),
    order_( other.order_ ),
    communicate_( other.communicate_ )
  {}


  template< class GridPart >
  template< class DiscreteFunctionType, class PartitionSet >
  inline typename DiscreteFunctionType::RangeFieldType
  DGNorm< GridPart >::norm ( const DiscreteFunctionType &u, const PartitionSet &partition ) const
  {
    typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > ReturnType ;

    ReturnType sum = BaseType :: forEach( u, ReturnType( 0 ), partition, order_ );

    // communicate_ indicates global norm
    if( communicate_ )
      sum[ 0 ] = comm().sum( sum[ 0 ] );

    // return result, e.g. sqrt of calculated sum
    return sqrt( sum[ 0 ] );
  }

  template< class GridPart >
  template< class UDiscreteFunctionType, class VDiscreteFunctionType, class PartitionSet >
  inline typename UDiscreteFunctionType::RangeFieldType
  DGNorm< GridPart >::distance ( const UDiscreteFunctionType &u,
                                 const VDiscreteFunctionType &v,
                                 const PartitionSet &partition ) const
  {
    typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
    typedef FieldVector< RangeFieldType, 1 > ReturnType ;

    ReturnType sum = BaseType :: forEach( u, v, ReturnType( 0 ), partition, order_ );

    // communicate_ indicates global norm
    if( communicate_ )
      sum[ 0 ] = comm().sum( sum[ 0 ] );

    // return result, e.g. sqrt of calculated sum
    return sqrt( sum[ 0 ] );
  }

  template< class GridPart >
  template< class ULocalFunction,
            class VLocalFunction,
            class ReturnType >
  inline void
  DGNorm< GridPart >::distanceLocal ( const EntityType& entity, const unsigned int order,
                                      const ULocalFunction &ulocal,
                                      const VLocalFunction &vlocal,
                                      ReturnType& sum ) const
  {
    typedef typename L2Norm< GridPart >::template FunctionDistance< ULocalFunction, VLocalFunction >
      LocalDistanceType;

    IntegratorType integrator( order );

    LocalDistanceType dist( ulocal, vlocal );
    FunctionJacobianSquare< LocalDistanceType > dist2( dist );

    integrator.integrateAdd( entity, dist2, sum );

    unsigned int enIdx = gridPart().indexSet().index(entity);
    const Geometry& geometry = entity.geometry();

    ULocalFunction ulocalNb( ulocal );
    VLocalFunction vlocalNb( vlocal );

    double jumpTerm = 0;
    {
      for (const auto& intersection : intersections(gridPart(), entity) )
      {
        if( intersection.neighbor() )
        {
          const EntityType& neighbor = intersection.outside();
          const Geometry& geometryNb = neighbor.geometry();

          unsigned int nbIdx = gridPart().indexSet().index(neighbor);
          if( (enIdx < nbIdx) || (neighbor.partitionType() != Dune::InteriorEntity) )
          {
            const auto intersectionGeometry = intersection.geometry();

            const double intersectionArea = intersectionGeometry.volume();
            const double heInverse = intersectionArea / std::min( geometry.volume(), geometryNb.volume() );
            vlocalNb.init( neighbor );
            ulocalNb.init( neighbor );
            LocalDistanceType distNb( ulocalNb, vlocalNb );

            unsigned int nbOrder = std::max( uint(2 * ulocalNb.order()) , order );
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
  template< class LocalFunction, class ReturnType >
  inline void
  DGNorm< GridPart >::normLocal ( const EntityType& entity, const unsigned int order,
                                  const LocalFunction &ulocal,
                                  ReturnType& sum ) const
  {
    // evaluate norm locally
    IntegratorType integrator( order );

    FunctionJacobianSquare< LocalFunction > ulocal2( ulocal );
    integrator.integrateAdd( entity, ulocal2, sum );

    unsigned int enIdx = gridPart().indexSet().index(entity);
    const Geometry& geometry = entity.geometry();

    // this should work as long as LocalFunction is either ConstLocalFunction< ... > or MutableLocalFunction< ... >
    LocalFunction ulocalNb( ulocal );

    double jumpTerm = 0;
    {
      for (const auto& intersection : intersections(gridPart(), entity) )
      {
        if( intersection.neighbor() )
        {
          const EntityType& neighbor = intersection.outside();
          const Geometry& geometryNb = neighbor.geometry();
          unsigned int nbIdx = gridPart().indexSet().index(neighbor);
          if( (enIdx < nbIdx) || (neighbor.partitionType() != Dune::InteriorEntity) )
          {
            const double intersectionArea = intersection.geometry().volume();
            const double heInverse = intersectionArea / std::min( geometry.volume(), geometryNb.volume() );
            ulocalNb.init( neighbor ); // local u on neighbor element

            unsigned int nbOrder = std::max( uint(2 * ulocalNb.order()) , order );
            FaceQuadratureType quadInside( gridPart(), intersection, nbOrder, FaceQuadratureType::INSIDE );
            FaceQuadratureType quadOutside( gridPart(), intersection, nbOrder, FaceQuadratureType::OUTSIDE );
            const size_t numQuadraturePoints = quadInside.nop();
            for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
            {
              const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
              typename LocalFunction::RangeType  distIn(0),distOut(0),jump(0);
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
