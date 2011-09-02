#ifndef PRIMALMATRIXASSEMBLY_HH
#define PRIMALMATRIXASSEMBLY_HH

#include <dune/fem-dg/operator/fluxes/dgprimalfluxes.hh>

namespace Dune {

template <class DestinationImp,
          class Model, 
          DGDiffusionFluxIdentifier = method_general >
class DGPrimalMatrixAssembly
{
  public:
  typedef DestinationImp DestinationType;
  typedef typename DestinationType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType; 
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType; 
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType; 
  typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename DestinationType::LocalFunctionType LocalFunctionType;

  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry IntersectionGeometryType;

  // need treatment of non conforming grids
  typedef Dune::ElementQuadrature< GridPartType, 1 > FaceQuadratureType; 
  typedef Dune::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef typename GridPartType :: GridType :: template Codim< 0 > :: EntityPointer EntityPointerType;

  typedef DGPrimalDiffusionFlux<DiscreteFunctionSpaceType,Model,method_general> FluxType;

  struct RangeValues
  {
    typedef std::vector< std::vector< RangeType > > VectorType;
    RangeValues(int col, const VectorType &vec) : col_(col), vec_(vec), zero_(0) 
    {}
    const RangeType &operator[](int row) const
    {
      if (col_==-1) return zero_;
      else return vec_[row][col_];
    }
    private:
    int col_;
    VectorType vec_;
    const RangeType zero_;
  };
  struct JacobianRangeValues
  {
    typedef std::vector< std::vector< JacobianRangeType > > VectorType;
    JacobianRangeValues(int col, const VectorType &vec) : col_(col), vec_(vec), zero_(0) 
    {}
    const JacobianRangeType &operator[](int row) const
    {
      if (col_==-1) return zero_;
      else return vec_[row][col_];
    }
    private:
    int col_;
    VectorType vec_;
    const JacobianRangeType zero_;
  };
  struct IntersectionStorage
  {
    IntersectionStorage( const EntityType &entity, const EntityType &neigbor, const double areaEn, const double areaNb )
    : en_(entity), nb_(neigbor), areaEn_(areaEn), areaNb_(areaNb)
    {}
    const EntityType &inside() const { return en_; }
    const EntityType &outside() const { return nb_; }
    const double enVolume() const { return areaEn_; }
    const double nbVolume() const { return areaNb_; }
    private:
    const EntityType &en_,&nb_;
    const double areaEn_, areaNb_;
  };

  public:

  DGPrimalMatrixAssembly( GridPartType& gridPart,
                          const Model& model )
    : model_(model),
      space_(gridPart),
      flux_(gridPart, model)
  {}

  const DiscreteFunctionSpaceType &space() const
  {
    return space_;
  }

  /* 
   * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
   * implementation.
   */
  template <class Matrix> 
  void assemble( const double time, 
                 Matrix& matrix, DestinationType& rhs ) const 
  {
    typedef typename Matrix::LocalMatrixType LocalMatrixType;
    // matrix.reserve();
    matrix.clear();
    rhs.clear();

    const DiscreteFunctionSpaceType &dfSpace = rhs.space();
    const size_t numBaseFunctions = dfSpace.mapper().maxNumDofs();

    std::vector< RangeType > phi( numBaseFunctions );
    std::vector< JacobianRangeType > dphi( numBaseFunctions );

    const RangeType uZero(0);
    const JacobianRangeType uJacZero(0);

    const IteratorType end = dfSpace.end();
    for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      LocalMatrixType localOpEn = matrix.localMatrix( entity, entity );
      LocalFunctionType rhsLocal = rhs.localFunction( entity );

      const BaseFunctionSetType &baseSet = localOpEn.domainBaseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
            
      QuadratureType quadrature( entity, 2*dfSpace.order() );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        const typename GeometryType::Jacobian &gjit = geometry.jacobianInverseTransposed( x );

        baseSet.evaluateAll( quadrature[ pt ], phi );
        baseSet.jacobianAll( quadrature[ pt ], gjit, dphi );

        RangeType arhs(0);
        // first assemble rhs (evaluate source with u=0)
        if ( model_.hasStiffSource() )
          model_.stiffSource( entity, time, x, uZero, uJacZero, arhs );
        if ( model_.hasNonStiffSource() )
        {
          RangeType sNonStiff (0);
          model_.nonStiffSource( entity, time, x, uZero, uJacZero, sNonStiff );
          arhs += sNonStiff;
        }

        for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
        {
          // now assemble source part depending on u (mass part)
          RangeType aphi(0);
          if ( model_.hasStiffSource() )
            model_.stiffSource( entity, time, x, phi[localCol], dphi[localCol], aphi );
          if ( model_.hasNonStiffSource() )
          {
            RangeType sNonStiff (0);
            model_.nonStiffSource( entity, time, x, phi[localCol], dphi[localCol], sNonStiff );
            aphi += sNonStiff;
          }
          // subtract affine part and move to left hand side
          aphi -= arhs;
          aphi *= -1;

          JacobianRangeType adphi;
          model_.diffusion(entity, time, x, phi[localCol], dphi[localCol], adphi);

          // get column object and call axpy method 
          localOpEn.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
        }
        // assemble rhs
        arhs *= weight;
        rhsLocal.axpy( quadrature[ pt ], arhs );
      }

      const IntersectionIteratorType endiit = dfSpace.gridPart().iend( entity );
      for ( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity );
            iit != endiit ; ++ iit )
      {
        const IntersectionType& intersection = *iit ;

        if( intersection.neighbor() )
        {
          EntityPointerType ep = intersection.outside();
          const EntityType& neighbor = *ep ;
          const GeometryType& nbGeometry = neighbor.geometry();

          // get local matrix for face entries 
          LocalMatrixType localOpNb = matrix.localMatrix( entity, neighbor );
          // get neighbor's base function set 
          const BaseFunctionSetType &baseSetNb = localOpNb.domainBaseFunctionSet();

          FaceQuadratureType faceQuadInside(dfSpace.gridPart(), intersection, 2*dfSpace.order()+1,
                                            FaceQuadratureType::INSIDE);
          FaceQuadratureType faceQuadOutside(dfSpace.gridPart(), intersection, 2*dfSpace.order()+1,
                                             FaceQuadratureType::OUTSIDE);

          const size_t numFaceQuadPoints = faceQuadInside.nop();

          // store all basis functions
          std::vector< std::vector< RangeType > > phiFaceEn( numFaceQuadPoints );
          std::vector< std::vector< JacobianRangeType > > dphiFaceEn( numFaceQuadPoints );
          std::vector< std::vector< RangeType > > phiFaceNb( numFaceQuadPoints );
          std::vector< std::vector< JacobianRangeType > > dphiFaceNb( numFaceQuadPoints );
          // evalute base functions
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            phiFaceEn[pt].resize( numBaseFunctions );
            dphiFaceEn[pt].resize( numBaseFunctions );
            phiFaceNb[pt].resize( numBaseFunctions );
            dphiFaceNb[pt].resize( numBaseFunctions );
            baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
            baseSet.jacobianAll( faceQuadInside[ pt ], 
                                 geometry.jacobianInverseTransposed( faceQuadInside.point( pt ) ), dphiFaceEn[pt] );
            baseSetNb.evaluateAll( faceQuadOutside[ pt ], phiFaceNb[pt] );
            baseSetNb.jacobianAll( faceQuadOutside[ pt ], 
                                   nbGeometry.jacobianInverseTransposed( faceQuadOutside.point( pt ) ), dphiFaceNb[pt] );
          }

          // storage for all flux values
          std::vector< RangeType > valueEn( numFaceQuadPoints );
          std::vector< JacobianRangeType > dvalueEn( numFaceQuadPoints );
          std::vector< RangeType > valueNb( numFaceQuadPoints );
          std::vector< JacobianRangeType > dvalueNb( numFaceQuadPoints );
          // compute fluxes and assemble matrix
          for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
          {
            // compute flux for one base function, i.e.,
            // - uLeft=phiFaceEn[.][localCol]
            // - uRight=0
            flux(intersection, entity, neighbor, time, faceQuadInside, faceQuadOutside,
                 RangeValues(localCol,phiFaceEn), JacobianRangeValues(localCol,dphiFaceEn),
                 RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
                 valueEn, dvalueEn, valueNb, dvalueNb );
            for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            {
              const double weight = faceQuadInside.weight( pt );
              localOpEn.column( localCol ).axpy( phiFaceEn[pt], dphiFaceEn[pt], 
                                                 valueEn[pt], dvalueEn[pt], weight );
              dvalueEn[pt] *= -1; // could also be taken from right value of numericalFlux...
              localOpNb.column( localCol ).axpy( phiFaceNb[pt], dphiFaceNb[pt], 
                                                 valueEn[pt], dvalueEn[pt], -weight );
            }
          }
        }
        else if ( intersection.boundary() )
        {
          IntersectionStorage intersectionStorage( entity, entity, geometry.volume(), geometry.volume() );
        
          FaceQuadratureType faceQuadInside(dfSpace.gridPart(), intersection, 2*dfSpace.order(),
                                            FaceQuadratureType::INSIDE);

          const size_t numFaceQuadPoints = faceQuadInside.nop();

          // store all basis functions
          std::vector< std::vector< RangeType > > phiFaceEn( numFaceQuadPoints );
          std::vector< std::vector< JacobianRangeType > > dphiFaceEn( numFaceQuadPoints );
          // evalute base functions
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            phiFaceEn[pt].resize( numBaseFunctions );
            dphiFaceEn[pt].resize( numBaseFunctions );
            baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
            baseSet.jacobianAll( faceQuadInside[ pt ], 
                                 geometry.jacobianInverseTransposed( faceQuadInside.point( pt ) ), dphiFaceEn[pt] );
          }

          // storage for all flux values
          std::vector< RangeType > valueEn( numFaceQuadPoints );
          std::vector< JacobianRangeType > dvalueEn( numFaceQuadPoints );
          std::vector< RangeType > boundaryValues( numFaceQuadPoints, RangeType(0) );

          // compute boundary fluxes
          for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
          {
            // compute flux for one base function, i.e.,
            // - uLeft=phiFaceEn[.][localCol]
            // - uRight=0
            flux_.initializeBoundary( intersection, entity, time, faceQuadInside, 
                                      RangeValues(localCol,phiFaceEn), boundaryValues ); // boundaryValues[i]=0
            for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            {
              flux_.boundaryFlux( intersection, intersectionStorage, 
                                  time, faceQuadInside, pt, 
                                  phiFaceEn[ pt ][ localCol ], uZero,
                                  dphiFaceEn[ pt ][ localCol ], 
                                  valueEn[ pt ], dvalueEn[ pt ]);
              const double weight = faceQuadInside.weight( pt );
              localOpEn.column( localCol ).axpy( phiFaceEn[pt], dphiFaceEn[pt], 
                                                 valueEn[pt], dvalueEn[pt], weight );
            }
          }

          // get boundary value
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            model_.boundaryValue(intersection, time,  faceQuadInside.localPoint( pt ), 
                                 uZero, boundaryValues[ pt ]);
          // rhs 
          RangeType rhsFlux;
          JacobianRangeType drhsFlux;
          flux_.initializeBoundary( intersection, entity, time, faceQuadInside, 
                                    RangeValues(-1,phiFaceEn), boundaryValues );
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            const double weight = faceQuadInside.weight( pt );
            flux_.boundaryFlux( intersection, intersectionStorage, 
                                time, faceQuadInside, pt, 
                                uZero, boundaryValues[ pt ], uJacZero, 
                                rhsFlux, drhsFlux);
            rhsFlux *= -weight;
            drhsFlux *= -weight;
            rhsLocal.axpy( faceQuadInside[ pt ], rhsFlux, drhsFlux );
          }
        }
      }
    }
    {
      double maxSymError_diagonal = 0;
      double maxSymError_offdiagonal = 0;
      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
      {
        const EntityType& entity = *it;
        LocalMatrixType localOpEn = matrix.localMatrix( entity, entity  );
        const BaseFunctionSetType &baseSet = localOpEn.domainBaseFunctionSet();
        const unsigned int numBaseFunctions = baseSet.numBaseFunctions();

        for( unsigned int r = 0; r < numBaseFunctions; ++r )
        {
          for( unsigned int c = 0; c < numBaseFunctions; ++c )
          {
            double v1 = localOpEn.get(r,c);
            double v2 = localOpEn.get(c,r);
            // if (std::abs(v1-v2)>1e-8)
            //  std::cout << "symmetry error on diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
            maxSymError_diagonal = std::max(maxSymError_diagonal,std::abs(v1-v2));
          }
        }
        const IntersectionIteratorType endiit = dfSpace.gridPart().iend( entity );
        for ( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity );
              iit != endiit ; ++ iit )
        {
          const IntersectionType& intersection = *iit ;
          if ( !intersection.neighbor() )
            continue;
          EntityPointerType ep = intersection.outside();
          const EntityType& neighbor = *ep ;

          LocalMatrixType localOpNb1 = matrix.localMatrix( entity, neighbor );
          LocalMatrixType localOpNb2 = matrix.localMatrix( neighbor, entity );
          for( unsigned int r = 0; r < numBaseFunctions; ++r )
          {
            for( unsigned int c = 0; c < numBaseFunctions; ++c )
            {
              double v1 = localOpNb1.get(r,c);
              double v2 = localOpNb2.get(c,r);
              // if (std::abs(v1-v2)>1e-8)
              //    std::cout << "symmetry error on off diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
              maxSymError_offdiagonal = std::max(maxSymError_offdiagonal,std::abs(v1-v2));
            }
          }
        }
      }
      if (maxSymError_diagonal > 1e-8)
        std::cout << "non symmetric diagonal: " << maxSymError_diagonal << std::endl;
      if (maxSymError_offdiagonal > 1e-8)
        std::cout << "non symmetric off diagonal: " << maxSymError_offdiagonal << std::endl;
    }
  }

  template <class Quadrature,class Value,class DValue,class RetType, class DRetType>
  void flux(const IntersectionType &intersection,
            const EntityType &entity, const EntityType &neighbor,
            const double time,
            const Quadrature &faceQuadInside, const Quadrature &faceQuadOutside,
            const Value &valueEn, const DValue &dvalueEn,
            const Value &valueNb, const DValue &dvalueNb,
            RetType &retEn, DRetType &dretEn, 
            RetType &retNb, DRetType &dretNb) const
  {
    IntersectionStorage intersectionStorage( entity, neighbor, entity.geometry().volume(), neighbor.geometry().volume() );
    flux_.initializeIntersection( intersection, entity, neighbor, time, faceQuadInside, faceQuadOutside,
                                  valueEn, valueNb );

    const size_t numFaceQuadPoints = faceQuadInside.nop();
    for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
    {
      flux_.numericalFlux( intersection, intersectionStorage, 
                           time, faceQuadInside, faceQuadOutside, pt, 
                           valueEn[pt], valueNb[pt], dvalueEn[pt], dvalueNb[pt],
                           retEn[ pt ], retNb[ pt ],
                           dretEn[ pt ], dretNb[ pt ]);
    }
  }


  private:
  const Model &model_;
  DiscreteFunctionSpaceType space_;
  mutable FluxType flux_;
};


};

#endif
