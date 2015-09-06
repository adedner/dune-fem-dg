#ifndef DUNE_FEM_FLUXPROJECTION_HH
#define DUNE_FEM_FLUXPROJECTION_HH

// #include <dune/grid/utility/twistutility.hh>

#include <dune/fem/misc/compatibility.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/lagrange.hh>

namespace Dune
{
  namespace Fem
  {
    struct FluxProjectionImpl
    {
      template <int dim>
      struct FakeQuad 
      {
        typedef Dune::FieldVector<double,dim-1> FaceDomainType;
        typedef Dune::FieldVector<double,dim> DomainType;
        FakeQuad( const FaceDomainType &xf, const DomainType &x) : xface_(xf), x_(x) {}
        const FaceDomainType& localPoint( int ) const
        {
          return xface_;
        }
        const DomainType& point( int ) const
        {
          return x_;
        }
        private:
        const FaceDomainType& xface_;
        const DomainType& x_;
      };
      template< class Function, class DiscreteFunction, class Weight, class Flux >
      static void project ( const Function &f, DiscreteFunction &u, Weight &weight, const Flux &flux )
      {
        typedef typename Function::FunctionSpaceType FunctionSpaceType;
        typedef typename Function::LocalFunctionType LocalFunctionType;

        typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef typename DiscreteFunction::LocalFunctionType LocalDiscreteFunctionType;

        typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
        typedef typename DiscreteFunctionSpaceType::LagrangePointSetType LagrangePointSetType;

        typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
        typedef typename GridPartType::IntersectionType IntersectionType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
        typedef typename GridPartType::template Codim< 0 >::GeometryType GeometryType;;

        typedef typename LagrangePointSetType::template Codim< 1 >::SubEntityIteratorType FaceDofIteratorType;

        typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
        typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
        typedef typename FunctionSpaceType::RangeType RangeType;
        typedef typename GeometryType::LocalCoordinate LocalCoordinateType;
        typedef typename IntersectionType::LocalGeometry LocalGeometryType;

        const int dim = GridPartType::dimension;
        const unsigned int dimRange = FunctionSpaceType::dimRange;
        const DiscreteFunctionSpaceType &space = u.space();

        u.clear();
        DiscreteFunction w ( "weight", space );
        w.clear();

        for ( const auto& entity : space )
        {
          weight.setEntity( entity );

          LocalDiscreteFunctionType lw = w.localFunction( entity );
          LocalDiscreteFunctionType lu = u.localFunction( entity );

          const LocalFunctionType lf = f.localFunction( entity );
          const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );

          const unsigned int numPoints = lagrangePointSet.nop();
          for( unsigned int pt = 0; pt < numPoints; ++pt )
          {
            // codim = 1 is handled in a special way
            if (lagrangePointSet.localKey(pt).codim() == 1)
              continue;
            RangeType val;
            lf.evaluate( lagrangePointSet[ pt ], val );

            RealType wght = weight( lagrangePointSet.point( pt ) );

            for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
            {
              lu[ dimRange*pt + coordinate ] += wght*val[ coordinate ];
              lw[ dimRange*pt + coordinate ] += wght;
            }
          }
          auto gridPart = space.gridPart();
          const IntersectionIteratorType iend = gridPart.iend( entity );
          for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != iend; ++iit )
          {
            const IntersectionType &intersection = *iit;
            const int indexInInside = intersection.indexInInside();
            const FaceDofIteratorType fdend = lagrangePointSet.template endSubEntity< 1 >( indexInInside );
            FaceDofIteratorType fdit = lagrangePointSet.template beginSubEntity< 1 >( indexInInside );

            if( intersection.neighbor() )
            {
              // get neighbor
              EntityType neighbor = make_entity( intersection.outside() );
              const LocalGeometryType &geoIn  = intersection.geometryInInside();
              const LocalGeometryType &geoOut = intersection.geometryInOutside();
              const LocalFunctionType lfNb = f.localFunction( neighbor );
              for( ; fdit != fdend; ++fdit )
              {
                unsigned int pt = *fdit;
                const LocalCoordinateType &xIn = lagrangePointSet.point( pt );
                auto xface = geoIn.local( xIn );
                const LocalCoordinateType xOut = geoOut.global( xface );
                RangeType valEn, valNb, ustar;
                lf.evaluate(  xIn, valEn );
                lfNb.evaluate( xOut, valNb );
                FakeQuad<dim> quadEn(xface,xOut);
                FakeQuad<dim> quadNb(xface,xIn);
                flux.uStar( intersection, entity, neighbor, 0.0 , quadEn, quadNb, 0,
                            valEn, valNb, ustar );
                for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
                {
                  lu[ dimRange*pt + coordinate ] += ustar[ coordinate ];
                  lw[ dimRange*pt + coordinate ] += 1.;
                }
              }
            }
            else
            {
              std::cout << "a boundary? Bad reconstruction..." << std::endl;
              assert(0);
            }
          }

        }

        u.communicate();
        w.communicate();

        typedef typename DiscreteFunction::DofIteratorType DofIteratorType;

        const DofIteratorType udend = u.dend();
        DofIteratorType udit = u.dbegin();
        DofIteratorType wdit = w.dbegin();
        for( ; udit != udend; ++udit, ++wdit )
        {
          // assert( (*wdit > 0.) || (*udit == 0.) );
          RealType weight = std::abs( *wdit );
          if ( weight > 1e-12 )
            *udit /= weight;
          else
            std::cout << "error in weight: " << weight << " [fluxprojection.hh]" << std::endl;
        }

        // make function continuous over hanging nodes

        if( !GridPartType::Traits::conforming && Fem::GridPartCapabilities::hasGrid< GridPartType >::v)
        {
          const GridPartType &gridPart =  space.gridPart();
          for( const auto& entity : space )
          {
            const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );

            const IntersectionIteratorType iend = gridPart.iend( entity );
            for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != iend; ++iit )
            {
              const IntersectionType &intersection = *iit;

              if( intersection.neighbor() )
              {
                // get neighbor
                EntityType neighbor = make_entity( intersection.outside() );

                // if non-conforming situation
                if( entity.level() > neighbor.level() )
                {
                  const int indexInInside = intersection.indexInInside();

                  typedef typename IntersectionType::LocalGeometry LocalGeometryType;
                  const LocalGeometryType &geoIn  = intersection.geometryInInside();
                  const LocalGeometryType &geoOut = intersection.geometryInOutside();

                  LocalDiscreteFunctionType uIn  = u.localFunction( entity );
                  LocalDiscreteFunctionType uOut = u.localFunction( neighbor );

                  const FaceDofIteratorType fdend = lagrangePointSet.template endSubEntity< 1 >( indexInInside );
                  FaceDofIteratorType fdit = lagrangePointSet.template beginSubEntity< 1 >( indexInInside );
                  for( ; fdit != fdend; ++fdit )
                  {
                    const LocalCoordinateType &xIn = lagrangePointSet.point( *fdit );
                    const LocalCoordinateType xOut = geoOut.global( geoIn.local( xIn ) );

                    RangeType val;
                    uOut.evaluate( xOut, val );

                    for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
                      uIn[ dimRange*(*fdit) + coordinate ] = val[ coordinate ];
                  }
                }
              }
            }
          }
        }
      }
      template< class Function, class DiscreteFunction, class Flux >
      static void project ( const Function &f, DiscreteFunction &u, const Flux &flux )
      {
        typedef typename DiscreteFunction::DiscreteFunctionSpaceType::GridPartType GridPartType;
        WeightDefault<GridPartType> weight;
        project(f,u,weight,flux);
      }
    };

    struct FullProjectionImpl
    {
      template <int dim>
      struct FakeQuad 
      {
        typedef Dune::FieldVector<double,dim-1> FaceDomainType;
        typedef Dune::FieldVector<double,dim> DomainType;
        FakeQuad( const FaceDomainType &xf, const DomainType &x) : xface_(xf), x_(x) {}
        const FaceDomainType& localPoint( int ) const
        {
          return xface_;
        }
        const DomainType& point( int ) const
        {
          return x_;
        }
        private:
        const FaceDomainType& xface_;
        const DomainType& x_;
      };
      template< class Function, class DiscreteFunction, class Flux >
      static void project ( const Function &f, DiscreteFunction &u, const Flux &flux )
      {
        typedef typename Function::FunctionSpaceType FunctionSpaceType;
        typedef typename Function::LocalFunctionType LocalFunctionType;

        typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef typename DiscreteFunction::LocalFunctionType LocalDiscreteFunctionType;

        typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
        typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
        typedef typename GridPartType::IntersectionType IntersectionType;
        typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
        typedef typename GridPartType::template Codim< 0 >::GeometryType GeometryType;;

        typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
        typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;
        typedef typename FunctionSpaceType::RangeType RangeType;
        typedef typename GeometryType::LocalCoordinate LocalCoordinateType;
        typedef typename IntersectionType::LocalGeometry LocalGeometryType;

        const int dim = GridPartType::dimension;
        const unsigned int dimRange = FunctionSpaceType::dimRange;
        const DiscreteFunctionSpaceType &space = u.space();

        u.clear();

        for ( const auto& entity : space )
        {
          const LocalFunctionType lf = f.localFunction( entity );
          LocalDiscreteFunctionType lu = u.localFunction( entity );
          assert(lu.size() == lf.size()+dimRange);
          for (int i=0;i<lf.size()-dimRange;++i)
            lu[i] = lf[i];

          // make continuous
          Dune::FieldMatrix<double,2*dimRange,2*dimRange> matrix;
          Dune::FieldVector<double,2*dimRange> rhs,sol;
          auto gridPart = space.gridPart();
          const IntersectionIteratorType iend = gridPart.iend( entity );
          for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != iend; ++iit )
          {
            const IntersectionType &intersection = *iit;
            const int indexInInside = intersection.indexInInside();
            if( intersection.neighbor() )
            {
              // U = u + u0*phi0(x) + u1*phi1(x)
              // u0r*phi0(x_i)[r] + u1r*phi1(x_i)[r] = ustar(x_i)[r] - u(x_i)[r]
              EntityType neighbor = make_entity( intersection.outside() );
              const LocalGeometryType &geoIn  = intersection.geometryInInside();
              const LocalGeometryType &geoOut = intersection.geometryInOutside();
              typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
              FaceQuadratureType faceQuadInside(space.gridPart(), intersection,1,FaceQuadratureType::INSIDE);
              auto &baseSet = lu.basisFunctionSet();
              std::vector< RangeType > phiFaceEn(lu.size());
              baseSet.evaluateAll( faceQuadInside[0], phiFaceEn );
              RangeType uVal;
              lu.evaluate( faceQuadInside[0],uVal);

              const LocalFunctionType lfNb = f.localFunction( neighbor );
              auto& xIn = faceQuadInside.point(0);
              auto xface = geoIn.local( xIn );
              const LocalCoordinateType xOut = geoOut.global( xface );
              RangeType valEn, valNb, ustar;
              lf.evaluate(  xIn, valEn );
              lfNb.evaluate( xOut, valNb );
              FakeQuad<dim> quadEn(xface,xOut);
              FakeQuad<dim> quadNb(xface,xIn);
              flux.uStar( intersection, entity, neighbor, 0.0 , quadEn, quadNb, 0,
                          valEn, valNb, ustar );
              int offset = intersection.indexInInside()*dimRange;
              for (int r=0;r<dimRange;++r)
              {
                rhs[offset+r] = ustar[r]-uVal[r];
                for (int j=0;j<2;++j)
                  for (int rr=0;rr<dimRange;++rr)
                  {
                    matrix[offset+r][j*dimRange+rr] = phiFaceEn[lf.size()/dimRange-1+j][rr];
                    std::cout << matrix[offset+r][j*dimRange+rr] << " ";
                  }
                std::cout << std::endl;
              }
              std::cout << std::endl;
              std::cout << std::endl;
            }
          }
          matrix.solve(sol,rhs);
          for (int i=0;i<2;++i)
          {
            int offset = lf.size()-dimRange+i*dimRange;
            for (int r=0;r<dimRange;++r)
              lu[offset+r] = sol[offset+r];
          }
        }
        u.communicate();
      }
    };

  } // namespace Fem

} // name space Dune

#endif // #ifndef DUNE_FEM_VTXPROJECTION_HH
