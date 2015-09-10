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
      template< class Function, class DiscreteFunction, class Flux >
      static void project ( const Function &f, DiscreteFunction &u, const Flux &flux )
      {
        typedef typename DiscreteFunction::DiscreteFunctionSpaceType::GridPartType GridPartType;
        WeightDefault<GridPartType> weight;

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
                FakeQuad<dim> quadEn(xface,xIn);
                FakeQuad<dim> quadNb(xface,xOut);
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
      }
      template< class Function, class DiscreteFunction, class Flux >
      static void project ( const Function &f, const Function &df, DiscreteFunction &u, const Flux &flux )
      {
        typedef typename DiscreteFunction::DiscreteFunctionSpaceType::GridPartType GridPartType;
        WeightDefault<GridPartType> weight;
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
        typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
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

          const LocalFunctionType ldf = df.localFunction( entity );
          const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );

          const unsigned int numPoints = lagrangePointSet.nop();
          for( unsigned int pt = 0; pt < numPoints; ++pt )
          {
            // codim = 1 is handled in a special way
            if (lagrangePointSet.localKey(pt).codim() == 1)
              continue;
            RangeType val;
            ldf.evaluate( lagrangePointSet[ pt ], val );

            RealType wght = weight( lagrangePointSet.point( pt ) );

            for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
            {
              lu[ dimRange*pt + coordinate ] += wght*val[ coordinate ];
              lw[ dimRange*pt + coordinate ] += wght;
            }
          }

          const LocalFunctionType lf = f.localFunction( entity );
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
              const LocalFunctionType ldfNb = df.localFunction( neighbor );
              for( ; fdit != fdend; ++fdit )
              {
                unsigned int pt = *fdit;
                const LocalCoordinateType &xIn = lagrangePointSet.point( pt );
                auto xface = geoIn.local( xIn );
                const LocalCoordinateType xOut = geoOut.global( xface );
                RangeType valEn, valNb, dvalEn, dvalNb, dustar;
                lf.evaluate(  xIn, valEn );
                lfNb.evaluate( xOut, valNb );
                ldf.evaluate(  xIn, dvalEn );
                ldfNb.evaluate( xOut, dvalNb );
                FakeQuad<dim> quadEn(xface,xIn);
                FakeQuad<dim> quadNb(xface,xOut);
                flux.duStar( intersection, entity, neighbor, 0.0 , quadEn, quadNb, 0,
                             valEn, valNb, dvalEn, dvalNb, dustar );
                for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
                {
                  lu[ dimRange*pt + coordinate ] += dustar[ coordinate ];
                  lw[ dimRange*pt + coordinate ] += 1.;
                }
              }
            }
            else
            {
              std::cout << "a boundary? Bad reconstruction..." << std::endl;
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
      }
    };

    struct L2FluxProjectionImpl
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
        typedef typename DiscreteFunction::DiscreteFunctionSpaceType::GridPartType GridPartType;

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

        typedef CachingQuadrature<GridPartType,0> QuadratureType;
      
        const int dim = GridPartType::dimension;
        const unsigned int dimRange = FunctionSpaceType::dimRange;
        const DiscreteFunctionSpaceType &space = u.space();

        u.clear();

        enum { localBlockSize = DiscreteFunctionSpaceType :: localBlockSize };
        enum { dgNumDofs = localBlockSize };

        const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;
        typedef Dune::Fem::DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polynomialOrder-2> TestSpaceType;
        TestSpaceType testSpace(space.gridPart());
        const int testDofs = TestSpaceType :: localBlockSize;
        const int numDofs = testDofs + 2*dimRange;

        typedef Dune::FieldMatrix< double, numDofs, numDofs > MatrixType;
        typedef Dune::FieldVector< double, numDofs >         VectorType;
        std::vector<RangeType> phi(numDofs);
        std::vector<RangeType> psi(testDofs);

        for ( const auto& entity : space )
        {
          MatrixType matrix(0);
          VectorType rhs(0),sol(0);
        
          const LocalFunctionType lf = f.localFunction( entity );

          // get geometry
          const auto& geo = entity.geometry();
          // get quadrature
          QuadratureType quad(entity, space.order()*2);
          const int quadNop = quad.nop();
          RangeType value ;
          const auto& set = space.basisFunctionSet(entity);
          const auto& test = testSpace.basisFunctionSet(entity);
          for(int qP = 0; qP < quadNop ; ++qP)
          {
            const double intel = quad.weight(qP) * geo.integrationElement( quad.point(qP) );
            // evaluate function
            lf.evaluate(quad[ qP ], value );
            set.evaluateAll(quad[qP], phi);
            test.evaluateAll(quad[qP],psi);
            for(int n=0; n<psi.size(); ++n)
            {
              const RangeType& psi_n = psi[n];
              const double val = intel * (value * psi_n);
              rhs[n] += val;
              for(int m=0; m<phi.size(); ++m)
              {
                const RangeType& phi_m = phi[m];
                const double val = intel * (phi_m * psi_n);
                matrix[n][m] += val;
              }
            }
          }
          int row = psi.size();
          const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );
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
                FakeQuad<dim> quadEn(xface,xIn);
                FakeQuad<dim> quadNb(xface,xOut);
                flux.uStar( intersection, entity, neighbor, 0.0 , quadEn, quadNb, 0,
                            valEn, valNb, ustar );
                // ustar = valEn; ustar += valNb; ustar *= 0.5;
                for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
                {
                  assert(row<rhs.size());
                  rhs[row] = ustar[ coordinate ];
                  matrix[row][dimRange*pt+coordinate] = 1.;
                  ++row;
                }
              }
            }
            else
            {
              std::cout << "a boundary? Bad reconstruction..." << std::endl;
            }
          }
          assert(row==matrix.rows);
          matrix.solve(sol,rhs);

          LocalDiscreteFunctionType lu = u.localFunction( entity );
          assert(sol.size()==lu.size());
          for (int i=0;i<sol.size();++i)
            lu[i] = sol[i];
        }
      }
      template< class Function, class DiscreteFunction, class Flux >
      static void project ( const Function &f, const Function &df, DiscreteFunction &u, const Flux &flux )
      {
        typedef typename DiscreteFunction::DiscreteFunctionSpaceType::GridPartType GridPartType;
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
        typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
        typedef typename GeometryType::LocalCoordinate LocalCoordinateType;
        typedef typename IntersectionType::LocalGeometry LocalGeometryType;

        typedef CachingQuadrature<GridPartType,0> QuadratureType;

        const int dim = GridPartType::dimension;
        const unsigned int dimRange = FunctionSpaceType::dimRange;
        const DiscreteFunctionSpaceType &space = u.space();

        u.clear();

        enum { localBlockSize = DiscreteFunctionSpaceType :: localBlockSize };
        enum { dgNumDofs = localBlockSize };

        const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;
        typedef Dune::Fem::DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, polynomialOrder-2> TestSpaceType;
        TestSpaceType testSpace(space.gridPart());
        const int testDofs = TestSpaceType :: localBlockSize;
        const int numDofs = testDofs + 2*dimRange;

        typedef Dune::FieldMatrix< double, numDofs, numDofs > MatrixType;
        typedef Dune::FieldVector< double, numDofs >         VectorType;
        std::vector<RangeType> phi(numDofs);
        std::vector<RangeType> psi(testDofs);

        for ( const auto& entity : space )
        {
          MatrixType matrix(0);
          VectorType rhs(0),sol(0);
        
          const LocalFunctionType ldf = df.localFunction( entity );

          // get geometry
          const auto& geo = entity.geometry();
          // get quadrature
          QuadratureType quad(entity, space.order()*2);
          const int quadNop = quad.nop();
          RangeType value ;
          const auto& set = space.basisFunctionSet(entity);
          const auto& test = testSpace.basisFunctionSet(entity);
          for(int qP = 0; qP < quadNop ; ++qP)
          {
            const double intel = quad.weight(qP) * geo.integrationElement( quad.point(qP) );
            // evaluate function
            ldf.evaluate(quad[ qP ], value );
            set.evaluateAll(quad[qP], phi);
            test.evaluateAll(quad[qP],psi);
            for(int n=0; n<psi.size(); ++n)
            {
              const RangeType& psi_n = psi[n];
              const double val = intel * (value * psi_n);
              rhs[n] += val;
              for(int m=0; m<phi.size(); ++m)
              {
                const RangeType& phi_m = phi[m];
                const double val = intel * (phi_m * psi_n);
                matrix[n][m] += val;
              }
            }
          }
          int row = psi.size();
          const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );
          const LocalFunctionType lf = f.localFunction( entity );
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
              const LocalFunctionType ldfNb = df.localFunction( neighbor );
              for( ; fdit != fdend; ++fdit )
              {
                unsigned int pt = *fdit;
                const LocalCoordinateType &xIn = lagrangePointSet.point( pt );
                auto xface = geoIn.local( xIn );
                const LocalCoordinateType xOut = geoOut.global( xface );
                RangeType valEn, valNb, dvalEn, dvalNb, dustar;
                lf.evaluate(  xIn, valEn );
                lfNb.evaluate( xOut, valNb );
                ldf.evaluate(  xIn, dvalEn );
                ldfNb.evaluate( xOut, dvalNb );
                FakeQuad<dim> quadEn(xface,xIn);
                FakeQuad<dim> quadNb(xface,xOut);
                flux.duStar( intersection, entity, neighbor, 0.0 , quadEn, quadNb, 0,
                             valEn, valNb, dvalEn, dvalNb, dustar );
                // dustar = dvalEn; dustar += dvalNb; dustar *= 0.5;
                for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
                {
                  assert(row<rhs.size());
                  rhs[row] = dustar[ coordinate ];
                  matrix[row][dimRange*pt+coordinate] = 1.;
                  ++row;
                }
              }
            }
            else
            {
              std::cout << "a boundary? Bad reconstruction..." << std::endl;
            }
          }
          assert(row==matrix.rows);
          matrix.solve(sol,rhs);

          LocalDiscreteFunctionType lu = u.localFunction( entity );
          assert(sol.size()==lu.size());
          for (int i=0;i<sol.size();++i)
            lu[i] = sol[i];
        }
      }
    };

  } // namespace Fem

} // name space Dune

#endif // #ifndef DUNE_FEM_VTXPROJECTION_HH
