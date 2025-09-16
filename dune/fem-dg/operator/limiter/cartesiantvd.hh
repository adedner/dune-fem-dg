#ifndef DUNE_FV_CARTESIANRECONSTRUCTION_HH
#define DUNE_FV_CARTESIANRECONSTRUCTION_HH

#include <cassert>
#include <cstddef>

#include <numeric>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/gridpart/common/capabilities.hh>

#include <dune/fem-dg/operator/limiter/limiterutility.hh>


namespace Dune
{

  namespace FV
  {

    // TVDReconstruction
    // ----------------

    /**
     * \class TVDReconstruction
     * \brief Minmod-type reconstruction based on dimensional splitting approach
     * \endcode
     **/
    template< class GP, class SV, class BV >
    class TVDReconstruction
    {
      typedef TVDReconstruction< GP, SV, BV > This;

    public:
      typedef GP GridPartType;
      typedef SV StateVector;
      typedef BV BoundaryValue;

      typedef StateVector RangeType;

      typedef FieldVector< typename GridPartType::ctype, GridPartType::dimensionworld > GlobalCoordinate;

      typedef typename GridPartType::Intersection Intersection;

      typedef typename FieldTraits< StateVector >::field_type Field;
      typedef typename FieldTraits< StateVector >::real_type Real;
      typedef FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > Jacobian;

      static const int dimension = GridPartType::dimension;
      static const bool isCartesian = Dune::Fem::GridPartCapabilities::isCartesian< GridPartType >::v;
      static_assert( isCartesian, "TVDReconstruction requires a Cartesian grid");

      static const int numFunc = 3;
    public:
      TVDReconstruction ( const GridPartType &gp, BoundaryValue boundaryValue, Real tolerance )
        : gridPart_( gp ),
          boundaryValue_( std::move( boundaryValue ) ),
          limiterFunction_()
      {
        const auto& mapper = gridPart().indexSet();
        neighbors_.resize( mapper.size(0) );

        const auto end = gridPart().template end< 0, Dune::InteriorBorder_Partition >();
        for( auto it = gridPart().template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it )
        {
          const auto& element = *it ;
          const std::size_t elIndex = mapper.index( element );
          auto& neighbors = neighbors_[ elIndex ];

          const GlobalCoordinate elCenter = element.geometry().center();

          const auto iend = gridPart().iend( element );
          for( auto iit = gridPart().ibegin( element ); iit != iend; ++iit )
          {
            const auto intersection = *iit;

            if( intersection.neighbor() )
            {
              const auto neighbor = intersection.outside();
              const std::size_t nbIndex = mapper.index( neighbor );
              neighbors[ intersection.indexInInside() ] = nbIndex;
              const GlobalCoordinate nbCenter = neighbor.geometry().center();
              centerDiff_[ intersection.indexInInside() ] = nbCenter - elCenter;
            }
            else if ( intersection.boundary() )
            {
              // for boundary faces set inside element index
              neighbors[ intersection.indexInInside() ] = elIndex;
            }
          }
        }

        //for( int i=0; i<2*dimension; ++i )
        //  std::cout << "diff[ " << i << " ] = " << centerDiff_[ i ] << std::endl;

        for( int d=0; d<dimension; ++d )
        {
          h_[ d ] = std::abs(centerDiff_[d*2][d]);
          //std::cout << "h[ " << d << " ] = " << h_[ d ] << std::endl;
        }

        //         inside
        //   1       0       2
        //   * ----- * ----- *
        //       h       h

        combos_[ 0 ][ 0 ] = 0;
        combos_[ 0 ][ 1 ] = 1;
        combos_[ 1 ][ 0 ] = 0;
        combos_[ 1 ][ 1 ] = 2;
        combos_[ 2 ][ 0 ] = 1;
        combos_[ 2 ][ 1 ] = 2;

        testset_[0] = std::vector< int8_t > (1, 1); // 1,0 test with 2,0 which is 1
        testset_[1] = std::vector< int8_t > (1, 0); // 2,0 test with 1,0 which is 0
        if( numFunc > 2 )
        {
          // 2,1
          testset_[2] = std::vector< int8_t > (2);
          testset_[2][0] = 0;
          testset_[2][1] = 1;
        }
      }

      template< class Entity, class Mapper, class Vector >
      void applyLocal ( const Entity& element,
                        const Mapper &mapper,
                        const Vector &u,
                        Jacobian& du ) const
      {
        static const int dim = dimension;
        static const int dimRange = Jacobian::rows;

        const std::size_t elIndex = mapper.index( element );
        const GlobalCoordinate elCenter = element.geometry().center();

        // TODO check boundary based on neighs

        const auto& neighs = neighbors_[ elIndex ];
        const RangeType &enVal = u[ elIndex ];

        // neighboring values
        std::array< RangeType, 2*dim > values_;
        bool hasBoundary = false ;
        for( int i=0; i<2*dim; ++i)
        {
          const std::size_t nbIndex = neighs[ i ];
          if( nbIndex != elIndex )
            values_[ i ] = u[ nbIndex ];
          else
            hasBoundary = true;
        }

        // if boundary present use intersection iterator to fill values
        if( hasBoundary )
        {
          const auto iend = gridPart().iend( element );
          for( auto iit = gridPart().ibegin( element ); iit != iend; ++iit )
          {
            const auto intersection = *iit;

            if ( intersection.boundary() )
            {
              const GlobalCoordinate iCenter = intersection.geometry().center();
              const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
              const StateVector uBnd = boundaryValue_( intersection, iCenter, iNormal, enVal );

              // TODO: Obtain boundary value
              values_[ intersection.indexInInside() ] = uBnd;
            }
          }
        }

        std::array< RangeType, 3 > vals_;
        std::array< RangeType, 3 > diffs_;

        // get element value (0 is the cell we are looking at)
        vals_[ 0 ] = enVal;

        std::array< RangeType, 3 > slopes;

        std::array< Field, 3 > dx_1;
        std::array< Field, 2 > dx;

        du = 0;

        // dimensional splitting
        for( int d = 0; d < dim; ++ d )
        {
          const Field h = h_[ d ];

          // 1/dx for computing slopes
          dx_1[ 0 ] = 1./h;
          dx_1[ 1 ] = 1./-h;
          dx_1[ 2 ] = 1./(-2.0 * h);

          // set u left (1) and u right (2)
          for( int i=0; i<2; ++i )
            vals_[ i+1 ] = values_[ 2*d + i ];

          // compute linear functions
          for( int i=0; i<numFunc; ++i )
          {
            for( int r=0; r<dimRange; ++r )
            {
              diffs_[ i ][ r ] = vals_[ combos_[ i ][0] ][ r ]  -  vals_[ combos_[ i ][1] ][ r ];
              slopes[ i ][ r ] = diffs_[ i ][ r ];
            }
            slopes[ i ] *= dx_1[ i ];
          }

          // recompute barycenter difference for limiting
          dx[ 0 ] =  h; // ( w_E,2 - w_E )
          dx[ 1 ] = -h; // ( w_E,1 - w_E )

          // limit slope
          for( int i=0; i<numFunc; ++i )
          {
            for(int r=0; r<dimRange; ++r)
            {
              Field minimalFactor = 1;

              for( const int c : testset_[ i ] )
              {
                // evaluate values for limiter function
                const Field d = diffs_[ c ][ r ];
                const Field g = slopes[ i ][ r ] * dx[  c ];

                // if the gradient in direction of the line
                // connecting the barycenters is very small
                // then neglect this direction since it does not give
                // valuable contribution to the linear function
                // call limiter function
                // g = grad L ( w_E,i - w_E ) ,  d = u_E,i - u_E
                Field localFactor = limiterFunction_( g, d );

                // take minimum
                minimalFactor = std::min( localFactor , minimalFactor );
                // if minimum is already zero stop computation here
                if( minimalFactor < 1e-12 )
                  break ;
              }
              //std::cout << minimalFactor << " min fac" << std::endl;

              // scale linear function
              slopes[ i ][ r ] *= minimalFactor;
            }
          }

          // select max slope
          RangeType& slope = slopes[0];
          for( int r=0; r<dimRange; ++r )
          {
            for( int i=1; i<numFunc; ++i )
              slope[ r ] = std::max( slope[r], slopes[ i ][ r ] );

            du[ r ][ d ] = slope[ r ];
          }
        }

        //std::cout << du << " computed gradient " << std::endl;
#if 0
        const auto iend = gridPart().iend( element );
        for( auto iit = gridPart().ibegin( element ); iit != iend; ++iit )
        {
          const auto intersection = *iit;

          std::cout << "Intersec : " << intersection.indexInInside() << std::endl;

          if( intersection.boundary() )
          {
            const GlobalCoordinate iCenter = intersection.geometry().center();
            const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
            const StateVector uBnd = boundaryValue_( intersection, iCenter, iNormal, u[ elIndex ] );
            differences.emplace_back( iCenter - elCenter, uBnd - u[ elIndex ] );
          }
          else if( intersection.neighbor() )
          {
            if constexpr ( isCartesian )
            {
              differences.emplace_back( centerDiff_[ intersection.indexInInside() ], u[ neighbors_[ elIndex ][ intersection.indexInInside() ] ] - u[ elIndex ] );
            }
            else
            {
              const auto neighbor = intersection.outside();
              const GlobalCoordinate nbCenter = neighbor.geometry().center();
              differences.emplace_back( nbCenter - elCenter, u[ mapper.index( neighbor ) ] - u[ elIndex ] );
            }
          }
        }
#endif
      }

      template< class Mapper, class Vector >
      void operator () ( const Mapper &mapper, const Vector &u, std::vector< Jacobian > &du ) const
      {
        du.resize( u.size() );

        const auto end = gridPart().template end< 0, Dune::InteriorBorder_Partition >();
        for( auto it = gridPart().template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it )
        {
          const auto element = *it;
          applyLocal( element, mapper, u, du[ mapper.index( element ) ] );
        }

        //auto handle = vectorCommDataHandle( mapper, du, [] ( Jacobian a, Jacobian b ) { return b; } );
        //gridPart().communicate( handle, InteriorBorder_All_Interface, ForwardCommunication );
      }

      const GridPartType &gridPart () const { return gridPart_; }

      [[deprecated]]
      const GridPartType &gridView () const { return gridPart_; }

    private:
      const GridPartType& gridPart_;
      BoundaryValue boundaryValue_;
      std::array< Field, dimension > h_;

      std::array< std::array< int8_t, 2>, numFunc > combos_;
      std::array< std::vector< int8_t >, numFunc > testset_;

      std::array< GlobalCoordinate, GridPartType::dimension*2 > centerDiff_;
      std::vector< std::array< int, GridPartType::dimension*2 > > neighbors_;

      //Dune::Fem::VanLeerLimiter< Field > limiterFunction_;
      //Dune::Fem::SuperBeeLimiter< Field > limiterFunction_;
      Dune::Fem::MinModLimiter< Field > limiterFunction_;
    };



    // lpReconstruction
    // ----------------

    template< class SV, class GP, class BV >
    inline static TVDReconstruction< GP, SV, BV > lpReconstruction ( const GP &gridPart, BV boundaryValue, typename FieldTraits< SV >::real_type tolerance )
    {
      return TVDReconstruction< GP, SV, BV >( gridPart, std::move( boundaryValue ), std::move( tolerance ) );
    }

  } // namespace FV

} // namespace Dune

#endif // #ifndef DUNE_FV_....
