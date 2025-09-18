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
#include <dune/fem-dg/operator/common/cartesianneighbors.hh>

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

      typedef Dune::Fem::CartesianNeighbors< GridPartType >  CartesianNeighborsType;

      static const int dimension = GridPartType::dimension;
      static const bool isCartesian = Dune::Fem::GridPartCapabilities::isCartesian< GridPartType >::v;
      static_assert( isCartesian, "TVDReconstruction requires a Cartesian grid");

      static const int numFunc = 3;
    public:
      TVDReconstruction ( const GridPartType &gp, BoundaryValue boundaryValue, Real tolerance )
        : neighbors_( gp ),
          boundaryValue_( std::move( boundaryValue ) ),
          limiterFunction_()
      {
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

        // update grid width
        update();
      }

      void update()
      {
        neighbors_.update();
      }

      template< class Entity, class Mapper, class Vector >
      void applyLocal ( const Entity& element,
                        const Mapper &mapper,
                        const Vector &u,
                        Jacobian& du ) const
      {
        static const int dim = dimension;
        static const int dimRange = Jacobian::rows;

        const auto& elIndex = mapper.index( element );

        const auto& neighs = neighbors_.indices()[ elIndex ];
        const auto& h_ = neighbors_.gridWidth();

        const RangeType &enVal = u[ elIndex ];

        // neighboring values
        std::array< RangeType, 2*dim > values_;
        bool hasBoundary = false ;
        for( int i=0; i<2*dim; ++i)
        {
          const auto& nbIndex = neighs[ i ];
          if( nbIndex != elIndex )
            values_[ i ] = u[ nbIndex ];
          else
            hasBoundary = true;
        }

        // if boundary present use intersection iterator to fill values
        if( hasBoundary )
        {
          const auto& gridPart = neighbors_.gridPart();
          const auto iend = gridPart.iend( element );
          for( auto iit = gridPart.ibegin( element ); iit != iend; ++iit )
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
                const Field g = slopes[ i ][ r ] * dx[ c ];

                // if the gradient in direction of the line
                // connecting the barycenters is very small
                // then neglect this direction since it does not give
                // valuable contribution to the linear function
                // call limiter function
                // g = grad L ( w_E,i - w_E ) ,  d = u_E,i - u_E
                Field localFactor = limiterFunction_( g, d );

                /*
                if( localFactor < 1.0 )
                {
                  const Field limitEps = 1e-8;
                  const Field length2 = dx[ c ]*dx[ c ];

                  // if length is to small then the grid is corrupted
                  assert( length2 > 1e-14 );

                  const Field factor = (g*g) / length2 ;
                  if( factor < limitEps )
                  {
                    //std::cout << "Using tanh" << std::endl;
                    //localFactor = 1.0 - std::tanh( factor / limitEps );
                    localFactor = 1.0 - ( factor / limitEps );
                  }
                }
                */

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
      }

      template< class Mapper, class Vector >
      void operator () ( const Mapper &mapper, const Vector &u, std::vector< Jacobian > &du ) const
      {
        du.resize( u.size() );

        const auto& gridPart = neighbors_.gridPart();
        const auto end = gridPart.template end< 0, Dune::InteriorBorder_Partition >();
        for( auto it = gridPart.template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it )
        {
          const auto element = *it;
          applyLocal( element, mapper, u, du[ mapper.index( element ) ] );
        }
      }

    protected:
      CartesianNeighborsType neighbors_;
      BoundaryValue boundaryValue_;

      std::array< std::array< int8_t, 2>, numFunc > combos_;
      std::array< std::vector< int8_t >, numFunc > testset_;

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
