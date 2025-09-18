#ifndef DUNE_FV_CARTESIANNEIGHBORS_HH
#define DUNE_FV_CARTESIANNEIGHBORS_HH

#include <cassert>
#include <cstddef>

#include <numeric>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/fem/gridpart/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    // CartesianNeighbors
    // ------------------

    /**
     * \class CartesianNeighbors
     * \brief Neighbor indices for Cartesian grids
     * \endcode
     **/
    template< class GP >
    class CartesianNeighbors
    {
      typedef CartesianNeighbors< GP > This;

    public:
      typedef GP GridPartType;

      typedef Dune::Fem::DofManager< typename GridPartType::GridType > DofManagerType;

      typedef typename GridPartType::ctype Field;
      typedef FieldVector< Field, GridPartType::dimensionworld > GlobalCoordinate;
      typedef FieldVector< Field, GridPartType::dimensionworld > LocalCoordinate;

      typedef typename GridPartType::Intersection Intersection;

      static const int dimension = GridPartType::dimension;
      static const bool isCartesian = Dune::Fem::GridPartCapabilities::isCartesian< GridPartType >::v;
      // static_assert( isCartesian, "CartesianNeighbors requires a Cartesian grid");

      typedef typename GridPartType::IndexSetType::IndexType IndexType;
      typedef std::vector< std::array< IndexType, dimension*2 > >  NeighborIndexVectorType;
      typedef std::array< GlobalCoordinate, dimension*2 > CenterDifferencesType;

    public:
      CartesianNeighbors ( const GridPartType &gp )
        : gridPart_( gp ),
          dofManager_( DofManagerType :: instance( gridPart_.grid() ) ),
          h_( 0 ),
          sequence_( -1 )
      {
        for( int i=0; i<2*dimension; ++i )
          centerDiff_[ i ] = 0;

        // update neighbors
        update();
      }

      void update()
      {
        // do nothing if up to date
        const int dmSequence = dofManager_.sequence();
        if( sequence_ == dmSequence )
          return ;

        const auto& mapper = gridPart().indexSet();
        neighbors_.resize( mapper.size(0) );

        h_ = 0;

        const auto end = gridPart().template end< 0, Dune::InteriorBorder_Partition >();
        for( auto it = gridPart().template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it )
        {
          const auto& element = *it ;
          const IndexType elIndex = mapper.index( element );
          auto& neighbors = neighbors_[ elIndex ];

          const GlobalCoordinate elCenter = element.geometry().center();

          const auto iend = gridPart().iend( element );
          for( auto iit = gridPart().ibegin( element ); iit != iend; ++iit )
          {
            const auto intersection = *iit;

            if( intersection.neighbor() )
            {
              const auto neighbor = intersection.outside();
              const IndexType nbIndex = mapper.index( neighbor );
              neighbors[ intersection.indexInInside() ] = nbIndex;

              // only real neighbors, no periodic boundaries
              if( ! intersection.boundary() )
              {
                const GlobalCoordinate nbCenter = neighbor.geometry().center();
                GlobalCoordinate centerDiff = nbCenter - elCenter;
                centerDiff_[ intersection.indexInInside() ] = centerDiff;

                const int d = intersection.indexInInside() / dimension;
                if( h_[ d ] > 0.0 )
                {
                  if( std::abs( h_[ d ] - std::abs(centerDiff[ d ]) )> 1e-10 )
                  {
                    DUNE_THROW(InvalidStateException, "TVDReconstruction::update: different dx detected in one direction");
                  }
                }
                else
                {
                  // compute h
                  h_[ d ] = std::abs(centerDiff[ d ]);
                }
              }
            }
            else if ( intersection.boundary() )
            {
              // for boundary faces set inside element index
              neighbors[ intersection.indexInInside() ] = elIndex;
            }
          }
        }

        // update sequence counter
        sequence_ = dmSequence;
      }

      const GridPartType &gridPart () const { return gridPart_; }

      const NeighborIndexVectorType& indices() const
      {
        // vector should be correct when grid is Cartesian, otherwise could be anything
        assert( isCartesian ? (neighbors_.size() == size_t(gridPart().indexSet().size(0))) : true );
        return neighbors_;
      }

      const LocalCoordinate& gridWidth() const { return h_; }
      const CenterDifferencesType& centerDifferences() const { return centerDiff_; }

    protected:
      const GridPartType& gridPart_;
      const DofManagerType& dofManager_;

      LocalCoordinate h_;
      CenterDifferencesType centerDiff_;
      NeighborIndexVectorType neighbors_;

      int sequence_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FV_....
