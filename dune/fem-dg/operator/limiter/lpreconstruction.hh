#ifndef DUNE_FV_LPRECONSTRUCTION_HH
#define DUNE_FV_LPRECONSTRUCTION_HH

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

#include <dune/optim/activeindexmapper.hh>
#include <dune/optim/common/axisalignedreferencefaces.hh>
#include <dune/optim/common/smallobject.hh>
#include <dune/optim/common/getreference.hh>
#include <dune/optim/constraint/linear.hh>
#include <dune/optim/lp.hh>
#include <dune/optim/solver/gaussjordan.hh>

namespace Dune
{

  namespace FV
  {

    // LPReconstruction
    // ----------------

    /**
     * \class LPReconstruction
     * \brief Minmod-type reconstruction based on linear programming problem
     *
     * The LPReconstruction was suggested in the following paper:
     * \code
     * @article{Chen:IntegratedLinearReconstruction,
     *   author  = {Chen, L. and Li, R.},
     *   title   = {An Integrated Linear Reconstruction for Finite Volume scheme
     *              on Unstructured Grids},
     *   journal = {J. Sci. Comput.},
     *   year    = {2016},
     *   volume  = {68},
     *   pages   = {1172--1197},
     *   doi     = {10.1007/s10915-016-0173-1}
     * }
     * \endcode
     **/
    template< class GV, class SV, class BV >
    class LPReconstruction
    {
      typedef LPReconstruction< GV, SV, BV > This;

    public:
      typedef GV GridView;
      typedef SV StateVector;
      typedef BV BoundaryValue;

      typedef FieldVector< typename GridView::ctype, GridView::dimensionworld > GlobalCoordinate;

      typedef typename GridView::Intersection Intersection;

      typedef typename FieldTraits< StateVector >::field_type Field;
      typedef typename FieldTraits< StateVector >::real_type Real;
      typedef FieldMatrix< Field, StateVector::dimension, GlobalCoordinate::dimension > Jacobian;

      static const int dimension = GridView::dimension;

    private:
      typedef Optim::LinearConstraint< GlobalCoordinate > Constraint;
      typedef std::vector< Constraint > Constraints;

      typedef Optim::GaussJordanSolver< FieldMatrix< Field, GlobalCoordinate::dimension, GlobalCoordinate::dimension > > Solver;
      typedef Optim::LinearProgramming< Solver, false > LP;

      typedef std::vector< std::pair< GlobalCoordinate, StateVector > > DifferencesVectorType;
    public:
      LPReconstruction ( const GridView &gridView, BoundaryValue boundaryValue, Real tolerance )
        : gridView_( gridView ),
          boundaryValue_( std::move( boundaryValue ) ),
          tolerance_( std::move( tolerance ) ),
          lp_( tolerance_ ),
          faceAxes_( LocalGeometryTypeIndex::size( dimension ) )
      {
        const unsigned int numTopo = Impl::numTopologies( dimension );
        std::unique_ptr< unsigned int[] > faceIndices( new unsigned int[ dimension * numTopo ] );
        std::unique_ptr< unsigned int[] > numFaces( new unsigned int[ numTopo ] );
        axisAlignedReferenceFaces( dimension, faceIndices.get(), numFaces.get() );
        for( unsigned int topologyId = 0; topologyId < numTopo; ++topologyId )
        {
          const GeometryType type( topologyId, dimension );
          std::vector< unsigned int > &faceAxes = faceAxes_[ LocalGeometryTypeIndex::index( type ) ];
          faceAxes.resize( numFaces[ topologyId ], dimension );
          for( int i = 0; i < dimension; ++i )
            faceAxes[ faceIndices[ topologyId*dimension + i ] ] = i;
        }
      }

      template< class Entity, class Mapper, class Vector >
      void applyLocal ( const Entity& element,
                        const Mapper &mapper,
                        const Vector &u,
                        Jacobian& du ) const
      {
        DifferencesVectorType& differences = differences_;
        differences.clear();

        Constraints constraints;

        const std::size_t elIndex = mapper.index( element );
        const GlobalCoordinate elCenter = element.geometry().center();

        std::array< unsigned int, dimension+1 > select;
        const std::vector< unsigned int > &faceAxes = faceAxes_[ LocalGeometryTypeIndex::index( element.type() ) ];
        if( !faceAxes.empty() )
        {
          const auto iend = gridView().iend( element );
          for( auto iit = gridView().ibegin( element ); iit != iend; ++iit )
          {
            const auto intersection = *iit;

            select[ faceAxes[ intersection.indexInInside() ] ] = differences.size();

            if( intersection.boundary() )
            {
              const GlobalCoordinate iCenter = intersection.geometry().center();
              const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
              const StateVector uBnd = boundaryValue_( intersection, iCenter, iNormal, u[ elIndex ] );
              differences.emplace_back( iCenter - elCenter, uBnd - u[ elIndex ] );
            }
            else if( intersection.neighbor() )
            {
              const auto neighbor = intersection.outside();
              const GlobalCoordinate nbCenter = neighbor.geometry().center();
              differences.emplace_back( nbCenter - elCenter, u[ mapper.index( neighbor ) ] - u[ elIndex ] );
            }
          }
        }
        else
        {
          Dune::ReservedVector< GlobalCoordinate, dimension > onb;

          const auto iend = gridView().iend( element );
          for( auto iit = gridView().ibegin( element ); iit != iend; ++iit )
          {
            const auto intersection = *iit;

            select[ onb.size() ] = differences.size();

            if( intersection.boundary() )
            {
              const GlobalCoordinate iCenter = intersection.geometry().center();
              const GlobalCoordinate iNormal = intersection.centerUnitOuterNormal();
              const StateVector uBnd = boundaryValue_( intersection, iCenter, iNormal, u[ elIndex ] );
              differences.emplace_back( iCenter - elCenter, uBnd - u[ elIndex ] );
            }
            else if( intersection.neighbor() )
            {
              const auto neighbor = intersection.outside();
              const GlobalCoordinate nbCenter = neighbor.geometry().center();
              differences.emplace_back( nbCenter - elCenter, u[ mapper.index( neighbor ) ] - u[ elIndex ] );
            }

            if( onb.size() < dimension )
            {
              GlobalCoordinate dx = differences.back().first;
              for( const GlobalCoordinate &v : onb )
                dx.axpy( -(dx*v), v );

              const Real dxNorm = dx.two_norm();
              if( dxNorm >= tolerance_ )
                onb.push_back( dx /= dxNorm );
            }
          }
        }

        // reserve memory for constraints
        const std::size_t numConstraints = differences.size();
        constraints.resize( 2u*numConstraints );
        Optim::ActiveIndexMapper< SmallObjectAllocator< unsigned int > > active( GlobalCoordinate::dimension, constraints.size() );
        for( int j = 0; j < StateVector::dimension; ++j )
        {
          GlobalCoordinate negGradient( 0 );
          for( std::size_t i = 0u; i < numConstraints; ++i )
          {
            const Field sign = (differences[ i ].second[ j ] >= 0 ? 1 : -1);

            negGradient.axpy( sign, differences[ i ].first );

            constraints[ 2*i ].normal() = differences[ i ].first;
            constraints[ 2*i ].normal() *= sign;
            constraints[ 2*i ].rhs() = sign*differences[ i ].second[ j ];

            constraints[ 2*i+1 ].normal() = constraints[ 2*i ].normal();
            constraints[ 2*i+1 ].normal() *= Field( -1 );
            constraints[ 2*i+1 ].rhs() = 0;
          }

          // activate GlobalCoordinate::dimension constraints active in the origin
          active.clear();
          for( int i = 0; i < dimension; ++i )
            active.activate( 2*select[ i ]+1 );

          // solve
          du[ j ] = 0;
          lp_( negGradient, constraints, du[ j ], active );
        }
      }

      template< class Mapper, class Vector >
      void operator () ( const Mapper &mapper, const Vector &u, std::vector< Jacobian > &du ) const
      {
        du.resize( u.size() );

        const auto end = gridView().template end< 0, Dune::InteriorBorder_Partition >();
        for( auto it = gridView().template begin< 0, Dune::InteriorBorder_Partition>(); it != end; ++it )
        {
          const auto element = *it;
          applyLocal( element, mapper, u, du[ mapper.index( element ) ] );
        }

        //auto handle = vectorCommDataHandle( mapper, du, [] ( Jacobian a, Jacobian b ) { return b; } );
        //gridView().communicate( handle, InteriorBorder_All_Interface, ForwardCommunication );
      }

      const GridView &gridView () const { return gridView_; }

    private:
      GridView gridView_;
      BoundaryValue boundaryValue_;
      Real tolerance_;
      LP lp_;
      std::vector< std::vector< unsigned int > > faceAxes_;
      mutable DifferencesVectorType differences_;
    };



    // lpReconstruction
    // ----------------

    template< class SV, class GV, class BV >
    inline static LPReconstruction< GV, SV, BV > lpReconstruction ( const GV &gridView, BV boundaryValue, typename FieldTraits< SV >::real_type tolerance )
    {
      return LPReconstruction< GV, SV, BV >( gridView, std::move( boundaryValue ), std::move( tolerance ) );
    }

  } // namespace FV

} // namespace Dune

#endif // #ifndef DUNE_FV_LPRECONSTRUCTION_HH
