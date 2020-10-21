#include <cstddef>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#if HAVE_DUNE_FEM_DG
#include <dune/fem-dg/operator/dg/defaultquadrature.hh>
#endif

#include <dune/fem/misc/linesegmentsampler.hh>

template <class DiscreteFunction>
int minMax( const DiscreteFunction& solution )
{
  typedef typename DiscreteFunction :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType ::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType ::GridPartType GridPartType;
  const DiscreteFunctionSpaceType& space = solution.space();
  const int dimRange = DiscreteFunctionSpaceType :: dimRange;

  typedef Dune::Fem::DefaultQuadrature< DiscreteFunctionSpaceType >
    DefaultQuadrature;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0, DefaultQuadrature::template DefaultQuadratureTraits > Quadrature;

  RangeType minVal( 1 );
  RangeType maxVal( -1 );

  bool isNan = false;

  std::vector< RangeType > values;
  for( const auto& element : space )
  {
    Quadrature quad( element, 2*space.order( element )+3 );
    auto lf = solution.localFunction( element );
    const int nop = quad.nop();
    values.resize( nop );
    lf.evaluateQuadrature( quad, values );
    for( int i=0; i<nop; ++i )
    {
      RangeType& val = values[ i ];
      for( int d=0; d<dimRange; ++d )
      {
        minVal[d] = std::min( minVal[d], val[ d ] );
        maxVal[d] = std::max( maxVal[d], val[ d ] );
        isNan = !(val[d]==val[d]);
      }
    }
  }

  const auto& comm = space.gridPart().grid().comm();

  comm.min( &minVal[0], dimRange );
  comm.max( &maxVal[0], dimRange );

  if( comm.rank() == 0 )
  {
    std::cout << "Min/max values: min = " << minVal[0] << "  max = " << maxVal[0] << std::endl;
  }

  if( minVal[ 0 ] < 0 )
    return -1;
  if( isNan )
    return 1e308;
  return 0;
}
