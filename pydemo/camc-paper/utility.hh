#include <cstddef>
#include <limits>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/linesegmentsampler.hh>

#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/bindguard.hh>

#include <dune/fem/io/parameter.hh>

template <class DiscreteFunction>
int minMax( const DiscreteFunction& solution )
{
  typedef typename DiscreteFunction :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType ::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType ::GridPartType GridPartType;
  const DiscreteFunctionSpaceType& space = solution.space();
  const int dimRange = DiscreteFunctionSpaceType :: dimRange;

  typedef Dune::Fem::Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >
    DefaultQuadrature;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0, DefaultQuadrature::template DefaultQuadratureTraits > Quadrature;

  RangeType minVal( 1 );
  RangeType maxVal( -1 );

  bool isNan = false;

  Dune::Fem::ConstLocalFunction< DiscreteFunction > lf( solution );

  std::vector< RangeType > values;
  for( const auto& element : space )
  {
    Quadrature quad( element, 2*space.order( element )+3 );
    auto guard = Dune::Fem::bindGuard( lf, element );
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

  if( Dune::Fem::Parameter::verbose() )
  {
    std::cout << "Min/max values: min = " << minVal[0] << "  max = " << maxVal[0] << std::endl;
  }

  if( minVal[ 0 ] < 0 )
    return -1;
  if( isNan )
    return std::numeric_limits< int >::max();
  return 0;
}
