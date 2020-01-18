#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/io/file/dataoutput.hh>

template <class DiscreteFunction>
void vtkout(const DiscreteFunction& u, int step )
{
  typedef std::tuple< const DiscreteFunction* > IOTuple;

  IOTuple tup( &u );

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType :: GridType GridType;
  typedef Dune::Fem::DataOutput< GridType, IOTuple > DataOutputType;

  const GridType& grid = u.space().grid();

  Dune::Fem::DataOutputParameters param;
  DUNE_EXPORT static DataOutputType dataOutput( grid, tup, param );

  dataOutput.writeData( 1.0 );
}
