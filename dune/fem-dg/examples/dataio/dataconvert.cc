// this programm only works without MPI and without threads
#undef ENABLE_PETSC
#undef ENABLE_MPI
#undef USE_PTHREADS
#undef _OPENMP

//************************************************************
//
//  (C) written and directed by Robert Kloefkorn
//
//************************************************************
#if defined YASPGRID
//&& HAVE_MPI == 1

#ifndef ENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID
#error "Put -DENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID to CPPFLAGS for simul and disp program"
#endif
#endif
#include <config.h>

#include <string>

#include <dune/fem-dg/misc/streams.hh>

///////////////////////////////////////////////////
//
// Include your header defining all necessary types
//
///////////////////////////////////////////////////
#include <dune/fem-dg/main/main_pol.cc>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/fem/operator/projection/vtxprojection.hh>

void appendUserParameter()
{
  Dune::Fem::Parameter :: append("parameter");
}
#define PARAMETER_APPEND_FUNCTION appendUserParameter

typedef Dune::GridSelector :: GridType GridType;
typedef AlgorithmCreator< GridType > ProblemTraits;

typedef CheckPointingAlgorithm<GridType,
                             ProblemTraits,
                             POLORDER> AlgorithmType ;

typedef AlgorithmType :: IOTupleType InTupleType ;

// type of discrete function tuple to restore
typedef InTupleType GR_InputType;


template < int dimRange, int probDimRange >
struct AdditionalVariables
{
  template <class DestinationType>
  static DestinationType* setup( const double, const DestinationType& )
  {
    return 0;
  }
};

template < int dimRange >
struct AdditionalVariables< dimRange, dimRange >
{
  template <class DestinationType>
  static DestinationType* setup( const double time, const DestinationType& Uh )
  {
    DestinationType* additionalVariables = 0;
    if( ParameterType :: getValue< bool >("femdg.additionalvariables", false) )
    {
      /*
      std::cout << "Setup additional variables" << std::endl;
      additionalVariables = new DestinationType("add", Uh.space() );
      typedef typename AlgorithmType::ProblemType ProblemType;
      ProblemType* problem = ProblemTraits :: problem() ;
      typedef typename AlgorithmType :: ModelType  ModelType;
      ModelType model( *problem );

      // create TimeProvider provider with given time
      Dune::Fem::TimeProviderBase tp( time );

      setupAdditionalVariables( tp, Uh, model, *additionalVariables );
      delete problem;
      */
    }
    return additionalVariables;
  }
};

template < int pos >
struct ProcessElement
{
template <class GR_GridType,
          class InTupleType>
  static void apply(const GR_GridType& grid,
                    const InTupleType& data,
                    const double time,
                    const int timestep,
                    const int myRank,
                    const int numProcs)
  {
    applyDF( grid, *(std::get< pos >( data )), time, timestep, myRank, numProcs );
  }

  template <class GR_GridType, class DestinationType>
  static void applyDF(const GR_GridType& grid,
                      const DestinationType& Uh,
                      const double time,
                      const int timestep,
                      const int myRank,
                      const int numProcs)
  {
    typedef typename DestinationType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename AlgorithmType :: ProblemTypeType ProblemType;

    DestinationType* additionalVariables =
      AdditionalVariables< DiscreteFunctionSpaceType::dimRange,ProblemType::dimRange >::setup( time, Uh );

    typedef Dune::Fem::AdaptiveLeafGridPart< typename GridPartType::GridType > VtxGridPartType;
    typedef Dune::Fem::LagrangeDiscreteFunctionSpace< typename DiscreteFunctionSpaceType::FunctionSpaceType,
                                           VtxGridPartType, 1 > VtxSpaceType;
                                          // DiscreteFunctionSpaceType::polynomialOrder >  VtxSpaceType;
    typedef Dune::Fem::AdaptiveDiscreteFunction< VtxSpaceType > VtxFunctionType;


    typedef Dune::Fem::Parameter  ParameterType ;

    const int subSamplingLevel = ParameterType :: getValue< int >("fem.io.subsamplinglevel", 0);

    // subsampling vtk output
    //Dune::Fem::VTKIO< GridPartType > vtkio( Uh.space().gridPart() );

    Dune::Fem::DataOutputParameters parameter;
    if( parameter.outputformat() == 1 ) // vtk-vertex
    {
      if( ParameterType :: verbose() )
        std::cout << "Writing vertex data" << std::endl;

      VtxGridPartType vxGridPart( Uh.space().gridPart().grid() );
      VtxSpaceType xvSpace( vxGridPart );
      Dune::Fem::SubsamplingVTKIO< VtxGridPartType > vtkio( vxGridPart, subSamplingLevel );

      VtxFunctionType vtxValues( "vx", xvSpace );
      VtxFunctionType* addValues = ( additionalVariables ) ? new VtxFunctionType( "add", xvSpace ) : 0 ;

      Dune::Fem::VtxProjection< DestinationType, VtxFunctionType > vxpro;
      vxpro( Uh, vtxValues );

      if( additionalVariables )
      {
        vxpro( *additionalVariables, *addValues );
      }

      vtkio.addVertexData( vtxValues );
      if( addValues )
        vtkio.addVertexData( *addValues );

      const bool verbose = ParameterType::verbose() ;
      // get data name
      std::string name( (Uh.name() == "") ? "sol" : Uh.name() );
      // get file name
      std::string filename = Dune :: Fem :: generateFilename(name,timestep,6);

      if( verbose )
        std::cout <<"Writing vtk output " <<filename <<"...";

      // write vtk output
      vtkio.write( filename,
                   Dune::VTK::appendedraw,
                   myRank, numProcs
                 );
      if( verbose )
        std::cout <<"[ok]" << std::endl;

      delete addValues;
      if( additionalVariables )
        delete additionalVariables ;
    }
    else // vtk-cell
    {
      Dune::Fem::SubsamplingVTKIO< GridPartType > vtkio( Uh.space().gridPart(), subSamplingLevel );
      if( Uh.space().order() > 0 )
      {
        vtkio.addVertexData(Uh);
        if( additionalVariables )
          vtkio.addVertexData( *additionalVariables );
      }
      else
      {
        vtkio.addCellData(Uh);
        if( additionalVariables )
          vtkio.addCellData( *additionalVariables );
      }

      const bool verbose = ParameterType::verbose() ;
      // get data name
      std::string name( (Uh.name() == "") ? "sol" : Uh.name() );
      // get file name
      std::string filename = Dune :: Fem :: generateFilename(name,timestep,6);

      if( verbose )
        std::cout <<"Writing vtk output " <<filename <<"...";

      // write vtk output
      vtkio.write( filename,
                   Dune::VTK::appendedraw,
                   myRank, numProcs
                 );
      if( verbose )
        std::cout <<"[ok]" << std::endl;

      if( additionalVariables )
        delete additionalVariables ;
    }
  }
};

#define WRITE_VTK
template <class GR_GridType,
          class InTupleType>
void process(const GR_GridType& grid,
             const InTupleType& data,
             const double time,
             const int timestep,
             const int myRank,
             const int numProcs)
{
#ifdef WRITE_VTK
  Dune :: ForLoop< ProcessElement, 0, 0 > // Dune::tuple_size< InTupleType > :: value-1 >
    :: apply( grid, data, time, timestep, myRank, numProcs );
#endif
}

// include main program
#include <dune/fem/io/visual/grape/datadisp/dataconvert.cc>
