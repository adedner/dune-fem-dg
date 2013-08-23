#undef ENABLE_MPI
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

#warning "Switching from YASPGRID to SGRID for parallel display!"
#undef YASPGRID 
#define SGRID 
#endif
#include <config.h>

#include <string>

///////////////////////////////////////////////////
//
// Include your header defining all necessary types 
//
///////////////////////////////////////////////////
#include "checkpointing.hh"

#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>


void appendUserParameter() 
{
  Dune::Fem::Parameter :: append("parameter");
}

typedef Dune::GridSelector :: GridType GridType;

// ProblemType is a Dune::Function that evaluates to $u_0$ and also has a
// method that gives you the exact solution.
typedef Dune::U0< GridType > ProblemType;

typedef  Stepper<GridType, ProblemType> StepperType ;

typedef StepperType :: IOTupleType InTupleType ;

// type of discrete function tuple to restore 
typedef InTupleType GR_InputType;

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
  typedef typename Dune :: TypeTraits< typename Dune :: tuple_element< 0, InTupleType >::type >::PointeeType DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  const DestinationType& Uh = *( Dune :: get<0>(data) );

  const int subSamplingLevel = Dune::Fem::Parameter :: 
        getValue< int >("fem.io.subsamplinglevel", 1);
  //std::cout <<"Got Sublevel = " << subSamplingLevel << std::endl;
  // vtk output 
  Dune::Fem::SubsamplingVTKIO< GridPartType > vtkio( Uh.space().gridPart(), subSamplingLevel );
  //VTKIO< GridPartType > vtkio(Uh.space().gridPart(),
  //                           (Uh.space().continuous()) ? 
  //                             VTKOptions::conforming :
  //                            VTKOptions::nonconforming);

  if( Uh.space().order() > 0 )
  {
    vtkio.addVertexData(Uh);
  }
  else 
  {
    vtkio.addCellData(Uh);
  }
  
  // get data name 
  std::string name( (Uh.name() == "") ? "rho" : Uh.name() ); 
  //std::cout << "Got discrete function " << name << std::endl;
  //std::string name( "rho" ); 
  // get file name
  std::string filename = Dune :: genFilename("",name,timestep,6);
  std::cout <<"Writing vtk output " << filename << " for time = " << time << " ...";
  // write vtk output 
  vtkio.write( filename, 
               Dune::VTK::appendedraw,
               myRank, numProcs
             );
  std::cout <<"[ok]" << std::endl;

#endif
}

// include main program 
#include <dune/fem/io/visual/grape/datadisp/dataconvert.cc>
