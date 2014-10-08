#undef ENABLE_MPI

#ifndef FEMHOWTO_POISSONSTEPPER_HH
#define FEMHOWTO_POISSONSTEPPER_HH
#include <config.h>

#include <dune/fem/main/codegen.hh>
#include "passtraits.hh"

// dune-grid includes
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// dune-fem includes
#include <dune/fem/io/parameter.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>


#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include "../common/inverseop.hh"

// local includes
#include "benchmarkproblem.hh"
#include "models.hh"
#include "estimator1.hh"

#define NS_ELLIPTIC_OPERATOR 
#ifndef TESTOPERATOR
#define TESTOPERATOR 
#endif

template <class GridType> 
struct ProblemGenerator 
{
  typedef Dune :: ProblemInterface<
             Dune::Fem::FunctionSpace< double, double, GridType::dimension, DIMRANGE> >  ProblemType;

  template <class GridPart>
  struct Traits
  {
    typedef ProblemType  InitialDataType;
    typedef PoissonModel< GridPart, InitialDataType > ModelType;
    typedef Dune::NoFlux< ModelType > FluxType;
    // choice of diffusion flux (see diffusionflux.hh for methods)
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId 
      = Dune :: method_general ;
  };


  static inline std::string diffusionFluxName()
  {
#if (defined PRIMALDG)
    return Dune::Fem::Parameter::getValue< std::string >("dgdiffusionflux.method");
#else
    return "LDG";
#endif
  }

  static inline Dune::GridPtr<GridType> initializeGrid(const std::string description)
  {
    // use default implementation 
    //GridPtr<GridType> gridptr = initialize< GridType >( description );
    std::string filename = Dune::Fem::Parameter::getValue< std::string >(Dune::Fem::IOInterface::defaultGridKey(GridType :: dimension, false)); 
  
    // initialize grid
    Dune::GridPtr< GridType > gridptr = initialize< GridType >( description );

    std::string::size_type position = std::string::npos;
    if( GridType::dimension == 2) 
      position = filename.find("mesh3");
    if( GridType::dimension == 3 )
      position = filename.find("_locraf.msh" );

    std::string::size_type locraf = filename.find("locrafgrid_");

    std::string::size_type checkerboard = filename.find("heckerboard" );;

    GridType& grid = *gridptr ;

    const int refineelement = 1 ;

    bool nonConformOrigin = Dune::Fem::Parameter::getValue< bool > ( "poisson.nonConformOrigin",false );

    if ( nonConformOrigin )
    {
      std::cout << "Create local refined grid" << std::endl;
      if( grid.comm().rank() == 0)
      {
        typedef typename GridType:: template Codim<0>::LeafIterator IteratorType;
        IteratorType endit = grid.template leafend<0>();
        for(IteratorType it = grid.template leafbegin<0>(); it != endit ; ++it)
        {
          const typename IteratorType :: Entity & entity = *it ;
          const typename IteratorType :: Entity :: Geometry& geo = entity.geometry();
          if (geo.center().two_norm() < 0.5)
            grid.mark(refineelement, entity );
        }
      }
    }
    else if ( position != std::string::npos || 
              locraf   != std::string::npos ) 
    {
      std::cout << "Create local refined grid" << std::endl;
      if( grid.comm().rank() == 0)
      {
        typedef typename GridType:: template Codim<0>::LeafIterator IteratorType;
        Dune::FieldVector<double,GridType::dimension> point(1.0);
        if( locraf != std::string::npos ) point[ 0 ] = 3.0;
        const int loops = ( GridType::dimension == 2 ) ? 2 : 1;
        for (int i=0; i < loops; ++i)
        {
          point /= 2.0 ;

          IteratorType endit = grid.template leafend<0>();
          for(IteratorType it = grid.template leafbegin<0>(); it != endit ; ++it)
          {
            const typename IteratorType :: Entity & entity = *it ;
            const typename IteratorType :: Entity :: Geometry& geo = entity.geometry();
            bool inside = true ; 
            const int corners = geo.corners();
            for( int c = 0; c<corners; ++ c) 
            {
              for( int i = 0; i<GridType::dimension; ++i )
              {
                if( geo.corner( c )[ i ] - point[ i ] > 1e-8 ) 
                {
                  inside = false ; 
                  break ;
                }
              }
            }
            if( inside ) grid.mark(refineelement, entity );
          }
        }
      }
    }
    else if ( checkerboard != std::string::npos ) 
    {
      std::cout << "Create Checkerboard mesh" << std::endl;
      int moduloNum = 32 ;
      for( int i = 2; i<32; i *= 2 ) 
      {
        std::stringstream key; 
        key << i << "x" << i << "x" << i;
        std::string::size_type counter = filename.find( key.str() );
        if( counter != std::string::npos ) 
        {
          moduloNum = i;
          break ;
        }
      }

      std::cout << "Found checkerboard cell number " << moduloNum << std::endl;

      if( grid.comm().rank() == 0)
      {
        typedef typename GridType:: template Codim<0>::LeafIterator IteratorType;
        IteratorType endit = grid.template leafend<0>();
        int count = 1;
        int modul = 1; // start with 1, such that the fine cube is at (0,0,0)
        for(IteratorType it = grid.template leafbegin<0>(); it != endit ; ++it, ++ count )
        {
          const typename IteratorType :: Entity & entity = *it ;
          if( count % 2 == modul ) 
            grid.mark(refineelement, entity );

          if( count % moduloNum == 0 ) 
            modul = 1 - modul ;

          if( count % (moduloNum * moduloNum) == 0 )
            modul = 1 - modul ;
        }
      }
    }

    // adapt the grid 
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    //Dune :: GrapeGridDisplay< GridType > grape( grid );
    //grape.display ();

    return gridptr ;
  }

  static ProblemType* problem()
  {
    // choice of benchmark problem 
    int probNr = Dune::Fem::Parameter::getValue< int > ( "femhowto.problem" );
    return new Dune :: BenchmarkProblems< GridType > ( probNr );
  }
};

#endif // FEMHOWTO_POISSONSTEPPER_HH
