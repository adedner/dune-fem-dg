#ifndef FEMDG_STOKESSTEPPER_HH
#define FEMDG_STOKESSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem/misc/l2norm.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/upwindflux.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/noflux.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/test/stokes/stokesalgorithm.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//--------- MODELS --------------------------
#include "models.hh"


enum AdvectionDiffusionOperatorType
{
  _unlimited = 0,
  _limited = 1,
};


template< class Op, class DiffusionOp, class AdvectionOp, bool advection, bool diffusion >
struct OperatorChooser
{
  typedef Op                   FullOperatorType;
  typedef DiffusionOp          ImplicitOperatorType;
  typedef AdvectionOp          ExplicitOperatorType;
};
template< class Op, class DiffusionOp, class AdvectionOp, bool advection >
struct OperatorChooser< Op, DiffusionOp, AdvectionOp, advection, false >
{
  typedef AdvectionOp          FullOperatorType;
  typedef FullOperatorType     ImplicitOperatorType;
  typedef FullOperatorType     ExplicitOperatorType;
};
template<class Op, class DiffusionOp, class AdvectionOp, bool diffusion >
struct OperatorChooser< Op, DiffusionOp, AdvectionOp, false, diffusion >
{
  typedef DiffusionOp          FullOperatorType;
  typedef FullOperatorType     ImplicitOperatorType;
  typedef FullOperatorType     ExplicitOperatorType;
};



template< class OperatorTraits, bool advection, bool diffusion, AdvectionDiffusionOperatorType op >
class AdvectionDiffusionOperators;


template< class OperatorTraits, bool advection, bool diffusion >
class AdvectionDiffusionOperators< OperatorTraits, advection, diffusion, _unlimited >
{
  typedef Dune::DGAdvectionDiffusionOperator< OperatorTraits >             DgType;
  typedef Dune::DGAdvectionOperator< OperatorTraits >                      DgAdvectionType;
  typedef Dune::DGDiffusionOperator< OperatorTraits >                      DgDiffusionType;
  typedef OperatorChooser< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                           OperatorChooserType;
public:
  typedef typename OperatorChooserType::FullOperatorType                            FullOperatorType;
  typedef typename OperatorChooserType::ImplicitOperatorType                        ImplicitOperatorType;
  typedef typename OperatorChooserType::ExplicitOperatorType                        ExplicitOperatorType;
};

template< class OperatorTraits, bool advection, bool diffusion >
class AdvectionDiffusionOperators< OperatorTraits, advection, diffusion, _limited >
{
  typedef Dune::DGLimitedAdvectionDiffusionOperator< OperatorTraits >      DgType;
  typedef Dune::DGLimitedAdvectionOperator< OperatorTraits >               DgAdvectionType;
  typedef Dune::DGDiffusionOperator< OperatorTraits >                      DgDiffusionType;
  typedef OperatorChooser< DgType, DgAdvectionType, DgDiffusionType, advection, diffusion >
                                                                           OperatorChooserType;
public:
  typedef typename OperatorChooserType::FullOperatorType                            FullOperatorType;
  typedef typename OperatorChooserType::ImplicitOperatorType                        ImplicitOperatorType;
  typedef typename OperatorChooserType::ExplicitOperatorType                        ExplicitOperatorType;
};


enum DiscreteFunctionSpacesType
{
  _lagrange = 0,
  _legendre = 1,
};


template< class FunctionSpaceImp, class GridPartImp, int polOrder, DiscreteFunctionSpacesType dfType, GalerkinType opType >
struct DiscreteFunctionSpaces;

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
struct DiscreteFunctionSpaces< FunctionSpaceImp, GridPartImp, polOrder, _lagrange, cg >
{
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, Dune::Fem::CachingStorage > type;
};

template< class FunctionSpaceImp, class GridPartImp, int polOrder >
struct DiscreteFunctionSpaces< FunctionSpaceImp, GridPartImp, polOrder, _legendre, dg >
{
  typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceImp, GridPartImp, polOrder, Dune::Fem::CachingStorage > type;
};

template< class DiscreteFunctionSpaceImp, SolverType solverType >
struct DiscreteFunctions;


template< class DiscreteFunctionSpaceImp >
struct DiscreteFunctions< DiscreteFunctionSpaceImp, fem >
{
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceImp > type;
  typedef Dune::Fem::SparseRowLinearOperator< type, type >                jacobian;
};

#if HAVE_DUNE_ISTL
template< class DiscreteFunctionSpaceImp >
struct DiscreteFunctions< DiscreteFunctionSpaceImp, istl >
{
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >  type;
  typedef Dune::Fem::ISTLLinearOperator< type, type >                             jacobian;
};
#endif

#if HAVE_PETSC
template< class DiscreteFunctionSpaceImp >
struct DiscreteFunctions< DiscreteFunctionSpaceImp, petsc >
{
  typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceImp > type;
  typedef Dune::Fem::PetscLinearOperator< type, type >                 jacobian;
};
#endif



template< class GridImp >
struct PoissonProblemCreator
{
  typedef GridImp                                         GridType;
  typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
  typedef HostGridPartType                                GridPartType;

  // define problem type here if interface should be avoided
  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE >
                                                                FunctionSpaceType;
  typedef Dune::ProblemInterface< FunctionSpaceType >           ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                ProblemType;
    typedef ProblemInterfaceType                                InitialDataType;
    typedef StokesModel< GridPart, InitialDataType >           ModelType;
    //typedef typename InitialDataType::TimeDependentFunctionType TimeDependentFunctionType;

    typedef std::vector< int >                                  EOCErrorIDs;

    static EOCErrorIDs initEoc()
    {
      EOCErrorIDs ids;
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-error" ) ) );
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "DG-error" ) ) );
      ids.push_back( Dune::Fem::FemEoc::addEntry( std::string( "sigma-norm" ) ) );
      return ids;
    }

    template< class SolutionImp, class Model, class ExactFunction, class SigmaEstimatorImp >
    static void addEOCErrors ( const EOCErrorIDs &ids, SolutionImp &u, const Model &model, const ExactFunction &f, SigmaEstimatorImp& sigma  )
    {
      // calculate L2 - Norm
      Dune::Fem::L2Norm< GridPartType > l2norm( u.space().gridPart() );
      const double l2error = l2norm.distance( f, u );

      Dune::Fem::DGNorm< GridPartType > dgnorm( u.space().gridPart() );
      const double dgerror = dgnorm.distance( f, u );

      Dune::Fem::H1Norm< GridPartType > sigmanorm( u.space().gridPart() );
      const double sigmaerror = sigmanorm.distance( f, sigma );

      Dune::Fem::FemEoc::setErrors( ids[ 0 ], l2error );
      Dune::Fem::FemEoc::setErrors( ids[ 1 ], dgerror );
      Dune::Fem::FemEoc::setErrors( ids[ 2 ], sigmaerror );
    }
  };

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
        // use default implementation
    std::string filename = Dune::Fem::Parameter::commonInputPath() + "/" +
               Dune::Fem::Parameter::getValue< std::string >(Dune::Fem::IOInterface::defaultGridKey(GridType :: dimension, false));

    // use default implementation
    Dune::GridPtr<GridType> gridptr = Dune::Fem::DefaultGridInitializer< GridType >::initialize();

    std::string::size_type position = std::string::npos;
    if( GridType::dimension == 2)
      position = filename.find("mesh3");
    if( GridType::dimension == 3 )
      position = filename.find("_locraf.msh" );

    std::string::size_type locraf = filename.find("locrafgrid_");

    std::string::size_type checkerboard = filename.find("heckerboard" );;

    GridType& grid = *gridptr ;

    const int refineelement = 1 ;

    bool nonConformOrigin = Dune::Fem::Parameter::getValue< bool > ( "nonConformOrigin",false );

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

    return gridptr ;
  }

  static ProblemInterfaceType* problem()
  {
    //int probNr = Dune::Fem::Parameter::getValue< int > ( "problem" );
    //return new Dune :: GeneralizedStokesProblem< GridType, DIMRANGE > ( probNr );
    return new Dune::GeneralizedStokesProblem< GridType > ();
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
public:
    typedef AnalyticalTraits< GridPartType >                              AnalyticalTraitsType;

    static const int polynomialOrder = polOrd;

    static inline std::string advectionFluxName()
    {
      return "LLF";
    }

    static inline std::string diffusionFluxName()
    {
      return Dune::Fem::Parameter::getValue< std::string >("dgdiffusionflux.method");
    }

    static const int quadOrder = polynomialOrder * 3 + 1;
    static const SolverType solverType = istl ;

    typedef typename DiscreteFunctionSpaces< FunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef std::tuple<> ExtraParameterTuple;

    typedef typename AnalyticalTraitsType::ProblemType::ExactSolutionType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;
    typedef std::tuple< DiscreteFunctionType*, GridExactSolutionType* > IOTupleType;

private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
public:
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    // type of restriction/prolongation projection for adaptive simulations
    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                           RestrictionProlongationType;
    // type of linear solver for implicit ode
    typedef Dune::AdaptationHandler< GridType, FunctionSpaceType >                              AdaptationHandlerType;


    //############################ poisson issues ############################
    static const bool symmetricSolver = true ;
    typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver> SolversType;

    //typedef Dune::UpwindFlux< typename AnalyticalTraitsType::ModelType >              FluxType;
    typedef Dune::NoFlux< typename AnalyticalTraitsType::ModelType >              FluxType;

    typedef typename SolversType::LinearOperatorType         FullOperatorType;
    typedef typename SolversType::LinearInverseOperatorType  BasicLinearSolverType;

    //----------- passes! ------------------------
    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType,  IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;
    typedef Dune::DGAdvectionDiffusionOperator< OperatorTraitsType >  AssemblyOperatorType;
    typedef Dune::DGPrimalMatrixAssembly< AssemblyOperatorType >            AssemblerType;
    //#########################################################################

    //------------- Limiter ---------------------------------------------
    typedef FullOperatorType                                                   LimiterOperatorType;
    //------------ Limiter ---------------------------------------------


  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::EllipticAlgorithm< GridType, PoissonProblemCreator<GridType>, polOrd > Type;
  };


};




template< class GridImp >
struct StokesProblemCreator
{
  typedef PoissonProblemCreator<GridImp>                           PoissonProblemCreatorType;

  typedef typename PoissonProblemCreatorType::GridType             GridType;
  typedef typename PoissonProblemCreatorType::HostGridPartType     HostGridPartType;
  typedef typename PoissonProblemCreatorType::GridPartType         GridPartType;


private:
  typedef typename PoissonProblemCreatorType::FunctionSpaceType  VelocityFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, 1 >
                                                                PressureFunctionSpaceType;
public:
  // define problem type here if interface should be avoided
  typedef typename PoissonProblemCreatorType::FunctionSpaceType FunctionSpaceType;


  typedef Dune::StokesProblemInterface< VelocityFunctionSpaceType, PressureFunctionSpaceType >       ProblemInterfaceType;

  template< class GridPart > // TODO: is this template parameter needed?
  struct AnalyticalTraits
  {
    typedef ProblemInterfaceType                                      ProblemType;
    typedef ProblemInterfaceType                                      InitialDataType;
    typedef StokesModel< GridPart, InitialDataType >                  ModelType;

    typedef typename PoissonProblemCreatorType::template AnalyticalTraits< GridPart >::EOCErrorIDs     EOCErrorIDs;

    static EOCErrorIDs initEoc ()
    {
      EOCErrorIDs ids;
      ids.push_back(Dune::Fem::FemEoc::addEntry( std::string( "$L^2$-p-error" ) ) );
      return ids;
    }

    template< class SolutionImp, class Model, class ExactFunction >
    static void addEOCErrors ( const EOCErrorIDs &ids, SolutionImp &p, const Model &model, const ExactFunction &g  )
    {
      // calculate L2 - p-Norm
      Dune::Fem::L2Norm< GridPartType > l2pnorm( p.space().gridPart() );
      const double l2perror = l2pnorm.distance( g, p );
      Dune::Fem::FemEoc::setErrors( ids[ 0 ], l2perror );
    }
  };

  static inline std::string moduleName()
  {
    return "";
  }

  static inline Dune::GridPtr<GridType>
  initializeGrid()
  {
    return PoissonProblemCreatorType::initializeGrid();
  }

  static ProblemInterfaceType* problem()
  {
    return new Dune::GeneralizedStokesProblem< GridType > ();
  }


  //Stepper Traits
  template< class GridPart, int polOrd > // TODO: is GridPart as a template parameter needed?
  struct DiscreteTraits
  {
    typedef typename PoissonProblemCreatorType::template DiscreteTraits< GridPart, polOrd >       PoissonDiscreteTraits;
public:
    typedef AnalyticalTraits<GridPart>                                                            AnalyticalTraitsType;

    static const int polynomialOrder = PoissonDiscreteTraits::polynomialOrder;

    static inline std::string advectionFluxName()
    {
      return PoissonDiscreteTraits::advectionFluxName();
    }

    static inline std::string diffusionFluxName()
    {
      return PoissonDiscreteTraits::diffusionFluxName();
    }

    static const int quadOrder = PoissonDiscreteTraits::quadOrder;
    static const SolverType solverType = PoissonDiscreteTraits::solverType;

    typedef typename DiscreteFunctionSpaces< PressureFunctionSpaceType, GridPartType, polynomialOrder, _legendre, dg >::type    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::type                                   DiscreteFunctionType;
    typedef typename DiscreteFunctions< DiscreteFunctionSpaceType, solverType >::jacobian                               JacobianOperatorType;

    typedef typename  std::tuple<> ExtraParameterTuple;


private:
    typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, AnalyticalTraitsType::ModelType::dimDomain, 3> FVFunctionSpaceType;
    typedef Dune::Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Dune::Fem::SimpleStorage> IndicatorSpaceType;
public:
    typedef Dune::Fem::AdaptiveDiscreteFunction<IndicatorSpaceType>                             IndicatorType;

    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >                           RestrictionProlongationType;

    typedef Dune::AdaptationHandler< GridType, VelocityFunctionSpaceType >                              AdaptationHandlerType;


    //############################ stokes ############################
    static const bool symmetricSolver = true ;
    //typedef Solvers<DiscreteFunctionSpaceType, solverType, symmetricSolver> SolversType;
    typedef Dune::NoFlux< typename AnalyticalTraitsType::ModelType >              FluxType;
    typedef typename PoissonDiscreteTraits::FullOperatorType                      FullOperatorType;

    //----------- passes! ------------------------
    typedef Dune::OperatorTraits< GridPartType, polynomialOrder, AnalyticalTraitsType,
                                  DiscreteFunctionType, FluxType, IndicatorType,
                                  AdaptationHandlerType, ExtraParameterTuple >                  OperatorTraitsType;
    typedef typename Dune::StokesAssembler< typename PoissonDiscreteTraits::DiscreteFunctionType,
                                            DiscreteFunctionType,
                                            OperatorTraitsType>                   AssemblerType;
    //------------------------------------------------------------

    static_assert( (int)DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == 1 , "pressure dimrange does not fit");

    typedef Dune::UzawaSolver< typename PoissonDiscreteTraits::DiscreteFunctionType, DiscreteFunctionType, AssemblerType,
                               typename PoissonDiscreteTraits::BasicLinearSolverType >    BasicLinearSolverType;

    //#########################################################################

    //------------- Limiter ---------------------------------------------
    typedef typename PoissonDiscreteTraits::LimiterOperatorType                                         LimiterOperatorType;
    //------------ Limiter ---------------------------------------------

    typedef typename AnalyticalTraitsType::ProblemType::ExactPressureType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;

    typedef std::tuple< typename std::tuple_element<0,typename PoissonDiscreteTraits::IOTupleType>::type, typename std::tuple_element<1,typename PoissonDiscreteTraits::IOTupleType>::type,
                        GridExactSolutionType*, DiscreteFunctionType* > IOTupleType;
  };

  template <int polOrd>
  struct Stepper
  {
    // this should be ok but could lead to a henn-egg problem
    typedef Dune::Fem::StokesAlgorithm< GridType, StokesProblemCreator<GridType>, PoissonProblemCreatorType, polOrd > Type;
  };


};

#define NEW_STEPPER_SELECTOR_USED
#endif // FEMHOWTO_HEATSTEPPER_HH
