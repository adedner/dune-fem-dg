#ifndef DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#define DUNE_FEM_DG_INCOMP_NAVIERSTOKES_HH
#include <config.h>

#ifndef GRIDDIM
#define GRIDDIM 2
#endif

#ifndef DIMRANGE
#define DIMRANGE GRIDDIM
#endif

#ifndef POLORDER
#define POLORDER 1
#endif

#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/fem/misc/l2norm.hh>
//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>
//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/algorithm/sub/advectiondiffusion.hh>
#include <dune/fem-dg/algorithm/sub/advection.hh>
#include <dune/fem-dg/examples/stokes/stokesalgorithm.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
#include "incompnavierstokesalgorithm.hh"
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/h1eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
//--------- PROBLEMS ------------------------
#include "problems.hh"
//#include "../stokes/problems.hh"
//--------- MODELS --------------------------
#include "stokesmodel.hh"
#include "models.hh"
//--------- PROBLEMCREATORSELECTOR ----------
#include <dune/fem-dg/misc/configurator.hh>
//--------- SUBALGORITHMCREATOR -------------
#include "substokesalgorithmcreator.hh"
#include "subnavierstokesalgorithmcreator.hh"

namespace Dune
{
namespace Fem
{

  template< class GridSelectorGridType >
  struct IncompressibleNavierStokesAlgorithmCreator
  {

    typedef AlgorithmConfigurator< GridSelectorGridType,
                                   Galerkin::Enum::dg,
                                   Adaptivity::Enum::yes,
                                   DiscreteFunctionSpaces::Enum::legendre,
                                   Solver::Enum::fem,
                                   AdvectionLimiter::Enum::unlimited,
                                   Matrix::Enum::assembled,
                                   AdvectionFlux::Enum::none,
                                   DiffusionFlux::Enum::primal > ACStokes;

    typedef AlgorithmConfigurator< GridSelectorGridType,
                                   Galerkin::Enum::dg,
                                   Adaptivity::Enum::yes,
                                   DiscreteFunctionSpaces::Enum::legendre,
                                   Solver::Enum::fem,
                                   AdvectionLimiter::Enum::unlimited,
                                   Matrix::Enum::matrixfree,
                                   AdvectionFlux::Enum::llf,
                                   DiffusionFlux::Enum::primal > ACNavier;



    template <int polOrd>
    using Algorithm = IncompNavierStokesAlgorithm< polOrd, UncoupledSubAlgorithms, SubStokesAlgorithmCreator<ACStokes>,
                                                                                   SubNavierStokesAlgorithmCreator<ACNavier>,
                                                                                   SubStokesAlgorithmCreator<ACStokes> >;

    typedef typename SubStokesAlgorithmCreator<ACStokes>::GridType                  GridType;

    static inline std::string moduleName() { return ""; }

    static inline GridPtr<GridType>
    initializeGrid() { return DefaultGridInitializer< GridType >::initialize(); }

    template< int polOrd >
    static decltype(auto) initContainer()
    {
      //Discrete Functions
      typedef typename SubStokesAlgorithmCreator<ACStokes>::SubCreatorType::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType1;
      typedef typename SubStokesAlgorithmCreator<ACStokes>::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType2;


      //Item1
      typedef std::tuple< _t< SubSteadyStateContainerItem >, _t< SubSteadyStateContainerItem >  >
                                                                      Steady;
      typedef std::tuple< _t< SubEvolutionContainerItem > >           Evol;
      typedef std::tuple< Steady, Evol, Steady >                      Item1TupleType;

      //Item2
      typedef _t< SubEllipticContainerItem, ISTLLinearOperator >      Istl;
      typedef _t< SubEllipticContainerItem, SparseRowLinearOperator > Sp;
      typedef _t< EmptyContainerItem >                                Empty;

      typedef std::tuple< std::tuple< Sp, Sp >,
                          std::tuple< Sp, Sp > >                      Stokes;
      typedef std::tuple< std::tuple< Empty > >                       AdvDiff;

      typedef std::tuple< Stokes, AdvDiff, Stokes >                   Item2TupleType;


      //Sub (discrete function argument ordering)
      typedef std::tuple<__0,__1 >                                    StokesOrder;
      typedef std::tuple<__0 >                                        AdvDiffOrder;

      typedef std::tuple< StokesOrder, AdvDiffOrder, StokesOrder >    SubOrderRowType;
      typedef SubOrderRowType                                         SubOrderColType;


      //Global container
      typedef CombinedGlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, DFType1, DFType2 > GlobalContainerType;


      ////TODO
      //typedef ThetaSchemeCouplingContainer< bla, bla, blubb >;

      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      typedef typename DFType1::DiscreteFunctionSpaceType SpaceType1;
      typedef typename DFType2::DiscreteFunctionSpaceType SpaceType2;

      typedef typename SpaceType1::GridPartType GridPartType1;
      typedef typename SpaceType2::GridPartType GridPartType2;

      static_assert( std::is_same< GridPartType1, GridPartType2 >::value, "GridParts does not match!");

      {
        std::shared_ptr< GridType > gridptr2 = gridptr;

        ////TODO: deref of gridptr dangerous
        GridPartType1 gridPart( *gridptr2 );
      }

      //SpaceType1 space1( gridPart );
      //SpaceType2 space2( gridPart );

      //auto shared_space1 = stackobject_to_shared_ptr(space1);
      //auto shared_space2 = stackobject_to_shared_ptr(space2);

      //std::make_shared< GlobalContainerType >( space1, space2 );

      //create container
      return std::make_shared< GlobalContainerType >( gridptr, "global" /*, gridptr, gridptr */ );

    }
  };

}
}



#endif
