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
//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- STEPPER -------------------------
#include "incompnavierstokesalgorithm.hh"
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
    using Algorithm = IncompNavierStokesAlgorithm< polOrd, SubStokesAlgorithmCreator<ACStokes>,
                                                           SubNavierStokesAlgorithmCreator<ACNavier>,
                                                           SubStokesAlgorithmCreator<ACStokes> >;


    static inline std::string moduleName() { return ""; }

    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef typename SubStokesAlgorithmCreator<ACStokes>::GridType  GridType;

      //Discrete Functions
      typedef typename SubStokesAlgorithmCreator<ACStokes>::SubCreatorType::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType1;
      typedef typename SubStokesAlgorithmCreator<ACStokes>::template DiscreteTraits<polOrd>::DiscreteFunctionType DFType2;


      //Item1
      typedef std::tuple< _t< SubSteadyStateContainerItem >, _t< SubSteadyStateContainerItem >  >
                                                                      Steady;
      typedef std::tuple< _t< SubEvolutionContainerItem > >           Evol;
      typedef std::tuple< Steady, Evol, Steady >                      Item1TupleType;

      //Item2
      // typedef _t< SubEllipticContainerItem, ISTLLinearOperator >      Istl;
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

      //create container
      return std::make_shared< GlobalContainerType >( gridptr, "global" /*, gridptr, gridptr */ );

    }
  };

}
}



#endif
