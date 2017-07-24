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

  private:
    //helper struct
    template< class DDF, class RDF >
    using DefaultContainer = typename ACStokes::template Containers< typename DDF::DiscreteFunctionSpaceType, typename RDF::DiscreteFunctionSpaceType >;

  public:
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
      typedef std::tuple< _t< SubEvolutionContainerItem > >               Evol;
      typedef std::tuple< Steady, Evol, Steady >                          Item1TupleType;

#if 0
      // solution()
      //  |
      //  V
      // solution(), solutionNTheta()
      // velocityN()
      typedef std::tuple< _t< TimeStepN, SubSteadyStateContainerItem >, _t< TimeStepN, SubSteadyStateContainerItem >  >
                                                                               Steady1;
      // solution()
      //  |
      //  V
      // solution(), solutionNPlus1MinusTheta(),
      // velocityN(), velocityNTheta() = velocityStar()
      typedef std::tuple< _t< ExternalVelocity, TimeStepNTheta, SubEvolutionContainerItem > >    Evol;

      //solution()
      //  |
      //  V
      // solution(), solutionNPlus1(),
      // velocityNPlus1MinusTheta(),
      // pressureNTheta()
      // velocityN(), velocityNTheta() = velocityStar()
      typedef std::tuple< _t< TimeStepNPlus1MinusTheta, SubSteadyStateContainerItem >, _t< TimeStepNPlus1MinusTheta, SubSteadyStateContainerItem >  >
                                                                               Steady2;

      typedef std::tuple< Steady1, Evol, Steady2 >                        Item1TupleType;
#endif

      //external parameters
      typedef _ee< __2, __0 >                                             Un;
      typedef _ee< __0, __0 >                                             UnTheta;
      typedef _ee< __1, __0 >                                             UnPlus1MinusTheta;

      typedef _ee< __2, __1 >                                             Pn;
      typedef _ee< __0, __1 >                                             PnTheta;
      //Does not exist:
      //typedef _ee< __1, __1 >                                             PnPlus1MinusTheta;


      //external params lists
      typedef ExtraArg< Un >                                              StokesExtra1Type;
      typedef ExtraArg< Un, UnTheta, PnTheta >                            AdvDiffExtraType;
      typedef ExtraArg< Un, UnTheta, UnPlus1MinusTheta >                  StokesExtra2Type;

      typedef std::tuple< StokesExtra1Type, AdvDiffExtraType, StokesExtra2Type > ExtraType;

#if 0
      //IOData (uses raw pointers...)
      //------
      typedef _e< IOSolutionSelect, __2, __0 >                            U;
      typedef _e< IOSolutionSelect, __2, __1 >                            P;
      typedef _e< IOExactSolutionSelect, __2, __0 >                       UExact;
      typedef _e< IOExactSolutionSelect, __2, __1 >                       PExact;

      //IO data list
      typedef ExtraArg< U, P, UExact, PExact >                            IODataTupleType;

      typedef std::tuple< IODataTupleType >                               IOExtraType;
#endif

      //Item2
      typedef _t< SubEllipticContainerItem, SparseRowLinearOperator >     Sp;
      typedef _t< EmptyContainerItem >                                    Empty;
      typedef _t< SubEllipticContainerItem, DefaultContainer >            Def;

      typedef std::tuple< std::tuple< Def, Sp >,
                          std::tuple< Sp,  Sp > >                         Stokes;
      typedef std::tuple< std::tuple< Empty > >                           AdvDiff;

      typedef std::tuple< Stokes, AdvDiff, Stokes >                       Item2TupleType;


      //Sub (discrete function argument ordering)
      typedef std::tuple<__0,__1 >                                        StokesOrder;
      typedef std::tuple<__0 >                                            AdvDiffOrder;

      typedef std::tuple< StokesOrder, AdvDiffOrder, StokesOrder >        SubOrderRowType;
      typedef SubOrderRowType                                             SubOrderColType;


      //Global container
      typedef CombinedGlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType, ExtraType, DFType1, DFType2 >
                                                                          GlobalContainerType;

      ////TODO
      //typedef ThetaSchemeCouplingContainer< bla, bla, blubb >;

      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( moduleName(), gridptr, gridptr, gridptr );

    }
  };

}
}



#endif
