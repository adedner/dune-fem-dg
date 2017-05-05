#ifndef FEMHOWTO_HEATSTEPPER_HH
#define FEMHOWTO_HEATSTEPPER_HH
#include <config.h>

#ifndef DIMRANGE
#define DIMRANGE 1
#endif

#ifndef POLORDER
#define POLORDER 2
#endif

#include <dune/fem/solver/timeprovider.hh>
#include <utility>
#include <dune/fem-dg/algorithm/createsubalgorithms.hh>

//--------- CALLER --------------------------------
#include <dune/fem-dg/algorithm/caller/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/sub/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/sub/adapt.hh>
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/algorithm/caller/eocwriter.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/examples/dataio/checkpointing.hh>
#include <dune/fem-dg/algorithm/evolution.hh>
#include <dune/fem-dg/algorithm/steadystate.hh>
//--------- EOCERROR ------------------------
#include <dune/fem-dg/misc/error/l2eocerror.hh>
//--------- PROBLEMS ------------------------
#include "problem.hh"
//--------- MODELS --------------------------
#include "models.hh"
//....
#include <dune/fem/space/discontinuousgalerkin/space.hh>


namespace Dune
{
namespace Fem
{
  // EvolutionAlgorithmTraits
  // -------------------------
  template< int polOrder, class ... ProblemTraits >
  struct CheckPointEvolutionAlgorithmTraits
  {
    // type of Grid
    typedef typename std::tuple_element_t<0, std::tuple< ProblemTraits... > >::GridType      GridType;

    typedef std::tuple< typename std::add_lvalue_reference<typename ProblemTraits::GridType>::type ... >
                                                                                             GridTypes;
    // wrap operator
    typedef Dune::Fem::GridTimeProvider< GridType >                                          TimeProviderType;

    typedef CreateSubAlgorithms< typename ProblemTraits::template Algorithm<polOrder>...  >  CreateSubAlgorithmsType;

    typedef typename CreateSubAlgorithmsType::SubAlgorithmTupleType                          SubAlgorithmTupleType;

    typedef typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >
                                                                                             IndexSequenceType;
    typedef std::index_sequence<>                                                            NoIndexSequenceType;

    typedef Dune::Fem::AdaptCaller< SubAlgorithmTupleType, NoIndexSequenceType >             AdaptCallerType;
    typedef Dune::Fem::DiagnosticsCaller < SubAlgorithmTupleType, NoIndexSequenceType >      DiagnosticsCallerType;
    typedef Dune::Fem::SolverMonitorCaller < SubAlgorithmTupleType, NoIndexSequenceType >    SolverMonitorCallerType;
    typedef Dune::Fem::CheckedCheckPointCaller < SubAlgorithmTupleType >                     CheckPointCallerType;
    typedef Dune::Fem::DataWriterCaller < SubAlgorithmTupleType >                            DataWriterCallerType;
    typedef Dune::Fem::EocWriterCaller < SubAlgorithmTupleType >                             EocWriterCallerType;
    typedef Dune::Fem::PostProcessingCaller < SubAlgorithmTupleType, NoIndexSequenceType >   PostProcessingCallerType;

  };


  /**
   *  \brief problem creator for an advection diffusion problem
   */
  template< class GridSelectorGridType >
  struct CheckPointingAlgorithmCreator
  {

    struct SubCheckPointingAlgorithmCreator
    {

      typedef GridSelectorGridType                            GridType;
      typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
      //typedef Dune::Fem::LeafGridPart< GridType >           HostGridPartType;
      //typedef AdaptiveLeafGridPart< GridType >              HostGridPartType;
      //typedef IdBasedLeafGridPart< GridType >               HostGridPartType;
      typedef HostGridPartType                                GridPartType;

      typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

      // define problem type here if interface should be avoided
      typedef U0< GridType >                                            ProblemInterfaceType;

      typedef ProblemInterfaceType                                    ProblemType;
      typedef NoModel< GridType, ProblemType >                        ModelType;

      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem()
      {
        return new ProblemInterfaceType();
      }

      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrd > DiscreteFunctionSpaceType;
      public:
        typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >                 DiscreteFunctionType;

        typedef void                                                                             AdaptIndicatorType;
        typedef void                                                                             SolverMonitorType;
        typedef void                                                                             DiagnosticsType;
        typedef SubDataWriter< SolutionOutput<DiscreteFunctionType> >                            DataWriterType;

      };

      template <int polOrd>
      using Algorithm = Dune::Fem::SubCheckPointingAlgorithm< GridType, SubCheckPointingAlgorithmCreator, polOrd >;

    };

    template <int polOrd>
    using Algorithm = Dune::Fem::EvolutionAlgorithmBase< CheckPointEvolutionAlgorithmTraits< polOrd, SubCheckPointingAlgorithmCreator >  >;


    static inline std::string moduleName() { return ""; }

    template< int polOrd >
    static decltype(auto) initContainer()
    {
      typedef typename SubCheckPointingAlgorithmCreator::GridType         GridType;
      //Discrete Functions
      typedef typename SubCheckPointingAlgorithmCreator::template DiscreteTraits<polOrd>::DiscreteFunctionType
                                                                          DFType;

      //Item1
      typedef _t< SubEvolutionContainerItem >                             Steady;
      typedef std::tuple< Steady >                                        Item1TupleType;

      //Item2
      typedef _t< EmptyContainerItem >                                    Empty;
      typedef std::tuple< std::tuple< Empty > >                           Item2TupleType;


      //Sub (discrete function argument ordering)
      typedef std::tuple<__0 >                                            AdvDiffOrder;

      typedef std::tuple< AdvDiffOrder >                                  SubOrderRowType;
      typedef SubOrderRowType                                             SubOrderColType;

      //external params lists
      typedef ExtraArg<>                                                  ExtraType;

      //Global container
      typedef GlobalContainer< Item2TupleType, Item1TupleType, SubOrderRowType, SubOrderColType,ExtraType,DFType >
                                                                          GlobalContainerType;

      //create grid
      std::shared_ptr< GridType > gridptr( DefaultGridInitializer< GridType >::initialize().release() );

      //create container
      return std::make_shared< GlobalContainerType >( moduleName(), gridptr );
    }
  };

}
}
#endif // FEMHOWTO_HEATSTEPPER_HH

