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
#include <dune/common/std/utility.hh>

//--------- HANDLER --------------------------------
#include <dune/fem-dg/algorithm/handler/sub/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/sub/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/sub/additionaloutput.hh>
#include <dune/fem-dg/algorithm/handler/sub/adapt.hh>
#include <dune/fem-dg/algorithm/monitor.hh>

//--------- GRID HELPER ---------------------
#include <dune/fem-dg/algorithm/gridinitializer.hh>
//--------- OPERATOR/SOLVER -----------------
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
//--------- FLUXES ---------------------------
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
//--------- STEPPER -------------------------
#include <dune/fem-dg/test/dataio/checkpointing.hh>
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
    typedef typename std::tuple_element<0, std::tuple< ProblemTraits... > >::type::GridType  GridType;

    // wrap operator
    typedef Dune::Fem::GridTimeProvider< GridType >                                   TimeProviderType;

    typedef std::tuple< typename std::add_pointer< typename ProblemTraits::template Stepper<polOrder>::Type >::type... > SubAlgorithmTupleType;

    typedef typename Dune::Std::make_index_sequence_impl< std::tuple_size< SubAlgorithmTupleType >::value >::type  IndexSequenceType;
    typedef Dune::Std::index_sequence<>                                                                       NoIndexSequenceType;

    typedef Dune::Fem::AdaptHandler< SubAlgorithmTupleType, NoIndexSequenceType >           AdaptHandlerType;
    typedef Dune::Fem::DiagnosticsHandler < SubAlgorithmTupleType, NoIndexSequenceType >    DiagnosticsHandlerType;
    typedef Dune::Fem::SolverMonitorHandler < SubAlgorithmTupleType, NoIndexSequenceType >  SolverMonitorHandlerType;
    typedef Dune::Fem::CheckedCheckPointHandler < SubAlgorithmTupleType >                   CheckPointHandlerType;
    typedef Dune::Fem::DataWriterHandler < SubAlgorithmTupleType >                          DataWriterHandlerType;
    typedef Dune::Fem::PostProcessingHandler < SubAlgorithmTupleType, NoIndexSequenceType > PostProcessingHandlerType;

    typedef typename DataWriterHandlerType::IOTupleType                                                                IOTupleType;

    template< std::size_t ...i >
    static SubAlgorithmTupleType createSubAlgorithm ( Dune::Std::index_sequence< i... >, GridType &grid, const std::string name = "" )
    {
      return std::make_tuple( new typename std::remove_pointer< typename std::tuple_element< i, SubAlgorithmTupleType >::type >::type( grid, name ) ... );
    }

    // create Tuple of contained sub algorithms
    static SubAlgorithmTupleType createSubAlgorithm( GridType &grid, const std::string name = "" )
    {
      return createSubAlgorithm( Dune::Std::index_sequence_for< ProblemTraits ... >(), grid, name );
    }
  };


  /**
   *  \brief problem creator for an advection diffusion problem
   */
  template< class GridImp >
  struct CheckPointingProblemCreator
  {

    struct SubCheckPointingProblemCreator
    {

      typedef GridImp                                         GridType;
      typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType >   HostGridPartType;
      //typedef Dune::Fem::LeafGridPart< GridType >           HostGridPartType;
      //typedef AdaptiveLeafGridPart< GridType >              HostGridPartType;
      //typedef IdBasedLeafGridPart< GridType >               HostGridPartType;
      typedef HostGridPartType                                GridPartType;

      typedef Dune::Fem::FunctionSpace< typename GridType::ctype, double, GridType::dimension, DIMRANGE> FunctionSpaceType;

      // define problem type here if interface should be avoided
      typedef U0< GridType >                                            ProblemInterfaceType;

      struct AnalyticalTraits
      {
        typedef ProblemInterfaceType                                    ProblemType;
        typedef ProblemInterfaceType                                    InitialDataType;
        typedef NoModel< GridPartType, ProblemType >                    ModelType;

        template< class Solution, class Model, class ExactFunction, class TimeProvider >
        static void addEOCErrors ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
        {
          //static L2EOCError l2EocError( "$L^2$-Error");
          //l2EocError.add( tp, u, model, f );
        }
      };

      static inline std::string moduleName() { return ""; }

      static ProblemInterfaceType* problem()
      {
        return new ProblemInterfaceType();
      }


      //Stepper Traits
      template< int polOrd >
      struct DiscreteTraits
      {
      private:
        typedef Dune::Fem::DiscontinuousGalerkinSpace < FunctionSpaceType, GridPartType, polOrd >     DiscreteFunctionSpaceType;
      public:
        typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >                      DiscreteFunctionType;

        typedef std::tuple< DiscreteFunctionType*, DiscreteFunctionType* >                            IOTupleType;

        typedef void                                                                                  AdaptIndicatorType;
        typedef void                                                                                  SolverMonitorHandlerType;
        typedef void                                                                                  DiagnosticsHandlerType;
        typedef void                                                                                  AdditionalOutputHandlerType;
      };


      template <int polOrd>
      struct Stepper
      {
       // this should be ok but could lead to a henn-egg problem
        typedef Dune::Fem::SubCheckPointingAlgorithm< GridType, SubCheckPointingProblemCreator, polOrd > Type;
      };

    };

    template <int polOrd>
    struct Stepper
    {
      typedef Dune::Fem::EvolutionAlgorithmBase< CheckPointEvolutionAlgorithmTraits< polOrd, SubCheckPointingProblemCreator >, DefaultSteadyStateCreator  > Type;
    };

    typedef GridImp                                         GridType;

    static inline std::string moduleName() { return ""; }

    static inline Dune::GridPtr<GridType>
    initializeGrid() { return Dune::Fem::DefaultGridInitializer< GridType >::initialize(); }


  };

}
}
#endif // FEMHOWTO_HEATSTEPPER_HH

