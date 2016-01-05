#ifndef DUNE_FEMDG_NAVIERSTOKES_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_NAVIERSTOKES_ALGORITHM_EVOLUTION_HH


#include "../../algorithm/evolution.hh"

namespace Dune
{
namespace Fem
{

  // ----------------------

  template< int polOrder, class ... ProblemTraits >
  class IncompNavierStokesAlgorithm
    : public EvolutionAlgorithm< polOrder, ProblemTraits... >
  {
    typedef EvolutionAlgorithm< polOrder, ProblemTraits... >     BaseType;
  public:
    typedef typename BaseType::GridType                          GridType;
    typedef typename BaseType::IOTupleType                       IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType          SolverMonitorHandlerType;

    typedef typename BaseType::StepperTupleType                  StepperTupleType;
    typedef typename BaseType::TimeProviderType                  TimeProviderType;

    typedef typename BaseType::DiagnosticsHandlerType            DiagnosticsHandlerType;
    typedef typename BaseType::CheckPointHandlerType             CheckPointHandlerType;
    typedef typename BaseType::DataWriterHandlerType             DataWriterHandlerType;
    typedef typename BaseType::SolutionLimiterHandlerType        SolutionLimiterHandlerType;
    typedef typename BaseType::AdaptHandlerType                  AdaptHandlerType;

    typedef typename BaseType::UInt64Type                        UInt64Type ;

    typedef typename BaseType::StepperParametersType             StepperParametersType;

    using BaseType::eocParams;
    using BaseType::grid;

    static const int numSteppers = std::tuple_size< StepperTupleType >::value;

    static_assert( numSteppers == 3, "This InCompNavierStokesAlgorithm needs three sub algorithms: 1. Stokes, 2. Oseen, 3. Stokes" );

  public:

    IncompNavierStokesAlgorithm ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  )
    {}

    //define your own time step
    virtual void step ( int loop, TimeProviderType &tp )
    {
      auto step1 = std::get<0>( BaseType::stepperTuple() );
      auto step2 = std::get<1>( BaseType::stepperTuple() );
      auto step3 = std::get<2>( BaseType::stepperTuple() );

      //stokes
      step1->solve( loop, &tp );
      //oseen
      step2->solve( loop, &tp );
      //stokes
      step3->solve( loop, &tp );
    }

  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
