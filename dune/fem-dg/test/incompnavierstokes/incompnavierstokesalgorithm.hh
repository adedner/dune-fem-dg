#ifndef DUNE_FEMDG_NAVIERSTOKES_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_NAVIERSTOKES_ALGORITHM_EVOLUTION_HH


#include "../../algorithm/evolution.hh"

namespace Dune
{
namespace Fem
{

  // ----------------------

  template< int polOrder, template<class, class > class EvolutionCreatorType, class ... ProblemTraits >
  class IncompNavierStokesAlgorithm
    : public EvolutionAlgorithm< polOrder, EvolutionCreatorType, ProblemTraits... >
  {
    typedef EvolutionAlgorithm< polOrder, EvolutionCreatorType, ProblemTraits... >     BaseType;
  public:
    typedef typename BaseType::GridType                          GridType;
    typedef typename BaseType::IOTupleType                       IOTupleType;
    typedef typename BaseType::SolverMonitorHandlerType          SolverMonitorHandlerType;

    typedef typename BaseType::SubAlgorithmTupleType             SubAlgorithmTupleType;
    typedef typename BaseType::TimeProviderType                  TimeProviderType;

    typedef typename BaseType::DiagnosticsHandlerType            DiagnosticsHandlerType;
    typedef typename BaseType::CheckPointHandlerType             CheckPointHandlerType;
    typedef typename BaseType::DataWriterHandlerType             DataWriterHandlerType;
    typedef typename BaseType::SolutionLimiterHandlerType        SolutionLimiterHandlerType;
    typedef typename BaseType::AdaptHandlerType                  AdaptHandlerType;

    typedef typename BaseType::UInt64Type                        UInt64Type ;

    typedef typename BaseType::TimeSteppingParametersType        TimeSteppingParametersType;

    using BaseType::eocParams;
    using BaseType::grid;

    static_assert( std::tuple_size< SubAlgorithmTupleType >::value == 3,
                   "This InCompNavierStokesAlgorithm needs three sub algorithms: 1. Stokes, 2. Oseen, 3. Stokes" );

  public:
    typedef EvolutionCreatorType< SubAlgorithmTupleType, GridType > CreatorType;

    IncompNavierStokesAlgorithm ( GridType &grid, const std::string name = "" )
    : BaseType( grid, name  )
    {}

    virtual void initialize ( int loop, TimeProviderType &tp )
    {
      auto step1 = std::get<0>( BaseType::subAlgorithmTuple() );
      auto step2 = std::get<1>( BaseType::subAlgorithmTuple() );
      auto step3 = std::get<2>( BaseType::subAlgorithmTuple() );


      //concate solutions()
      //step3->setSolution( std::make_shared< typename std::remove_reference< decltype( step1->solution() ) >::type >( step1->solution() ) );

      BaseType::initialize( loop, tp );
    }

    virtual void preStep ( int loop, TimeProviderType &tp )
    {
      auto step1 = std::get<0>( BaseType::subAlgorithmTuple() );
      auto step2 = std::get<1>( BaseType::subAlgorithmTuple() );
      auto step3 = std::get<2>( BaseType::subAlgorithmTuple() );

      const double theta = step1->problem().get<0>().theta();
      const double time  = tp.time();
      const double dt    = tp.deltaT();

      //stokes
      step1->problem().get<0>().setTime( time );
      step1->problem().get<1>().setTime( time );
      step1->problem().get<0>().setDeltaT( dt );
      step1->problem().get<1>().setDeltaT( dt );

      //oseen
      TimeProviderType subTimeProvider( 0.0, 1.0, grid() );
      subTimeProvider.init();
      subTimeProvider.provideTimeStepEstimate( time + dt * theta );
      subTimeProvider.next();
      subTimeProvider.next( (1.0 - 2.0 *theta) * dt );

      //stokes
      step3->problem().get<0>().setTime( time + dt * ( 1.0 - theta ));
      step3->problem().get<1>().setTime( time + dt * ( 1.0 - theta ));
      step3->problem().get<0>().setDeltaT( dt );
      step3->problem().get<1>().setDeltaT( dt );

      BaseType::preStep( loop, tp );
    }

    //define your own time step
    virtual void step ( int loop, TimeProviderType &tp )
    {
      auto step1 = std::get<0>( BaseType::subAlgorithmTuple() );
      auto step2 = std::get<1>( BaseType::subAlgorithmTuple() );
      auto step3 = std::get<2>( BaseType::subAlgorithmTuple() );

      const double theta = step1->problem().get<0>().theta();
      const double time  = tp.time();
      const double dt    = tp.deltaT();

      //stokes
      step1->solve( loop, &tp );

      //oseen
      TimeProviderType subTimeProvider( 0.0, 1.0, grid() );
      subTimeProvider.init();
      subTimeProvider.provideTimeStepEstimate( time + dt * theta );
      subTimeProvider.next();
      subTimeProvider.next( (1.0 - 2.0 *theta) * dt );

      step2->solve( loop, &subTimeProvider );

      //stokes
      step3->solve( loop, &tp );



      const double dtEstimate = subTimeProvider.timeStepEstimate();
      tp.provideTimeStepEstimate( dtEstimate / (1.0 - 2.0 * theta) );
    }

  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
