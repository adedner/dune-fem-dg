#ifndef DUNE_FEMDG_NAVIERSTOKES_ALGORITHM_EVOLUTION_HH
#define DUNE_FEMDG_NAVIERSTOKES_ALGORITHM_EVOLUTION_HH


#include <dune/fem-dg/algorithm/evolution.hh>

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

    typedef typename BaseType::SubAlgorithmTupleType             SubAlgorithmTupleType;
    typedef typename BaseType::TimeProviderType                  TimeProviderType;

    typedef typename BaseType::SolverMonitorCallerType           SolverMonitorCallerType;
    typedef typename BaseType::DiagnosticsCallerType             DiagnosticsCallerType;
    typedef typename BaseType::CheckPointCallerType              CheckPointCallerType;
    typedef typename BaseType::DataWriterCallerType              DataWriterCallerType;
    typedef typename BaseType::PostProcessingCallerType          PostProcessingCallerType;
    typedef typename BaseType::AdaptCallerType                   AdaptCallerType;

    typedef typename BaseType::UInt64Type                        UInt64Type ;

    typedef typename BaseType::TimeSteppingParametersType        TimeSteppingParametersType;

    using BaseType::eocParams;
    using BaseType::grid;

    static_assert( std::tuple_size< SubAlgorithmTupleType >::value == 3,
                   "This InCompNavierStokesAlgorithm needs three sub algorithms: 1. Stokes, 2. Oseen, 3. Stokes" );

  public:
    template< class GlobalContainerImp >
    IncompNavierStokesAlgorithm ( const std::string name, const std::shared_ptr<GlobalContainerImp>& cont )
    : BaseType( name, cont  )
    {}


    virtual void initialize ( int loop, TimeProviderType &tp )
    {
      auto step1 = std::get<0>( BaseType::subAlgorithmTuple() );
      auto step2 = std::get<1>( BaseType::subAlgorithmTuple() );
      auto step3 = std::get<2>( BaseType::subAlgorithmTuple() );

      BaseType::initialize( loop, tp );
    }

    virtual void preSolve ( int loop, TimeProviderType &tp )
    {
      auto step1 = std::get<0>( BaseType::subAlgorithmTuple() );
      auto step2 = std::get<1>( BaseType::subAlgorithmTuple() );
      auto step3 = std::get<2>( BaseType::subAlgorithmTuple() );

      const double theta = step1->model().problem().template get<0>().theta();
      const double time  = tp.time();
      const double dt    = tp.deltaT();

      //stokes
      step1->model().problem().template get<0>().setTime( time );
      step1->model().problem().template get<1>().setTime( time );
      step1->model().problem().template get<0>().setDeltaT( dt );
      step1->model().problem().template get<1>().setDeltaT( dt );
      //step1->setTime( time );
      //step1->setDeltaT( dt );

      //oseen
      TimeProviderType subTimeProvider( 0.0, 1.0, grid() );
      subTimeProvider.init();
      subTimeProvider.provideTimeStepEstimate( time + dt * theta );
      subTimeProvider.next();
      subTimeProvider.next( (1.0 - 2.0 *theta) * dt );

      //stokes
      step3->model().problem().template get<0>().setTime( time + dt * ( 1.0 - theta ));
      step3->model().problem().template get<1>().setTime( time + dt * ( 1.0 - theta ));
      step3->model().problem().template get<0>().setDeltaT( dt );
      step3->model().problem().template get<1>().setDeltaT( dt );
      //step3->setTime( time + dt * ( 1.0 - theta ));
      //step3->setDeltaT( dt );

      BaseType::preSolve( loop, tp );
    }

    //define your own time step
    virtual void solve ( int loop, TimeProviderType &tp )
    {
      auto step1 = std::get<0>( BaseType::subAlgorithmTuple() );
      auto step2 = std::get<1>( BaseType::subAlgorithmTuple() );
      auto step3 = std::get<2>( BaseType::subAlgorithmTuple() );

      const double theta = step1->model().problem().template get<0>().theta();
      const double time  = tp.time();
      const double dt    = tp.deltaT();


      //stokes
      std::cout << "###### STOKES I ###### time: " << step1->model().problem().template get<0>().time() << " dt: " << step1->model().problem().template get<0>().deltaT() << std::endl;
      step1->solve( loop, &tp );

      //oseen
      TimeProviderType subTimeProvider( time, 1.0, grid() );
      subTimeProvider.init( time + dt * theta );
      subTimeProvider.next( (1.0 - 2.0 *theta) * dt );
      //subTimeProvider.init();
      //subTimeProvider.provideTimeStepEstimate( time + dt * theta );
      //subTimeProvider.next();
      //subTimeProvider.next( (1.0 - 2.0 *theta) * dt );

      std::cout << "###### OSEEN I ###### time: " << subTimeProvider.time() << " dt: " << dt << std::endl;
      step2->solve( loop, &subTimeProvider );


      //stokes
      std::cout << "###### STOKES II ###### time: " << step3->model().problem().template get<0>().time() << " dt: " << step3->model().problem().template get<0>().deltaT() << std::endl;
      step3->solve( loop, &tp );



      const double dtEstimate = subTimeProvider.timeStepEstimate();
      tp.provideTimeStepEstimate( dtEstimate / (1.0 - 2.0 * theta) );
    }

  };



} // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_ALGORITHM_EVOLUTION_HH
