#ifndef DUNE_FEMDG_ALGORITHM_COMBINEDSTEADYSTATEALGORITHM_HH
#define DUNE_FEMDG_ALGORITHM_COMBINEDSTEADYSTATEALGORITHM_HH

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/createsubalgorithms.hh>

#include <dune/fem-dg/algorithm/caller/solvermonitor.hh>

#include <dune/fem-dg/algorithm/caller/diagnostics.hh>
#include <dune/fem-dg/algorithm/caller/solvermonitor.hh>
#include <dune/fem-dg/algorithm/caller/checkpoint.hh>
#include <dune/fem-dg/algorithm/caller/datawriter.hh>
#include <dune/fem-dg/algorithm/caller/eocwriter.hh>
#include <dune/fem-dg/algorithm/caller/postprocessing.hh>
#include <dune/fem-dg/algorithm/caller/adapt.hh>

#include "containers.hh"

namespace Dune
{
namespace Fem
{


  //
  template< int polOrder, class ... ProblemTraits >
  struct SteadyStateTraits
  {
    // type of Grid
    typedef typename std::tuple_element_t<0, std::tuple< ProblemTraits... > >::GridType
                                                                           GridType;

    typedef std::tuple< typename std::add_lvalue_reference<typename ProblemTraits::GridType>::type ... >
                                                                           GridTypes;

    // wrap operator
    typedef GridTimeProvider< GridType >                                   TimeProviderType;

    typedef CreateSubAlgorithms< typename ProblemTraits::template Algorithm<polOrder>...  >
                                                                           CreateSubAlgorithmsType;

    typedef typename CreateSubAlgorithmsType::SubAlgorithmTupleType        SubAlgorithmTupleType;

    //typedef typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >
    //                                                                     IndexSequenceType;

    typedef Dune::Fem::SolverMonitorCaller< SubAlgorithmTupleType >        SolverMonitorCallerType;
    typedef Dune::Fem::DataWriterCaller< SubAlgorithmTupleType >           DataWriterCallerType;
    typedef Dune::Fem::EocWriterCaller< SubAlgorithmTupleType >            EocWriterCallerType;
    typedef Dune::Fem::AdaptCaller< SubAlgorithmTupleType >                AdaptCallerType;
  };


  /**
   * \brief An algorithm modelling a stationary partial differential equation.
   *
   * \ingroup Algorithms
   */
  template< int polOrder, class... ProblemTraits>
  class SteadyStateAlgorithm
    : public EOCAlgorithm< SteadyStateTraits< polOrder, ProblemTraits... > >
  {

    typedef SteadyStateTraits< polOrder, ProblemTraits... >             Traits;
    typedef EOCAlgorithm< Traits >                                      BaseType;
  public:
    typedef typename BaseType::GridType                                 GridType;
    typedef typename BaseType::SolverMonitorCallerType                  SolverMonitorCallerType;
    typedef typename Traits::DataWriterCallerType                       DataWriterCallerType;
    typedef typename Traits::EocWriterCallerType                        EocWriterCallerType;
    typedef typename Traits::AdaptCallerType                            AdaptCallerType;
    typedef typename Traits::SubAlgorithmTupleType                      SubAlgorithmTupleType;

    typedef typename Traits::CreateSubAlgorithmsType                    CreateSubAlgorithmsType;

    typedef typename BaseType::UInt64Type                               UInt64Type;

  private:
    struct Initialize {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->initialize( std::forward<Args>(a)... ); }
    };
    struct PreSolve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->preSolve( std::forward<Args>(a)... ); }
    };
    struct Solve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->solve( std::forward<Args>(a)... ); }
    };
    struct PostSolve {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->postSolve( std::forward<Args>(a)... ); }
    };
    struct Finalize {
      template< class T, class ... Args > static void apply ( T& e, Args && ... a )
      { e->finalize( std::forward<Args>(a)... ); }
    };
    struct GridWidth {
      template< class T, class ... Args > static void apply ( T& e, double& res, Args && ... a )
      { res = std::max( res, e->gridWidth(std::forward<Args>(a)... ) ); }
    };
    struct GridSize {
      template<class T, class... Args > static void apply ( T& e, UInt64Type& res, Args && ... a )
      { res = std::max( res, e->gridSize(std::forward<Args>(a)... ) ); }
    };
    struct CheckSolutionValid {
      template<class T, class... Args > static void apply( T e, bool& res, Args&& ... a )
      { res &= e->checkSolutionValid( std::forward<Args>(a)... ); }
    };
      template< class Caller >
    class LoopCallee
    {
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... Args >
        static void apply ( Tuple &tuple, Args&& ... a )
        {
          Caller::apply( std::get<i>( tuple ), std::forward<Args>(a)... );
        }
      };
    };

    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0,  sizeof ... ( ProblemTraits )-1 >;

    //static const int size = std::tuple_size< std::tuple< ProblemTraits... > >::value;
    //typedef std::make_index_sequence<size> sequence;

    using BaseType::tuple_;
  public:
    using BaseType::grid;

    template< class GlobalContainerImp >
    SteadyStateAlgorithm ( const std::string name, const std::shared_ptr<GlobalContainerImp>& cont )
    : BaseType( name, cont ),
      adaptCaller_( tuple_ )
    {}



    // return grid width of grid (overload in derived classes)
    virtual double gridWidth () const
    {
      double res=0.0;
      //Hybrid::forEach( sequence{}, [&](auto i){ res = std::max( std::get<i>( tuple_ )->gridWidth(), res ); } );
      ForLoopType< GridWidth >::apply( tuple_, res );
      return res;
    }

    // return size of grid
    virtual UInt64Type gridSize () const
    {
      UInt64Type res=0;
      ForLoopType< GridSize >::apply( tuple_, res );
      return res;
    }

    virtual void initialize ( const int loop )
    {
      ForLoopType< Initialize >::apply( tuple_, loop );
      adaptCaller_.initializeEnd( this, loop );
    }

    virtual void preSolve( const int loop )
    {
      ForLoopType< PreSolve >::apply( tuple_, loop );
    }

    virtual void eocSolve ( const int loop )
    {
      initialize( loop );
      preSolve( loop );
      //refine until grid is fine enough...
      //TODO improve?
      while( adaptCaller_.needsAdaptation( this, loop ) )
      {
        adaptCaller_.solveStart( this, loop );
      }
      ForLoopType< Solve >::apply( tuple_, loop );
      postSolve( loop );
      finalize( loop );
    }

    virtual void postSolve( const int loop )
    {
      ForLoopType< PostSolve >::apply( tuple_, loop );
    }

    void finalize ( const int loop )
    {
      adaptCaller_.finalizeStart( this, loop );
      ForLoopType< Finalize >::apply( tuple_, loop );
    }

  protected:
    AdaptCallerType                adaptCaller_;

  };

}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
