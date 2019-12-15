#ifndef FEMDG_ALGORITHM_SOLVERMONITORCALLER_HH
#define FEMDG_ALGORITHM_SOLVERMONITORCALLER_HH

#include <string>
#include <dune/fem/common/forloop.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Caller class to collect all data from solver monitor from sub-algorithms.
   *
   * \ingroup Caller
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename std::make_index_sequence< std::tuple_size< AlgTupleImp >::value > >
  class SolverMonitorCaller;


  /**
   * \brief Specialization of a caller class collecting all data from solver monitor from sub-algorithms.
   *
   * \ingroup Callers
   *
   * This class manages data collecting from solver monitor for a tuple of sub-algorithms.
   * For each sub-algorithm data collecting can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef DataWriterCaller< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                            std::index_sequence< 0, 2 > >
   *                                           MyCaller;
   * \endcode
   * This would enable data collecting for `Alg1` and `Alg3`;
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam std::index_sequence< Ints... > Index sequence for enabling the data collecting feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class SolverMonitorCaller< AlgTupleImp, std::index_sequence< Ints... > >
    : public CallerInterface
  {
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    enum CombinationType { max, min, sum, avg };

    template< class Caller >
    class LoopCallee
    {
      template<class C, class T, class... A >
      static typename std::enable_if< std::is_void< typename T::element_type::SolverMonitorType >::value >::type
      getMonitor( T&, A&& ... ){}
      template<class C, class T, class... A >
      static typename std::enable_if< !std::is_void< typename T::element_type::SolverMonitorType >::value >::type
      getMonitor( T& elem, A &&... a )
      {
        if( elem->monitor() )
          C::apply(elem->monitor(), std::forward<A>(a)... );
      }
    public:
      template< int i >
      struct Apply
      {
        template< class Tuple, class ... A >
        static void apply ( Tuple &tuple, A&& ... a )
        { getMonitor< Caller >( std::get<i>( tuple ), std::forward<A>(a)... ); }
      };
    };

    struct Step {
      template<class T, class... A > static void apply( T e, A&& ... a )
      { e->step( std::forward<A>(a)... ); }
    };

    struct Finalize {
      template<class T, class... A > static void apply( T e, A&& ... a )
      { e->finalize( std::forward<A>(a)... ); }
    };

    struct GetData {
      template<class T, class... A > static void apply( T e, double& res, const std::string name, CombinationType comb, A&& ... a )
      {
        switch( comb )
        {
          case CombinationType::max:
            res = std::max( res, e->getData( name ) ); break;
          case CombinationType::min:
            res = std::min( res, e->getData( name ) ); break;
          case CombinationType::sum:
            res += e->getData( name ); break;
          case CombinationType::avg:
            res += e->getData( name ); break;
        }
      }
    };

    template< class Caller >
    using ForLoopType = ForLoop< LoopCallee<Caller>::template Apply, 0, numAlgs - 1 >;

  public:

    /**
     * \brief Constructor.
     *
     * \param[in] tuple Tuple of all sub-algorithms.
     */
    SolverMonitorCaller( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) )
    {}

    /**
     * \brief Collects solver monitor information.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void postSolveEnd( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Step >::apply( tuple_, tp );
    }

    /**
     * \brief Finalizes the collection of solver monitor information.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void finalizeStart( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      ForLoopType< Finalize >::apply( tuple_, alg->gridWidth(), alg->gridSize() );
    }

    double getData( const std::string name, CombinationType comb = CombinationType::max ) const
    {
      double res = (comb == CombinationType::min) ? std::numeric_limits<double>::max() : 0.0;
      ForLoopType< GetData >::apply( tuple_, res, name, comb );
      return res;
    }

  private:
    TupleType  tuple_;
  };

  /**
   * \brief Specizalization of hander class doing no solver monitor handling.
   */
  template< class AlgTupleImp >
  class SolverMonitorCaller< AlgTupleImp, std::index_sequence<> >
    : public CallerInterface
  {
  public:
    template <class ... Args>
    SolverMonitorCaller(Args&& ... )
    {}

    template <class ... Args>
    double getData( Args&& ... ) const
    { return 0.0; }

  };
}
}

#endif
