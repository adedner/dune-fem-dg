#ifndef FEMDG_DATAWRITERCALLER_HH
#define FEMDG_DATAWRITERCALLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{
  /**
   * \brief Caller class managing the data writing.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename std::make_index_sequence< std::tuple_size< AlgTupleImp >::value > >
  class DataWriterCaller;

  /**
   * \brief Specialization of a caller class managing the data writing.
   *
   * \ingroup Callers
   *
   * This class manages data writing for a tuple of sub-algorithms.
   * For each sub-algorithm data writing can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef DataWriterCaller< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                            std::index_sequence< 0, 2 > >
   *                                           MyCaller;
   * \endcode
   * This would enable data writing for `Alg1` and `Alg3`;
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam std::index_sequence< Ints... > Index sequence for enabling the data writing feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class DataWriterCaller< AlgTupleImp, std::index_sequence< Ints... > >
    : public CallerInterface
  {
    template< class TupleType > struct IOTupleExtractor;
    template< class ... Args > struct IOTupleExtractor< std::tuple< Args... > >
    { typedef typename tuple_concat< typename Args::element_type::IOTupleType::type... >::type type; };

    typedef AlgTupleImp                                                            AlgTupleType;

    typedef std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef typename std::tuple_element< 0, TupleType >::type::element_type::GridType
                                                                                    GridType;

  public:
    typedef typename IOTupleExtractor< TupleType >::type                            IOTupleType;

    typedef DataWriter< GridType, IOTupleType >                                     DataWriterType;

    template <class Grid, class DataTuple>
    struct DataOutput
    {
      typedef Dune::Fem::DataWriter< Grid, DataTuple > Type;
    };


  private:
    template< int i >
    struct AdditionalOutput
    {
      template<class T, class... Args >
      static typename std::enable_if< std::is_void< typename T::element_type::AdditionalOutputType >::value >::type
      additionalOutput( T&, Args&& ... ){}
      template<class T, class TimeProviderImp, class... Args >
      static typename std::enable_if< !std::is_void< typename T::element_type::AdditionalOutputType >::value >::type
      additionalOutput( T& elem, TimeProviderImp& tp, Args && ... args )
      {
        if( elem->additionalOutput() )
          elem->additionalOutput()->step( tp, *elem, args... );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        additionalOutput( std::get<i>( tuple ), std::forward<Args>(args)... );
      }
    };

    template< template< int > class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;

  public:

    /**
     * \brief Constructor.
     *
     * \param[in] tuple Tuple of all sub-algorithms.
     */
    DataWriterCaller( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) ),
        dataWriter_(),
        dataTuple_( dataTuple( tuple_, IndexSequenceType() ) )
    {
    }

    /**
     * \brief Creates the data writer.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void initializeEnd( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      dataWriter_.reset( new DataWriterType( std::get<0>(tuple_)->solution().space().grid(), dataTuple_, tp,
                                             alg->eocParams().dataOutputParameters( loop, alg->dataPrefix() ) ) );
    }

    /**
     * \brief Write current solution and additional output data to disk.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void preSolveStart( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      finalizeStart( alg, loop, tp );
    }

    /**
     * \brief Write current solution and additional output data to disk.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void postSolveEnd( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      // Check that no NAN have been generated
      if( !alg->checkSolutionValid( loop, tp ) )
      {
        finalizeStart( alg, loop, tp );
        std::cerr << "Solution is not valid. Aborting." << std::endl;
        std::abort();
      }
    }

    /**
     * \brief Write current, final solution and additional output data to disk.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void finalizeStart( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        //update all additional Output
        ForLoopType< AdditionalOutput >::apply( tuple_, tp );
        //writeData
        dataWriter_->write( tp );
      }
    }

    /**
     * \brief Returns a tuple of pointer to all discrete functions that should be
     * written to disk.
     */
    IOTupleType dataTuple()
    {
      return dataTuple_;
    }

  private:
    template< std::size_t ... i >
    IOTupleType dataTuple ( const TupleType &tuple, std::index_sequence< i ... > )
    {
      return std::tuple_cat( (*std::get< i >( tuple )->dataTuple() )... );
    }

    TupleType                         tuple_;
    std::unique_ptr< DataWriterType > dataWriter_;
    IOTupleType                       dataTuple_;
  };


  /**
   * \brief Caller class doing no the data writing.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp >
  class DataWriterCaller< AlgTupleImp, std::index_sequence<> >
    : public CallerInterface
  {
  public:
    template <class Grid, class DataTuple>
    struct DataOutput
    {
      typedef Dune::Fem::DataOutput< Grid, DataTuple > Type;
    };

    template< class ... Args >
    DataWriterCaller( Args&& ... ) {}

  };

}
}
#endif
