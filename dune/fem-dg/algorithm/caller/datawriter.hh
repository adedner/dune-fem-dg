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
    { typedef tuple_concat_t< typename Args::element_type::DataWriterType::IOTupleType... > type; };

    typedef AlgTupleImp                                                            AlgTupleType;

    typedef std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef typename std::tuple_element_t< 0, TupleType >::element_type::GridType  GridType;

  public:
    typedef typename IOTupleExtractor< TupleType >::type                            IOTupleType;

    typedef DataWriter< GridType, IOTupleType >                                     DataWriterType;

    typedef Dune::Fem::DataOutput< GridType, IOTupleType >                         DataOutputType;

  private:
    template< int i >
    struct DataWriterOutput
    {
      template<class T, class... Args >
      static typename std::enable_if< std::is_void< typename T::element_type::DataWriterType >::value >::type
      dataWriter( T&, Args&& ... ){}
      template<class T, class TimeProviderImp, class... Args >
      static typename std::enable_if< !std::is_void< typename T::element_type::DataWriterType >::value >::type
      dataWriter( T& elem, TimeProviderImp& tp, Args && ... args )
      {
        if( elem->dataWriter() )
          elem->dataWriter()->prepare( tp, elem /*args...*/ );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        dataWriter( std::get<i>( tuple ), std::forward<Args>(args)... );
      }
    };

    template< int i >
    struct Init
    {
      template<class T, class... Args >
      static typename std::enable_if< std::is_void< typename T::element_type::DataWriterType >::value >::type
      dataWriter( T&, Args&& ... ){}
      template<class T, class... Args >
      static typename std::enable_if< !std::is_void< typename T::element_type::DataWriterType >::value >::type
      dataWriter( T& elem, Args && ... args )
      {
        if( elem->dataWriter() )
          elem->dataWriter()->init( elem /*std::forward<Args>(args)...*/ );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        dataWriter( std::get<i>( tuple ), std::forward<Args>(args)... );
      }
    };

    template< int i >
    struct Write
    {
      template<class T, class... Args >
      static typename std::enable_if< std::is_void< typename T::element_type::DataWriterType >::value >::type
      dataWriter( T&, Args&& ... ){}
      template<class T, class TimeProviderImp, class... Args >
      static typename std::enable_if< !std::is_void< typename T::element_type::DataWriterType >::value >::type
      dataWriter( T& elem, TimeProviderImp& tp, Args && ... args )
      {
        if( elem->dataWriter() )
          elem->dataWriter()->write( tp, elem /*args...*/ );
      }

      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        dataWriter( std::get<i>( tuple ), std::forward<Args>(args)... );
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
      dataOutput_(),
        dataWriter_(),
      ioTuple_()
    {}

    template< class AlgImp >
    void eocInitializeStart( AlgImp* alg )
    {
      ForLoopType< Init >::apply( tuple_, alg );

      ioTuple_ = dataTuple( tuple_, IndexSequenceType() );

      dataOutput_ = std::make_unique<DataWriterType>( std::get<0>(tuple_)->solution().space().grid(), ioTuple_ );
    }

    template< class AlgImp >
    void eocPostSolveEnd( AlgImp* alg, int loop )
    {
      dataOutput_->writeData( loop );
    }


    /**
     * \brief Creates the data writer.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void initializeEnd( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      dataWriter_ = std::make_unique<DataWriterType>( std::get<0>(tuple_)->solution().space().grid(), ioTuple_, tp,
                                                      alg->eocParams().dataOutputParameters( loop ) );
    }

    /**
     * \brief Write current solution and additional output data to disk.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void preSolveStart( AlgImp* alg, int loop, TimeProviderImp& tp )
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
    template< class AlgImp, class TimeProviderImp >
    void postSolveEnd( AlgImp* alg, int loop, TimeProviderImp& tp )
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
    template< class AlgImp, class TimeProviderImp >
    void finalizeStart( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( dataWriter_ && dataWriter_->willWrite( tp ) )
      {
        //update all additional Output
        ForLoopType< DataWriterOutput >::apply( tuple_, tp, alg );
        //writeData
        dataWriter_->write( tp );

        //write everthing a sub data writer wants to write...
        ForLoopType< Write >::apply( tuple_, tp, alg );
      }
    }

  private:
    template< std::size_t ... i >
    IOTupleType dataTuple ( const TupleType &tuple, std::index_sequence< i ... > )
    {
      return std::tuple_cat( (std::get< i >( tuple )->dataWriter()->dataTuple() )... );
    }

    TupleType                         tuple_;
    std::unique_ptr< DataOutputType > dataOutput_;
    std::unique_ptr< DataWriterType > dataWriter_;
    IOTupleType                       ioTuple_;
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
    template< class ... Args >
    DataWriterCaller( Args&& ... ) {}

  };

}
}
#endif
