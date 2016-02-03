#ifndef FEMDG_CHECKPOINTHANDLER_HH
#define FEMDG_CHECKPOINTHANDLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Helper class for the CheckPointHandler containing some static methods.
   *
   * This class is needed by the CheckPointHandler. All methods in this class
   * have to be static because the grid is constructed from the checkpoint, here
   * which (obviuously) has to be done before the construction of algorithms.
   */
  template< class GridImp >
  class GridCheckPointHandler
  {
  public:
    typedef GridImp                                            GridType;
    typedef Dune::Fem::CheckPointer< GridType, std::tuple<> >  CheckPointerType;
    typedef Dune::Fem::CheckPointerParameters                  CheckPointerParametersType;

    /**
     * \brief Returns true if checkpoint file exists.
     */
    static bool checkPointExists( const std::string keyPrefix = "" )
    {
      return (checkPointFile(keyPrefix).size() > 0 );
    }

    /**
     * \brief Restores the grid.
     */
    static Dune::GridPtr< GridType > restoreGrid( const std::string keyPrefix = "" )
    {
      if( 0 == Dune::Fem::MPIManager::rank () )
      {
        std::cout << std::endl;
        std::cout << "********************************************************" << std::endl;
        std::cout << "**   Restart from checkpoint `" << checkPointFile() << "' " << std::endl;
        std::cout << "********************************************************" << std::endl;
        std::cout << std::endl;
      }
      Dune::GridPtr< GridType > gridptr ;
      gridptr = CheckPointerType::restoreGrid( checkPointFile( keyPrefix ), -1,
                                               CheckPointerParametersType( ParameterKey::generate( keyPrefix, "fem.io." ) ) );
      return gridptr;
    }

    /**
     * \brief Returns the name of the checkpoint file.
     *
     * \param[in] keyPrefix an additional key prefix
     */
    static inline std::string checkPointFile ( const std::string keyPrefix = "" )
    {
      static std::string checkFileName ;
      const std::string key( ParameterKey::generate( keyPrefix, "fem.io.checkpointrestartfile" ) );
      if( Fem::Parameter::exists( key ) )
        checkFileName = Fem::Parameter::getValue<std::string> ( key );
      else
        checkFileName = "";
      return checkFileName;
    }


  };


  /**
   * \brief Handler class managing the checkpointing.
   *
   * \ingroup Handlers
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename Std::make_index_sequence_impl< std::tuple_size< AlgTupleImp >::value >::type >
  class CheckPointHandler;

  /**
   * \brief Specialization of a handler class managing the checkpointing.
   *
   * \ingroup Handlers
   *
   * This class manages checkpointing for a tuple of sub-algorithms.
   * For each sub-algorithm checkpointing can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef CheckPointHandler< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                            Std::index_sequence< 0, 2 > >
   *                                           MyHandler;
   * \endcode
   * This would enable checkpointing for `Alg1` and `Alg3`;
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam Std::index_sequence< Ints... > Index sequence for enabling the checkpointing feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class CheckPointHandler< AlgTupleImp, Std::index_sequence< Ints... > >
    : public HandlerInterface
  {

    template< int i >
    struct RegisterData
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        if( std::get< i >( tuple )->checkPointSolution() )
          Dune::Fem::persistenceManager << *(std::get< i >( tuple )->checkPointSolution());
      }
    };

  public:
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef Std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef GridCheckPointHandler< typename std::remove_pointer< typename std::tuple_element< 0, TupleType >::type >::type::GridType >
                                                                                BaseType;

    typedef typename BaseType::GridType                                         GridType;
    typedef typename BaseType::CheckPointerType                                 CheckPointerType;
    typedef typename BaseType::CheckPointerParametersType                       CheckPointerParametersType;

    template< template<int> class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;

  public:

    /**
     * \brief Constructor.
     *
     * \param[in] tuple The whole tuple of all sub-algorithms
     */
    CheckPointHandler( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) ),
        checkPointer_(),
        keyPrefix_( std::get<0>( tuple_ )->name() ),
        checkParam_( ParameterKey::generate( keyPrefix_, "fem.io." ) )
    {}

    /**
     * \brief Restore data from a checkpoint file.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    bool initialize_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( BaseType::checkPointExists(keyPrefix_) )
      {
        // restore data
        checkPointer( tp ).restoreData( std::get<0>(tuple_)->solution().space().grid(), BaseType::checkPointFile( keyPrefix_ ) );
        return false;
      }
      return true;
    }

    /**
     * \brief Register data.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void initialize_post( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( BaseType::checkPointExists(keyPrefix_) )
        ForLoopType< RegisterData >::apply( tuple_ );
    }

    /**
     * \brief Write checkpoint data to disk
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class SubAlgImp, class TimeProviderImp >
    void preSolve_pre( SubAlgImp* alg, int loop, TimeProviderImp& tp )
    {
      checkPointer( tp ).write( tp );
    }

    /**
     * \brief Returns the checkpointer.
     *
     * \param[in] tp the time provider
     */
    template< class TimeProviderImp >
    CheckPointerType& checkPointer( const TimeProviderImp& tp  ) const
    {
      // create check point if not existent
      if( ! checkPointer_ )
        checkPointer_.reset( new CheckPointerType( std::get<0>(tuple_)->solution().space().grid(), tp, checkParam_) );
      return *checkPointer_;
    }

  protected:
    TupleType               tuple_;
    mutable std::unique_ptr< CheckPointerType > checkPointer_;
    const std::string keyPrefix_;
    CheckPointerParametersType checkParam_;
  };


  /**
   * \brief Specialization of a handler class without checkpointing.
   *
   * \ingroup Handlers
   */
  template< class AlgTupleImp >
  class CheckPointHandler< AlgTupleImp, Std::index_sequence<> >
    : public HandlerInterface
  {
    public:

    template< class ... Args >
    CheckPointHandler( Args&& ... ) {}

    template< class ... Args>
    static bool checkPointExists( Args&& ... ) {return false;}

    template< class GridImp, class ... Args>
    static Dune::GridPtr< GridImp > restoreGrid( Args&& ... )
    {
      Dune::GridPtr< GridImp > gridptr;
      return gridptr;
    }

  };

}
}
#endif
