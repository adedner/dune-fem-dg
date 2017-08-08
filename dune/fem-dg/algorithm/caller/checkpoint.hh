#ifndef FEMDG_CHECKPOINTCALLER_HH
#define FEMDG_CHECKPOINTCALLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/misc/gridsolution.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include "interface.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Helper class for the CheckPointCaller containing some static methods.
   *
   * This class is needed by the CheckPointCaller. All methods in this class
   * have to be static because the grid is constructed from the checkpoint, here
   * which (obviuously) has to be done before the construction of algorithms.
   */
  template< class GridImp >
  class GridCheckPointCaller
  {
  public:
    typedef GridImp                                GridType;
    typedef Dune::Fem::CheckPointer< GridType >    CheckPointerType;
    typedef Dune::Fem::CheckPointerParameters      CheckPointerParametersType;

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
   * \brief Caller class managing the checkpointing.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp,
            class IndexSequenceImp=typename std::make_index_sequence< std::tuple_size< AlgTupleImp >::value > >
  class CheckPointCaller;

  /**
   * \brief Specialization of a caller class managing the checkpointing.
   *
   * \ingroup Callers
   *
   * This class manages checkpointing for a tuple of sub-algorithms.
   * For each sub-algorithm checkpointing can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef CheckPointCaller< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                            std::index_sequence< 0, 2 > >
   *                                           MyCaller;
   * \endcode
   * This would enable checkpointing for `Alg1` and `Alg3`;
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam std::index_sequence< Ints... > Index sequence for enabling the checkpointing feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class CheckPointCaller< AlgTupleImp, std::index_sequence< Ints... > >
    : public CallerInterface
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

    template< int i >
    struct WriteGridSolution
    {
      template< class Tuple, class TimeProvider, class ... Args >
      static void apply ( Tuple &tuple, const TimeProvider& tp, int counter, Args && ... args )
      {
        typedef Dune::Fem::GridSolutionVector< GridType, std::decay_t<decltype( std::get<i>(tuple_)->solution() )> > ExSolGrid;
        ExSolGrid::writeDiscreteFunction( std::get<i>(tuple)->solution().space().grid(), std::get<i>(tuple)->solution(), tp.time(), counter );
      }
    };
  public:
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );

    typedef GridCheckPointCaller< typename std::tuple_element_t< 0, TupleType >::element_type::GridType >
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
    CheckPointCaller( const AlgTupleType& tuple )
      : tuple_( TupleReducerType::apply( tuple ) ),
        checkPointer_(),
        checkPointRestored_( false ),
        keyPrefix_( std::get<0>( tuple_ )->name() ),
        checkParam_( ParameterKey::generate( keyPrefix_, "fem.io." ) ),
        writeGridSolution_( Fem::Parameter::getValue<bool>("gridsol.writesolution", false) ),
        saveStep_( Fem::Parameter::getValue< double >("gridsol.firstwrite", 0.0 ) ),
        saveInterval_( Fem::Parameter::getValue< double >("gridsol.savestep", 0.0 ) ),
        writeCounter_( 0 )
    {}

    /**
     * \brief Restore data from a checkpoint file.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void initializeStart( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      if( alg->eocParams().steps() == 1 && BaseType::checkPointExists(keyPrefix_) )
      {
        // restore data
        checkPointer( tp ).restoreData( std::get<0>(tuple_)->solution().space().grid(), BaseType::checkPointFile( keyPrefix_ ) );
        checkPointRestored_ = true;
      }
      checkPointRestored_ = false;
    }

    /**
     * \brief Register data.
     *
     * \param[in] alg pointer to the calling sub-algorithm
     * \param[in] loop number of eoc loop
     * \param[in] tp the time provider
     */
    template< class AlgImp, class TimeProviderImp >
    void initializeEnd( AlgImp* alg, int loop, TimeProviderImp& tp )
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
    template< class AlgImp, class TimeProviderImp >
    void preSolveStart( AlgImp* alg, int loop, TimeProviderImp& tp )
    {
      checkPointer( tp ).write( tp );
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
      if( writeGridSolution_ && tp.time() > saveStep_ )
      {
        ForLoopType< WriteGridSolution >::apply( tuple_, tp, writeCounter_ );
        //increase counter
        ++writeCounter_;
        saveStep_ += saveInterval_;

        ////later usage for eoc calculation against discrete solution
        //static const int rank = 0;
        //ExSolGrid gridSol;
        //gridsol.discreteFunction(rank);
      }
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

    bool checkPointRestored()
    {
      return checkPointRestored_;
    }

  protected:
    TupleType               tuple_;
    mutable std::unique_ptr< CheckPointerType > checkPointer_;
    bool checkPointRestored_;
    const std::string keyPrefix_;
    CheckPointerParametersType checkParam_;
    //grid solution writer...
    const bool                        writeGridSolution_;
    mutable double                    saveStep_ ;
    const double                      saveInterval_ ;
    mutable int                       writeCounter_ ;
  };


  /**
   * \brief Specialization of a caller class without checkpointing.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp >
  class CheckPointCaller< AlgTupleImp, std::index_sequence<> >
    : public CallerInterface
  {
    public:

    template< class ... Args >
    CheckPointCaller( Args&& ... ) {}

    template< class ... Args>
    static bool checkPointRestored( Args&& ... ) {return false;}

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
