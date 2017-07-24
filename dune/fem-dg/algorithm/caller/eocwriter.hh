#ifndef FEMDG_EOCWRITERCALLER_HH
#define FEMDG_EOCWRITERCALLER_HH

#include <memory>
#include <tuple>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/femeoc.hh>
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
  class EocWriterCaller;

  /**
   * \brief Specialization of a caller class managing the data writing.
   *
   * \ingroup Callers
   *
   * This class manages eoc writing for a tuple of sub-algorithms.
   * For each sub-algorithm eoc data writing can be disabled using an `index_sequence`.
   *
   * Example:
   * \code
   * typedef EocWriterCaller< std::tuple< Alg1, Alg2, Alg3, Alg4 >,
   *                            std::index_sequence< 0, 2 > >
   *                                           MyCaller;
   * \endcode
   * This would enable data writing for `Alg1` and `Alg3`;
   *
   * \tparam AlgTupleImp A tuple of all known sub-algorithms.
   * \tparam std::index_sequence< Ints... > Index sequence for enabling the data writing feature.
   */
  template< class AlgTupleImp, std::size_t... Ints >
  class EocWriterCaller< AlgTupleImp, std::index_sequence< Ints... > >
    : public CallerInterface
  {
    typedef AlgTupleImp                                                            AlgTupleType;

    typedef std::index_sequence< Ints... >                                         IndexSequenceType;
    static const int numAlgs = IndexSequenceType::size();
    typedef tuple_reducer<AlgTupleType, IndexSequenceType >                        TupleReducerType;
    typedef typename TupleReducerType::type                                        TupleType;

    static_assert( std::tuple_size< TupleType >::value>=1, "Empty Tuples not allowed..." );
  public:

    /**
     * \brief Constructor.
     *
     * \param[in] tuple Tuple of all sub-algorithms.
     */
    EocWriterCaller( const AlgTupleType& tuple )
    {}

    template< class SubAlgImp >
    void eocInitializeStart( SubAlgImp* alg )
    {
      FemEoc::clear();
      FemEoc::initialize( alg->eocParams().outputPath(), alg->eocParams().fileName(), "results" );
    }

    template< class SubAlgImp >
    void eocPreSolveEnd( SubAlgImp* alg, int loop )
    {
      timer_.reset();
    }


    template< class SubAlgImp >
    void eocPostSolveEnd( SubAlgImp* alg, int loop )
    {
      const double runTime = timer_.elapsed();

      std::stringstream eocInfo;

      // generate EOC information
      writeFemEoc( alg->monitor(), runTime, eocInfo );

      // in verbose mode write EOC info to std::cout
      if( Fem::Parameter::verbose() )
      {
        std::cout << std::endl << "EOC info: " << std::endl << eocInfo.str() << std::endl;
      }
    }

  private:
    //! add solver Monitor data to Fem Eoc
    template< class Monitor >
    void writeFemEoc ( const Monitor &m, const double runTime, std::stringstream &out )
    {
      Fem::FemEoc::write( m.getData( "GridWidth" ), m.getData( "Elements" ), runTime, m.getData( "TimeSteps" ),
                          m.getData( "AvgTimeStep" ), m.getData( "MinTimeStep" ), m.getData( "MaxTimeStep" ),
                          m.getData( "Newton" ), m.getData( "ILS" ), m.getData( "MaxNewton" ), m.getData( "MaxILS" ), out );
    }
    Dune::Timer                    timer_;

  };


  /**
   * \brief Caller class doing no the data writing.
   *
   * \ingroup Callers
   */
  template< class AlgTupleImp >
  class EocWriterCaller< AlgTupleImp, std::index_sequence<> >
    : public CallerInterface
  {
  public:
    template< class ... Args >
    EocWriterCaller( Args&& ... ) {}

  };

}
}
#endif
