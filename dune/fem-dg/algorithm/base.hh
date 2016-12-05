#ifndef DUNE_FEMDG_ALGORITHM_BASE_HH
#define DUNE_FEMDG_ALGORITHM_BASE_HH

// system includes
#include <sstream>
#include <sys/resource.h>

// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// local includes
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/misc/femeoc.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Parameter class describing where to write eoc data.
   *
   * \ingroup ParameterClass
   */
  struct EocDataOutputParameters :
         public Fem::LocalParameter<Fem::DataWriterParameters,EocDataOutputParameters>
  {
    protected:
    std::string loop_;
    public:
    EocDataOutputParameters(int loop, const std::string& name)
    {
      std::stringstream ss;
      ss << name << loop;
      loop_ = ss.str();
    }
    EocDataOutputParameters(const EocDataOutputParameters& other)
    : loop_(other.loop_)
    {}

    std::string path() const
    {
      return loop_;
    }
  };

  /**
   * \brief Parameter class describing information about eoc steps.
   *
   * \ingroup ParameterClass
   */
  struct EocParameters :
         public Fem::LocalParameter<EocParameters,EocParameters>
  {
    protected:
    std::string keyPrefix_;

    public:
    EocParameters( const std::string keyPrefix = "fem.eoc.")
      : keyPrefix_( keyPrefix )
    {}

    /**
     * \brief returns the EocDataOutputParameters
     *
     * \param[in] loop the current number of the eoc loop
     * \param[in] name an additional name
     */
    const EocDataOutputParameters dataOutputParameters( int loop, const std::string& name ) const
    {
      return EocDataOutputParameters( loop, name );
    }

    /**
     * \brief returns the number of eoc steps
     */
    virtual int steps() const
    {
      return Fem::Parameter::getValue<int>( keyPrefix_ + "steps", 1);
    }

    /**
     * \brief returns the output path where eoc data should be written to
     */
    virtual std::string outputPath() const
    {
      return Fem::Parameter::getValue<std::string>( keyPrefix_ + "outputpath", "./");
    }

    /**
     * \brief returns the name of the file for the eoc data.
     */
    virtual std::string fileName() const
    {
      return Fem::Parameter::getValue<std::string>( keyPrefix_ + "filename", "eoc" );
    }

    /**
     * \brief returns the quadrature order
     *
     * \deprecated{This method is not needed here}
     */
    virtual int quadOrder() const
    {
      return Fem::Parameter::getValue< int >( keyPrefix_ + "quadorder", -1 );
    }

  };


  //! add solver Monitor data to Fem Eoc
  template< class Monitor >
  void writeFemEoc ( const Monitor &m, const double runTime, std::stringstream &out )
  {
    Fem::FemEoc::write( m.getData( "GridWidth" ), m.getData( "Elements" ), runTime, m.getData( "TimeSteps" ),
                        m.getData( "AvgTimeStep" ), m.getData( "MinTimeStep" ), m.getData( "MaxTimeStep" ),
                        m.getData( "Newton" ), m.getData( "ILS" ), m.getData( "MaxNewton" ), m.getData( "MaxILS" ), out );
  }


  //! get memory in MB
  inline double getMemoryUsage()
  {
    struct rusage info;
    getrusage( RUSAGE_SELF, &info );
    return (info.ru_maxrss / 1024.0);
  }


  // AlgorithmInterface
  // -------------

  /**
   *  \brief Interface class for all algorithms.
   *
   *  \ingroup Algorithms
   */
  template< class AlgorithmTraits >
  class AlgorithmInterface
  {
    typedef AlgorithmTraits Traits;

  public:
    // type of Grid
    typedef typename Traits::GridType                 GridType;
    typedef typename Traits::GridTypes                GridTypes;

    // type of IOTuple
    typedef typename Traits::IOTupleType              IOTupleType;

    // type of statistics monitor
    typedef typename Traits::SolverMonitorCallerType  SolverMonitorCallerType;

    typedef EocParameters                             EocParametersType;

    //! constructor
    AlgorithmInterface ( const std::string algorithmName, const GridTypes& grids )
      : grids_( grids ),
        algorithmName_( algorithmName ),
        eocParam_( ParameterKey::generate( "", "fem.eoc." ) )
    {}

    //! destructor
    virtual~ AlgorithmInterface()
    {}

    //! return reference to hierarchical grid
    GridType &grid () const { return std::get<0>( grids_ ); }

    // return size of grid
    virtual std::size_t gridSize () const
    {
      std::size_t grSize = std::get<0>(grids_).size( 0 );
      return std::get<0>(grids_).comm().sum( grSize );
    }

    //! return default data tuple for data output
    virtual IOTupleType dataTuple() = 0;

    // solve the problem for eoc loop 'loop'
    virtual void solve ( const int loop ) = 0;

    virtual SolverMonitorCallerType& monitor() = 0;

    EocParametersType& eocParams() { return eocParam_; }

    virtual const std::string name() { return algorithmName_; }

  private:
    GridTypes grids_;
    const std::string algorithmName_;
    EocParameters eocParam_;
  };



  /////////////////////////////////////////////////////////////////////////////
  //
  //  compute algorithm
  //
  /////////////////////////////////////////////////////////////////////////////
  template <class Algorithm>
  void compute(Algorithm& algorithm)
  {
    typedef typename Algorithm::GridType             GridType;
    typedef typename Algorithm::EocParametersType    EocParametersType;

    auto& grid = algorithm.grid();
    auto dataTup = algorithm.dataTuple() ;

    typedef typename Algorithm::DataWriterCallerType::template DataOutput< decltype( grid ), decltype( dataTup ) > ::Type DataOutputType;
    DataOutputType dataOutput( grid, dataTup );

    // initialize FemEoc if eocSteps > 1
    EocParametersType& eocParam( algorithm.eocParams() );
    FemEoc::clear();
    FemEoc::initialize( eocParam.outputPath(), eocParam.fileName(), "results" );

    const unsigned int femTimerId = FemTimer::addTo("timestep");
    for(int eocloop=0; eocloop < eocParam.steps(); ++eocloop )
    {
      // start fem timer (enable with -DFEMTIMER)
      FemTimer :: start(femTimerId);

      // additional timer
      Dune::Timer timer ;

      // call algorithm and return solver statistics and some info
      algorithm.solve( eocloop );

      // get run time
      const double runTime = timer.elapsed();

      // also get times for FemTimer if enabled
      FemTimer::stop(femTimerId);
      FemTimer::printFile("./timer.out");
      FemTimer::reset(femTimerId);

      std::stringstream eocInfo ;

      // generate EOC information
      writeFemEoc( algorithm.monitor(), runTime, eocInfo );

      // in verbose mode write EOC info to std::cout
      if( Fem::Parameter::verbose() )
      {
        std::cout << std::endl << "EOC info: " << std::endl << eocInfo.str() << std::endl;
      }

      // write eoc step
      dataOutput.writeData( eocloop );

      // Refine the grid for the next EOC Step. If the scheme uses adaptation,
      // the refinement level needs to be set in the algorithms' initialize method.
      if( eocParam.steps() > 1 )
        Fem::GlobalRefine::apply(grid,Dune::DGFGridInfo<GridType>::refineStepsForHalf());

    } /***** END of EOC Loop *****/

    // FemTimer cleanup
    FemTimer::removeAll();
  }




} // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_BASE_HH
