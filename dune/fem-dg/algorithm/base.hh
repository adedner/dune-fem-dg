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
  struct EocDataOutputParameters :
         public Fem::LocalParameter<Fem::DataWriterParameters,EocDataOutputParameters>
  {
    protected:
    std::string loop_;
    public:
    EocDataOutputParameters(int loop, const std::string& name) {
      std::stringstream ss;
      ss << name << loop;
      loop_ = ss.str();
    }
    EocDataOutputParameters(const EocDataOutputParameters& other)
    : loop_(other.loop_) {}

    std::string path() const {
      return loop_;
    }
  };

  struct EocParameters :
         public Fem::LocalParameter<EocParameters,EocParameters>
  {
    protected:
    std::string keyPrefix_;

    public:
    EocParameters( const std::string keyPrefix = "fem.eoc.")
      : keyPrefix_( keyPrefix )
    {}

    const EocDataOutputParameters dataOutputParameters( int loop, const std::string& name ) const
    {
      return EocDataOutputParameters( loop, name );
    }

    virtual int steps() const
    {
      return Fem::Parameter::getValue<int>( keyPrefix_ + "steps", 1);
    }

    virtual std::string outputPath() const
    {
      return Fem::Parameter::getValue<std::string>( keyPrefix_ + "outputpath", "./");
    }

    virtual std::string fileName() const
    {
      return Fem::Parameter::getValue<std::string>( keyPrefix_ + "filename", "eoc" );
    }

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


  // AlgorithmBase
  // -------------

  template< class AlgorithmTraits >
  class AlgorithmBase
  {
    typedef AlgorithmTraits Traits;

  public:
    // type of Grid
    typedef typename Traits::GridType                 GridType;

    // type of IOTuple
    typedef typename Traits::IOTupleType              IOTupleType;

    // type of statistics monitor
    typedef typename Traits::SolverMonitorHandlerType SolverMonitorHandlerType;

    typedef EocParameters                             EocParametersType;

    //! constructor
    AlgorithmBase ( GridType &grid, const std::string algorithmName = "" )
      : grid_( grid ),
        algorithmName_( algorithmName ),
        eocParam_( Dune::ParameterKey::generate( "", "fem.eoc." ) )
    {}

    virtual ~AlgorithmBase () {}

    //! return reference to hierarchical grid
    GridType &grid () const { return grid_; }

    // return size of grid
    virtual std::size_t gridSize () const
    {
      std::size_t grSize = grid_.size( 0 );
      return grid_.comm().sum( grSize );
    }

    //! return default data tuple for data output
    virtual IOTupleType dataTuple() = 0;

    // solve the problem for eoc loop 'loop'
    virtual void solve ( const int loop ) = 0;

    virtual SolverMonitorHandlerType& monitor() = 0;

    EocParametersType& eocParams() { return eocParam_; }

    //! name of the scheme
    virtual std::string description () const = 0;

    virtual const std::string name() { return algorithmName_; }

  protected:
    GridType &grid_;
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
    typedef typename Algorithm::IOTupleType          IOTupleType;

    GridType& grid = algorithm.grid();
    IOTupleType dataTup = algorithm.dataTuple() ;
    Fem::DataOutput<GridType, IOTupleType > dataOutput( grid, dataTup );

    // initialize FemEoc if eocSteps > 1
    EocParametersType& eocParam( algorithm.eocParams() );
    FemEoc::clear();
    FemEoc::initialize( eocParam.outputPath(), eocParam.fileName(), algorithm.description() );

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
