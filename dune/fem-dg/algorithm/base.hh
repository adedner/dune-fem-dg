#ifndef DUNE_FEMDG_ALGORITHM_BASE_HH
#define DUNE_FEMDG_ALGORITHM_BASE_HH

// system includes
#include <sstream>
#include <sys/resource.h>

// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/adaptationmanager.hh>

// local includes
#include <dune/fem-dg/algorithm/monitor.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem/common/utility.hh>
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
    EocDataOutputParameters(int loop, const std::string& name = "" )
    {
      std::stringstream ss;
      ss << (name == "" ? "loop" : name) << loop;
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
    const EocDataOutputParameters dataOutputParameters( int loop, const std::string& name = "" ) const
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



  //! get memory in MB
  inline double getMemoryUsage()
  {
    struct rusage info;
    getrusage( RUSAGE_SELF, &info );
    return (info.ru_maxrss / 1024.0);
  }

  /**
   *  \brief Interface class for all algorithms.
   *
   *  \ingroup Algorithms
   */
  class EOCAlgorithmInterface
  {
  public:
    //! constructor
    EOCAlgorithmInterface ( )
    {}

    //! destructor
    virtual~ EOCAlgorithmInterface()
    {}

    virtual void eocInitialize ()
    {}

    virtual void eocPreSolve ( const int loop )
    {}

    // solve the problem for eoc loop 'loop'
    virtual void eocSolve ( const int loop )
    {}

    virtual void eocPostSolve ( const int loop )
    {}

    virtual void eocFinalize()
    {}

    int eocSteps() const
    {
      return 0;
    }
  };



  // EOCAlgorithm
  // -------------

  /**
   *  \brief Interface class for all algorithms.
   *
   *  \ingroup Algorithms
   */
  template< class AlgorithmTraits >
  class EOCAlgorithm
    : public EOCAlgorithmInterface
  {
    typedef AlgorithmTraits Traits;

  public:
    // type of Grid
    typedef typename Traits::GridType                 GridType;
    typedef typename Traits::GridTypes                GridTypes;

    // type of statistics monitor
    typedef typename Traits::SolverMonitorCallerType  SolverMonitorCallerType;
    typedef typename Traits::DataWriterCallerType     DataWriterCallerType;
    typedef typename Traits::EocWriterCallerType      EocWriterCallerType;

    typedef typename Traits::CreateSubAlgorithmsType  CreateSubAlgorithmsType;

    typedef typename Traits::SubAlgorithmTupleType    SubAlgorithmTupleType;

    static const int size = std::tuple_size< SubAlgorithmTupleType >::value;

    typedef typename std::make_index_sequence< size > IndexSequenceType;

    typedef EocParameters                             EocParametersType;

    typedef uint64_t                                  UInt64Type;

    //! constructor
    template< class GlobalContainerImp >
    EOCAlgorithm ( const std::string algorithmName, const std::shared_ptr<GlobalContainerImp>& cont )
      : grids_( CreateSubAlgorithmsType::grids( cont ) ),
        tuple_( CreateSubAlgorithmsType::apply( cont ) ),
        solverMonitorCaller_( tuple_ ),
        dataWriterCaller_( tuple_ ),
        eocWriterCaller_( tuple_ ),
        algorithmName_( algorithmName ),
        eocParam_( ParameterKey::generate( "", "fem.eoc." ) )
    {}

    //! destructor
    virtual~ EOCAlgorithm()
    {}

    //! return default data tuple for data output
    virtual DataWriterCallerType& dataWriter ()
    {
      return dataWriterCaller_;
    }

    virtual SolverMonitorCallerType& monitor()
    {
      return solverMonitorCaller_;
    }

    //! return default data tuple for data output
    virtual EocWriterCallerType& eocWriter ()
    {
      return eocWriterCaller_;
    }

    //! return reference to hierarchical grid
    GridType& grid () const
    {
      return std::get<0>( grids_ );
    }

    //! return reference to hierarchical grids
    const GridTypes& grids () const
    {
      return grids_;
    }


    // return size of grid
    virtual UInt64Type gridSize () const
    {
      return gridSize( grids_, IndexSequenceType() );
    }

    // return size of grid
    virtual typename tuple_copy< size, UInt64Type >::type gridSizes () const
    {
      return gridSizes( grids_, IndexSequenceType() );
    }

    virtual void eocInitialize ()
    {
      //CALLER
      dataWriter().eocInitializeStart( this );
      eocWriter().eocInitializeStart( this );
    }

    virtual void eocPreSolve ( const int loop )
    {
      eocWriter().eocPreSolveEnd( this, loop );
    }

    // solve the problem for eoc loop 'loop'
    virtual void eocSolve ( const int loop )
    {}

    virtual void eocPostSolve ( const int loop )
    {
      //CALLER
      eocWriter().eocPostSolveEnd( this, loop );
      dataWriter().eocPostSolveEnd( this, loop );

      // Refine the grid for the next EOC Step.
      refineGrid( IndexSequenceType() );
    }

    virtual void eocFinalize()
    {}

    int eocSteps() const
    {
      return eocParam_.steps();
    }

    const EocParametersType& eocParams() const { return eocParam_; }

    virtual const std::string name() const { return algorithmName_; }

    SubAlgorithmTupleType &subAlgorithmTuple () { return tuple_; }
    const SubAlgorithmTupleType &subAlgorithmTuple () const { return tuple_; }

  protected:
    template< int i >
    struct RefineGrid
    {
      template< class T, class GridSizes >
      static void apply ( T &t, const GridSizes& oldGridSizes, GridSizes& newGridSizes  )
      {
        typedef std::remove_reference_t< std::tuple_element_t<i,GridTypes> > Element;
        const int step = std::get<i>(gridSizes( t, IndexSequenceType() ))==(int)std::get<i>(oldGridSizes) ?
                         Dune::DGFGridInfo<Element>::refineStepsForHalf() : 0;
        Fem::GlobalRefine::apply( std::get<i>(t), step );

        //manipulate newGridSizes to prevent double refinement
        //std::get<i>(newGridSizes) = std::get<0>(gridSizes( t, std::index_sequence<i>() ) );

      }
    };
    template< template< int > class Caller >
    using ForLoopType = ForLoop< Caller, 0, size-1 >;

    template< unsigned long int... i >
    static decltype(auto) gridSizes ( const GridTypes& grids, std::index_sequence<i...> seq )
    {
      auto grSizes = std::make_tuple( std::get<i>(grids).size( 0 )... );
      return std::make_tuple( std::get<i>(grids).comm().sum( std::get<i>( grSizes ) )... );
    }

    template< unsigned long int... i >
    static decltype(auto) gridSize ( const GridTypes& grids, std::index_sequence<i...> seq )
    {
      const auto& grdSizes = gridSizes( grids, seq );
      return Std::max( std::get<i>(grdSizes)...);
    }

    template< unsigned long int... i >
    void refineGrid( std::index_sequence<i...> seq )
    {
      if( eocSteps() > 1 )
      {
        auto newGridSizes = gridSizes();
        ForLoopType< RefineGrid >::apply( grids(), gridSizes(), newGridSizes );
      }
    }

  protected:
    GridTypes                      grids_;
    SubAlgorithmTupleType          tuple_;
    SolverMonitorCallerType        solverMonitorCaller_;
    DataWriterCallerType           dataWriterCaller_;
    EocWriterCallerType            eocWriterCaller_;
    const std::string              algorithmName_;
    EocParameters                  eocParam_;
  };




  /////////////////////////////////////////////////////////////////////////////
  //
  //  compute algorithm
  //
  /////////////////////////////////////////////////////////////////////////////
  template <class Algorithm>
  void compute(Algorithm& algorithm)
  {
    algorithm.eocInitialize();

    const unsigned int femTimerId = FemTimer::addTo("timestep");
    for(int eocloop=0; eocloop < algorithm.eocSteps(); ++eocloop )
    {
      // start fem timer (enable with -DFEMTIMER)
      FemTimer::start(femTimerId);

      algorithm.eocPreSolve( eocloop );

      // call algorithm and return solver statistics and some info
      algorithm.eocSolve( eocloop );

      // also get times for FemTimer if enabled
      FemTimer::stop(femTimerId);
      FemTimer::printFile("./timer.out");
      FemTimer::reset(femTimerId);

      algorithm.eocPostSolve( eocloop );
    } /***** END of EOC Loop *****/

    // FemTimer cleanup
    FemTimer::removeAll();

    algorithm.eocFinalize();
  }




} // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_BASE_HH
