#ifndef DUNE_FEMDG_ALGORITHM_INTERFACE_HH
#define DUNE_FEMDG_ALGORITHM_INTERFACE_HH

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem-dg/misc/parameterkey.hh>
#include <dune/fem-dg/misc/typedefcheck.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/misc/covarianttuple.hh>

#include <dune/fem-dg/algorithm/base.hh>
#include <dune/fem-dg/algorithm/sub/evolution.hh>

#include <dune/fem-dg/algorithm/handler/solvermonitor.hh>
#include <dune/fem-dg/algorithm/handler/diagnostics.hh>
#include <dune/fem-dg/algorithm/handler/additionaloutput.hh>


namespace Dune
{

namespace Fem
{
  /**
   *  \addtogroup Algorithms
   *
   *  Algorithms are mainly partial differential equations solvers.
   *  They provide all important information for an efficient solve of the PDE.
   *
   *
   */

  /**
   *  \addtogroup SubAlgorithms
   *
   *  \ingroup Algorithms
   *
   *  Structure of SubAlgorithms
   */

  /**
   *  \addtogroup Algorithm Handlers
   *
   *  \ingroup Algorithms
   *
   *  Handlers...
   */

  template< class Grid, class ProblemTraits, int polOrder >
  struct SubAlgorithmInterfaceTraits
  {

  private:
    CHECK_TYPEDEF_EXISTS( AdaptIndicatorType )
    CHECK_TYPEDEF_EXISTS( AdditionalOutputHandlerType )
    CHECK_TYPEDEF_EXISTS( SolverMonitorHandlerType )
    CHECK_TYPEDEF_EXISTS( DiagnosticsHandlerType )

    typedef ProblemTraits                                         Traits;

    typedef typename Traits::AnalyticalTraits                     AnalyticalTraits;
    typedef typename Traits::template DiscreteTraits< polOrder >  DiscreteTraits;

  public:

    //! type of the grid all discrete functions are belonging to
    typedef typename Traits::GridType                             GridType;
    //! type of the grid part
    typedef typename Traits::GridPartType                         GridPartType;
    //! type of the host grid part
    typedef typename Traits::HostGridPartType                     HostGridPartType;

    //! type of the time provider
    typedef GridTimeProvider< GridType >                          TimeProviderType;

    //! type of the model
    typedef typename AnalyticalTraits::ModelType                  ModelType;
    //! type of the problem
    typedef typename AnalyticalTraits::ProblemType                ProblemType;

    //! type of the discrete function type
    typedef typename DiscreteTraits::DiscreteFunctionType         DiscreteFunctionType;

    // type of the unsigned 64 bit integer
    typedef uint64_t                                              UInt64Type ;

    typedef DiscreteFunctionType                                  CheckPointDiscreteFunctionType;
    typedef DiscreteFunctionType                                  LimitDiscreteFunctionType;
    typedef DiscreteFunctionType                                  AdaptationDiscreteFunctionType;

    typedef CovariantTuple< typename DiscreteTraits::IOTupleType >          IOTupleType;
    typedef typename AdaptIndicatorTypes< DiscreteTraits >::type            AdaptIndicatorType;
    typedef typename AdditionalOutputHandlerTypes< DiscreteTraits >::type   AdditionalOutputHandlerType;
    typedef typename SolverMonitorHandlerTypes< DiscreteTraits >::type      SolverMonitorHandlerType;
    typedef typename DiagnosticsHandlerTypes< DiscreteTraits >::type        DiagnosticsHandlerType;
  };


  /**
   *  \brief This is the interface all sub algorithms (stationary and instationary) depends on.
   *
   *  \ingroup Algorithms
   *
   *  The class SubAlgorithmInterface is the base class for all sub algorithms. We suggest to derive
   *  your own Algorithm from this interface.
   *
   *  \note For this interface we de not really care about time dependency or the PDE.
   *  Both interfaces (stationary and instationary) are provided through this interface.
   *
   *  \note The SubEvolutionAlgorithm and the SubSteadyStateAlgorithm classes encapsulate the
   *  stationary and instationary behaviour.
   *
   *  \tparam Grid grid type
   *  \tparam ProblemTraits Traits class defining all important types of the algorithm
   *  \tparam polOrder general polynomial order of the scheme
   */
  template< class Grid, class ProblemTraits, int polOrder >
  class SubAlgorithmInterface
  {

  public:
    typedef SubAlgorithmInterfaceTraits< Grid, ProblemTraits, polOrder > Traits;

    //! type of the grid
    typedef typename Traits::GridType                             GridType;
    //! type of the grid part
    typedef typename Traits::GridPartType                         GridPartType;
    //! type of the host grid part
    typedef typename Traits::HostGridPartType                     HostGridPartType;

    //! type of the time provider
    typedef typename Traits::TimeProviderType                     TimeProviderType;

    //! type of the model
    typedef typename Traits::ModelType                            ModelType;
    //! type of the problem
    typedef typename Traits::ProblemType                          ProblemType;

    //! type of the discrete solution
    typedef typename Traits::DiscreteFunctionType                 DiscreteFunctionType;

    //! type of the unsigned 64 bit integer
    typedef typename Traits::UInt64Type                           UInt64Type;

    typedef typename Traits::CheckPointDiscreteFunctionType       CheckPointDiscreteFunctionType;
    typedef typename Traits::LimitDiscreteFunctionType            LimitDiscreteFunctionType;
    typedef typename Traits::AdaptationDiscreteFunctionType       AdaptationDiscreteFunctionType;

    typedef typename Traits::IOTupleType                          IOTupleType;
    typedef typename Traits::AdaptIndicatorType                   AdaptIndicatorType;
    typedef typename Traits::DiagnosticsHandlerType               DiagnosticsHandlerType;
    typedef typename Traits::SolverMonitorHandlerType             SolverMonitorHandlerType;
    typedef typename Traits::AdditionalOutputHandlerType          AdditionalOutputHandlerType;


    /**
     * \brief Constructor
     *
     * \param grid A grid
     * \param name A unique name for the sub algorithm
     */
    SubAlgorithmInterface ( GridType &grid, const std::string name = "" )
    : grid_( grid ),
      algorithmName_( name ),
      problem_( ProblemTraits::problem() ),
      model_( problem() )
    {}

    /**
     * \brief return the name of the algorithm.
     *
     * The idea of is to cleary identify an algorithm. This could be used for
     * reading sub algorithm specific parameters in the parameter file.
     * A sub algorithm name "foo" will use (nearly everywhere)
     * \code{.txt}
     *   foo.parameter: boo
     * \endcode
     *
     * instead of
     *
     * \code{.txt}
     *   parameter: boo
     * \endcode
     *
     * Since this is an interface, it is up to the user whether he wants to
     * use this feature or not.
     *
     * \return name of the algorithm.
     */
    const std::string name () const { return algorithmName_; }

    /**
     * \brief returns the grid
     *
     * \return the grid \f$ G \f$ all discrete solutions \f$ u_h \f$ depend on.
     */
    GridType& grid () const { return grid_; }

    /**
     *  \brief return minimal grid width \f$ h \f$ of grid
     *
     *  This method is mainly needed for documentation of numerical results.
     *
     *  \note override this method to for example calculate a grid width which
     *  is dependent on a grid part only.
     */
    virtual double gridWidth () const
    {
      GridWidthProvider< GridType > gwp( &grid_ );
      return gwp.gridWidth();
    }

     /**
     *  \brief return number of (codim 0 ) elements of a grid
     *
     *  This method is mainly needed for documentation of numerical results.
     *
     *  \note override this method to for example calculate a grid size which
     *  is dependent on a grid part only.
     */
    virtual UInt64Type gridSize () const
    {
      UInt64Type grSize = grid().size(0);
      return grid().comm().sum( grSize );
    }


    // return reference to discrete function holding solution
    virtual DiscreteFunctionType& solution () = 0;

    //DATAWRITING
    virtual IOTupleType& dataTuple () = 0;

    //SOLVERMONITOR
    virtual SolverMonitorHandlerType* monitor() { return nullptr; }

    //DIAGNOSTICS
    virtual DiagnosticsHandlerType* diagnostics() { return nullptr; }

    //ADDITIONALOUTPUT
    virtual AdditionalOutputHandlerType* additionalOutput() { return nullptr; }

    //LIMITING
    virtual void limit(){}
    virtual LimitDiscreteFunctionType* limitSolution () { return nullptr; }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() { return nullptr; }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () { return nullptr; }

    //CHECKPOINTING
    virtual CheckPointDiscreteFunctionType* checkPointSolution () { return nullptr; }

    /**
     * \brief Checks whether a discrete solution is valid (contains NANs or other unphysical solutions) or not.
     *
     * \note Overide the methods doCheckSolutionValid(int, TimeProviderType&) const or doCheckSolutionValid(int) const to implement a check.
     *
     * \param loop number of eoc loop
     * \param tp time provider if provided
     * \return true if valid, false otherwise
     */
    bool checkSolutionValid ( int loop, TimeProviderType* tp = nullptr ) const
    {
      if( tp == nullptr ) return doCheckSolutionValid( loop );
      else return doCheckSolutionValid( loop, *tp );
    }

    void initialize ( int loop, TimeProviderType* tp = nullptr )
    {
      if( tp == nullptr ) doInitialize( loop );
      else doInitialize( loop, *tp );
    }

    void preSolve ( int loop, TimeProviderType* tp = nullptr )
    {
      if( tp == nullptr ) doPreSolve( loop );
      else doPreSolve( loop, *tp );
    }

    void solve ( int loop, TimeProviderType* tp = nullptr )
    {
      if( tp == nullptr ) doSolve( loop );
      else doSolve( loop, *tp );
    }

    void postSolve ( int loop, TimeProviderType* tp = nullptr )
    {
      if( tp == nullptr ) doPostSolve( loop );
      else doPostSolve( loop, *tp );
    }

    void finalize ( int loop, TimeProviderType* tp = nullptr )
    {
      if( tp == nullptr ) doFinalize( loop );
      else doFinalize( loop, *tp );
    }

    virtual ProblemType& problem ()
    {
      assert( problem_ );
      return *problem_;
    }

    virtual const ProblemType& problem () const
    {
      assert( problem_ );
      return *problem_;
    }

    virtual const ModelType& model () const { return model_; }
    virtual ModelType& model () { return model_; }


  protected:
    virtual bool doCheckSolutionValid ( const int loop, TimeProviderType& tp ) const { return true; }
    /**
     * \brief instationary version
     */
    virtual bool doCheckSolutionValid ( const int loop ) const { return true; }

    virtual void doInitialize ( int loop, TimeProviderType& tp ){}
    virtual void doInitialize ( int loop ){}

    virtual void doPreSolve( int loop, TimeProviderType& tp ){}
    virtual void doPreSolve( int loop ){}

    virtual void doSolve ( int loop, TimeProviderType& tp ){}
    virtual void doSolve ( int loop ){}

    virtual void doPostSolve( int loop, TimeProviderType& tp ){}
    virtual void doPostSolve( int loop ){}

    virtual void doFinalize ( int loop, TimeProviderType& tp ){}
    virtual void doFinalize ( int loop ){}

    GridType&    grid_;
    std::string  algorithmName_;
    ProblemType* problem_;
    ModelType    model_;

  };



}  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_INTERFACE_HH