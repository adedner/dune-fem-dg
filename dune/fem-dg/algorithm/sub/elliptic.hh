#ifndef DUNE_FEMDG_ALGORITHM_ELLIPTIC_ALGORITHM_HH
#define DUNE_FEMDG_ALGORITHM_ELLIPTIC_ALGORITHM_HH
#include <config.h>

// include std libs
#include <iostream>
#include <string>
#include <memory>

// dune-fem includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem-dg/misc/dgnorm.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/space/discontinuousgalerkin/legendre.hh>
#include <dune/fem/space/discontinuousgalerkin/space.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>

// dune-fem-dg includes
#include <dune/fem/operator/projection/dgl2projection.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/misc/gridname.hh>

// include local header files
#include "steadystate.hh"




namespace Dune
{
namespace Fem
{
  template <class SteadyStateContainerImp, class MatrixContainerImp >
  struct SubEllipticContainer
  {

    typedef typename SteadyStateContainerImp::DiscreteFunctionType        DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridType                  GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType              GridPartType;

    typedef MatrixContainerImp                                            MatrixType;

    using DiscreteFunction = DiscreteFunctionType;
    using DiscreteFunctionSpace = DiscreteFunctionSpaceType;
    using Matrix = MatrixType;

    typedef Fem::SubSteadyStateContainer< DiscreteFunctionType >          ContainerType;

  public:

    SubEllipticContainer( GridType& grid, const std::string name = "" )
    : stringId_( FunctionIDGenerator::instance().nextId() ),
      container_( grid, name ),
      matrix_( new MatrixType( name + "matrix",space(), space() ) )
    {}

    //container list
    const ContainerType& adapter() const
    {
      return container_;
    }
    ContainerType& adapter()
    {
      return container_;
    }

    //grid
    const GridType& grid() const
    {
      return adapter().grid();
    }
    GridType& grid()
    {
      return adapter().grid();
    }

    //grid part
    const GridPartType& gridPart() const
    {
      return adapter().gridPart();
    }
    GridPartType& gridPart()
    {
      return adapter().gridPart();
    }

    //spaces
    const DiscreteFunctionSpaceType& space() const
    {
      return adapter().space();
    }
    DiscreteFunctionSpaceType& space()
    {
      return adapter().space();
    }

    //solution
    std::shared_ptr< DiscreteFunction > solution() const
    {
      return adapter().solution();
    }
    void setSolution( std::shared_ptr< DiscreteFunction > solution )
    {
      adapter().solution() = solution;
    }

    //exact solution
    std::shared_ptr< DiscreteFunction > exactSolution() const
    {
      return adapter().exactSolution();
    }
    void setExactSolution( std::shared_ptr< DiscreteFunction > exactSolution )
    {
      adapter().exactSolution() = exactSolution;
    }

    //rhs
    std::shared_ptr< DiscreteFunction > rhs() const
    {
      return adapter().rhs();
    }
    void setRhs( std::shared_ptr< DiscreteFunction > rhs )
    {
      adapter().rhs() = rhs;
    }

    //matrix for assembly
    std::shared_ptr< MatrixType > matrix() const
    {
      assert( matrix_ );
      return matrix_;
    }
    void setMatrix( std::shared_ptr< MatrixType > matrix )
    {
      matrix_ = matrix;
    }


  private:
    const std::string                   stringId_;
    ContainerType                       container_;
    std::shared_ptr< Matrix >           matrix_;
  };


  template< class DiscreteFunctionImp, class AssemblerImp, int polOrder>
  class PoissonSigmaEstimator
  {
  public:

    typedef DiscreteFunctionImp                                      DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::GridPartType              GridPartType;
    typedef typename GridPartType::GridType                          GridType;
    typedef AssemblerImp                                             AssemblerType;
    static const int polynomialOrder = polOrder;

    typedef typename DiscreteFunctionSpaceType ::
      template ToNewDimRange< GridType::dimension * DiscreteFunctionSpaceType::FunctionSpaceType::dimRange >::NewFunctionSpaceType
                                                                     SigmaFunctionSpaceType;

    PoissonSigmaEstimator( GridPartType& gridPart,
                           const DiscreteFunctionType& solution,
                           const AssemblerType& assembler,
                           const std::string name = "" )
    : gridPart_( gridPart ),
      solution_( solution ),
      assembler_( assembler ),
      sigmaSpace_( gridPart_ ),
      sigmaDiscreteFunction_( "sigma-"+name, sigmaSpace_ ),
      sigmaLocalEstimate_( solution_, assembler_ ),
      sigmaLocalFunction_( solution_, sigmaDiscreteFunction_, sigmaLocalEstimate_ ),
      sigma_( "sigma function", sigmaLocalFunction_, gridPart_, solution_.space().order() ),
      sigmaEstimateFunction_( "function 4 estimate-"+name, sigmaLocalEstimate_, gridPart_, solution_.space().order() )
    {}


    //- hybrid spaces use PAdaptiveDGSpace
    template <class Grid, int topoId>
    struct SigmaSpaceChooser
    {
      typedef Fem::PAdaptiveDGSpace< SigmaFunctionSpaceType, GridPartType, polOrder > Type ;
    };

    //- cubes use LegendreDGSpace
    template <class Grid>
    struct SigmaSpaceChooser< Grid, 1 >
    {
      typedef Dune::Fem::LegendreDiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polOrder > Type ;
    };

    //- cubes use OrthonormalDGSpace
    template <class Grid>
    struct SigmaSpaceChooser< Grid, 0 >
    {
      typedef Dune::Fem::DiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polOrder > Type ;
    };

    // work arround internal compiler error
    enum { simplexTopoId =  GenericGeometry::SimplexTopology< GridType::dimension >::type::id,
           cubeTopoId    =  GenericGeometry::CubeTopology< GridType::dimension >::type::id,
           myTopo        =  Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
         };

    enum { topoId = (simplexTopoId == myTopo) ? 0 :  // 0 = simplex, 1 = cube, -1 = hybrid
                      (myTopo == cubeTopoId) ? 1 : -1 };

    typedef typename SigmaSpaceChooser< GridType, topoId >::Type  SigmaDiscreteFunctionSpaceType ;

    //typedef Dune::Fem::AdaptiveDiscreteFunction< SigmaDiscreteFunctionSpaceType >
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< SigmaDiscreteFunctionSpaceType >
                                                                  SigmaDiscreteFunctionType;



    // compute the function sigma = grad u + sum_e r_e
    template <class DF, class Operator>
    struct SigmaLocal : public Fem::LocalFunctionAdapterHasInitialize
    {
      typedef typename DF::DiscreteFunctionSpaceType                     UDFS;
      typedef typename UDFS::GridPartType                                GridPartType;
      typedef typename GridPartType::GridType::template Codim<0>::Entity EntityType;
      typedef typename GridPartType::IntersectionIteratorType            IntersectionIteratorType;
      typedef typename GridPartType::IntersectionType                    IntersectionType;

      typedef typename Operator::DiffusionFluxType::LiftingFunctionType  LiftingFunctionType;
      typedef typename LiftingFunctionType::RangeType                    RangeType;
      typedef typename LiftingFunctionType::DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType      FunctionSpaceType;

      typedef typename DF::RangeType                                     URangeType;
      typedef typename DF::JacobianRangeType                             UJacobianRangeType;

      SigmaLocal( const DF &df, const Operator &oper )
      : df_(df), oper_(oper), localdf_(df_), reSpace_( oper.gradientSpace() ), localre_( reSpace_ )
      {}
      SigmaLocal(const SigmaLocal &other)
      : df_(other.df_), oper_(other.oper_), localdf_(df_), reSpace_(other.reSpace_), localre_( reSpace_ )
      {}
      ~SigmaLocal()
      {
      }
      template <class PointType>
      void evaluate(const PointType& x, RangeType& val) const
      {
        typename DF::JacobianRangeType jac;
        localdf_.jacobian(x,jac);
        localre_.evaluate(x,val);
        Dune::Fem::FieldMatrixConverter< RangeType, typename DF::JacobianRangeType> val1( val );
        val1 += jac;
      }
      void init(const EntityType& entity)
      {
        localdf_.init(entity);

        localre_.init(entity);
        localre_.clear();
        IntersectionIteratorType end = df_.space().gridPart().iend( entity );
        for( IntersectionIteratorType it = df_.space().gridPart().ibegin( entity ); it != end; ++it )
        {
          const IntersectionType &intersection = *it;
          if ( intersection.neighbor() && df_.space().continuous(intersection) )
          {
            if( ! intersection.conforming() )
              getLifting< false > ( intersection, entity ) ;
            else
              getLifting< true > ( intersection, entity );
          }
        }
      }
      private:
      template <bool conforming>
      void getLifting( const IntersectionType &intersection, const EntityType &entity)
      {
        // CACHING
        typedef typename Operator::FaceQuadratureType                               FaceQuadratureType ;
        typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
        typedef typename IntersectionQuadratureType::FaceQuadratureType             QuadratureImp;

        const EntityType &outside = intersection.outside();

        typename DF::LocalFunctionType uOutside = df_.localFunction(outside);

        const int enOrder = df_.space().order( entity );
        const int nbOrder = df_.space().order( outside );

        const int quadOrder = 2 * std::max( enOrder, nbOrder ) + 1;

        IntersectionQuadratureType interQuad( df_.space().gridPart(), intersection, quadOrder );
        const QuadratureImp &quadInside  = interQuad.inside();
        const QuadratureImp &quadOutside = interQuad.outside();
        const int numQuadraturePoints = quadInside.nop();

        // obtain all required function values on intersection
        std::vector< URangeType > uValuesEn( numQuadraturePoints );
        std::vector< URangeType > uValuesNb( numQuadraturePoints );
        localdf_.evaluateQuadrature( quadInside, uValuesEn );
        uOutside.evaluateQuadrature( quadOutside, uValuesNb );
        oper_.lifting(df_.space().gridPart(),
                      intersection, entity, outside, 0, quadInside, quadOutside,
                      uValuesEn, uValuesNb,
                      localre_
                     );
      }
      const DF&                        df_;
      const Operator&                  oper_;
      typename DF::LocalFunctionType   localdf_;
      const DiscreteFunctionSpaceType& reSpace_;
      LiftingFunctionType              localre_;
    };

    template <class SigmaLocalType>
    struct SigmaLocalFunction : public Fem::LocalFunctionAdapterHasInitialize
    {
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionType::RangeType                 RangeType;
      typedef typename DiscreteFunctionType::JacobianRangeType         JacobianRangeType;
      typedef typename DiscreteFunctionType::EntityType                EntityType;
      typedef typename DiscreteFunctionSpaceType::FunctionSpaceType    FunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridPartType         GridPartType;

      SigmaLocalFunction( const DiscreteFunctionType &u,
                          const SigmaDiscreteFunctionType &q,
                          const SigmaLocalType &sigmaLocal )
      : u_(u), uLocal_(u), q_(q), qLocal_(q), sigmaLocal_(sigmaLocal)
      {}
      SigmaLocalFunction(const SigmaLocalFunction &other)
      : u_(other.u_), uLocal_(u_), q_(other.q_), qLocal_(q_), sigmaLocal_(other.sigmaLocal_)
      {}
      ~SigmaLocalFunction()
      {
      }
      template <class PointType>
      void evaluate(const PointType& x, RangeType& val) const
      {
        uLocal_.evaluate(x,val);
      }
      template <class PointType>
      void jacobian(const PointType& x, JacobianRangeType& val) const
      {
        typename SigmaLocalType::RangeType qval;
        qLocal_.evaluate(x,qval);
        // sigmaLocal_.evaluate(x,qval);
        Dune::Fem::FieldMatrixConverter< typename SigmaLocalType::RangeType, JacobianRangeType> val1( qval );
        val = val1;
        // uLocal_.jacobian(x,val);
      }
      void init(const EntityType& entity)
      {
        uLocal_.init(entity);
        qLocal_.init(entity);
        sigmaLocal_.init(entity);
      }
      private:
      const DiscreteFunctionType&                           u_;
      typename DiscreteFunctionType::LocalFunctionType      uLocal_;
      const SigmaDiscreteFunctionType&                      q_;
      typename SigmaDiscreteFunctionType::LocalFunctionType qLocal_;
      SigmaLocalType                                        sigmaLocal_;
    };

    typedef Dune::Fem::LocalFunctionAdapter< SigmaLocal<DiscreteFunctionType, AssemblerType> >
                                                                    SigmaEstimateFunctionType;
    typedef SigmaLocal<DiscreteFunctionType, AssemblerType>         SigmaLocalType;
    typedef SigmaLocalFunction<SigmaLocalType >                     SigmaLocalFunctionType;
    typedef Dune::Fem::LocalFunctionAdapter<SigmaLocalFunctionType> SigmaLocalFunctionAdapterType;

    void update ()
    {
      Dune::Fem::DGL2ProjectionImpl::project( sigmaEstimateFunction_, sigmaDiscreteFunction_ );
    }

    //const SigmaLocalFunctionAdapterType& sigma () const
    //{
    //  return sigma_;
    //}

    const SigmaDiscreteFunctionType& sigma () const
    {
      return sigmaDiscreteFunction_;
    }


  public:

    GridPartType&                   gridPart_;
    const DiscreteFunctionType&     solution_;
    const AssemblerType&            assembler_;
    SigmaDiscreteFunctionSpaceType  sigmaSpace_;
    SigmaDiscreteFunctionType       sigmaDiscreteFunction_;

    SigmaLocal<DiscreteFunctionType, AssemblerType>
                                    sigmaLocalEstimate_;
    SigmaLocalFunctionType          sigmaLocalFunction_;
    SigmaLocalFunctionAdapterType   sigma_;
    SigmaEstimateFunctionType       sigmaEstimateFunction_;

  };



  template< int polOrder >
  class NoEstimator
  {
    public:

    NoEstimator(){}


    int minimalOrder()
    {
      return polOrder;
    }

    bool isPadaptive()
    {
      return false;
    }

    template< class ProblemImp>
    double estimate( ProblemImp problem )
    { return 0.0; }

    template< class EntityImp >
    int newOrder( double tolerance, EntityImp& entity )
    {
      return minimalOrder();
    }

    bool mark( double tolerance ){ return false; }

    void closure(){}
  };


  template< int polOrder, class DiscreteFunctionSpaceImp, class EstimatorImp = NoEstimator<polOrder> >
  class PAdaptivity
  {
    struct PolOrderStructure
    {
      // set polynomial order to 2 by default
      PolOrderStructure() : val_( -1 ) {}
      explicit PolOrderStructure(int init) : val_( init ) {}
      const int value() const
      {
        assert( val_ > 0 );
        return val_;
      }
      int &value() { return val_; }
      int val_;
    };

  public:
    typedef DiscreteFunctionSpaceImp                        DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridType    GridType;
    typedef EstimatorImp                                    EstimatorType;

    typedef PersistentContainer<GridType,PolOrderStructure> PolOrderContainer;

    PAdaptivity( GridType& grid, const DiscreteFunctionSpaceType& space, EstimatorType& estimator, const std::string name = ""  )
      : polOrderContainer_( grid, 0 ),
        space_( space ),
        estimator_( estimator ),
        param_( AdaptationParameters() )
    {
#ifdef PADAPTSPACE
      // we start with max order
      typedef typename PolOrderContainer::Iterator Iterator ;
      const Iterator end = polOrderContainer_.end();
      const int minimalOrder = estimator_.minimalOrder() ;
      for( Iterator it = polOrderContainer_.begin(); it != end; ++it )
      {
        (*it).value() = minimalOrder ;
      }
#endif
    }

    bool adaptive() const
    {
      return estimator_.isPadaptive();
    }

    void prepare()
    {
#ifdef PADAPTSPACE
      const int minimalOrder = estimator_.minimalOrder();
      // only implemented for PAdaptiveSpace
      std::vector<int> polOrderVec( space_.gridPart().indexSet().size(0) );
      std::fill( polOrderVec.begin(), polOrderVec.end(), polOrder );

      polOrderContainer_.resize();
      if ( estimator_.isPadaptive() )
      {
        typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
        typedef typename IteratorType::Entity EntityType ;
        const IteratorType end = space_.end();
        for( IteratorType it = space_.begin(); it != end; ++it )
        {
          const EntityType& entity = *it;
          int order = polOrderContainer_[ entity ].value();
          while (order == -1) // is a new element
          {
            if ( entity.level() == 0)
              order = minimalOrder;
            else
            {
              // don't call father twice
              order = polOrderContainer_[ entity.father() ].value();
              assert(order > 0);
            }
          }
          polOrderVec[space_.gridPart().indexSet().index(entity)] = order;
        }
      }
      space_.adapt( polOrderVec );
#endif
    }


    template< class ProblemImp >
    bool estimateMark( const ProblemImp& problem )
    {
      double tolerance = param_.refinementTolerance();
#ifdef PADAPTSPACE
      // resize container
      polOrderContainer_.resize();

      const double error = estimator_.estimate( problem );
      std::cout << "ESTIMATE: " << error << std::endl;
      typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
      const IteratorType end = space_.end();
      for( IteratorType it = space_.begin(); it != end; ++it )
      {
        const typename IteratorType::Entity &entity = *it;
        polOrderContainer_[entity].value() =
          estimator_.newOrder( 0.98*tolerance, entity );
      }
      return (error < std::abs(tolerance) ? false : estimator_.mark( 0.98 * tolerance));
#else
      return false;
#endif
    }
    void closure()
    {
#ifdef PADAPTSPACE
      estimator_.closure();
#endif
    }

    private:

    PolOrderContainer                polOrderContainer_;
    const DiscreteFunctionSpaceType& space_;
    EstimatorType&                   estimator_;
    AdaptationParameters             param_;

  };


  /**
   * \brief Adaptation indicator doing no indication and marking of the entities.
   */
  template< class EstimatorImp, class SigmaEstimatorImp, class ProblemImp >
  class PAdaptIndicator
  {
    typedef ProblemImp                                                   ProblemType;
    typedef SigmaEstimatorImp                                            SigmaEstimatorType;
    typedef EstimatorImp                                                 EstimatorType;
    typedef typename SigmaEstimatorType::DiscreteFunctionType            DiscreteFunctionType;
    typedef typename SigmaEstimatorType::AssemblerType                   AssemblerType;
    typedef typename SigmaEstimatorType::SigmaFunctionSpaceType          SigmaFunctionSpaceType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::GridPartType                  GridPartType;
    typedef typename GridPartType::GridType                              GridType;
    static const int polOrder = SigmaEstimatorType::polynomialOrder;

    typedef PAdaptivity< polOrder, DiscreteFunctionSpaceType, EstimatorImp > PAdaptivityType;

    enum { dimension = GridType::dimension  };

  public:
   typedef uint64_t                          UInt64Type;

    PAdaptIndicator( GridPartType& gridPart,
                     DiscreteFunctionType& solution,
                     const ProblemType& problem,
                     AssemblerType& assembler,
                     const std::string name = "" )
      : sigmaEstimator_( gridPart, solution, assembler, name ),
        estimator_( solution, sigmaEstimator_.sigma(), assembler, gridPart.grid(), name /*, AdaptationParameters( param ) */),
        pAdapt_( gridPart.grid(), solution.space(), estimator_ ),
        problem_( problem )
    {}

    bool adaptive() const
    {
      return pAdapt_.adaptive();
    }

    size_t numberOfElements() const
    {
      return 0;
    }

    UInt64Type globalNumberOfElements() const { return 0; }

    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp )
    {
    }

    void preAdapt()
    {
      pAdapt_.prepare();
    }

    void estimateMark( const bool initialAdapt = false )
    {
      // calculate new sigma
      sigmaEstimator_.update();

      const bool marked = pAdapt_.estimateMark( problem_ );
      if( marked )
        pAdapt_.closure();
    }

    void postAdapt()
    {
    }

    void finalize()
    {
    }

    int minNumberOfElements() const { return 0; }

    int maxNumberOfElements() const { return 0; }

    const int finestLevel() const { return 0; }

    // return some info
    const typename SigmaEstimatorType::SigmaDiscreteFunctionType& sigma()
    {
      return sigmaEstimator_.sigma();
    }

  private:

    SigmaEstimatorType sigmaEstimator_;
    EstimatorType      estimator_;
    PAdaptivityType    pAdapt_;
    const ProblemType& problem_;
  };




  /**
   *  \brief Algorithm for solving an elliptic PDE.
   *
   *  \ingroup SubAlgorithms
   */
  template <class GridImp,
            class ProblemTraits,
            int polOrder>
  class SubEllipticAlgorithm
    : public SubSteadyStateAlgorithm< GridImp, ProblemTraits, polOrder >
  {
  public:
    typedef SubSteadyStateAlgorithm< GridImp, ProblemTraits, polOrder > BaseType;

    // type of Grid
    typedef typename BaseType::GridType                    GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                GridPartType;

    // initial data type
    typedef typename BaseType::ProblemType                 ProblemType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                   ModelType;

    // type of linear operator (i.e. matrix implementation)
    typedef typename BaseType::OperatorType::type          OperatorType;

    // The DG space operator
    typedef typename BaseType::OperatorType::AssemblerType AssemblerType;


    typedef typename AssemblerType::ContainerType          ContainerType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType        DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;

    // type of inverse operator (i.e. linear solver implementation)
    typedef typename BaseType::SolverType::type            SolverType;

    enum { dimension = GridType::dimension  };

    typedef typename  BaseType::AnalyticalTraits           AnalyticalTraits;

    typedef typename BaseType::IOTupleType                 IOTupleType;

    typedef typename BaseType::TimeProviderType            TimeProviderType;

    typedef typename BaseType::AdaptIndicatorType          AdaptIndicatorType;

    typedef typename BaseType::AdaptationDiscreteFunctionType AdaptationDiscreteFunctionType;

    using BaseType::problem;
    using BaseType::model;
    using BaseType::name;
    using BaseType::grid;
    using BaseType::gridWidth;
    using BaseType::gridSize;
    using BaseType::solution;
    using BaseType::rhs;
    using BaseType::exactSolution;
    using BaseType::solver;

  public:
    SubEllipticAlgorithm( GridType& grid, ContainerType& container )
    : BaseType( grid, container.adapter() ),
      container_( container ),
      gridPart_( container_.gridPart() ),
      space_( container_.space() ),
      assembler_( container_, model() ),
      matrix_( container_.matrix() ),
      adaptIndicator_( std::make_unique<AdaptIndicatorType>( gridPart_, *container_.solution().get(), problem(), assembler_, name() ) ),
      step_( 0 ),
      time_( 0 )
    {
      std::string gridName = Fem::gridName( grid );
      if( gridName == "ALUGrid" || gridName == "ALUConformGrid" || gridName == "ALUSimplexGrid" )
      {
        if( space_.begin() != space_.end() )
        {
          if( space_.begin()->type().isSimplex() && space_.order() > 2 && space_.continuous() && GridType::dimension > 2 )
          {
            std::cerr << std::endl<< "ERROR: Lagrange spaces for p>2 do not work on simplex grids due to the twist problem !!!" << std::endl << std::endl;
          }
        }
      }
    }

    void virtual setTime ( const double time ) override
    {
      time_ = time;
    }

    const AssemblerType& assembler () const
    {
      return assembler_;
    }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator()
    {
      return adaptIndicator_.get();
    }
    virtual AdaptationDiscreteFunctionType* adaptationSolution ()
    {
      return &solution();
    }

  protected:
    virtual std::shared_ptr< SolverType > doCreateSolver() override
    {
      Dune::Timer timer;
      timer.start();

      if( space_.continuous() ) // Lagrange case
      {
        typedef Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType ;
        StencilType stencil( space_, space_ );
        matrix_->reserve( stencil );
      }
      else // DG case
      {
        typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType ;
        StencilType stencil( space_, space_ );
        matrix_->reserve( stencil );
      }

      assembler_.assemble( time_ );
      std::cout << "Solver (Poisson) assemble time: " << timer.elapsed() << std::endl;

      assembler_.testSymmetrie();

      double absLimit   = Dune::Fem:: Parameter::getValue<double>("istl.absLimit",1.e-10);
      double reduction  = Dune::Fem:: Parameter::getValue<double>("istl.reduction",1.e-6);

      return std::make_shared< SolverType >(*matrix_, reduction, absLimit );
    }

    //! default time loop implementation, overload for changes
    virtual void doPreSolve ( const int loop ) override
    {
      BaseType::doPreSolve( loop );
    }

    //! finalize computation by calculating errors and EOCs
    virtual void doFinalize ( const int loop ) override
    {
      AnalyticalTraits::addEOCErrors( solution(), model(), problem().exactSolution(), adaptIndicator()->sigma() );
    }


  protected:

    ContainerType&                               container_;
    GridPartType&                                gridPart_;       // reference to grid part, i.e. the leaf grid
    const DiscreteFunctionSpaceType&             space_;

    AssemblerType                                assembler_;

    std::shared_ptr< OperatorType >              matrix_;
    std::unique_ptr< AdaptIndicatorType >        adaptIndicator_;
    int                                          step_;
    double                                       time_;

  };

}
}
#endif // FEMHOWTO_STEPPER_HH
