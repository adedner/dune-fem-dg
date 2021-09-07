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
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>
#include <dune/fem/function/localfunction/const.hh>

// dune-fem-dg includes
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem-dg/pass/context.hh>

#include <dune/fem/misc/gridname.hh>

// include local header files
#include "steadystate.hh"
#include "containers.hh"



namespace Dune
{
namespace Fem
{

  template< class MatrixContainerImp >
  struct SubEllipticContainerItem
  {
  public:
    using Matrix = MatrixContainerImp;

    template< class ContainerItem1, class ContainerItem2 >
    SubEllipticContainerItem ( const ContainerItem1& row, const ContainerItem2& col, const std::string name = "" )
    : matrix_( std::make_shared< Matrix >( name + "matrix", col->solution()->space(), row->solution()->space() ) )
    {}

    //matrix for assembly
    std::shared_ptr< Matrix > matrix()
    {
      return matrix_;
    }

  private:
    std::shared_ptr< Matrix >            matrix_;
  };


  template< template<class,class> class MatrixImp, class... DiscreteFunctions >
  struct SubEllipticContainer
    : public TwoArgContainer< ArgContainerArgWrapperUnique< SubEllipticContainerItem, MatrixImp >::template _t2Inv,
                              ArgContainerArgWrapperUnique< SubSteadyStateContainerItem >::template _t1,
                              ArgContainerArgWrapperUnique< SubSteadyStateContainerItem >::template _t1,
                              std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > >
  {
    typedef TwoArgContainer< ArgContainerArgWrapperUnique< SubEllipticContainerItem, MatrixImp >::template _t2Inv,
                             ArgContainerArgWrapperUnique< SubSteadyStateContainerItem >::template _t1,
                             ArgContainerArgWrapperUnique< SubSteadyStateContainerItem >::template _t1,
                             std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > > BaseType;

  public:
    using BaseType::operator();

    // constructor: do not touch/delegate everything
    template< class ... Args>
    SubEllipticContainer( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

  };


  template< class ErrorEstimatorImp >
  class PoissonSigmaEstimator
  {
  public:
    typedef ErrorEstimatorImp                                           ErrorEstimatorType;

    typedef typename ErrorEstimatorType::DiscreteFunctionType           DiscreteFunctionType;
    typedef typename ErrorEstimatorType::SigmaDiscreteFunctionType      SigmaDiscreteFunctionType;
    typedef typename ErrorEstimatorType::SigmaDiscreteFunctionSpaceType SigmaDiscreteFunctionSpaceType;
    typedef typename ErrorEstimatorType::DGOperatorType                 DGOperatorType;

    //typedef typename DGOperatorType::ContainerType                      ContainerType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;
    static const int polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder;
    typedef typename DiscreteFunctionType::GridPartType                 GridPartType;
    typedef typename GridPartType::GridType                             GridType;

    template< class ContainerImp >
    PoissonSigmaEstimator( const std::shared_ptr< ContainerImp >& cont,
                           const DGOperatorType& assembler,
                           const std::string name = "" )
    : gridPart_( (*cont)(_0)->solution()->space().gridPart() ),
      solution_( *(*cont)(_0)->solution() ),
      assembler_( assembler ),
      sigmaSpace_( gridPart_ ),
      sigmaDiscreteFunction_( "sigma-"+name, sigmaSpace_ ),
      sigmaLocalEstimate_( solution_, assembler_ ),
      sigmaLocalFunction_( solution_, sigmaDiscreteFunction_, sigmaLocalEstimate_ ),
      sigma_( "sigma function", sigmaLocalFunction_, gridPart_, solution_.space().order() ),
      sigmaEstimateFunction_( "function 4 estimate-"+name, sigmaLocalEstimate_, gridPart_, solution_.space().order() )
    {}

    // compute the function sigma = grad u + sum_e r_e
    template <class DF, class Operator>
    struct SigmaLocal /*: public Fem::LocalFunctionAdapterHasInitialize*/
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
      : df_(df), oper_(oper), localdf_(df_), uOutside_( df_ ),
        reSpace_( oper.gradientSpace() ), localre_( reSpace_ )
      {}
      SigmaLocal(const SigmaLocal &other)
      : df_(other.df_), oper_(other.oper_), localdf_(df_), uOutside_( df_ ),
        reSpace_(other.reSpace_), localre_( reSpace_ )
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
        for (const auto& intersection : intersections( df_.space().gridPart(), entity ) )
        {
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

        const EntityType &neighbor = intersection.outside();

        auto uGuard = Dune::Fem::bindGuard( uOutside_, neighbor );
        const auto& uOutside = uOutside_;

        const int enOrder = df_.space().order( entity );
        const int nbOrder = df_.space().order( neighbor );

        const int quadOrder = 2 * std::max( enOrder, nbOrder ) + 1;

        IntersectionQuadratureType interQuad( df_.space().gridPart(), intersection, quadOrder );
        const auto& quadInside  = interQuad.inside();
        const auto& quadOutside = interQuad.outside();
        const int numQuadraturePoints = quadInside.nop();

        // obtain all required function values on intersection
        std::vector< URangeType > uValuesEn( numQuadraturePoints );
        std::vector< URangeType > uValuesNb( numQuadraturePoints );
        localdf_.evaluateQuadrature( quadInside, uValuesEn );
        uOutside.evaluateQuadrature( quadOutside, uValuesNb );

        typedef QuadratureContext< EntityType, IntersectionType, typename IntersectionQuadratureType::FaceQuadratureType > ContextType;
        typedef LocalEvaluation< ContextType, std::vector< URangeType >, std::vector< URangeType >  > LocalEvaluationType;

        ContextType cLeft( entity, intersection, quadInside, entity.geometry().volume() );
        ContextType cRight( neighbor, intersection, quadOutside, neighbor.geometry().volume() );

        LocalEvaluationType left( cLeft, uValuesEn, uValuesEn );
        LocalEvaluationType right( cRight, uValuesNb, uValuesNb );

        oper_.lifting( left, right, uValuesEn, uValuesNb, localre_ );
      }

      const DF&                        df_;
      const Operator&                  oper_;
      ConstLocalFunction< DF >         localdf_;
      ConstLocalFunction< DF >         uOutside_;
      const DiscreteFunctionSpaceType& reSpace_;
      LiftingFunctionType              localre_;
    };

    template <class SigmaLocalType>
    struct SigmaLocalFunction /*: public Fem::LocalFunctionAdapterHasInitialize*/
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

    typedef Dune::Fem::LocalFunctionAdapter< SigmaLocal<DiscreteFunctionType, DGOperatorType> >
                                                                    SigmaEstimateFunctionType;
    typedef SigmaLocal<DiscreteFunctionType, DGOperatorType>        SigmaLocalType;
    typedef SigmaLocalFunction<SigmaLocalType >                     SigmaLocalFunctionType;
    typedef Dune::Fem::LocalFunctionAdapter<SigmaLocalFunctionType> SigmaLocalFunctionAdapterType;

    void update ()
    {
      interpolate( sigmaEstimateFunction_, sigmaDiscreteFunction_ );
      // sigmaDiscreteFunction_.communicate();
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

    GridPartType&                    gridPart_;
    const DiscreteFunctionType&      solution_;
    const DGOperatorType&            assembler_;
    SigmaDiscreteFunctionSpaceType   sigmaSpace_;
    SigmaDiscreteFunctionType        sigmaDiscreteFunction_;

    SigmaLocal<DiscreteFunctionType, DGOperatorType>
                                     sigmaLocalEstimate_;
    SigmaLocalFunctionType           sigmaLocalFunction_;
    SigmaLocalFunctionAdapterType    sigma_;
    SigmaEstimateFunctionType        sigmaEstimateFunction_;
  };



  template< class DiscreteFunctionSpaceImp, int polOrder, class SigmaEstimatorImp >
  class PAdaptivity
  {
    struct PolOrderStructure
    {
      // set polynomial order to 2 by default
      PolOrderStructure() : val_( -1 ) {}
      explicit PolOrderStructure(int init) : val_( init ) {}
      int value() const
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
    typedef SigmaEstimatorImp                               SigmaEstimatorType;
    typedef typename SigmaEstimatorType::ErrorEstimatorType ErrorEstimatorType;

    typedef typename ErrorEstimatorType::DGOperatorType     DGOperatorType;

    //typedef typename SigmaEstimatorType::ContainerType      ContainerType;

    typedef PersistentContainer<GridType,PolOrderStructure> PolOrderContainer;

    template< class ContainerImp >
    PAdaptivity( const std::shared_ptr< ContainerImp >& cont, DGOperatorType& assembler, const std::string name = ""  )
      : polOrderContainer_( (*cont)(_0)->solution()->space().gridPart().grid(), 0 ),
        space_( (*cont)(_0)->solution()->space() ),
        sigmaEstimator_( cont, assembler, name ),
        errorEstimator_( *(*cont)(_0)->solution(), assembler, sigmaEstimator_.sigma() ),
        param_( AdaptationParameters() )
    {
#ifdef PADAPTSPACE
      // we start with max order
      typedef typename PolOrderContainer::Iterator Iterator ;
      const Iterator end = polOrderContainer_.end();
      const int minimalOrder = errorEstimator_.minimalOrder() ;
      for( Iterator it = polOrderContainer_.begin(); it != end; ++it )
      {
        (*it).value() = minimalOrder ;
      }
#endif
    }

    bool adaptive() const
    {
      return errorEstimator_.isPadaptive();
    }

    void prepare()
    {
#ifdef PADAPTSPACE
      const int minimalOrder = errorEstimator_.minimalOrder();
      // only implemented for PAdaptiveSpace
      std::vector<int> polOrderVec( space_.gridPart().indexSet().size(0) );
      std::fill( polOrderVec.begin(), polOrderVec.end(), polOrder );

      polOrderContainer_.resize();
      if ( errorEstimator_.isPadaptive() )
      {
        for( const auto& entity : elements( space_.gridPart() ) )
        {
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
#ifdef PADAPTSPACE
      double tolerance = param_.refinementTolerance();
      // resize container
      polOrderContainer_.resize();

      const double error = errorEstimator_.estimate( problem );
      std::cout << "ESTIMATE: " << error << std::endl;

      for( const auto& entity : elements( space_.gridPart() ) )
        polOrderContainer_[entity].value() = errorEstimator_.newOrder( 0.98*tolerance, entity );

      return (error < std::abs(tolerance) ? false : errorEstimator_.mark( 0.98 * tolerance));
#else
      return false;
#endif
    }
    void closure()
    {
#ifdef PADAPTSPACE
      errorEstimator_.closure();
#endif
    }

    ErrorEstimatorType& errorEstimator()
    {
      return errorEstimator_;
    }

    SigmaEstimatorType& sigmaEstimator()
    {
      return sigmaEstimator_;
    }

  private:

    PolOrderContainer                polOrderContainer_;
    const DiscreteFunctionSpaceType& space_;
    SigmaEstimatorType               sigmaEstimator_;
    ErrorEstimatorType               errorEstimator_;
    AdaptationParameters             param_;
  };

  /**
   * \brief Adaptation indicator doing no indication and marking of the entities.
   */
  class NoPAdaptIndicator
  {
  public:
   typedef uint64_t                          UInt64Type;

    template< class... Args >
    NoPAdaptIndicator( Args&& ...a )
    {}

    bool adaptive() const
    {
      return false;
    }

    size_t numberOfElements() const
    {
      return 0;
    }

    UInt64Type globalNumberOfElements() const { return 0; }

    template< class TimeProviderImp >
    void setAdaptation( TimeProviderImp& tp ) {}

    void preAdapt() {}

    void estimateMark( const bool initialAdapt = false ) {}

    void postAdapt(){}

    void finalize(){}

    int minNumberOfElements() const { return 0; }

    int maxNumberOfElements() const { return 0; }

    int finestLevel() const { return 0; }

    // return some info
   // const typename SigmaEstimatorType::SigmaDiscreteFunctionType& sigma()
   // {
   //   return pAdapt_.sigmaEstimator().sigma();
   // }

  };


  /**
   * \brief Adaptation indicator doing no indication and marking of the entities.
   */
  template< class PAdaptivityImp, class ModelImp >
  class PAdaptIndicator
  {
    protected:
    typedef PAdaptivityImp                                               PAdaptivityType;
    typedef ModelImp                                                     ModelType;
    typedef typename ModelType::ProblemType                              ProblemType;
    typedef typename PAdaptivityType::SigmaEstimatorType                 SigmaEstimatorType;
    typedef typename PAdaptivityType::ErrorEstimatorType                 ErrorEstimatorType;
    typedef typename SigmaEstimatorType::DiscreteFunctionType            DiscreteFunctionType;
    typedef typename SigmaEstimatorType::DGOperatorType                  DGOperatorType;
    //typedef typename DGOperatorType::ContainerType                       ContainerType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::GridPartType                  GridPartType;
    typedef typename GridPartType::GridType                              GridType;
    static const int polOrder = SigmaEstimatorType::polynomialOrder;

    enum { dimension = GridType::dimension  };

  public:
   typedef uint64_t                          UInt64Type;

    template< class ContainerImp >
    PAdaptIndicator( const std::shared_ptr< ContainerImp >& cont,
                     DGOperatorType& assembler,
                     const ModelType& model,
                     const std::string name = "" )
      : pAdapt_( cont, assembler, name ),
        problem_( model.problem() )
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
      pAdapt_.sigmaEstimator().update();

      const bool marked = pAdapt_.estimateMark( problem_ );
      if( marked )
        pAdapt_.closure();
    }

    void postAdapt(){}

    void finalize(){}

    int minNumberOfElements() const { return 0; }

    int maxNumberOfElements() const { return 0; }

    int finestLevel() const { return 0; }

    // return some info
    const typename SigmaEstimatorType::SigmaDiscreteFunctionType& sigma()
    {
      return pAdapt_.sigmaEstimator().sigma();
    }

  private:

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

    // type of linear operator (i.e. matrix implementation)
    typedef typename BaseType::OperatorType::type          OperatorType;

    // The DG space operator
    typedef typename BaseType::OperatorType::AssemblerType AssemblerType;

    //type for a standalone container
    typedef typename AssemblerType::ContainerType          ContainerType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType        DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;

    // type of inverse operator (i.e. linear solver implementation)
    typedef typename BaseType::SolverType::type            SolverType;

    enum { dimension = GridType::dimension  };

    typedef typename BaseType::TimeProviderType            TimeProviderType;

    typedef typename BaseType::AdaptIndicatorType          AdaptIndicatorType;

    typedef typename BaseType::AdaptationDiscreteFunctionType AdaptationDiscreteFunctionType;

    using BaseType::name;
    using BaseType::grid;
    using BaseType::gridWidth;
    using BaseType::gridSize;
    using BaseType::solution;
    using BaseType::rhs;
    using BaseType::exactSolution;
    using BaseType::solver;

  protected:
    using BaseType::model_;

  public:
    template< class ContainerImp, class ExtraArgsImp >
    SubEllipticAlgorithm( const std::shared_ptr< ContainerImp >& cont,
                          const std::shared_ptr< ExtraArgsImp >& extra )
    : BaseType( cont, extra ),
      assembler_( cont, extra, model_ ),
      matrix_( (*cont)(_0,_0)->matrix() ),
      adaptIndicator_( std::make_unique<AdaptIndicatorType>( cont, assembler_, model_, name() ) ),
      step_( 0 )
    {
      std::string gridName = Fem::gridName( grid() );
      if( gridName == "ALUGrid" || gridName == "ALUConformGrid" || gridName == "ALUSimplexGrid" )
      {
        if( solution().space().begin() != solution().space().end() )
        {
          if( solution().space().begin()->type().isSimplex() && solution().space().order() > 2 && solution().space().continuous() && GridType::dimension > 2 )
          {
            std::cerr << std::endl<< "ERROR: Lagrange spaces for p>2 do not work on simplex grids due to the twist problem !!!" << std::endl << std::endl;
          }
        }
      }
    }

    void virtual setTime ( const double time ) override
    {
      model_.setTime( time );
    }

    const AssemblerType& assembler () const
    {
      return assembler_;
    }

    AssemblerType& assembler ()
    {
      return assembler_;
    }

    //ADAPTATION
    virtual AdaptIndicatorType* adaptIndicator() override
    {
      return adaptIndicator_.get();
    }
    virtual AdaptationDiscreteFunctionType* adaptationSolution () override
    {
      return &solution();
    }

  protected:
    virtual std::shared_ptr< SolverType > doCreateSolver() override
    {
      Dune::Timer timer;
      timer.start();

      if( solution().space().continuous() ) // Lagrange case
      {
        typedef Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType ;
        StencilType stencil( solution().space(), solution().space() );
        matrix().reserve( stencil );
      }
      else // DG case
      {
        typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType ;
        StencilType stencil( solution().space(), solution().space() );
        matrix().reserve( stencil );
      }

      assembler_.assemble();
      std::cout << "Solver (Poisson) assemble time: " << timer.elapsed() << std::endl;

      assembler_.testSymmetry();

      SolverParameter solverParameter;
      return std::make_shared< SolverType >( matrix(), solverParameter );
    }

    //! default time loop implementation, overload for changes
    virtual void doPreSolve ( const int loop ) override
    {
      BaseType::doPreSolve( loop );
    }

    //! finalize computation by calculating errors and EOCs
    virtual void doFinalize ( const int loop ) override
    {
      model_.eocErrors( solution(), adaptIndicator()->sigma() );
    }


  protected:
    const OperatorType& matrix () const
    {
      assert( matrix_ );
      return *matrix_;
    }

    OperatorType& matrix ()
    {
      assert( matrix_ );
      return *matrix_;
    }

    AssemblerType                                assembler_;

    std::shared_ptr< OperatorType >              matrix_;
    std::unique_ptr< AdaptIndicatorType >        adaptIndicator_;
    int                                          step_;
  };

}
}
#endif // FEMHOWTO_STEPPER_HH
