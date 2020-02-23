#ifndef STOKES_ALGORITHM_HH
#define STOKES_ALGORITHM_HH

#ifndef NDEBUG
// enable fvector and fmatrix checking
#define DUNE_ISTL_WITH_CHECKING
#endif

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>
#include <dune/fem/operator/common/stencil.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/adaptation/stokesestimator.hh>
#include <dune/fem-dg/algorithm/sub/elliptic.hh>
#include <dune/fem-dg/solver/uzawa.hh>
#include <dune/fem-dg/misc/tupleutility.hh>

// include local header files
#include <dune/fem-dg/examples/stokes/stokesassembler.hh>

#include <dune/fem-dg/algorithm/sub/containers.hh>


namespace Dune
{
namespace Fem
{


  template< class PoissonErrorEstimatorImp, class AssemblerImp >
  class StokesSigmaEstimator
    : public PoissonSigmaEstimator< PoissonErrorEstimatorImp >
  {
    typedef PoissonSigmaEstimator< PoissonErrorEstimatorImp >        BaseType;

  public:
    typedef typename BaseType::DiscreteFunctionType                  DiscreteVelocityFunctionType;
    typedef typename BaseType::DGOperatorType                        BaseAssemblerType;
    typedef AssemblerImp                                             AssemblerType;
    typedef typename AssemblerType::DiscretePressureFunctionType     DiscretePressureFunctionType;
    typedef typename BaseType::DiscreteFunctionSpaceType             DiscreteFunctionSpaceType;
    typedef typename BaseType::GridPartType                          GridPartType;
    typedef typename BaseType::GridType                              GridType;
    typedef typename BaseType::DGOperatorType                        DGOperatorType;
    //typedef typename BaseType::AssemblerType                         AssemblerType;
    static const int polynomialOrder = BaseType::polynomialOrder;


    typedef typename AssemblerType::ModelType                        ModelType;

    using BaseType::sigmaSpace_;
    using BaseType::sigmaDiscreteFunction_;
    using BaseType::sigmaLocalFunction_;
    using BaseType::sigmaLocalEstimate_;
    using BaseType::assembler_;
    using BaseType::gridPart_;
    using BaseType::solution_;

    template<class DiscreteFunction,class DiscretePressureFunction,class Operator>
    class SigmaEval: public BaseType::SigmaLocalType
    {
      typedef typename BaseType::SigmaLocalType SigmaBaseType;
    public:
      typedef DiscretePressureFunction DiscretePressureFunctionType;
      typedef typename DiscretePressureFunctionType::RangeType PRangeType;

    private:
      const DiscretePressureFunctionType& ph_;
      typename DiscretePressureFunctionType::LocalFunctionType localp_;

    public:
      SigmaEval(const DiscreteFunction &uh,
                const DiscretePressureFunctionType &ph,
                const Operator& oper)
        : SigmaBaseType(uh,oper),
          ph_(ph),
          localp_(ph_)
      {}
      SigmaEval(const SigmaEval &other)
      : SigmaBaseType(other), ph_(other.ph_), localp_(ph_)
      {}
      void init(const typename SigmaBaseType::EntityType &en)
      {
        SigmaBaseType::init(en);
        localp_.init(en);
      }
      template< class PointType >
      void evaluate(const PointType& x,typename SigmaBaseType::RangeType& res) const
      {
        SigmaBaseType::evaluate(x,res);
        PRangeType p;
        localp_.evaluate(x,p);
        Dune::Fem::FieldMatrixConverter< typename SigmaBaseType::RangeType, typename DiscreteFunction::JacobianRangeType> res1( res );
        for(int i=0;i<res1.rows;++i)
          res1[i][i] -= p;
      }
    };


    struct StokesFlux
    {
      typedef typename GridPartType::GridType::template Codim<0>::Entity EntityType;
      typedef typename GridPartType::IntersectionIteratorType            IntersectionIteratorType;
      typedef typename GridPartType::IntersectionType                    IntersectionType;

      typedef typename AssemblerType::FaceQuadratureType                 FaceQuadratureType;
      typedef typename AssemblerType::VolumeQuadratureType               VolumeQuadratureType;

      typedef typename AssemblerType::ContainerType                    ContainerType;


      StokesFlux(const DiscretePressureFunctionType &p, const BaseAssemblerType &oper)
      : p_(p),
        oper_(oper)
        {}
      template <class Quadrature,class Value,class DValue,class RetType, class DRetType>
      void numericalFlux(const GridPartType &gridPart,
                         const IntersectionType &intersection,
                         const EntityType &entity, const EntityType &neighbor,
                         const double time,
                         const Quadrature &faceQuadInside, const Quadrature &faceQuadOutside,
                         const Value &valueEn, const DValue &dvalueEn,
                         const Value &valueNb, const DValue &dvalueNb,
                         RetType &retEn, DRetType &dretEn,
                         RetType &retNb, DRetType &dretNb) const
      {
        //\hat{K}
        oper_.numericalFlux(gridPart,intersection,entity,neighbor,time,faceQuadInside,faceQuadOutside,
                            valueEn,dvalueEn,valueNb,dvalueNb,retEn,dretEn,retNb,dretNb);
        typename DiscretePressureFunctionType::LocalFunctionType pEn = p_.localFunction(entity);
        typename DiscretePressureFunctionType::LocalFunctionType pNb = p_.localFunction(neighbor);
        std::vector< typename DiscretePressureFunctionType::RangeType > pValuesEn( faceQuadInside.nop() );
        std::vector< typename DiscretePressureFunctionType::RangeType > pValuesNb( faceQuadOutside.nop() );

        pEn.evaluateQuadrature( faceQuadInside, pValuesEn );
        pNb.evaluateQuadrature( faceQuadOutside, pValuesNb );

        assert(retEn.size() == faceQuadInside.nop() );

        for (unsigned int i=0;i<retEn.size();++i)
        {
          typename IntersectionType::GlobalCoordinate normal = intersection.integrationOuterNormal( faceQuadInside.localPoint(i) );
          double value=0.5*(pValuesEn[i]+pValuesNb[i]);
          normal*=value;

          retEn[i]+=normal;
          retNb[i]+=normal;

        }

      }

      template<class Quadrature,class RetType>
      void boundaryValues(const GridPartType &gridPart,
                          const IntersectionType &intersection,
                          const EntityType &entity,
                          const double time,
                          const Quadrature &faceQuadInside,
                          RetType &retEn) const
      {
        oper_.boundaryValues(gridPart,
                             intersection, entity, time, faceQuadInside,
                             retEn);
      }


      template<class Quadrature,class Value,class DValue,class RetType, class DRetType>
      void boundaryFlux(const GridPartType &gridPart,
                const IntersectionType &intersection,
                const EntityType &entity,
                const double time,
                const Quadrature &faceQuadInside,
                const Value &valueEn, const DValue &dvalueEn,
                const Value &valueNb,
                RetType &retEn, DRetType &dretEn) const

      {
        oper_.boundaryFlux(gridPart,
                           intersection, entity, time, faceQuadInside,
                           valueEn, dvalueEn, valueNb,
                           retEn, dretEn );
        typename DiscretePressureFunctionType::LocalFunctionType pEn = p_.localFunction(entity);
        std::vector< typename DiscretePressureFunctionType::RangeType > pValuesEn( faceQuadInside.nop() );

        pEn.evaluateQuadrature( faceQuadInside, pValuesEn );
        assert(retEn.size() == faceQuadInside.nop() );

        for (unsigned int i=0;i<retEn.size();++i)
        {
          typename IntersectionType::GlobalCoordinate normal = intersection.integrationOuterNormal( faceQuadInside.localPoint(i) );
          normal*=pValuesEn[i];

          retEn[i]+=normal;

        }
      }
      const ModelType &model() const
      {
        return oper_.model();
      }
      private:
      const DiscretePressureFunctionType &p_;
      const BaseAssemblerType &oper_;
    };

    typedef StokesFlux                                                    StokesFluxType;

    typedef Dune::Fem::LocalFunctionAdapter< SigmaEval<DiscreteVelocityFunctionType,DiscretePressureFunctionType,BaseAssemblerType> > StokesEstimateFunction;
    typedef StokesErrorEstimator< DiscreteVelocityFunctionType, StokesEstimateFunction, StokesFluxType > StokesErrorEstimatorType;

    typedef Dune::Fem::LocalFunctionAdapter< StokesErrorEstimatorType >                                  StokesEstimateDataType;

    typedef typename AssemblerType::ContainerType                    ContainerType;

    template< class ContainerImp >
    StokesSigmaEstimator( const std::shared_ptr< ContainerImp >& cont,
                          const BaseAssemblerType& ass,
                          const AssemblerType& assembler,
                          const std::string keyPrefix = "" )
    : BaseType( (*cont)(std::make_tuple(_0),std::make_tuple(_0) ), ass, keyPrefix ),
      stkFlux_( *(*cont)(_1)->solution(), ass ),
      stkLocalEstimate_( *(*cont)(_0)->solution(), *(*cont)(_1)->solution(), ass ),
      stkEstimateFunction_("stokes estimate", stkLocalEstimate_, (*cont)(_0)->solution()->gridPart(), (*cont)(_0)->solution()->space().order() ),
      stkEstimator_( *(*cont)(_0)->solution(), stkFlux_, stkEstimateFunction_ ),
      stkEstimateData_("stokesEstimator", stkEstimator_, (*cont)(_0)->solution()->gridPart(), (*cont)(_0)->solution()->space().order() )
    {}


    StokesErrorEstimatorType& estimator() const
    {
      return stkEstimator_;
    };

    StokesEstimateDataType& data() const
    {
      return stkEstimateData_;
    }


  private:
    StokesFluxType stkFlux_;
    SigmaEval<DiscreteVelocityFunctionType,DiscretePressureFunctionType,BaseAssemblerType> stkLocalEstimate_;

    StokesEstimateFunction stkEstimateFunction_;
    StokesErrorEstimatorType  stkEstimator_;
    StokesEstimateDataType  stkEstimateData_;


  };


  template< class PAdaptivityImp, class DiscreteFunctionSpaceImp, int polOrder, class SigmaEstimatorImp >
  class StokesPAdaptivity
  {
  public:
    typedef PAdaptivityImp                                      PAdaptivityType;

    typedef DiscreteFunctionSpaceImp                            DiscreteFunctionSpaceType;
    typedef SigmaEstimatorImp                                   SigmaEstimatorType;
    typedef typename SigmaEstimatorType::ErrorEstimatorType     ErrorEstimatorType;

    typedef typename DiscreteFunctionSpaceType::GridType        GridType;

    typedef typename SigmaEstimatorType::DGOperatorType         DGOperatorType;

    typedef typename SigmaEstimatorType::AssemblerType          AssemblerType;
    typedef typename SigmaEstimatorType::BaseAssemblerType      BaseAssemblerType;

    typedef typename SigmaEstimatorType::ContainerType          ContainerType;

    typedef typename PAdaptivityType::PolOrderContainer         PolOrderContainer;

    template< class ContainerImp >
    StokesPAdaptivity( const std::shared_ptr< ContainerImp >& cont, BaseAssemblerType& ass, AssemblerType& assembler, const std::string name = ""  )
      : pAdapt_( (*cont)(std::make_tuple(_0),std::make_tuple(_0) ), ass, name ),
        polOrderContainer_( (*cont)(_0)->solution()->gridPart().grid(), 0 ),
        space_( (*cont)(_1)->solution()->space() ),
        sigmaEstimator_( cont, ass, assembler, name ),
        errorEstimator_( *(*cont)(_0)->solution(), ass, sigmaEstimator_.sigma() ),
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
      return errorEstimator_.isPadaptive() && pAdapt_.adaptive();
    }

    void prepare()
    {
      pAdapt_.prepare();
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
      return pAdapt_.estimateMark( problem );
      //note: no h-Adaptation regarding pressure space!?
    }
    void closure()
    {
      pAdapt_.closure();
      //note: no closure regarding pressure space!?
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

    PAdaptivityType                  pAdapt_;
    PolOrderContainer                polOrderContainer_;
    const DiscreteFunctionSpaceType& space_;
    SigmaEstimatorType               sigmaEstimator_;
    ErrorEstimatorType               errorEstimator_;
    AdaptationParameters             param_;

  };



  /**
   * \brief Adaptation indicator doing no indication and marking of the entities.
   */
  template< class PAdaptivityImp, class ModelImp >
  class StokesPAdaptIndicator
  {

  protected:
    typedef PAdaptivityImp                                               PAdaptivityType;
    typedef ModelImp                                                     ModelType;
    typedef typename ModelType::ProblemType                              ProblemType;
    typedef typename PAdaptivityType::SigmaEstimatorType                 SigmaEstimatorType;
    typedef typename PAdaptivityType::ErrorEstimatorType                 ErrorEstimatorType;
    typedef typename SigmaEstimatorType::DiscreteFunctionType            DiscreteFunctionType;
    typedef typename SigmaEstimatorType::AssemblerType                   AssemblerType;
    typedef typename SigmaEstimatorType::BaseAssemblerType               BaseAssemblerType;
    typedef typename SigmaEstimatorType::DGOperatorType                  DGOperatorType;
    typedef typename AssemblerType::ContainerType                        ContainerType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::GridPartType                  GridPartType;
    typedef typename GridPartType::GridType                              GridType;
    static const int polOrder = SigmaEstimatorType::polynomialOrder;

    enum { dimension = GridType::dimension  };

  public:
   typedef uint64_t                          UInt64Type;

    template< class ContainerImp >
    StokesPAdaptIndicator( const std::shared_ptr< ContainerImp >& cont,
                           BaseAssemblerType& ass,
                           AssemblerType& assembler,
                           const ModelType& model,
                           const std::string name = "" )
      : pAdapt_( cont, ass, assembler, name ),
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

    const int finestLevel() const { return 0; }

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
   *  \brief Algorithm for solving the Stokes equation.
   *
   *  \ingroup SubAlgorithms
   */
  template <class GridImp, class ProblemTraits, class ElliptProblemTraits, int polOrd >
  class SubStokesAlgorithm : public SubSteadyStateAlgorithm<GridImp,ProblemTraits,polOrd>
  {

    typedef SubSteadyStateAlgorithm<GridImp,ProblemTraits,polOrd> BaseType;

  public:
    typedef typename ElliptProblemTraits::template Algorithm<polOrd> EllipticalAlgorithmType;

    // type of Grid
    typedef typename BaseType::GridType                             GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                         GridPartType;

    typedef typename BaseType::OperatorType::AssemblerType          AssemblerType;

    typedef typename EllipticalAlgorithmType::DiscreteTraits::DiscreteFunctionType
                                                                    DiscreteVelocityFunctionType;
    typedef typename BaseType::DiscreteFunctionType                 DiscreteFunctionType;
    // ... as well as the Space type
    typedef typename DiscreteVelocityFunctionType::DiscreteFunctionSpaceType
                                                                    DiscreteVelocityFunctionSpaceType;
    typedef typename BaseType::DiscreteFunctionSpaceType            DiscreteFunctionSpaceType;

    typedef typename BaseType::SolverMonitorType                    SolverMonitorType;

    // type of inverse operator (i.e. linear solver implementation)
    typedef typename BaseType::SolverType::type                     SolverType;

    enum { dimension = GridType::dimension  };

    typedef typename BaseType::AdaptIndicatorType                   AdaptIndicatorType;

    typedef typename BaseType::AdaptationDiscreteFunctionType       AdaptationDiscreteFunctionType;

    typedef typename BaseType::TimeProviderType                     TimeProviderType;

    using BaseType::model_;
    using BaseType::name;
    using BaseType::grid;
    using BaseType::rhs;
    using BaseType::rhs_;
    using BaseType::gridWidth;
    using BaseType::gridSize;
    using BaseType::solution;
    using BaseType::solver_;
    using BaseType::exactSolution_;

    typedef typename AssemblerType::ContainerType                   ContainerType;

  public:

    template< class ContainerImp, class ExtraArgsImp >
    explicit SubStokesAlgorithm( const std::shared_ptr< ContainerImp >& cont,
                                 const std::shared_ptr< ExtraArgsImp >& extra )
    : BaseType( (*cont)(std::make_tuple(_1),std::make_tuple(_1 ) ), extra ),
      space_( (*cont)(_1)->solution()->space() ),
      assembler_( cont, extra, model_ ),
      ellAlg_( (*cont)(std::make_tuple(_0),std::make_tuple(_0)), extra ),
      stokesSolver_( std::make_shared< SolverType >( *cont, ellAlg_ ) ),
      adaptIndicator_( std::make_unique<AdaptIndicatorType>( cont, ellAlg_.assembler(), assembler_, model_.problem(), BaseType::name() ) )
    {
    }

    void virtual setTime ( const double time ) override
    {
      model_.setTime( time );
      ellAlg_.setTime( time );
    }


  private:
    virtual std::shared_ptr< SolverType > doCreateSolver() override
    {
      Dune::Timer timer;
      timer.start();

      //if( rhsOperator_ ) //rhs by external rhs operator
      //  rhsOperator_( solution(), rhs() );
      rhs().clear();
      assembler_.assemble();

      std::cout << "Solver (Stokes) assemble time: " << timer.elapsed() << std::endl;

      stokesSolver_->maxIterations( 3*ellAlg_.solution().space().size() );
      return stokesSolver_;
      //return std::make_shared< SolverType >( *container_, *ellAlg_.solver(), absLimit, );
    }

    virtual void doInitialize ( const int loop ) override
    {
      ellAlg_.initialize( loop );
      BaseType::doInitialize( loop );
    }

    virtual void doPreSolve( const int loop ) override
    {
      ellAlg_.preSolve( loop );
      BaseType::doPreSolve( loop );
    }

    virtual void doSolve( const int loop ) override
    {
      BaseType::doSolve( loop );
    }

    virtual void doPostSolve( const int loop ) override
    {
      ellAlg_.postSolve( loop );
      BaseType::doPostSolve( loop );
    }

    //! finalize computation by calculating errors and EOCs
    virtual void doFinalize( const int loop ) override
    {
      ellAlg_.finalize( loop );

      //add error
      model_.eocErrors( solution() );

      BaseType::doFinalize( loop );
    }

  protected:
    //std::shared_ptr< ContainerType >       container_;
    const DiscreteFunctionSpaceType&       space_;
    AssemblerType                          assembler_;

    EllipticalAlgorithmType                ellAlg_;
    std::shared_ptr< SolverType >          stokesSolver_;
    std::unique_ptr< AdaptIndicatorType >  adaptIndicator_;
  };


}
}
#endif // FEMHOWTO_STEPPER_HH
