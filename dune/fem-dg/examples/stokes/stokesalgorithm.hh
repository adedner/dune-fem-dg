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
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
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

    StokesSigmaEstimator( ContainerType& container,
                          const BaseAssemblerType& ass,
                          const AssemblerType& assembler,
                          const std::string keyPrefix = "" )
    : BaseType( container.template adapter<0>(), ass, keyPrefix ),
      container_( container ),
      stkFlux_( *container_.template solution<1>(), ass ),
      stkLocalEstimate_( *container_.template solution<0>(), *container_.template solution<1>(), ass ),
      stkEstimateFunction_("stokes estimate", stkLocalEstimate_, container_.template solution<0>()->gridPart(), container_.template space<0>().order() ),
      stkEstimator_( *container_.template solution<0>(), stkFlux_, stkEstimateFunction_ ),
      stkEstimateData_("stokesEstimator", stkEstimator_, container_.template solution<0>()->gridPart(), container_.template space<0>().order() )
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
    ContainerType& container_;
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

    StokesPAdaptivity( ContainerType& container, BaseAssemblerType& ass, AssemblerType& assembler, const std::string name = ""  )
      : pAdapt_( container.template adapter<0>(), ass, name ),
        polOrderContainer_( container.template solution<0>()->gridPart().grid(), 0 ),
        space_( container.template space<1>() ),
        sigmaEstimator_( container, ass, assembler, name ),
        errorEstimator_( *container.template solution<0>(), ass, sigmaEstimator_.sigma() ),
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
  template< class PAdaptivityImp, class ProblemImp >
  class StokesPAdaptIndicator
  {

  protected:
    typedef PAdaptivityImp                                               PAdaptivityType;
    typedef ProblemImp                                                   ProblemType;
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

    StokesPAdaptIndicator( ContainerType& container,
                           BaseAssemblerType& ass,
                           AssemblerType& assembler,
                           const ProblemType& problem,
                           const std::string name = "" )
      : pAdapt_( container, ass, assembler, name ),
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

    // initial data type
    typedef typename BaseType::ProblemType                          ProblemType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                            ModelType;

    typedef typename BaseType::OperatorType::AssemblerType          AssemblerType;

    typedef typename EllipticalAlgorithmType::DiscreteTraits::DiscreteFunctionType
                                                                    DiscreteVelocityFunctionType;
    typedef typename BaseType::DiscreteFunctionType                 DiscreteFunctionType;
    // ... as well as the Space type
    typedef typename EllipticalAlgorithmType::DiscreteTraits::DiscreteFunctionType::DiscreteFunctionSpaceType
                                                                    DiscreteVelocityFunctionSpaceType;
    typedef typename BaseType::DiscreteFunctionSpaceType            DiscreteFunctionSpaceType;

    typedef typename BaseType::SolverMonitorType                    SolverMonitorType;

    // type of inverse operator (i.e. linear solver implementation)
    typedef typename BaseType::SolverType::type                     SolverType;

    enum { dimension = GridType::dimension  };

    typedef typename BaseType::AdaptIndicatorType                   AdaptIndicatorType;

    typedef typename BaseType::AdaptationDiscreteFunctionType       AdaptationDiscreteFunctionType;

    typedef CovariantTuple< typename BaseType::IOTupleType::type, typename EllipticalAlgorithmType::IOTupleType::type >
                                                                    IOTupleType;

    typedef typename BaseType::TimeProviderType                     TimeProviderType;

    typedef typename BaseType::AnalyticalTraits                     AnalyticalTraits;

    using BaseType::problem;
    using BaseType::model;
    using BaseType::name;
    using BaseType::grid;
    using BaseType::rhs;
    using BaseType::rhs_;
    using BaseType::gridWidth;
    using BaseType::gridSize;
    using BaseType::solution;
    using BaseType::solver_;
    using BaseType::exactSolution_;

    typedef typename EllipticalAlgorithmType::ContainerType         EllipticContainerType;

    typedef typename AssemblerType::ContainerType                   ContainerType;

  public:

    explicit SubStokesAlgorithm( GridType& grid, ContainerType& container ) :
      BaseType( grid, container.template adapter<1>() ),
      container_( container ),
      space_( container_.template space<1>() ),
      assembler_( container_, model() ),
      ellAlg_( grid, container_.template adapter<0>() ),
      adaptIndicator_( Std::make_unique<AdaptIndicatorType>( container_, ellAlg_.assembler(), assembler_, problem(), name() ) ),
      ioTuple_( new IOTupleType( *BaseType::dataTuple(), *ellAlg_.dataTuple() ) ),
      time_(0)
    {
    }

    virtual IOTupleType& dataTuple () override
    {
      return *ioTuple_;
    }

    void virtual setTime ( const double time ) override
    {
      time_ = time;
      assembler_.setTime( time_ );
      ellAlg_.setTime( time_ );
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

      double absLimit = Dune::Fem::Parameter::getValue<double>("istl.absLimit",1.e-10);
      return std::make_shared< SolverType >( container_, *ellAlg_.solver(), absLimit, 3*ellAlg_.solution().space().size() );
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
      AnalyticalTraits::addEOCErrors( solution(), ellAlg_.model(), problem().template get<1>().exactSolution() );
    }

  protected:
    ContainerType&                         container_;
    const DiscreteFunctionSpaceType&       space_;
    AssemblerType                          assembler_;

    EllipticalAlgorithmType                ellAlg_;
    std::unique_ptr< AdaptIndicatorType >  adaptIndicator_;
    std::unique_ptr< IOTupleType >         ioTuple_;
    double                                 time_;
  };


}
}
#endif // FEMHOWTO_STEPPER_HH
