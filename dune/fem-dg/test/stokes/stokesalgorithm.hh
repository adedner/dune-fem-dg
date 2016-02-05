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
#include <dune/fem-dg/test/stokes/stokesassembler.hh>

namespace Dune
{

namespace Fem
{

  template< class GridPartImp, class DiscreteVelocityFunctionImp, class DiscretePressureFunctionImp, class SigmaFunctionSpaceImp, class AssemblerImp, class ModelImp, int polOrder>
  class StokesSigmaEstimator: public PoissonSigmaEstimator<GridPartImp, DiscreteVelocityFunctionImp, SigmaFunctionSpaceImp, AssemblerImp, polOrder >
  {
    typedef PoissonSigmaEstimator<GridPartImp, DiscreteVelocityFunctionImp, SigmaFunctionSpaceImp, AssemblerImp, polOrder > BaseType;

  public:
    typedef DiscreteVelocityFunctionImp     DiscreteVelocityFunctionType;
    typedef DiscretePressureFunctionImp     DiscretePressureFunctionType;
    typedef SigmaFunctionSpaceImp           SigmaFunctionSpaceType;
    typedef GridPartImp                     GridPartType;
    typedef typename GridPartType::GridType GridType;
    typedef AssemblerImp                    AssemblerType;
    typedef ModelImp                        ModelType;

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
      typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename GridPartType::IntersectionType IntersectionType;

      typedef typename AssemblerType::FaceQuadratureType   FaceQuadratureType;
      typedef typename AssemblerType::VolumeQuadratureType VolumeQuadratureType;

      StokesFlux(const DiscretePressureFunctionType &p, const AssemblerType &oper)
      : p_(p),
        oper_(oper)
        {}
      template <class Quadrature,class Value,class DValue,class RetType, class DRetType>
      void flux(const GridPartType &gridPart,
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
        oper_.flux(gridPart,intersection,entity,neighbor,time,faceQuadInside,faceQuadOutside,
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
      const AssemblerType &oper_;
    };



    typedef Dune::Fem::LocalFunctionAdapter< SigmaEval<DiscreteVelocityFunctionType,DiscretePressureFunctionType,AssemblerType> > StokesEstimateFunction;
    typedef StokesErrorEstimator< DiscreteVelocityFunctionType, StokesEstimateFunction, StokesFlux > StokesEstimatorType;
    typedef Dune::Fem::LocalFunctionAdapter< StokesEstimatorType >          StokesEstimateDataType;


    StokesSigmaEstimator( GridPartType& gridPart,
                          const DiscreteVelocityFunctionType& solution,
                          const DiscretePressureFunctionType& pressureSolution,
                          const AssemblerType& assembler,
                          const std::string keyPrefix = "" )
    : BaseType( gridPart, solution, assembler, keyPrefix ),
      stkFlux_( pressureSolution, assembler_),
      stkLocalEstimate_( solution_, pressureSolution, assembler_ ),
      stkEstimateFunction_("stokes estimate", stkLocalEstimate_, gridPart_, solution_.space().order() ),
      stkEstimator_( solution_, stkEstimateFunction_, stkFlux_, gridPart_.grid() ),
      stkEstimateData_("stokesEstimator", stkEstimator_, gridPart_, solution_.space().order() )
    {}


    StokesEstimatorType& estimator() const
    {
      return stkEstimator_;
    };

    StokesEstimateDataType& data() const
    {
      return stkEstimateData_;
    }


  private:
    StokesFlux stkFlux_;
    SigmaEval<DiscreteVelocityFunctionType,DiscretePressureFunctionType,AssemblerType> stkLocalEstimate_;

    StokesEstimateFunction stkEstimateFunction_;
    StokesEstimatorType  stkEstimator_;
    StokesEstimateDataType  stkEstimateData_;


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

    typedef typename DiscreteVelocityFunctionSpaceType ::
      template ToNewDimRange< dimension * ElliptProblemTraits::AnalyticalTraits::ModelType::dimRange >::NewFunctionSpaceType
                                                                    SigmaFunctionSpaceType;

    typedef StokesSigmaEstimator< GridPartType, DiscreteVelocityFunctionType, DiscreteFunctionType,
                                  SigmaFunctionSpaceType, typename EllipticalAlgorithmType::AssemblerType,
                                  typename ElliptProblemTraits::AnalyticalTraits::ModelType, polOrd >
                                                                    StokesSigmaEstimatorType;

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
      space_( container_.template adapter<1>().space() ),
      ellAlg_( grid, container_.template adapter<0>() ),
      assembler_( container_, model() ),
      ioTuple_( new IOTupleType( *BaseType::dataTuple(), *ellAlg_.dataTuple() ) ),
      stokesSigmaEstimator_( new StokesSigmaEstimatorType( container_.template adapter<1>().gridPart(), ellAlg_.solution(), solution(), ellAlg_.assembler(), name() ) )
    {
    }

    virtual IOTupleType& dataTuple ()
    {
      return *ioTuple_;
    }

    bool adaptation(const double tolerance)
    {
#if 0
#if PADAPTSPACE
      // update to current
      polOrderContainer_.resize();

      const double error = stokesSigmaEstimator_.estimator().estimate( problem() );
      typedef typename DiscreteVelocityFunctionSpaceType::IteratorType IteratorType;
      const IteratorType end = ellAlg_.space().end();
      for( IteratorType it = ellAlg_.space().begin(); it != end; ++it )
      {
        const typename IteratorType::Entity &entity = *it;
        polOrderContainer_[entity].value() =
          stokesSigmaEstimator_.estimator().newOrder( 0.98*tolerance, entity );
      }
#endif
#endif
      return (error < std::abs(tolerance) ? false : stokesSigmaEstimator_->estimator().mark( 0.98 * tolerance));
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
#if 0
#ifdef PADAPTSPACE
      int polOrder = Dune::Fem::Parameter::getValue<double>("femdg.polOrd",1);
      // only implemented for PAdaptiveSpace
      std::vector<int> polOrderVec( ellAlg_.space().gridPart().indexSet().size(0) );
      std::vector<int> polOrderVecPressure( space_.gridPart().indexSet().size(0) );
      std::fill( polOrderVec.begin(), polOrderVec.end(), polOrder );
      std::fill( polOrderVecPressure.begin(), polOrderVecPressure.end(), polOrder-1 );
      ellAlg_.space().adapt( polOrderVec );
      space_.adapt( polOrderVecPressure);
#endif
#endif
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

      // TODO check wheather we need the following line
      stokesSigmaEstimator_->update();
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
      AnalyticalTraits::addEOCErrors( solution(), ellAlg_.model(), problem().get<1>().exactSolution() );
    }

  protected:
    ContainerType&                              container_;
    const DiscreteFunctionSpaceType&            space_;
    AssemblerType                               assembler_;

    EllipticalAlgorithmType                     ellAlg_;
    std::unique_ptr< IOTupleType >              ioTuple_;
    std::unique_ptr< StokesSigmaEstimatorType > stokesSigmaEstimator_;
  };


}
}
#endif // FEMHOWTO_STEPPER_HH
