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
#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>
#include <dune/fem-dg/operator/adaptation/stokesestimator.hh>
#include <dune/fem-dg/algorithm/ellipticalgorithm.hh>
#include <dune/fem-dg/solver/uzawa.hh>
#include <dune/fem-dg/solver/linearsolvers.hh>

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
    typedef Dune::StokesErrorEstimator< DiscreteVelocityFunctionType, StokesEstimateFunction, StokesFlux > StokesEstimatorType;
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


  template <class GridImp, class ProblemTraits, class ElliptProblemTraits, int polynomialOrder >
  class StokesAlgorithm : public SteadyStateAlgorithm<GridImp,ProblemTraits,polynomialOrder>
  {

    typedef SteadyStateAlgorithm<GridImp,ProblemTraits,polynomialOrder> BaseType;

  public:
    typedef typename ElliptProblemTraits::template Stepper<polynomialOrder>::Type  EllipticalAlgorithmType;

    // type of Grid
    typedef typename BaseType::GridType                             GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                         GridPartType;

    // initial data type
    typedef typename BaseType::ProblemType                          ProblemType;

    typedef typename BaseType::DiscreteTraits::ExactPressureType    ExactSolutionType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                            ModelType;

    typedef typename BaseType::AssemblerType                        AssemblerType;
    // The discrete function for the unknown solution is defined in the DgOperator

    typedef typename EllipticalAlgorithmType::DiscreteTraits::DiscreteFunctionType
                                                                    DiscreteVelocityFunctionType;
    typedef typename BaseType::DiscreteFunctionType                 DiscreteFunctionType;
    // ... as well as the Space type
    typedef typename EllipticalAlgorithmType::DiscreteTraits::DiscreteFunctionSpaceType
                                                                    DiscreteVelocityFunctionSpaceType;
    typedef typename BaseType::DiscreteFunctionSpaceType            DiscreteFunctionSpaceType;

    typedef typename BaseType::SolverMonitorType                    SolverMonitorType;

    // type of inverse operator (i.e. linear solver implementation)
    typedef typename BaseType::BasicLinearSolverType                BasicLinearSolverType;

    enum { dimension = GridType::dimension  };

    typedef typename DiscreteVelocityFunctionSpaceType ::
      template ToNewDimRange< dimension * ModelType::dimRange >::NewFunctionSpaceType SigmaFunctionSpaceType;

    typedef StokesSigmaEstimator< GridPartType, DiscreteVelocityFunctionType, DiscreteFunctionType,
                                  SigmaFunctionSpaceType, typename EllipticalAlgorithmType::DiscreteTraits::AssemblerType, ModelType, polynomialOrder > StokesSigmaEstimatorType;

    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridExactSolutionType;

    typedef typename BaseType::IOTupleType                                         IOTupleType;
    typedef Dune::Fem::DataOutput<GridType,IOTupleType>         DataWriterType;

    typedef typename  BaseType::AnalyticalTraits    AnalyticalTraits;

    using BaseType::problem;
    using BaseType::model;
    using BaseType::grid;
    using BaseType::gridWidth;
    using BaseType::gridSize;
    using BaseType::eocIds_;
    using BaseType::space;
    using BaseType::solution;
    using BaseType::solver_;

  public:
    explicit StokesAlgorithm(GridType& grid, std::string moduleName = "" ) :
      BaseType( grid, moduleName ),
      ellAlg_( grid, moduleName ),
      gridPart_( grid ),
      space_( gridPart_ ),
      solution_("pressuresolution-" + moduleName, space_ ),
      stokesSigmaEstimator_( gridPart_, ellAlg_.solution(), solution_, ellAlg_.assembler(), moduleName ),
      exact_( "exact pressure-" + moduleName, problem().exactPressure(), gridPart_, 2 ),
      assembler_( ellAlg_.space() , space_, problem() )
    {}

    //! return reference to discrete space
    const DiscreteFunctionSpaceType & space() const { return space_; }

    virtual BasicLinearSolverType* createSolver( DiscreteFunctionType* rhs )
    {
      assembler_.assemble( problem() );

      double absLimit   = Dune::Fem:: Parameter::getValue<double>("istl.absLimit",1.e-10);
      solution_.clear();
#if 0
#ifdef PADAPTSPACE
      int polOrder = Dune::Fem::Parameter::getValue<double>("femdg.polynomialOrder",1);
      // only implemented for PAdaptiveSpace
      std::vector<int> polOrderVec( ellAlg_.space().gridPart().indexSet().size(0) );
      std::vector<int> polOrderVecPressure( space_.gridPart().indexSet().size(0) );
      std::fill( polOrderVec.begin(), polOrderVec.end(), polOrder );
      std::fill( polOrderVecPressure.begin(), polOrderVecPressure.end(), polOrder-1 );
      ellAlg_.space().adapt( polOrderVec );
      space_.adapt( polOrderVecPressure);
#endif
#endif

      return new BasicLinearSolverType( assembler_, ellAlg_.solver(), ellAlg_.rhs(), absLimit, 3*ellAlg_.space().size() );
    }

    virtual DiscreteFunctionType& rhs()
    {
      return assembler_.pressureRhs();
    }

    void solve( const int loop )
    {
      BaseType::solve( loop );

      ellAlg_.solution().assign( solver_->velocity() );

      // TODO check wheather we need the following line
      stokesSigmaEstimator_.update();
    }

    //! finalize computation by calculating errors and EOCs
    void finalize( const int eocloop )
    {
      ellAlg_.finalize( eocloop );
      AnalyticalTraits::addEOCErrors( eocIds_, solution_, ellAlg_.model(), exact_ );
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
      return (error < std::abs(tolerance) ? false : stokesSigmaEstimator_.estimator().mark( 0.98 * tolerance));
    }

    DiscreteFunctionType& solution()
    {
      return solution_;
    }

  private:
    //dummy
    const GridExactSolutionType& exact() const
    {
      return exact_;
    }

  private:
    EllipticalAlgorithmType           ellAlg_;
    GridPartType                      gridPart_;
    DiscreteFunctionSpaceType         space_;
    DiscreteFunctionType              solution_;

    StokesSigmaEstimatorType          stokesSigmaEstimator_;
    GridExactSolutionType             exact_;
    AssemblerType                     assembler_;

  public:

    virtual auto dataTuple () -> decltype( std::tuple_cat( this->ellAlg_.dataTuple(), std::make_tuple( &this->solution(), &this->exact() ) ) )
    {
      return std::tuple_cat( ellAlg_.dataTuple(), std::make_tuple( &this->solution(), &this->exact() ) );
      //return make_tuple( &solution_, &solution_,&stokesSigmaEstimator_.data(), &ugrid_, &pgrid_ );
    }

  };


}
}
#endif // FEMHOWTO_STEPPER_HH
