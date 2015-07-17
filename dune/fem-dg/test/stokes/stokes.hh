#ifndef STOKES_ALGORITHM_HH
#define STOKES_ALGORITHM_HH
#include <config.h>

#ifdef LOCALDEBUG
static double sum_ = 0.;
static double sum2_ = 0.;
static double localMaxRatio_ = 0.;
static double localMinRatio_ = 1e+100;
static double maxRatioOfSums = 0.;
static double minRatioOfSums = 1e+100;
#endif

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
// include local header files
#include "stokesassembler.hh"
#include <dune/fem-dg/solver/uzawa.hh>
#include <dune/fem-dg/stepper/base.hh>

// #include "../base/baseevolution.hh"
//#include "ellipt.hh"
#include <dune/fem-dg/stepper/ellipticalgorithm.hh>
#include <dune/fem-dg/operator/adaptation/stokesestimator.hh>

using namespace Dune;


template <class GridImp,
          class ProblemTraits,
          int polOrd>
struct StokesTraits : public ElliptTraits<GridImp, ProblemTraits,polOrd>
{

public:
  typedef ElliptTraits<GridImp, ProblemTraits,polOrd> BaseType;
  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::ModelType ModelType;
  typedef typename ProblemTraits :: ProblemType ProblemType;
  typedef typename ProblemType::ExactPressureType ExactPressureType;
  typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, 1 >   PressureFunctionSpaceType;


#if DGSCHEME
#if ONB
    #warning using DG space with ONB
    typedef Dune::Fem::DiscontinuousGalerkinSpace
#elif LAG
    #warning using DG space with Lagrange base functions
    typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace
#elif LEG
    #warning using DG space with Legendre base functions
    typedef Dune::Fem::LegendreDiscontinuousGalerkinSpace
#else
    #warning using p-adaptive DG space
    typedef Dune::Fem::PAdaptiveDGSpace
#define PADAPTSPACE
#endif
#else
#warning using p-adaptive Lagrange space
    typedef Fem::PAdaptiveLagrangeSpace
#define PADAPTSPACE
#endif
    < PressureFunctionSpaceType,GridPartType, ( polOrd > 0 ) ? polOrd-1 : 0, Dune::Fem::CachingStorage>          DiscretePressureSpaceType;

#if WANT_ISTL
	typedef  Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscretePressureSpaceType >   DiscretePressureFunctionType;
#else
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscretePressureSpaceType >   DiscretePressureFunctionType;
#endif
	typedef StokesAssembler< DiscreteFunctionType,
			   DiscretePressureFunctionType, typename BaseType::OperatorTraits>
  StokesAssemblerType;


};


template <class GridImp,
          class ProblemTraits,
					int polynomialOrder>
class StokesAlgorithm : public EllipticAlgorithm<GridImp,ProblemTraits,polynomialOrder>
{

  typedef EllipticAlgorithm<GridImp,ProblemTraits,polynomialOrder> BaseType;


  // my traits class
  typedef StokesTraits< GridImp, ProblemTraits, polynomialOrder> Traits ;
public:
  // type of Grid
  typedef typename Traits :: GridType                 GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType             GridPartType;

  // initial data type
  typedef typename Traits :: InitialDataType          InitialDataType;

	typedef typename Traits::ExactPressureType           ExactPressureType;


  // An analytical version of our model
  typedef typename Traits :: ModelType                 ModelType;

  // The flux for the discretization of advection terms
  typedef typename Traits :: FluxType                  FluxType;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  //typedef typename Traits :: DgOperatorType            DgOperatorType;
  typedef typename Traits :: DgAssembledOperatorType   DgAssembledOperatorType;
  typedef typename Traits :: StokesAssemblerType       StokesAssemblerType;

  // The discrete function for the unknown solution is defined in the DgOperator

  typedef typename Traits :: DiscreteFunctionType              DiscreteFunctionType;
  typedef typename Traits :: DiscretePressureFunctionType      DiscretePressureFunctionType;
  // ... as well as the Space type
  typedef typename Traits :: DiscreteSpaceType                 DiscreteSpaceType;
  typedef typename Traits :: DiscretePressureSpaceType         DiscretePressureSpaceType;

	typedef typename BaseType::LinearOperatorType         LinearOperatorType;
	typedef typename BaseType::LinearInverseOperatorType  LinearInverseOperatorType;

  typedef SolverMonitor  SolverMonitorType;


  template<class DiscreteFunction,class DiscretePressureFunction,class Operator>
  class SigmaEval: public BaseType::template SigmaLocal<DiscreteFunction,Operator>
  {
  public:
    typedef typename BaseType::template SigmaLocal<DiscreteFunction,Operator> SigmaBaseType;
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

    typedef typename DgAssembledOperatorType::FaceQuadratureType   FaceQuadratureType;
    typedef typename DgAssembledOperatorType::VolumeQuadratureType VolumeQuadratureType;

    StokesFlux(const DiscretePressureFunctionType &p, const DgAssembledOperatorType &oper)
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
    const DgAssembledOperatorType &oper_;
  };

  typedef Dune::Fem::LocalFunctionAdapter< SigmaEval<DiscreteFunctionType,DiscretePressureFunctionType,DgAssembledOperatorType> > StokesEstimateFunction;
  typedef Dune::StokesErrorEstimator< DiscreteFunctionType, StokesEstimateFunction, StokesFlux > StokesEstimatorType;

  typedef Dune::Fem::LocalFunctionAdapter< StokesEstimatorType >          StokesEstimateDataType;
  typedef tuple< const DiscreteFunctionType*, const DiscretePressureFunctionType*, const StokesEstimateDataType* >  IOTupleType;
  typedef Dune::Fem::DataOutput<GridType,IOTupleType>         DataWriterType;

  //---- Local Restriction and Prolongation Operator -------------------------
  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
  //---- Adaptation Manager --------------------------------------------------
  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

public:
  explicit StokesAlgorithm(GridType& grid, std::string moduleName = "" ) :
    BaseType(grid),
    pressurespace_( gridPart_ ),
    pressuresolution_("pressuresolution", pressurespace_ ),
    stkFlux_(pressuresolution_,dgAssembledOperator_),
    stkLocalEstimate_(solution_,pressuresolution_,dgAssembledOperator_),
    stkEstimateFunction_("stokes estimate",stkLocalEstimate_,gridPart_,space_.order()),
    stkEstimator_( solution_, stkEstimateFunction_,stkFlux_, grid),
    stkEstimateData_("stokesEstimator",stkEstimator_,gridPart_,space_.order()),
    ioTuple_( &solution_, &pressuresolution_,&stkEstimateData_ ),
    stokesAssembler_(space_ , pressurespace_, problem()),
		averageIter_(0),
    eocpId_(-1)
  {
    std::string filename(Dune::Fem:: Parameter::commonOutputPath() );
    filename += "/run.gnu";
    runFile_.open( filename.c_str() );
    if( ! runFile_.is_open() )
      {
        std::cerr << filename << "runfile not open" << std::endl;
	      abort();
      }
    runFile_ << "# h | elements | CPU time | iter | l_min | l_max | cond  | L2 error" << std::endl;
    const std::string eocpDescription[]={"$L^2$-p-error","Average Iterations"};
    eocpId_ = Dune::Fem::FemEoc::addEntry(eocpDescription,2);


    std::string name = Fem :: gridName( grid_ );
    if( name == "ALUGrid" || name == "ALUConformGrid" || name == "ALUSimplexGrid" )
    {
      if( space_.begin() != space_.end() )
      {
        if( space_.begin()->type().isSimplex() && space_.order() > 2 && space_.continuous() )
        {
          std::cerr << std::endl<< "ERROR: Lagrange spaces for p>2 do not work on simplex grids due to the twist problem !!!";
        }
      }
    }
  }                                                                        /*@LST1E@*/

  //! return reference to discrete space
  DiscreteSpaceType & space() { return space_; }                    /*@LST0E@*/
  DiscretePressureSpaceType & pressurespace() { return space_; }
  //! returns data prefix for EOC loops ( default is loop )
  virtual std::string dataPrefix() const
  {
    return problem_->dataPrefix();
  }

  // gather information from the space operator, the time integratior
  // and the problem to output before each table in tex file
  std::string description() const
  {
    std::string latexInfo;

    latexInfo = dgAssembledOperator_.description();

    // latexInfo = dgAdvectionOperator_.description()
    //            + dgDiffusionOperator_.description();

    std::stringstream odeInfo;

    latexInfo += odeInfo.str()
      + "\n"
      + problem_->description()
      + "\n\n";

    return latexInfo;
  }


  //! default time loop implementation, overload for changes
  SolverMonitorType solve( int loop )
  {
    numbers_.resize( 0 );

    // calculate grid width
    const double h = Dune::Fem::GridWidth::calcGridWidth(gridPart_);
    numbers_.push_back( h );


    const double size = grid_.size(0);
    numbers_.push_back( size );

#ifdef PADAPTSPACE
    int polOrder = Dune::Fem::Parameter::getValue<double>("femdg.polynomialOrder",1);
		// only implemented for PAdaptiveSpace
    std::vector<int> polOrderVec( space_.gridPart().indexSet().size(0) );
    std::vector<int> polOrderVecPressure( pressurespace_.gridPart().indexSet().size(0) );
    std::fill( polOrderVec.begin(), polOrderVec.end(), polOrder );
		std::fill( polOrderVecPressure.begin(), polOrderVecPressure.end(), polOrder-1 );
		space_.adapt( polOrderVec );
		pressurespace_.adapt( polOrderVecPressure);
#endif


		typedef Dune::UzawaSolver< DiscreteFunctionType,DiscretePressureFunctionType,StokesAssemblerType,LinearInverseOperatorType >UzawaType;

		DiscretePressureFunctionType pressure("press",pressurespace_);
		DiscreteFunctionType veloprojection("veloprojection",space_);
		pressure.clear();
		veloprojection.clear();
    double absLimit   = Dune::Fem:: Parameter::getValue<double>("istl.absLimit",1.e-10);
    double reduction  = Dune::Fem:: Parameter::getValue<double>("istl.reduction",1.e-10);
    double uzawareduction  = Dune::Fem:: Parameter::getValue<double>("uzawareduction",reduction*100.);

#if WANT_ISTL
    linDgOperator_.reset( new LinearOperatorType("dg operator", space_, space_ ) );

#if DGSCHEME // for all dg schemes including pdg (later not working)
      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> StencilType ;
#else
      typedef Dune::Fem::DiagonalStencil<DiscreteSpaceType,DiscreteSpaceType> StencilType ;
#endif
    StencilType stencil( space_, space_);

    SolverMonitorType monitor;
    monitor.gridWidth = h; // space_.size();

    linDgOperator_->reserve(stencil);
		linDgOperator_->clear();
		dgAssembledOperator_.assemble(0, *linDgOperator_, rhs_);
		invDgOperator_.reset( new LinearInverseOperatorType(*linDgOperator_, reduction, absLimit ) );
#else
// 			{
abort();
				invDgOperator_.reset( new LinearInverseOperatorType(dgAssembledOperator_, 1e-12, 1e-12, step_++ ) );

// 			}
#endif

		stokesAssembler_.assemble( *problem_);

#if WANT_ISTL
				UzawaType uzawa(stokesAssembler_,*invDgOperator_,rhs_,uzawareduction,uzawareduction,100000);
#else
				UzawaType uzawa(stokesAssembler_,*invDgOperator_,uzawareduction, uzawareduction,100000);
#endif
		pressuresolution_.clear();
		uzawa(stokesAssembler_.pressureRhs(),pressuresolution_);

		solution_.assign(uzawa.velocity());
		monitor.ils_iterations = uzawa.iterations();
		averageIter_=uzawa. averageLinIter();
	  // calculate new sigma
    Dune::Fem:: DGL2ProjectionImpl :: project( sigmaEstimateFunction_, sigmaDiscreteFunction_ );

	  return monitor;
  }



  //! finalize computation by calculating errors and EOCs
  void finalize( const int eocloop )
  {
    //---- Adapter for exact solution ------------------------------------------
    typedef typename InitialDataType :: ExactSolutionType ExactSolutionType;
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >
      GridExactSolutionType;

    typedef Dune::Fem::GridFunctionAdapter< ExactPressureType, GridPartType >  GridExactPressureType;

		// create grid function adapter
    GridExactSolutionType ugrid( "exact solution", problem().exactSolution(), gridPart_, 1 );
		GridExactPressureType pgrid( "exact pressure", problem().exactPressure_, gridPart_,1);

		  // calculate L2 - Norm
		Dune::Fem::L2Norm< GridPartType > l2norm( gridPart_ );
		const double l2error = l2norm.distance( ugrid, solution_ );
    numbers_.push_back( l2error );
		const double l2press = l2norm.distance( pgrid,pressuresolution_);
    Dune::Fem::DGNorm< GridPartType > dgnorm( gridPart_ );
    double dgerror = dgnorm.distance( ugrid, solution_ );
    //double dgerror = dgnorm.norm( ugrid );

    //		dgerror+=l2press;
    Dune::Fem::H1Norm< GridPartType > sigmanorm( gridPart_ );

   typedef typename BaseType::template SigmaLocalFunction<typename BaseType::template SigmaLocal<DiscreteFunctionType,DgAssembledOperatorType> >
      SigmaLocalFunctionType;
    SigmaLocalFunctionType sigmaLocalFunction( solution_, sigmaDiscreteFunction_, sigmaLocalEstimate_ );
    Dune::Fem::LocalFunctionAdapter<SigmaLocalFunctionType> sigma( "sigma function", sigmaLocalFunction, gridPart_, space_.order() );
   const double sigmaerror = sigmanorm.distance( ugrid, sigma );




    for(size_t i=0; i<numbers_.size(); ++i)
      runFile_ << numbers_[ i ] << " ";

    runFile_ << std::endl;


    // store values
    std::vector<double> errors;
    errors.push_back( l2error );
    errors.push_back( dgerror );
    errors.push_back( sigmaerror );
    std::vector<double> perr;
    perr.push_back(l2press);
		perr.push_back(averageIter_);
		// submit error to the FEM EOC calculator
    Dune::Fem::FemEoc :: setErrors(eocId_, errors);
    Dune::Fem::FemEoc :: setErrors(eocpId_,perr);

		invDgOperator_.reset();
    linDgOperator_.reset();
  }


  bool adaptation(const double tolerance)
  {
    // update to current
    polOrderContainer_.resize();

    const double error = stkEstimator_.estimate( problem() );
    typedef typename DiscreteSpaceType::IteratorType IteratorType;
    const IteratorType end = space_.end();
    for( IteratorType it = space_.begin(); it != end; ++it )
    {
      const typename IteratorType::Entity &entity = *it;
      polOrderContainer_[entity].value() =
        stkEstimator_.newOrder( 0.98*tolerance, entity );
    }
    return (error < std::abs(tolerance) ? false : stkEstimator_.mark( 0.98 * tolerance));
  }




//   bool adaptation(const double tolerance)
//   {
//     const double error = stkEstimator_.estimate( problem() );
//         return (error < std::abs(tolerance) ? false : stkEstimator_.mark( 0.98 * tolerance));
//   }

  const InitialDataType& problem() const { assert( problem_ ); return *problem_; }
  const DiscreteFunctionType& solution() const { return solution_; }
  DiscreteFunctionType& solution() { return solution_; }
  IOTupleType& ioTuple() { return ioTuple_; }
  IOTupleType dataTuple() { return IOTupleType( ioTuple_ ); }

private:
  using BaseType::grid_;
  using BaseType::gridPart_;       // reference to grid part, i.e. the leaf grid
  using BaseType::problem_;
	using BaseType::model;
	using BaseType::dgAssembledOperator_;
	using BaseType::invDgOperator_;
#if WANT_ISTL
  using BaseType::linDgOperator_;
#endif
	using BaseType::space_;    // the discrete function space
  using BaseType::rhs_;    // the discrete function space
	using BaseType::sigmaSpace_;
  using BaseType::sigmaDiscreteFunction_;
  using BaseType::sigmaLocalEstimate_;
  using BaseType::sigmaEstimateFunction_;
  using BaseType::solution_;
	using BaseType::polOrderContainer_;
  using BaseType::eocId_;
	using BaseType::step_;
	using BaseType::numbers_;
  using BaseType::runFile_;

	DiscretePressureSpaceType pressurespace_;
  DiscretePressureFunctionType pressuresolution_;

  StokesFlux stkFlux_;
  SigmaEval<DiscreteFunctionType,DiscretePressureFunctionType,DgAssembledOperatorType> stkLocalEstimate_;

  StokesEstimateFunction stkEstimateFunction_;
  StokesEstimatorType  stkEstimator_;
  StokesEstimateDataType  stkEstimateData_;
  IOTupleType ioTuple_;
	StokesAssemblerType         stokesAssembler_;

  // Initial flux for advection discretization (UpwindFlux)

	mutable double averageIter_;
	int eocpId_;
};
#endif // FEMHOWTO_STEPPER_HH
