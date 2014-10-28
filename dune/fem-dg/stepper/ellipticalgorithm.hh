#ifndef DUNE_FEMDG_ELLIPTIC_ALGORITHM_HH
#define DUNE_FEMDG_ELLIPTIC_ALGORITHM_HH
#include <config.h>

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem-dg/misc/dgnorm.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem/operator/projection/dgl2projection.hh>

#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/misc/gridname.hh>

// include local header files
#include "base.hh"
#include "baseevolution.hh"


using namespace Dune;                                        


template <class GridImp,
          class ProblemTraits, 
          int polOrd>
struct ElliptTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                        GridType;

  // Choose a suitable GridView
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType >    GridPartType;

  // problem dependent types 
  typedef typename ProblemTraits :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: FluxType         FluxType;
  static const Dune :: DGDiffusionFluxIdentifier DiffusionFluxId = 
    ProblemTraits :: template Traits< GridPartType > :: PrimalDiffusionFluxId ;
  static const int dimRange = InitialDataType :: dimRange ;

  typedef PassTraits<ModelType,dimRange,polynomialOrder>      PassTraitsType;
  typedef typename PassTraitsType::DestinationType            DiscreteFunctionType;
  typedef typename PassTraitsType::LinearOperatorType         LinearOperatorType;
  typedef typename PassTraitsType::LinearInverseOperatorType  LinearInverseOperatorType;
  
  typedef DGPrimalMatrixAssembly<DiscreteFunctionType,ModelType,FluxType > DgType;

  // ... as well as the Space type
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType   DiscreteSpaceType;

  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};


template <class GridImp,
          class ProblemTraits, 
          int polynomialOrder>             
class EllipticAlgorithm
{
public:

  // my traits class 
  typedef ElliptTraits< GridImp, ProblemTraits, polynomialOrder> Traits ;

  // type of Grid
  typedef typename Traits :: GridType                    GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType                GridPartType;

  // initial data type 
  typedef typename Traits :: InitialDataType             InitialDataType;

  // An analytical version of our model
  typedef typename Traits :: ModelType                   ModelType;

  // The flux for the discretization of advection terms
  typedef typename Traits :: FluxType                    FluxType;

  // The DG space operator
  typedef typename Traits :: DgType                      DgType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename Traits :: DiscreteFunctionType        DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename Traits :: DiscreteSpaceType           DiscreteSpaceType;

  // type of linear operator (i.e. matrix implementation)
  typedef typename Traits :: LinearOperatorType          LinearOperatorType;

  // type of inverse operator (i.e. linear solver implementation)
  typedef typename Traits :: LinearInverseOperatorType   LinearInverseOperatorType;

  enum { dimension = GridType :: dimension  };
#if STOKES
  typedef typename DiscreteSpaceType ::  
   template ToNewDimRange< dimension*dimension > :: NewFunctionSpaceType SigmaFunctionSpaceType ;
#else
  typedef typename DiscreteSpaceType ::  
   template ToNewDimRange< dimension > :: NewFunctionSpaceType SigmaFunctionSpaceType ;
#endif

  //- hybrid spaces use PAdaptiveDGSpace
  template <class Grid, int topoId> 
  struct SigmaSpaceChooser 
  {
    typedef Fem :: PAdaptiveDGSpace< SigmaFunctionSpaceType, GridPartType, polynomialOrder > Type ;
  };

  //- cubes use LegendreDGSpace
  template <class Grid> 
  struct SigmaSpaceChooser< Grid, 1 >
  {
    typedef Dune::Fem:: LegendreDiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polynomialOrder > Type ;
  };

  //- cubes use OrthonormalDGSpace
  template <class Grid> 
  struct SigmaSpaceChooser< Grid, 0 >
  {
    typedef Dune::Fem:: DiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polynomialOrder > Type ;
  };

  // work arround internal compiler error 
  enum { simplexTopoId =  GenericGeometry :: SimplexTopology< dimension > :: type :: id,
         cubeTopoId    =  GenericGeometry :: CubeTopology< dimension > :: type :: id,
         myTopo        =  Capabilities::hasSingleGeometryType< GridType > :: topologyId  
       };
    
  enum { topoId = (simplexTopoId == myTopo) ? 0 :  // 0 = simplex, 1 = cube, -1 = hybrid  
                    (myTopo == cubeTopoId) ? 1 : -1 };

  typedef typename SigmaSpaceChooser< GridType, topoId > :: Type  SigmaDiscreteFunctionSpaceType ;

  typedef Dune::Fem::AdaptiveDiscreteFunction< SigmaDiscreteFunctionSpaceType >
    SigmaDiscreteFunctionType;

  // ... as well as the Space type
  //typedef typename PassTraits<ModelType, GridType :: dimension, polynomialOrder>::DestinationType
  //  SigmaDiscreteFunctionType;

  //typedef typename SigmaDiscreteFunctionType :: DiscreteFunctionSpaceType
  //  SigmaDiscreteFunctionSpaceType;

  // compute the function sigma = grad u + sum_e r_e
  template <class DF, class Operator>
  struct SigmaLocal : public Fem :: LocalFunctionAdapterHasInitialize
  {
    typedef typename DF::DiscreteFunctionSpaceType UDFS;
    typedef typename UDFS::GridPartType GridPartType;
    typedef typename GridPartType::GridType::template Codim<0>::Entity EntityType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType IntersectionType;

    typedef typename Operator::FluxType::LiftingFunctionType LiftingFunctionType;
    typedef typename LiftingFunctionType::RangeType RangeType;
    typedef typename LiftingFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

    typedef typename DF::RangeType URangeType;
    typedef typename DF::JacobianRangeType UJacobianRangeType;

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
      typedef typename Operator :: FaceQuadratureType  FaceQuadratureType ;
      typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
      typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;
      const typename EntityType::EntityPointer pOutside = intersection.outside();
      const EntityType &outside = *pOutside;
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
    const DF &df_;
    const Operator &oper_;
    typename DF::LocalFunctionType localdf_;
    const DiscreteFunctionSpaceType &reSpace_;
    LiftingFunctionType localre_;
  };
  
  typedef Dune::Fem::LocalFunctionAdapter< SigmaLocal<DiscreteFunctionType, DgType> > SigmaEstimateFunction;
  //typedef Estimator1< DiscreteFunctionType, SigmaDiscreteFunctionType, DgType > EstimatorType;

  //typedef Dune::Fem::LocalFunctionAdapter< EstimatorType > EstimateDataType;
  //typedef tuple< const DiscreteFunctionType*, const SigmaEstimateFunction*, const EstimateDataType* >  IOTupleType;
  typedef tuple< const DiscreteFunctionType* > IOTupleType;
  typedef Dune::Fem::DataWriter<GridType,IOTupleType> DataWriterType;

  template <class SigmaLocalType>
  struct SigmaLocalFunction : public Fem :: LocalFunctionAdapterHasInitialize
  {
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
                     DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::RangeType RangeType;
    typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionType::EntityType EntityType;
    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

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
    const DiscreteFunctionType &u_;
    typename DiscreteFunctionType::LocalFunctionType uLocal_;
    const SigmaDiscreteFunctionType &q_;
    typename SigmaDiscreteFunctionType::LocalFunctionType qLocal_;
    SigmaLocalType sigmaLocal_;
  };


  //---- Local Restriction and Prolongation Operator -------------------------
  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
  //---- Adaptation Manager --------------------------------------------------
  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;
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
  typedef PersistentContainer<GridType,PolOrderStructure> PolOrderContainer;

  // type of statistics monitor 
  typedef SolverMonitor  SolverMonitorType ;

public:
  EllipticAlgorithm(GridType& grid) :
    grid_( grid ),
    gridPart_( grid_ ),
    problem_( ProblemTraits::problem() ),
    model_( new ModelType( problem() ) ),
    convectionFlux_( *model_ ),
    dgOperator_(gridPart_, *model_),
    invDgOperator_(),
    linDgOperator_(),
    space_( const_cast<DiscreteSpaceType &> (dgOperator_.space()) ),
    solution_("solution", space_ ),
    rhs_("rhs", space_ ),
    sigmaSpace_( gridPart_ ),
    sigmaDiscreteFunction_( "sigma", sigmaSpace_ ),
    sigmaLocalEstimate_( solution_, dgOperator_ ),
    sigmaEstimateFunction_( "function 4 estimate", sigmaLocalEstimate_, gridPart_, space_.order() ),
    //estimator_( solution_, sigmaDiscreteFunction_, dgOperator_, grid ),
    //estimateData_( "estimator", estimator_, gridPart_, space_.order() ),
    ioTuple_( &solution_ ), // &sigmaEstimateFunction_), &estimateData_ ),
    polOrderContainer_(grid_ , 0),
    eocId_( -1 ),
    step_( 0 )
  {
    std::string filename( Dune::Fem::Parameter::commonOutputPath() );
    filename += "/run.gnu";
    runFile_.open( filename.c_str() );
    if( ! runFile_.is_open() ) 
    {
      std::cerr << filename << "runfile not open" << std::endl;
    }

    runFile_ << "# h | elements | CPU time | iter | l_min | l_max | cond  | L2 error" << std::endl;

    std::string name = Fem :: gridName( grid_ );
    if( name == "ALUGrid" || name == "ALUConformGrid" || name == "ALUSimplexGrid" )
    {
      if( space_.begin() != space_.end() )
      {
        if( space_.begin()->type().isSimplex() && space_.order() > 2 && space_.continuous() && GridType :: dimension > 2 )
        {
          std::cerr << std::endl<< "ERROR: Lagrange spaces for p>2 do not work on simplex grids due to the twist problem !!!" << std::endl << std::endl;
        }
      }
    }

    solution_.clear();

    // we start with max order
    typedef typename PolOrderContainer :: Iterator Iterator ;
    const Iterator end = polOrderContainer_.end();
    const int minimalOrder = solution_.space().order(); // estimator_.minimalOrder() ;
    for( Iterator it = polOrderContainer_.begin(); it != end; ++it ) 
    {
      (*it).value() = minimalOrder ;
    }

    const std::string eocDescription[] = { "$L^2$-error", "DG-error", "sigma-norm" };
    eocId_ = Dune::Fem::FemEoc::addEntry( eocDescription, 3);
  }                                                                        /*@LST1E@*/

  GridType& grid () { return grid_; }

  IOTupleType dataTuple()
  {
    // tuple with additionalVariables
    return IOTupleType( &solution_ );
  }

  //! return reference to discrete space 
  DiscreteSpaceType & space() { return space_; }                    /*@LST0E@*/

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

    latexInfo = dgOperator_.description();

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
  SolverMonitorType solve( const int loop ) 
  {
    SolverMonitorType monitor;
    numbers_.resize( 0 );

    // calculate grid width
    const double h = Dune::Fem::GridWidth::calcGridWidth(gridPart_);
    numbers_.push_back( h );

    const double size = grid_.size(0);
    numbers_.push_back( size );
    
    //assert( solution_.space().size() > 0 );


    //Dune::Fem::DGL2ProjectionImpl :: project( problem(), solution_ );
    //return ;
#if 0
#ifdef PADAPTSPACE
    int polOrder = Dune::Fem:: Parameter::getValue<double>("femhowto.polynomialOrder",1);
    const int minimalOrder = estimator_.minimalOrder() ;

    // only implemented for PAdaptiveSpace
    std::vector<int> polOrderVec( space_.gridPart().indexSet().size(0) );
    std::fill( polOrderVec.begin(), polOrderVec.end(), polOrder );

    polOrderContainer_.update();
    if ( estimator_.isPadaptive() )
    {
      typedef typename DiscreteSpaceType::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType ;
      const IteratorType end = space_.end();
      for( IteratorType it = space_.begin(); it != end; ++it )
      {
        const EntityType& entity = *it;
        int order = polOrderContainer_[ entity ].value(); 

        // use constrcutor here, operator = is deprecated 
        typename EntityType::EntityPointer hit ( it );
        while (order == -1) // is a new element 
        {
          if ( entity.level() == 0) 
          {
            order = minimalOrder;
          }
          else
          {
            hit = hit->father();
            // don't call father twice 
            order = polOrderContainer_[ *hit ].value();
            assert(order > 0);
          }
        }
        polOrderVec[space_.gridPart().indexSet().index(entity)] = order;
      }
    }
    space_.adapt( polOrderVec );
#endif
#endif
    
    if (!invDgOperator_)
    {
      linDgOperator_.reset( new LinearOperatorType("dg operator", space_, space_ ) );

#if DGSCHEME // for all dg schemes including pdg (later not working)
      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> StencilType ;
#else
      typedef Dune::Fem::DiagonalStencil<DiscreteSpaceType,DiscreteSpaceType> StencilType ;
#endif  
      StencilType stencil( space_, space_ );

      linDgOperator_->reserve( stencil );
      linDgOperator_->clear();
      dgOperator_.assemble(0, *linDgOperator_, rhs_);
			dgOperator_.testSymmetrie(*linDgOperator_);

      double absLimit   = Dune::Fem:: Parameter::getValue<double>("istl.absLimit",1.e-10);
      double reduction  = Dune::Fem:: Parameter::getValue<double>("istl.reduction",1.e-10);
      // this describes the factor for max iterations in terms of spaces size
      int maxIterFactor = Dune::Fem:: Parameter::getValue<double>("istl.maxiterfactor", int(-1) );

      // if max iterations is negative set it to twice the spaces size
      if( maxIterFactor < 0 ) 
        maxIterFactor = 2; 

      invDgOperator_.reset( new LinearInverseOperatorType(*linDgOperator_, reduction, absLimit ));
      // invDgOperator_.reset( new InverseOperatorType(*linDgOperator_, reduction, absLimit, maxIterFactor * space_.size() ));
      (*invDgOperator_)(rhs_,solution_);
      monitor.ils_iterations = invDgOperator_->iterations();
    }

    // calculate new sigma 
    Dune::Fem:: DGL2ProjectionImpl :: project( sigmaEstimateFunction_, sigmaDiscreteFunction_ );

    //Dune::Fem::DGL2ProjectionImpl :: project( problem(), solution_ );
    return monitor;
  }


  //! finalize computation by calculating errors and EOCs 
  void finalize( const int eocloop )
  {
    typedef typename InitialDataType :: ExactSolutionType  ExactSolutionType ;
    //---- Adapter for exact solution ------------------------------------------
    typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >
        GridExactSolutionType;

    // create grid function adapter 
    GridExactSolutionType ugrid( "exact solution", problem().exactSolution(), gridPart_, 1 );

    // calculate L2 - Norm 
    Dune::Fem::L2Norm< GridPartType > l2norm( gridPart_ );
    const double l2error = l2norm.distance( ugrid, solution_ );

    numbers_.push_back( l2error );

    for(size_t i=0; i<numbers_.size(); ++i)
      runFile_ << numbers_[ i ] << " ";

    runFile_ << std::endl;

    Dune::Fem::DGNorm< GridPartType > dgnorm( gridPart_ );
    const double dgerror = dgnorm.distance( ugrid, solution_ );

    Dune::Fem::H1Norm< GridPartType > sigmanorm( gridPart_ );
    typedef SigmaLocalFunction<SigmaLocal<DiscreteFunctionType, DgType> >
      SigmaLocalFunctionType;
    SigmaLocalFunctionType sigmaLocalFunction( solution_, sigmaDiscreteFunction_, sigmaLocalEstimate_ );
    Dune::Fem::LocalFunctionAdapter<SigmaLocalFunctionType> sigma( "sigma function", sigmaLocalFunction, gridPart_, space_.order() );
    const double sigmaerror = sigmanorm.distance( ugrid, sigma );

    // store values 
    std::vector<double> errors;
    errors.push_back( l2error );
    errors.push_back( dgerror );
    errors.push_back( sigmaerror );

    // submit error to the FEM EOC calculator 
    Dune::Fem::FemEoc :: setErrors(eocId_, errors);

    // delete solver and linear operator for next step
    delete invDgOperator_.release();
    delete linDgOperator_.release();
  }

  bool adaptation(const double tolerance)
  {
    // resize container 
    polOrderContainer_.resize();

    const double error = 0.0;//estimator_.estimate( problem() );
    std::cout << "ESTIMATE: " << error << std::endl;
    typedef typename DiscreteSpaceType::IteratorType IteratorType;
    const IteratorType end = space_.end();
    for( IteratorType it = space_.begin(); it != end; ++it )
    {
      const typename IteratorType::Entity &entity = *it;
      //polOrderContainer_[entity].value() = 
      //  estimator_.newOrder( 0.98*tolerance, entity );
    }
    return false ;
    //return (error < std::abs(tolerance) ? false : estimator_.mark( 0.98 * tolerance));
  }
  void closure() 
  {
    //estimator_.closure();
  }

  const InitialDataType& problem() const { assert( problem_ ); return *problem_; }
  const ModelType& model() const { assert( model_ ); return *model_; }
  const DiscreteFunctionType& solution() const { return solution_; }
  DiscreteFunctionType& solution() { return solution_; }
  IOTupleType& ioTuple() { return ioTuple_; }
  const DgType &oper() const
  {
    return dgOperator_;
  }
  protected:

  GridType& grid_;
  GridPartType gridPart_;       // reference to grid part, i.e. the leaf grid 
  // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  std::unique_ptr< const InitialDataType > problem_;
  std::unique_ptr< ModelType >             model_;
  FluxType                convectionFlux_;
  DgType                  dgOperator_;

  std::unique_ptr< LinearInverseOperatorType > invDgOperator_;
  std::unique_ptr< LinearOperatorType        > linDgOperator_;

  DiscreteSpaceType &space_;    // the discrete function space
  // the solution 
  DiscreteFunctionType   solution_, rhs_;

  SigmaDiscreteFunctionSpaceType sigmaSpace_;
  SigmaDiscreteFunctionType sigmaDiscreteFunction_;

  SigmaLocal<DiscreteFunctionType, DgType> sigmaLocalEstimate_;
  SigmaEstimateFunction sigmaEstimateFunction_;
  //EstimatorType estimator_;
  //EstimateDataType estimateData_;
  IOTupleType ioTuple_;
  PolOrderContainer polOrderContainer_;

  // Initial flux for advection discretization (UpwindFlux)
  int  eocId_;
  int  step_;

  std::vector<double> numbers_;
  std::ofstream runFile_;
};
#endif // FEMHOWTO_STEPPER_HH
