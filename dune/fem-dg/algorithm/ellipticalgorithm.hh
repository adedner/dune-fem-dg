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
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>
#include <dune/fem-dg/assemble/primalmatrix.hh>
#include <dune/fem/operator/projection/dgl2projection.hh>

#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/misc/gridname.hh>

// include local header files
#include "steadystate.hh"




namespace Dune
{
namespace Fem
{


  template< class GridPartImp, class DiscreteFunctionImp, class SigmaFunctionSpaceImp, class AssembledOperatorImp, int polOrder>
  class PoissonSigmaEstimator
  {
  public:

    typedef DiscreteFunctionImp             DiscreteFunctionType;
    typedef SigmaFunctionSpaceImp           SigmaFunctionSpaceType;
    typedef GridPartImp                     GridPartType;
    typedef typename GridPartType::GridType GridType;
    typedef AssembledOperatorImp            AssembledOperatorType;


    PoissonSigmaEstimator( GridPartType& gridPart,
                           const DiscreteFunctionType& solution,
                           const AssembledOperatorType& assembledOperator,
                           const std::string name = "" )
    : gridPart_( gridPart ),
      solution_( solution ),
      assembledOperator_( assembledOperator ),
      sigmaSpace_( gridPart_ ),
      sigmaDiscreteFunction_( "sigma-"+name, sigmaSpace_ ),
      sigmaLocalEstimate_( solution_, assembledOperator_ ),
      sigmaLocalFunction_( solution_, sigmaDiscreteFunction_, sigmaLocalEstimate_ ),
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
      typedef Dune::Fem:: LegendreDiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polOrder > Type ;
    };

    //- cubes use OrthonormalDGSpace
    template <class Grid>
    struct SigmaSpaceChooser< Grid, 0 >
    {
      typedef Dune::Fem:: DiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polOrder > Type ;
    };

    // work arround internal compiler error
    enum { simplexTopoId =  GenericGeometry::SimplexTopology< GridType::dimension >::type::id,
           cubeTopoId    =  GenericGeometry::CubeTopology< GridType::dimension >::type::id,
           myTopo        =  Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
         };

    enum { topoId = (simplexTopoId == myTopo) ? 0 :  // 0 = simplex, 1 = cube, -1 = hybrid
                      (myTopo == cubeTopoId) ? 1 : -1 };

    typedef typename SigmaSpaceChooser< GridType, topoId >::Type  SigmaDiscreteFunctionSpaceType ;

    typedef Dune::Fem::AdaptiveDiscreteFunction< SigmaDiscreteFunctionSpaceType >
      SigmaDiscreteFunctionType;

    // compute the function sigma = grad u + sum_e r_e
    template <class DF, class Operator>
    struct SigmaLocal : public Fem::LocalFunctionAdapterHasInitialize
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
        typedef typename Operator::FaceQuadratureType  FaceQuadratureType ;
        typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
        typedef typename IntersectionQuadratureType::FaceQuadratureType QuadratureImp;

        const EntityType &outside = Dune::Fem::make_entity( intersection.outside() );

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

    typedef Dune::Fem::LocalFunctionAdapter< SigmaLocal<DiscreteFunctionType, AssembledOperatorType> > SigmaEstimateFunction;

    template <class SigmaLocalType>
    struct SigmaLocalFunction : public Fem::LocalFunctionAdapterHasInitialize
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

    typedef SigmaLocal<DiscreteFunctionType, AssembledOperatorType> SigmaLocalEstimateType;

    typedef SigmaLocalFunction<SigmaLocal<DiscreteFunctionType, AssembledOperatorType> >
      SigmaLocalFunctionType;

    typedef Dune::Fem::LocalFunctionAdapter<SigmaLocalFunctionType> SigmaLocalFunctionAdapterType;
                                            //StokesEstimatorType   StokesEstimateDataType;

    void update ()
    {
      Dune::Fem::DGL2ProjectionImpl::project( sigmaEstimateFunction_, sigmaDiscreteFunction_ );
    }

    SigmaLocalFunctionAdapterType* sigma ()
    {
      return new SigmaLocalFunctionAdapterType( "sigma function", sigmaLocalFunction_, gridPart_, solution_.space().order() );
    }


  public:

    GridPartType& gridPart_;
    const DiscreteFunctionType& solution_;
    const AssembledOperatorType& assembledOperator_;
    SigmaDiscreteFunctionSpaceType sigmaSpace_;
    SigmaDiscreteFunctionType sigmaDiscreteFunction_;

    SigmaLocal<DiscreteFunctionType, AssembledOperatorType> sigmaLocalEstimate_;
    SigmaLocalFunctionType sigmaLocalFunction_;
    SigmaEstimateFunction sigmaEstimateFunction_;

  };




  template <class GridImp,
            class ProblemTraits,
            int polOrder>
  class EllipticAlgorithm
    : public SteadyStateAlgorithm< GridImp, ProblemTraits, polOrder >
  {
  public:
    typedef SteadyStateAlgorithm< GridImp, ProblemTraits, polOrder > BaseType;

    // type of Grid
    typedef typename BaseType::GridType                    GridType;

    // Choose a suitable GridView
    typedef typename BaseType::GridPartType                GridPartType;

    // initial data type
    typedef typename BaseType::ProblemType                 ProblemType;

    // An analytical version of our model
    typedef typename BaseType::ModelType                   ModelType;

    // The flux for the discretization of advection terms
    typedef typename BaseType::DiscreteTraits::FluxType    FluxType;

    // type of linear operator (i.e. matrix implementation)
    typedef typename BaseType::FullOperatorType            FullOperatorType;

    // The DG space operator
    typedef typename BaseType::AssemblyOperatorType        AssemblyOperatorType;
    typedef typename BaseType::AssembledOperatorType       AssembledOperatorType;

    // The discrete function for the unknown solution is defined in the DgOperator
    typedef typename BaseType::DiscreteFunctionType        DiscreteFunctionType;

    // ... as well as the Space type
    typedef typename BaseType::DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;

    // type of inverse operator (i.e. linear solver implementation)
    typedef typename BaseType::BasicLinearSolverType       BasicLinearSolverType;

    typedef typename BaseType::RestrictionProlongationType RestrictionProlongationType;

    typedef typename BaseType::EOCErrorIDs                 EOCErrorIDs;

    enum { dimension = GridType::dimension  };

    typedef typename DiscreteFunctionSpaceType ::
      template ToNewDimRange< dimension * ModelType::dimRange >::NewFunctionSpaceType SigmaFunctionSpaceType;

    typedef PoissonSigmaEstimator< GridPartType, DiscreteFunctionType, SigmaFunctionSpaceType, AssembledOperatorType, polOrder > PoissonSigmaEstimatorType;

    typedef typename  BaseType::AnalyticalTraits    AnalyticalTraits;

    using BaseType::problem;
    using BaseType::model;
    using BaseType::grid;
    using BaseType::gridWidth;
    using BaseType::gridSize;
    using BaseType::eocIds_;
    using BaseType::space;
    using BaseType::solution;
    using BaseType::rhs_;


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

  public:
    EllipticAlgorithm(GridType& grid, const std::string name = "" )
    : BaseType( grid, name ),
      grid_( grid ),
      gridPart_( grid_ ),
      dgOperator_( gridPart_, problem() ),
      dgAssembledOperator_( gridPart_, dgOperator_ ),
      invDgOperator_(),
      linDgOperator_(),
      space_( const_cast<DiscreteFunctionSpaceType &> (dgAssembledOperator_.space()) ),
      poissonSigmaEstimator_( gridPart_, solution(), dgAssembledOperator_, name ),
      polOrderContainer_(grid_ , 0),
      step_( 0 )
    {
      std::string gridName = Fem::gridName( grid_ );
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

      // we start with max order
      typedef typename PolOrderContainer::Iterator Iterator ;
      const Iterator end = polOrderContainer_.end();
      const int minimalOrder = solution().space().order(); // estimator_.minimalOrder() ;
      for( Iterator it = polOrderContainer_.begin(); it != end; ++it )
      {
        (*it).value() = minimalOrder ;
      }
    }

    virtual BasicLinearSolverType* createSolver()
    {
      linDgOperator_.reset( new FullOperatorType("dg operator", space_, space_ ) );

      if( space_.continuous() ) // Lagrange case
      {
        typedef Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType ;
        StencilType stencil( space_, space_ );
        linDgOperator_->reserve( stencil );
      }
      else // DG case
      {
        typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType ;
        StencilType stencil( space_, space_ );
        linDgOperator_->reserve( stencil );
      }

      linDgOperator_->clear();
      dgAssembledOperator_.assemble(0, *linDgOperator_, rhs_);
  		dgAssembledOperator_.testSymmetrie(*linDgOperator_);

      double absLimit   = Dune::Fem:: Parameter::getValue<double>("istl.absLimit",1.e-6);
      double reduction  = Dune::Fem:: Parameter::getValue<double>("istl.reduction",1.e-6);
      // this describes the factor for max iterations in terms of spaces size
      int maxIterFactor = Dune::Fem:: Parameter::getValue<double>("istl.maxiterfactor", int(-1) );

      // if max iterations is negative set it to twice the spaces size
      if( maxIterFactor < 0 )
        maxIterFactor = 2;

      return new BasicLinearSolverType(*linDgOperator_, reduction, absLimit );
    }

    //! default time loop implementation, overload for changes
    void solve ( const int loop )
    {
  #if 0
  #ifdef PADAPTSPACE
      int polOrder = Dune::Fem:: Parameter::getValue<double>("femdg.polynomialOrder",1);
      const int minimalOrder = estimator_.minimalmemory
      // only implemented for PAdaptiveSpace
      std::vector<int> polOrderVec( space_.gridPart().indexSet().size(0) );
      std::fill( polOrderVec.begin(), polOrderVec.end(), polOrder );

      polOrderContainer_.update();
      if ( estimator_.isPadaptive() )
      {
        typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
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
              order = minimalOrder;
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
      BaseType::solve( loop );

      // calculate new sigma
      poissonSigmaEstimator_.update();
    }


    //! finalize computation by calculating errors and EOCs
    void finalize( const int eocloop )
    {
      typedef typename ProblemType::ExactSolutionType  ExactSolutionType ;
      typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType > GridExactSolutionType;
      GridExactSolutionType ugrid( "exact solution", problem().exactSolution(), gridPart_, 1 );

      AnalyticalTraits::addEOCErrors( eocIds_, solution(), dgOperator_.model(), ugrid, *poissonSigmaEstimator_.sigma() );

      // delete solver and linear operator for next step
      invDgOperator_.reset();
      linDgOperator_.reset();
    }

    bool adaptation(const double tolerance)
    {
      // resize container
      polOrderContainer_.resize();

      const double error = 0.0;//estimator_.estimate( problem() );
      std::cout << "ESTIMATE: " << error << std::endl;
      typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
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

    const AssembledOperatorType &oper() const
    {
      return dgAssembledOperator_;
    }

    protected:

    GridType& grid_;
    GridPartType gridPart_;       // reference to grid part, i.e. the leaf grid

    AssemblyOperatorType           dgOperator_;
    AssembledOperatorType          dgAssembledOperator_;

    std::unique_ptr< BasicLinearSolverType > invDgOperator_;
    std::unique_ptr< FullOperatorType      > linDgOperator_;

    DiscreteFunctionSpaceType &space_;    // the discrete function space

    PoissonSigmaEstimatorType poissonSigmaEstimator_;

    PolOrderContainer polOrderContainer_;

    int  step_;

  };

}
}
#endif // FEMHOWTO_STEPPER_HH
