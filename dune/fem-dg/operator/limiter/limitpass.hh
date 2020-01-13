#ifndef DUNE_LIMITERPASS_HH
#define DUNE_LIMITERPASS_HH

#include <vector>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/quadrature/cornerpointset.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/operator/1order/localmassmatrix.hh>

#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>

#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange/lagrangepoints.hh>

#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/misc/compatibility.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

#include <dune/fem-dg/pass/pass.hh>
#include <dune/fem-dg/pass/context.hh>
#include <dune/fem-dg/pass/discretemodel.hh>
#include <dune/fem-dg/pass/modelcaller.hh>

#include <dune/fem-dg/operator/limiter/limiterutility.hh>

#if HAVE_DUNE_OPTIM
//#define WANT_DUNE_OPTIM 1
#define WANT_DUNE_OPTIM 0
#endif

#if WANT_DUNE_OPTIM
#include <dune/fv/lpreconstruction.hh>
#endif


//*************************************************************
namespace Dune
{
namespace Fem
{

  template <class DiscreteModel, class Argument, class PassIds >
  class LimiterDiscreteModelCaller
  : public CDGDiscreteModelCaller< DiscreteModel, Argument, PassIds >
  {
    typedef CDGDiscreteModelCaller< DiscreteModel, Argument, PassIds > BaseType;

  public:
    typedef typename BaseType::ArgumentType ArgumentType;
    typedef typename BaseType::DiscreteModelType DiscreteModelType;

    typedef typename BaseType::IntersectionType IntersectionType;

    typedef typename BaseType::RangeTupleType RangeTupleType;

    using BaseType::time;

    LimiterDiscreteModelCaller( ArgumentType &argument, DiscreteModelType &discreteModel )
    : BaseType( argument, discreteModel )
#ifndef NDEBUG
      , quadId_(size_t(-1))
      , quadPoint_(-1)
#endif
    {}

    // check whether we have inflow or outflow direction
    template< class QuadratureType >
    bool checkPhysical( const IntersectionType &intersection,
                        QuadratureType &quadrature,
                        const int qp )
    {
#ifndef NDEBUG
      // store quadature info
      quadId_    = quadrature.id();
      quadPoint_ = qp;
#endif
      // evaluate data
      localFunctionsInside_.evaluate( quadrature[ qp ], ranges_ );
      // BaseType :: evaluateQuad( quadrature, qp, localFunctionsInside_, values_ );

      // call problem checkDirection
      const typename BaseType::EntityType &entity = intersection.inside();
      return discreteModel().checkPhysical( entity, entity.geometry().local( intersection.geometry().global( quadrature.localPoint( qp ) ) ), ranges_ );
    }

    // check whether we have inflow or outflow direction
    template< class QuadratureType>
    bool checkDirection( const IntersectionType &intersection,
                         QuadratureType &quadrature,
                         const int qp )
    {
      // check quadature info
      assert( quadId_    == quadrature.id() );
      assert( quadPoint_ == qp );

      // call checkDirection() on discrete model
      return discreteModel().checkDirection( intersection, time(), quadrature.localPoint( qp ), ranges_ );
    }
  protected:
    using BaseType::discreteModel;
    using BaseType::localFunctionsInside_;

  private:
    Dune::TypeIndexedTuple< RangeTupleType, typename DiscreteModelType::Selector > ranges_;

#ifndef NDEBUG
    size_t quadId_ ;
    int quadPoint_ ;
#endif
  };

  template < class Model, class DomainFieldType, class dummy = double >
  struct ExistsLimiterFunction
  {
    typedef MinModLimiter< DomainFieldType > LimiterFunctionType ;
  };

  template < class Model >
  struct ExistsLimiterFunction< Model, typename Model::Traits::LimiterFunctionType >
  {
    typedef typename  Model::Traits::LimiterFunctionType  LimiterFunctionType;
  };

  template <class GlobalTraitsImp, class Model, int passId = -1 >
  class LimiterDefaultDiscreteModel;

  template <class GlobalTraitsImp, class Model, int passId >
  struct LimiterDefaultTraits
  {
    typedef Model ModelType;
    typedef typename ModelType::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };
    enum { dimGrid = GridType::dimension };

    typedef GlobalTraitsImp Traits;
    typedef typename ModelTraits::FunctionSpaceType FunctionSpaceType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainFieldType DomainFieldType;
    typedef typename DestinationType::DomainType DomainType;
    typedef FieldVector<DomainFieldType, dimGrid >  LocalDomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;
    typedef FieldVector<DomainFieldType, dimGrid - 1> FaceLocalDomainType;

    // Indicator function type for Limiter (for output mainly)
    typedef Fem::FunctionSpace< DomainFieldType, double, dimDomain, 3> FVFunctionSpaceType;
    typedef Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Fem::SimpleStorage> IndicatorSpaceType;
    typedef Fem::AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;

    typedef LimiterDefaultDiscreteModel <GlobalTraitsImp,Model,passId> DGDiscreteModelType;

    //typedef typename ExistsLimiterFunction< Model, DomainFieldType > :: LimiterFunctionType  LimiterFunctionType;
    //typedef MinModLimiter< DomainFieldType > LimiterFunctionType;
    typedef typename Model :: Traits :: LimiterFunctionType  LimiterFunctionType;
  };


  /**
   * \brief Default discrete model for limiter
   *
   * \ingroup PassBased
   * \ingroup DiscreteModel
   */
  template <class GlobalTraitsImp, class Model, int passId >
  class LimiterDefaultDiscreteModel :
    public Fem::DGDiscreteModelDefaultWithInsideOutside< LimiterDefaultTraits<GlobalTraitsImp,Model,passId >, passId >
  {
    std::integral_constant< int, passId > uVar;

  public:
    typedef LimiterDefaultTraits<GlobalTraitsImp,Model,passId> Traits;
    typedef Fem::DGDiscreteModelDefaultWithInsideOutside< Traits,passId > BaseType;

    typedef typename Traits::DomainType                     DomainType;
    typedef typename Traits::LocalDomainType                LocalDomainType;
    typedef typename Traits::FaceLocalDomainType            FaceLocalDomainType;
    typedef typename Traits::RangeType                      RangeType;
    typedef typename Traits::GridType                       GridType;
    typedef typename Traits::GridPartType                   GridPartType;
    typedef typename Traits::JacobianRangeType              JacobianRangeType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType         IntersectionType;
    typedef typename GridPartType::template Codim<0>::EntityType EntityType;
    typedef typename DomainType::field_type                 DomainFieldType;

    typedef typename Traits::LimiterFunctionType            LimiterFunctionType;

    enum { dimRange = RangeType::dimension };

  public:
    using BaseType::inside;
    using BaseType::outside;

    /** \brief default limiter discrete model */
    LimiterDefaultDiscreteModel(const Model& mod,
                                const DomainFieldType veloEps = 1e-8,
                                const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : model_(mod)
      , limiterFunction_(parameter)
      , velocity_(0)
      , veloEps_( veloEps )
    {
    }

    void setEntity(const EntityType& en)
    {
      BaseType :: setEntity ( en );
      model().setEntity( en );
    }
    //! \brief returns false
    bool hasSource() const { return false; }
    //! \brief returns true
    bool hasFlux() const   { return true;  }

    template <class LocalEvaluationVec >
    void initializeIntersection( const LocalEvaluationVec& left,
                                 const LocalEvaluationVec& right )
    {}

    template <class LocalEvaluationVec >
    void initializeBoundary(const LocalEvaluationVec& local )
    {}

    // old version
    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeIntersection(const IntersectionType& it,
                                const QuadratureImp& quadInner,
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec)
    {
    }

    // old version
    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeBoundary(const IntersectionType& it,
                            const QuadratureImp& quadInner,
                            const ArgumentTupleVector& uLeftVec)
    {
    }

    /** \brief numericalFlux of for limiter evaluateing the difference
         of the solution in the current integration point if we are a an
         inflow intersection.
         This is needed for the shock detection.
    */
    template <class FaceQuadratureImp,
              class ArgumentTuple,
              class JacobianTuple>
    double numericalFlux(const IntersectionType& it,
                         const double time,
                         const FaceQuadratureImp& innerQuad,
                         const FaceQuadratureImp& outerQuad,
                         const int quadPoint,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         const JacobianTuple& jacLeft,
                         const JacobianTuple& jacRight,
                         RangeType& shockIndicator,
                         RangeType& adaptIndicator,
                         JacobianRangeType&,
                         JacobianRangeType& ) const
    {
      const FaceLocalDomainType& x = innerQuad.localPoint( quadPoint );

      if( checkDirection(it,time, x, uLeft) )
      {
        shockIndicator  = uLeft[ uVar ];
        shockIndicator -= uRight[ uVar ];
        adaptIndicator = shockIndicator;
        return it.integrationOuterNormal( x ).two_norm();
      }
      else
      {
        adaptIndicator = shockIndicator = 0;
        return 0.0;
      }
    }

    /** \brief boundaryFlux to evaluate the difference of the interior solution
        with the boundary value. This is needed for the limiter.
        The default returns 0, meaning that we use the interior value as ghost value.
    */
    //! returns difference between internal value and boundary
    //! value
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& adaptIndicator,
                        JacobianRangeType& gDiffLeft ) const
    {
      adaptIndicator = 0 ;
      return 0.0;
    }

    /** \brief returns true if model provides boundary values for this
        intersection */
    template <class LocalEvaluation>
    inline bool hasRobinBoundaryValue(const LocalEvaluation& local) const
    {
      return false;
    }

    /** \brief returns true if model provides boundary values for this
        intersection */
    template <class LocalEvaluation>
    inline bool hasBoundaryValue(const LocalEvaluation& local) const
    {
      return true;
    }

    //! returns difference between internal value and boundary
    //! value
    template <class LocalEvaluation>
    inline void boundaryValue(const LocalEvaluation& local,
                              const RangeType& uLeft,
                              RangeType& uRight) const
    {
      uRight = uLeft;
    }

    //! \brief returns true if physical check does something useful */
    bool hasPhysical () const { return false; }

    /** \brief check physical values (called from LimiterDiscreteModelCaller) */
    template <class ArgumentTuple>
    bool checkPhysical( const EntityType& en,
                        const LocalDomainType& x,
                        const ArgumentTuple& u ) const
    {
      // take first component
      return physical( en, x, u[ uVar ] );
    }

    /** \brief check physical values
     *  \param[in] xglBary Point in the local coordinates
     */
    bool physical( const EntityType& en,
                   const LocalDomainType& x,
                   const RangeType& u) const
    {
      return true;
    }

    void indicatorMax() const
    {
    }

    /** \brief adaptation method */
    void adaptation(GridType& grid, const EntityType& en,
                    const RangeType& shockIndicator,
                    const RangeType& adaptindicator) const
    {
    }

    //! return reference to model
    const Model& model() const {
      return model_;
    }

    //! returns true, if we have an inflow boundary
    template <class ArgumentTuple>
    bool checkDirection(const IntersectionType& it,
                        const double time, const FaceLocalDomainType& x,
                        const ArgumentTuple& uLeft) const
    {
      // evaluate velocity
      model_.velocity(this->inside(), this->inside().geometry().local( it.geometry().global(x) ), uLeft[ uVar ], velocity_);
      return checkDirection(it, x, velocity_);
    }

    // g = grad L ( w_E,i - w_E ) ,  d = u_E,i - u_E
    // default limiter function is minmod
    const LimiterFunctionType& limiterFunction() const
    {
      return limiterFunction_;
    }

  protected:
    //! returns true, if we have an inflow boundary
    bool checkDirection(const IntersectionType& it,
                        const FaceLocalDomainType& x,
                        const DomainType& velocity) const
    {
      // calculate scalar product of normal and velocity
      const double scalarProduct = it.outerNormal(x) * velocity;

      // if scalar product is larger than veloEps it's considered to be
      // an outflow intersection, otherwise an inflow intersection
      return (scalarProduct < veloEps_);
    }

  protected:
    const Model& model_;
    const LimiterFunctionType limiterFunction_;
    mutable DomainType velocity_;
    const DomainFieldType veloEps_;
  };

  /** \brief Concrete implementation of Pass for Limiting.
   *
   *  \ingroup Pass
   *
   *  \note: A detailed description can be found in:
   *
   *   A. Dedner and R. Klöfkorn.
   *   \b A Generic Stabilization Approach for Higher Order
   *   Discontinuous Galerkin Methods for Convection Dominated Problems. \b
   *   J. Sci. Comput., 47(3):365-388, 2011. http://link.springer.com/article/10.1007%2Fs10915-010-9448-0
   */
  template <class DiscreteModelImp, class PreviousPassImp, int passId >
  class LimitDGPass
  : public LocalPass<DiscreteModelImp, PreviousPassImp , passId >
  {
    typedef LimitDGPass< DiscreteModelImp, PreviousPassImp, passId > ThisType;
    typedef LocalPass< DiscreteModelImp, PreviousPassImp, passId >   BaseType;

  public:
    //- Typedefs and enums

    //! Repetition of template arguments
    typedef DiscreteModelImp                                             DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp                                              PreviousPassType;

    typedef typename BaseType::PassIds                                   PassIds;

    // Types from the base class
    typedef typename BaseType::EntityType                                EntityType;

    typedef typename BaseType::ArgumentType                              ArgumentType;

  private:
   typedef typename DiscreteModelType::Selector                          Selector;
   typedef std::tuple_element_t< 0, Selector >                           ArgumentIdType;
   static const std::size_t argumentPosition
     = Dune::FirstTypeIndex< PassIds, ArgumentIdType >::type::value;
   typedef std::tuple_element_t< argumentPosition, ArgumentType >         ArgumentFunctionPtrType;

  public:
    typedef typename PreviousPassType::GlobalArgumentType                 ArgumentFunctionType;
    typedef typename ArgumentFunctionType::LocalFunctionType              LocalFunctionType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType           DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType      VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType        FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::IteratorType              IteratorType;
    typedef TemporaryLocalFunction< DiscreteFunctionSpaceType >           TemporaryLocalFunctionType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType                  GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType              GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType                DomainType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType           DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType            RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType         JacobianRangeType;

    typedef typename GridType::Traits::LocalIdSet                         LocalIdSetType;
    typedef typename LocalIdSetType::IdType                               IdType;

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType               IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType                       IntersectionType;
    typedef typename GridPartType::template Codim<0>::GeometryType        Geometry;
    typedef typename Geometry::LocalCoordinate                            LocalDomainType;

    // define limiter utility class
    typedef LimiterUtility< typename DiscreteFunctionSpaceType::FunctionSpaceType, GridType::dimension > LimiterUtilityType;

    // Various other types
    typedef typename DestinationType::LocalFunctionType                   DestLocalFunctionType;

    typedef LimiterDiscreteModelCaller< DiscreteModelType, ArgumentType, PassIds > DiscreteModelCallerType;

    // type of Communication Manager
    typedef CommunicationManager< DiscreteFunctionSpaceType >             CommunicationManagerType;

    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::dimRange,
           dimDomain = DiscreteFunctionSpaceType::dimDomain};
    enum { dimGrid = GridType :: dimension };
    typedef typename GridType::ctype                                       ctype;
    typedef FieldVector<ctype, dimGrid-1>                                  FaceLocalDomainType;

    typedef PointBasedDofConversionUtility< dimRange >                     DofConversionUtilityType;

    static const bool StructuredGrid     = GridPartCapabilities::isCartesian< GridPartType >::v;
    static const bool conformingGridPart = GridPartCapabilities::isConforming< GridPartType >::v;

    typedef typename LimiterUtilityType::GradientType  GradientType;
    typedef typename LimiterUtilityType::MatrixType    MatrixType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef AllGeomTypes< IndexSetType, GridType> GeometryInformationType;

    typedef GeometryInformation< GridType, 1 > FaceGeometryInformationType;

    // get LagrangePointSet of pol order 1
    typedef CornerPointSet< GridPartType >                                 CornerPointSetType;
    // get lagrange point set of order 1
    typedef std::map< Dune::GeometryType, CornerPointSetType >             CornerPointSetContainerType;

    typedef typename LimiterUtilityType::KeyType         KeyType;
    typedef typename LimiterUtilityType::CheckType       CheckType;
    typedef typename LimiterUtilityType::VectorCompType  VectorCompType;
    typedef typename LimiterUtilityType::ComboSetType    ComboSetType;

    typedef std::map< int, ComboSetType > ComboSetMapType ;

    typedef typename LimiterUtilityType::MatrixStorage MatrixCacheEntry;
    typedef std::map< KeyType, MatrixCacheEntry > MatrixCacheType;

    //! type of local mass matrix
    typedef LocalMassMatrix< DiscreteFunctionSpaceType,
                  VolumeQuadratureType > LocalMassMatrixType;

    //! type of used adaptation method
    typedef AdaptationMethod<GridType> AdaptationMethodType;

    //! id for choosing admissible linear functions
    enum AdmissibleFunctions { DGFunctions = 0, ReconstructedFunctions = 1 , BothFunctions = 2 };

    //! returns true if pass is currently active in the pass tree
    using BaseType :: active ;

    //! type of cartesian grid checker
    typedef CheckCartesian< GridPartType >  CheckCartesianType;

  protected:
    template <class DiscreteSpace>
    struct HierarchicalBasis
    {
      static const bool v = false ;
    };

    template < class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct HierarchicalBasis< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
    {
      static const bool v = true ;
    };

    template < class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct HierarchicalBasis< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
    {
      static const bool v = true ;
    };

#if WANT_DUNE_OPTIM
    typedef typename GridPartType :: GridViewType  GridViewType ;

    struct BoundaryValue
    {
      const ThisType& op_;
      BoundaryValue( const ThisType& op ) : op_( op ) {}

      RangeType operator () ( const typename GridViewType::Intersection &i,
                              const DomainType &x,
                              const DomainType &n,
                              const RangeType &uIn ) const
      {
        return uIn;
      }
    };

    typedef Dune::FV::LPReconstruction< GridViewType, RangeType, BoundaryValue > LinearProgramming;
#endif


  public:
    //- Public methods
    /** \brief constructor
     *
     *  \param  problem    Actual problem definition (see problem.hh)
     *  \param  pass       Previous pass
     *  \param  spc        Space belonging to the discrete function local to this pass
     *  \param  vQ         order of volume quadrature
     *  \param  fQ         order of face quadrature
     */
    LimitDGPass(DiscreteModelType& problem,
                PreviousPassType& pass,
                const DiscreteFunctionSpaceType& spc,
                const int vQ = -1,
                const int fQ = -1 ) :
      BaseType(pass, spc),
      caller_( 0 ),
      discreteModel_(problem),
      currentTime_(0.0),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
#if WANT_DUNE_OPTIM
      linProg_( static_cast< GridViewType > (gridPart_), BoundaryValue( *this ),
                Dune::Fem::Parameter::getValue<double>("finitevolume.linearprogramming.tol", 1e-8 )
              ),
#endif
      indexSet_( gridPart_.indexSet() ),
      localIdSet_( gridPart_.grid().localIdSet()),
      cornerPointSetContainer_(),
      uTmpLocal_( spc_ ),
      orderPower_( -((spc_.order()+1.0) * 0.25)),
      dofConversion_(dimRange),
      faceQuadOrd_( fQ ),
      volumeQuadOrd_( vQ ),
      argOrder_( spc_.order() ),
      storedComboSets_(),
      tolFactor_( getTolFactor() ),
      tol_1_( 1.0/getTol() ),
      geoInfo_( gridPart_.indexSet() ),
      faceGeoInfo_( geoInfo_.geomTypes(1) ),
      phi0_( 0 ),
      matrixCacheVec_( gridPart_.grid().maxLevel() + 1 ),
      factors_(),
      numbers_(),
      localMassMatrix_( spc_ , 2*spc_.order() ),
      adaptive_((AdaptationMethodType(gridPart_.grid())).adaptive()),
      cartesianGrid_( CheckCartesianType::check( gridPart_ ) ),
      stepTime_(3, 0.0),
      calcIndicator_(discreteModel_.calculateIndicator()),
      reconstruct_(false),
      admissibleFunctions_( getAdmissibleFunctions() ),
      usedAdmissibleFunctions_( admissibleFunctions_ ),
      counter_( 0 ),
      computeTime_( 0 )
    {
      if( Parameter :: verbose () )
      {
        std::cout << "LimitPass: Grid is ";
        if( cartesianGrid_ )
          std::cout << "cartesian";
        else
          std::cout << "unstructured";
        //std::cout << "! LimitEps: "<< limitEps_ << ", LimitTol: "<< 1./tol_1_ << std::endl;
        std::cout << "! Limiter.tolerance: "<< 1./tol_1_ << std::endl;
      }

      // we need the flux here
      assert(problem.hasFlux());
    }

    //! Destructor
    virtual ~LimitDGPass() {
      std::cout << "~LimitDGPass: op calls " << counter_ << " T_l = " << computeTime_ << std::endl;
    }

    //! return default face quadrature order
    static int defaultVolumeQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return (2 * space.order( entity ));
    }

    //! return default face quadrature order
    static int defaultFaceQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return (2 * space.order( entity )) + 1;
    }

  protected:
    //! return appropriate quadrature order, default is 2 * order(entity)
    int volumeQuadratureOrder( const EntityType& entity ) const
    {
      return ( volumeQuadOrd_ < 0 ) ? ( defaultVolumeQuadratureOrder( spc_, entity ) ) : volumeQuadOrd_ ;
    }

    //! return appropriate quadrature order, default is 2 * order( entity ) + 1
    int faceQuadratureOrder( const EntityType& entity ) const
    {
      return ( faceQuadOrd_ < 0 ) ? defaultFaceQuadratureOrder( spc_, entity ) : faceQuadOrd_ ;
    }


    //! get tolerance factor for shock detector
    double getTolFactor() const
    {
      const double dim = dimGrid;
      const double dimFactor = 1./dim;
      const double order = spc_.order();
      return dimFactor * 0.016 * std::pow(5.0, order);
    }

    //! get tolerance for shock detector
    double getTol() const
    {
      double tol = 1.0;
      tol = Parameter::getValue("femdg.limiter.tolerance", tol );
      return tol;
    }

    //! get tolerance for shock detector
    AdmissibleFunctions getAdmissibleFunctions() const
    {
      // default value
      int val = 1;
      val = Parameter::getValue("femdg.limiter.admissiblefunctions", val);
      assert( val >= DGFunctions || val <= BothFunctions );
      return (AdmissibleFunctions) val;
    }

    template <class S1, class S2>
    struct AssignFunction
    {
      template <class ArgImp, class DestImp>
      static bool assign(const ArgImp& arg, DestImp& dest, const bool firstThread)
      {
        // reconstruct if this combination of orders has been given
        return (arg.space().order() == 0) && (dest.space().order() == 1);
      }
    };

    template <class S1>
    struct AssignFunction<S1,S1>
    {
      template <class ArgImp, class DestImp>
      static bool assign(const ArgImp& arg, DestImp& dest, const bool firstThread )
      {
        if( firstThread )
        {
          dest.assign(arg);
        }
        return false;
      }
    };

  public:
    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      compute( arg, dest, std::numeric_limits<size_t>::max() );
    }

    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest, const size_t breakAfter) const
    {
      // get stopwatch
      Dune::Timer timer;

      // if polOrder of destination is > 0 then we have to do something
      if( spc_.order() > 0 && active() )
      {
        //std::cout << "LimitPass::compute is active" << std::endl;
        //std::cout << " is active";
        // prepare, i.e. set argument and destination
        prepare(arg, dest);

        elementCounter_ = 0;
        // do limitation
        const auto endit = spc_.end();
        for( auto it = spc_.begin(); (it != endit); ++it )
        {
          // for initialization of thread passes for only a few iterations
          //if( elementCounter_ > breakAfter) break;
          const auto& entity = *it;
          Dune::Timer localTime;
          applyLocalImp( entity );
          stepTime_[2] += localTime.elapsed();
          ++elementCounter_;
        }

        // finalize
        finalize(arg, dest);
        ++counter_;
        computeTime_ += timer.elapsed();
      }
      else
      {
        // get reference to U and pass on to dest
        const ArgumentFunctionType &U = *(std::get< argumentPosition >( arg ));
        dest.assign( U );
      }

      //std::cout << std::endl;

      // accumulate time
      this->computeTime_ += timer.elapsed();
    }

  protected:
    struct EvalAverage
    {
      const ThisType& op_;
      const GridPartType& gridPart_;
      const ArgumentFunctionType &U_;
      const DiscreteModelType& discreteModel_;
      const double volume_;

      typedef typename IntersectionType::Geometry IntersectionGeometry;

      EvalAverage( const ThisType& op, const ArgumentFunctionType& U, const DiscreteModelType& model, const double volume = -1.0 )
        : op_( op ), gridPart_( U.space().gridPart() ), U_( U ), discreteModel_( model ), volume_( volume )
      {}

      // return true is average value is non-physical
      bool evaluate( const EntityType& entity, RangeType& value ) const
      {
        // get U on entity
        const LocalFunctionType uEn = U_.localFunction(entity);
        return op_.evalAverage( entity, uEn, value );
      }

      bool boundaryValue( const EntityType& entity,
                          const IntersectionType& intersection,
                          const IntersectionGeometry& interGeo,
                          const DomainType& globalPoint,
                          const RangeType& entityValue,
                          RangeType& neighborValue ) const
      {
        FaceQuadratureType faceQuadInner( gridPart_, intersection, 0, FaceQuadratureType::INSIDE );
        typedef QuadratureContext< EntityType, IntersectionType, FaceQuadratureType > ContextType;
        typedef LocalEvaluation< ContextType, RangeType, RangeType > EvalType;

        ContextType cLeft( entity, intersection, faceQuadInner, volume_ );
        // create quadrature of low order
        EvalType local( cLeft, entityValue, entityValue, 0 );

        // check for boundary Value
        if( discreteModel_.hasBoundaryValue( local ) )
        {
          discreteModel_.boundaryValue( local, entityValue, neighborValue );
          return true ;
        }
        return false ;
      }
    };

  public:
    virtual std::vector<double> computeTimeSteps () const
    {
      //std::cout << stepTime_[1] << " step time limit \n";
      std::vector<double> tmp( stepTime_ );
      stepTime_[0] = stepTime_[1] = stepTime_[2] = 0.0;
      return tmp;
    }

    size_t leafElements() const
    {
      return elementCounter_;
    }

    //! In the preparations, store pointers to the actual arguments and
    //! destinations. Filter out the "right" arguments for this pass.
    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      prepare( arg, dest, true );
    }

    //! In the preparations, store pointers to the actual arguments and
    //! destinations. Filter out the "right" arguments for this pass.
    void prepare(const ArgumentType& arg, DestinationType& dest, const bool firstThread ) const
    {
      // get reference to U
      const ArgumentFunctionType &U = *(std::get< argumentPosition >( arg ));

      // initialize dest as copy of U
      // if reconstruct_ false then only reconstruct in some cases
      reconstruct_ =
          AssignFunction<typename ArgumentFunctionType ::
          DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>::
               assign( U , dest, firstThread );

      // if case of finite volume scheme set admissible functions to reconstructions
      usedAdmissibleFunctions_ = reconstruct_ ? ReconstructedFunctions : admissibleFunctions_;

      // in case of reconstruction
      if( reconstruct_ )
      {
        // in case of non-adaptive scheme indicator not needed
        calcIndicator_ = adaptive_;

        // adjust quadrature orders
        argOrder_ = U.space().order();
        faceQuadOrd_ = 2 * argOrder_ + 1;
        volumeQuadOrd_ = 2 * argOrder_;
      }

      limitedElements_ = 0;
      notPhysicalElements_ = 0;

      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      // time initialisation
      currentTime_ = this->time();

      // initialize caller
      caller_ = new DiscreteModelCallerType( *arg_, discreteModel_ );
      caller_->setTime(currentTime_);

      // calculate maximal indicator (if necessary)
      discreteModel_.indicatorMax();

      const size_t size = indexSet_.size( 0 ) ;
      // reset visited vector
      visited_.resize( size );
      std::fill( visited_.begin(), visited_.end(), false );

      factors_.resize( size );
      numbers_.clear();
      numbers_.resize( size );

      const int numLevels = gridPart_.grid().maxLevel() + 1;
      // check size of matrix cache vec
      if( (int) matrixCacheVec_.size() < numLevels )
      {
        matrixCacheVec_.resize( numLevels );
      }

#if WANT_DUNE_OPTIM
      values_.resize( size );
      gradients_.resize( size );

      // helper class for evaluation of average value of discrete function
      EvalAverage average( *this, U, discreteModel_);

      // do limitation
      const auto endit = spc_.end();
      for( auto it = spc_.begin(); (it != endit); ++it )
      {
        const auto& en = *it;
        average.evaluate( en, values_[ gridPart_.indexSet().index( en ) ] );
      }

      // get reconstructions
      linProg_( gridPart_.indexSet(), values_, gradients_ );
#endif
    }

    //! Some management (interface version)
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      finalize( arg, dest, true );
    }

    //! Some management (thread parallel version)
    void finalize(const ArgumentType& arg, DestinationType& dest, const bool doCommunicate) const
    {
      /*
      if( limitedElements_ > 0 )
      {
        std::cout << " Time: " << currentTime_
                  << " Elements limited: " << limitedElements_
                  << " due to side effects: " << notPhysicalElements_
                  << std::endl;
      }
      */

      if( doCommunicate )
      {
        // communicate dest
        dest.communicate();
      }

      // finalize caller
      if( caller_ )
      {
        delete caller_;
        caller_ = 0;
      }
    }

    //! apply local is virtual
    void applyLocal( const EntityType& entity ) const
    {
      applyLocalImp( entity );
    }

    //! apply local with neighbor checker (does nothing here)
    template <class NeighborChecker>
    void applyLocal( const EntityType& entity,
                     const NeighborChecker& ) const
    {
      // neighbor checking not needed in this case
      applyLocalImp( entity );
    }

    //! apply limiter only to elements without neighboring process boundary
    template <class NeighborChecker>
    void applyLocalInterior( const EntityType& entity,
                             const NeighborChecker& nbChecker ) const
    {
      if( nbChecker.isActive() )
      {
        // check whether on of the intersections is with ghost element
        // and if so, skip the computation of the limited solution for now
        for (const auto& intersection : intersections(gridPart_, entity) )
        {
          if( intersection.neighbor() )
          {
            // get neighbor
            const EntityType& nb = intersection.outside();

            // check whether we have to skip this intersection
            if( nbChecker.skipIntersection( nb ) )
            {
              return ;
            }
          }
        }
      }

      // otherwise apply limiting process
      applyLocalImp( entity );
    }

    //! apply limiter only to elements with neighboring process boundary
    template <class NeighborChecker>
    void applyLocalProcessBoundary( const EntityType& entity,
                                    const NeighborChecker& nbChecker ) const
    {
      assert( nbChecker.isActive() );
      assert( indexSet_.index( entity ) < int(visited_.size()) );
      // if entity was already visited, do nothing in this turn
      if( visited_[ indexSet_.index( entity ) ] ) return ;

      // apply limiter otherwise
      applyLocalImp( entity );
    }

  protected:
    //! Perform the limitation on all elements.
    void applyLocalImp(const EntityType& en) const
    {
      // timer for shock detection
      Dune::Timer indiTime;

      // extract types
      enum { dim = EntityType :: dimension };

      // check argument is not zero
      assert( arg_ );

      //- statements
      // set entity to caller
      caller().setEntity( en );

      // get function to limit
      const ArgumentFunctionType &U = *(std::get< argumentPosition >( *arg_ ));

      // get U on entity
      const LocalFunctionType uEn = U.localFunction(en);

      // get geometry
      const Geometry& geo = en.geometry();

      // cache geometry type
      const GeometryType geomType = geo.type();
      // get bary center of element
      const LocalDomainType& x = geoInfo_.localCenter( geomType );
      const DomainType enBary = geo.global( x );

      const IntersectionIteratorType endnit = gridPart_.iend( en );
      IntersectionIteratorType niter = gridPart_.ibegin(en);
      if( niter == endnit ) return ;

      // if a component is true, then this component has to be limited
      FieldVector<bool,dimRange> limit(false);
      // determ whether limitation is necessary
      // total jump vector
      RangeType shockIndicator(0);
      RangeType adaptIndicator(0);

      RangeType enVal;

      // if limiter is true then limitation is done
      // when we want to reconstruct in any case then
      // limiter is true but indicator is calculated
      // because of adaptation
      bool limiter = reconstruct_;

      // check physicality of data
      // evaluate average returns true if not physical
      if( evalAverage( en, uEn, enVal ) )
      {
        limiter = true;
        // enable adaptation because of not physical
        adaptIndicator = -1;
      }
      else if ( calcIndicator_ )
      {
        // check shock indicator
        limiter = calculateIndicator(en, uEn, enVal, geo, limiter, limit, shockIndicator, adaptIndicator);
      }
      else if( !reconstruct_ )
      {
        // check physical values for quadrature
        VolumeQuadratureType quad( en, spc_.order( en ) + 1 );
        if( ! checkPhysicalQuad( quad, uEn ) )
        {
          limiter = true;
          shockIndicator = 1.5;
        }
      }

      {
        // get barycenter of entity
        const DomainType& enBaryLocal = (int(dimGrid) == int(dimDomain)) ?
          geoInfo_.localCenter( geomType ) :
          geo.local( enBary ) ;

        // check average value
        if( discreteModel_.hasPhysical() && !discreteModel_.physical( en, enBaryLocal, enVal ) )
        {
          std::cerr << "Average Value "
                    << enVal
                    << " in point "
                    << enBary
                    << " is Unphysical!"
                    << std::endl << "ABORTED" << std::endl;
          assert( false );
          abort();
        }
      }

      stepTime_[0] += indiTime.elapsed();
      indiTime.reset();

      // boundary is true if boundary segment was found
      // nonconforming is true if entity has at least one non-conforming intersections
      // cartesian is true if the grid is cartesian and no nonconforming refinement present
      typename LimiterUtilityType::Flags flags( cartesianGrid_, limiter );

      // evaluate function
      if( limiter )
      {
        // helper class for evaluation of average value of discrete function
        EvalAverage average( *this, U, discreteModel_, geo.volume() );

        // setup neighbors barycenter and mean value for all neighbors
        LimiterUtilityType::setupNeighborValues( gridPart_, en, average, enBary, enVal,
                                                 StructuredGrid, flags, barys_, nbVals_ );
      }

      // if limit, then limit all components
      limit = limiter;
      {
        // check whether not physical occurred
        if (limiter && shockIndicator[0] < 1.)
        {
          shockIndicator = -1;
          ++notPhysicalElements_;
        }

        // if limiter also adapt
        if( limiter && ! reconstruct_ )
        {
          adaptIndicator = -1;
        }

        // call problem adaptation for setting refinement marker
        discreteModel_.adaptation( gridPart_.grid() , en, shockIndicator, adaptIndicator );
      }

      // mark entity as finished, even if not limited everything necessary was done
      assert( indexSet_.index( en ) < int(visited_.size()) );
      visited_[ indexSet_.index( en ) ] = true ;

      // if nothing to limit then just return here
      if ( ! limiter ) return ;

      // increase number of limited elements
      ++limitedElements_;

      const unsigned int enIndex = indexSet_.index( en );
#if WANT_DUNE_OPTIM
      bool useLinProg = true ;
      if( useLinProg )
      {
        deoMod_ = gradients_[ enIndex ];
      }
      else
#endif
      {
        // obtain combination set
        ComboSetType& comboSet = storedComboSets_[ nbVals_.size() ];
        if( comboSet.empty() )
        {
          // create combination set
          LimiterUtilityType::buildComboSet( nbVals_.size(), comboSet );
        }

        // reset values
        deoMods_.clear();
        comboVec_.clear();

        if( usedAdmissibleFunctions_ >= ReconstructedFunctions )
        {
          // level is only needed for Cartesian grids to access the matrix caches
          const int matrixCacheLevel = ( flags.cartesian ) ? en.level() : 0 ;
          assert( matrixCacheLevel < (int) matrixCacheVec_.size() );
          MatrixCacheType& matrixCache = matrixCacheVec_[ matrixCacheLevel ];

          // calculate linear functions, stored in deoMods_ and comboVec_
          LimiterUtilityType::calculateLinearFunctions( comboSet, geomType, flags,
                                    barys_, nbVals_,
                                    matrixCache,
                                    deoMods_,
                                    comboVec_ );
        }

        // add DG Function
        if( (usedAdmissibleFunctions_ % 2) == DGFunctions )
        {
          addDGFunction( en, geo, uEn, enVal, enBary );
        }

        // Limiting
        std::vector< RangeType > factors;
        LimiterUtilityType::limitFunctions(
            discreteModel_.limiterFunction(), comboVec_, barys_, nbVals_, deoMods_, factors );

        // take maximum of limited functions
        LimiterUtilityType::getMaxFunction(deoMods_, deoMod_, factors_[ enIndex ], numbers_[ enIndex ], factors );
      } // end if linProg

      // get local funnction for limited values
      DestLocalFunctionType limitEn = dest_->localFunction(en);

      // project deoMod_ to limitEn
      L2project(en, geo, enBary, enVal, limit, deoMod_, limitEn );

      assert( checkPhysical(en, geo, limitEn) );

      stepTime_[1] += indiTime.elapsed();

      //end limiting process
    }

  protected:
    // add linear components of the DG function
    void addDGFunction(const EntityType& en,
                       const Geometry& geo,
                       const LocalFunctionType& uEn,
                       const RangeType& enVal,
                       const DomainType& enBary) const
    {
      GradientType D;
      FieldMatrix<double,dimDomain,dimDomain> A;
      RangeType b[dimDomain];
      uTmpLocal_.init( en );
      uTmpLocal_.clear();

      // assume that basis functions are hierarchical
      assert( HierarchicalBasis< DiscreteFunctionSpaceType > :: v );

      for (int r=0;r<dimRange;++r)
      {
        for (int i=0; i<dimDomain+1; ++i)
        {
          const int idx = dofConversion_.combinedDof(i,r);
          uTmpLocal_[ idx ] = uEn[ idx ];
        }
      }

      // use LagrangePointSet to evaluate on cornners of the
      // geometry and also use caching
      const CornerPointSetType quad( en );
      for(int i=0; i<dimDomain; ++i)
      {
        uTmpLocal_.evaluate( quad[ i ], b[ i ]);
        b[i] -= enVal;
        A[i]  = geo.corner( i );
        A[i] -= enBary;
      }

      A.invert();
      DomainType rhs;

      // setup matrix
      for (int r=0; r<dimRange;++r)
      {
        for (int l=0; l<dimDomain; ++l)
        {
          rhs[ l ] = b[ l ][ r ];
        }
        A.mv( rhs, D[ r ]);
      }

      deoMods_.push_back( D );
      std::vector<int> comb( nbVals_.size() );
      const size_t combSize = comb.size();
      for (size_t i=0;i<combSize; ++i) comb[ i ] = i;
      comboVec_.push_back(comb);
    }

    // check physicality on given quadrature
    template <class QuadratureType, class LocalFunctionImp>
    bool checkPhysicalQuad(const QuadratureType& quad,
                           const LocalFunctionImp& uEn) const
    {
      // use LagrangePointSet to evaluate on corners of the
      // geometry and also use caching
      RangeType u;
      const int quadNop = quad.nop();
      const EntityType& en = uEn.entity();
      //const Geometry& geo = en.geometry();
      for(int l=0; l<quadNop; ++l)
      {
        uEn.evaluate( quad[l] , u );
        // warning!!! caching of geo can be more efficient
        //const DomainType& xgl = geo.global( quad[l] );
        // check data
        if ( ! discreteModel_.physical( en, quad.point(l), u ) )
        {
          // return notPhysical
          return false;
        }
      }
      // solution is physical
      return true;
    }

    //! check physicallity of data
    template <class LocalFunctionImp>
    bool checkPhysical(const EntityType& en,
                       const Geometry& geo,
                       const LocalFunctionImp& uEn) const
    {
      enum { dim = dimGrid };
      if( discreteModel_.hasPhysical() )
      {
#if 1
        // use LagrangePointSet to evaluate on corners of the
        // geometry and also use caching
        return checkPhysicalQuad( CornerPointSetType( en ), uEn );
#else
        {
          VolumeQuadratureType volQuad(en, volumeQuadratureOrder( en ) );
          if( ! checkPhysicalQuad(volQuad, uEn) ) return false;
        }

        const IntersectionIteratorType endnit = gridPart_.iend(en);
        for (IntersectionIteratorType nit = gridPart_.ibegin(en);
             nit != endnit; ++nit)
        {
          const IntersectionType& intersection = *nit;
          if( intersection.neighbor() && ! intersection.conforming() )
          {
            typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
            NonConformingQuadratureType faceQuadInner(gridPart_,intersection, faceQuadratureOrder( faceQuadratureOrder(en ), FaceQuadratureType::INSIDE);
            if( ! checkPhysicalQuad( faceQuadInner, uEn ) ) return false;
          }
          else
          {
            // conforming case
            FaceQuadratureType faceQuadInner(gridPart_,intersection, faceQuadratureOrder( faceQuadratureOrder(en ), FaceQuadratureType::INSIDE);
            if( ! checkPhysicalQuad( faceQuadInner, uEn ) ) return false;
          }
        }
#endif
      } // end physical
      return true;
    }

    template <class LocalFunctionImp, class SpaceImp>
    struct NumLinearBasis
    {
      inline static int numBasis(const LocalFunctionImp& lf)
      {
        return lf.numDofs()/dimRange;
      }
    };

    template <class LocalFunctionImp, class FunctionSpaceImp, class
      GridPartImp, int polOrd, template <class> class StrorageImp >
    struct NumLinearBasis<LocalFunctionImp,
              DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,
                                         StrorageImp> >
    {
      inline static int numBasis(const LocalFunctionImp& lf)
      {
        return dimGrid + 1;
      }
    };

    /*
    template <class LocalFunctionImp, class FunctionSpaceImp, class
      GridPartImp, int polOrd, template <class> class StrorageImp,
          int N , DofStoragePolicy policy >
    struct NumLinearBasis<LocalFunctionImp,
              CombinedSpace<
              DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,
                                         StrorageImp> , N , policy > >
    {
      inline static int numBasis(const LocalFunctionImp& lf)
      {
        return dimGrid + 1;
      }
    };
    */

    // L2 projection
    template <class LocalFunctionImp>
    void L2project(const EntityType& en,
       const Geometry& geo,
       const DomainType& enBary,
       const RangeType& enVal,
       const FieldVector<bool,dimRange>& limit,
       const GradientType& deoMod,
       LocalFunctionImp& limitEn,
       const bool constantValue = false ) const
    {
      enum { dim = dimGrid };

      // true if geometry mapping is affine
      const bool affineMapping = localMassMatrix_.affine();

      // set zero dof to zero
      uTmpLocal_.init( en );
      uTmpLocal_.clear();

      // get quadrature for destination space order
      VolumeQuadratureType quad( en, spc_.order() + 1 );

      const int quadNop = quad.nop();
      for(int qP = 0; qP < quadNop ; ++qP)
      {
        // get quadrature weight
        const double intel = (affineMapping) ?
          quad.weight(qP) : // affine case
          quad.weight(qP) * geo.integrationElement( quad.point(qP) ); // general case

        RangeType retVal( enVal );
        // if we don't have only constant value then evaluate function
        if( !constantValue )
        {
          // get global coordinate
          DomainType point = geo.global( quad.point( qP ) );
          point -= enBary;

          // evaluate linear function
          for( int r = 0; r < dimRange; ++r )
            retVal[ r ] += (deoMod[ r ] * point);
        }

        retVal *= intel;
        uTmpLocal_.axpy( quad[ qP ], retVal );
      }

      // apply local inverse mass matrix for non-linear mappings
      if( !affineMapping )
        localMassMatrix_.applyInverse( en, uTmpLocal_ );

      // check physicality of projected data
      if ( (! constantValue) && (! checkPhysical(en, geo, uTmpLocal_)) )
      {
        // for affine mapping we only need to set higher moments to zero
        if( affineMapping )
        {
          const int numBasis = uTmpLocal_.numDofs()/dimRange;
          for(int i=1; i<numBasis; ++i)
          {
            for(int r=0; r<dimRange; ++r)
            {
              const int dofIdx = dofConversion_.combinedDof(i,r);
              uTmpLocal_[dofIdx] = 0;
            }
          }
        }
        else
        {
          // if not physical project to mean value
          L2project(en,geo,enBary,enVal,limit,deoMod_, limitEn, true);
          return ;
        }
      }

      // in case of higher order FV the whole local functions needs to be assigned
      // since dest is not previously initialized with the current solution
      if( reconstruct_ )
      {
        limitEn.assign( uTmpLocal_ );
      }
      else
      {
        // copy to limitEn skipping components that should not be limited
        const int numBasis = uTmpLocal_.numDofs()/dimRange;
        for(int i=1; i<numBasis; ++i)
        {
          for( const auto& r : discreteModel_.model().limitedRange() )
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[ dofIdx ] = uTmpLocal_[ dofIdx ];
          }
        }
      }
    }

    template <class BasisFunctionSetType, class PointType>
    const RangeType& evaluateConstantBasis( const BasisFunctionSetType& basisSet,
                                           const PointType& x ) const
    {
      // calculate constant part of the basis functions
      if( ! (phi0_[ 0 ] > 0 ) )
      {
        std::vector< RangeType > phi( basisSet.size() );
        basisSet.evaluateAll( x, phi );
        phi0_ = phi[ 0 ];
      }

#ifndef NDEBUG
      // check that phi0 is valid
      {
        std::vector< RangeType > phi( basisSet.size() );
        basisSet.evaluateAll( x, phi );
        assert( (phi0_ - phi[ 0 ]).infinity_norm() < 1e-8 );
      }
#endif

      // return constant part of basis functions
      return phi0_ ;
    }

    // evaluate average of local function lf on entity en
    bool evalAverage(const EntityType& en,
                     const LocalFunctionType& lf,
                     RangeType& val) const
    {
      bool notphysical = false;
      if( HierarchicalBasis< DiscreteFunctionSpaceType > :: v
          && localMassMatrix_.affine() )
      {
        // get point quadrature
        VolumeQuadratureType quad( en, 0 );

        const RangeType& phi0 = evaluateConstantBasis( lf.basisFunctionSet(), quad[ 0 ] );
        for(int r=0; r<dimRange; ++r)
        {
          const int dofIdx = dofConversion_.combinedDof(0, r);
          // here evaluateScalar could be used
          val[r] = lf[dofIdx] * phi0 [ 0 ];
        }

        // possibly adjust average value, e.g. calculate primitive vairables and so on
        discreteModel_.adjustAverageValue( en, quad.point( 0 ), val );

        // return whether value is physical
        notphysical = (discreteModel_.hasPhysical() && !discreteModel_.physical( en, quad.point( 0 ), val ) );
      }
      else
      {
        const Geometry& geo = en.geometry();

        // get quadrature
        VolumeQuadratureType quad( en, volumeQuadratureOrder( en ) );

        // set value to zero
        val = 0;

        const int quadNop = quad.nop();
        /*
        if( int(aver_.size()) < quadNop )
        {
          // resize value vector
          aver_.resize( quadNop );
        }
        */

        // evaluate quadrature at once (does not qork correctly yet)
        // lf.evaluateQuadrature( quad, aver_ );
        RangeType aver;

        for(int qp=0; qp<quadNop; ++qp)
        {
          lf.evaluate( quad[ qp ], aver );

          // check whether value is physical
          notphysical |= (discreteModel_.hasPhysical() && !discreteModel_.physical( en, quad.point( qp ), aver ) );

          // possibly adjust average value, e.g. calculate primitive vairables and so on
          discreteModel_.adjustAverageValue( en, quad.point( qp ), aver );

          // apply integration weight
          aver *= quad.weight(qp) * geo.integrationElement( quad.point(qp) );
          // sum up
          val += aver;
        }

        // mean value, i.e. devide by volume
        val *= 1.0/geo.volume();
      }

      return notphysical;
    }

    template <bool conforming>
    bool applyLocalNeighbor(const IntersectionType & intersection,
                            const EntityType& nb,
                            RangeType& shockIndicator,
                            RangeType& adaptIndicator) const
    {
      // make sure we got the right conforming statement
      assert( intersection.conforming() == conforming );

      // use IntersectionQuadrature to create appropriate face quadratures
      typedef IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
      typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;

      // create intersection quadrature (no neighbor check here)
      IntersectionQuadratureType interQuad( gridPart_, intersection, faceQuadratureOrder( nb ), true );

      // get appropriate references
      const QuadratureImp &faceQuadInner = interQuad.inside();
      const QuadratureImp &faceQuadOuter = interQuad.outside();

      // set neighbor and initialize intersection
      caller().initializeIntersection( nb, intersection, faceQuadInner, faceQuadOuter );

      typedef typename IntersectionType :: Geometry LocalGeometryType;
      const LocalGeometryType& interGeo = intersection.geometry();

      JacobianRangeType dummy ;
      RangeType jump, adapt;
      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l)
      {
        // calculate jump
        const double val = caller().numericalFlux( intersection,
                                                   faceQuadInner, faceQuadOuter, l,
                                                   jump , adapt, dummy, dummy );

        // non-physical solution
        if(val < 0.0)
        {
          return true;
        }

        // get integration factor
        const double intel =
             interGeo.integrationElement(faceQuadInner.localPoint(l))
           * faceQuadInner.weight(l);

        // shock indicator
        jump *= intel;
        shockIndicator += jump;
        // adapt indicator
        adapt *= intel;
        adaptIndicator += adapt;
      }
      return false;
    }

    template <class QuadratureImp>
    bool applyBoundary(const EntityType& entity,
                       const IntersectionType & intersection,
                       const QuadratureImp & faceQuadInner,
                       const RangeType& entityValue,
                       RangeType& shockIndicator,
                       RangeType& adaptIndicator) const
    {
      typedef typename IntersectionType :: Geometry LocalGeometryType;
      const LocalGeometryType& interGeo = intersection.geometry();

      RangeType jump;
      JacobianRangeType dummy ;

      typedef QuadratureContext< EntityType, IntersectionType, QuadratureImp > ContextType;
      typedef LocalEvaluation< ContextType, RangeType, RangeType > EvalType;

      ContextType cLeft( entity, intersection, faceQuadInner, 0.0 );
      // create quadrature of low order
      EvalType local( cLeft, entityValue, entityValue );

      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l)
      {
        // calculate jump
        const double val = caller().boundaryFlux( intersection, faceQuadInner, l, jump, dummy);

        // non-physical solution
        if (val < 0.0)
        {
          return true;
        }

        const double intel =
              interGeo.integrationElement(faceQuadInner.localPoint(l))
            * faceQuadInner.weight(l);

        jump *= intel;
        shockIndicator += jump;
        adaptIndicator += jump;
      }
      return false;
    }

    // calculate shock detector
    bool calculateIndicator(const EntityType& en,
                            const LocalFunctionType& uEn,
                            const RangeType& enVal,
                            const Geometry& geo,
                            const bool initLimiter,
                            FieldVector<bool,dimRange>& limit,
                            RangeType& shockIndicator,
                            RangeType& adaptIndicator) const
    {
      enum { dim = EntityType :: dimension };

      // calculate circume during neighbor check
      double circume = 0.0;

      //   for min faceVol
      double faceVol = 1e10;
      int numberInSelf = -1;
      double currVol = -1e10;
      //double refFaceVol = -1e10;

      bool limiter = initLimiter;
      limit = false;

      shockIndicator = 0;
      adaptIndicator = 0;

      for (const auto& intersection : intersections(gridPart_, en) )
      {
        typedef typename IntersectionType :: Geometry LocalGeometryType;
        const LocalGeometryType& interGeo = intersection.geometry();

        // calculate max face volume
        if( numberInSelf != ( int ) intersection.indexInInside() )
        {
          if (numberInSelf >= 0) {
            //    min faceVol
            //faceVol = std::min( faceVol, currVol * refFaceVol );
            // omit ref vol in newest version
            faceVol = std::min( faceVol, currVol );
          }
          //refFaceVol = faceGeoInfo_.referenceVolume( interGeo.type() );
          numberInSelf = intersection.indexInInside();
          currVol = 0.0;
        }

        // add face volume to sum of local volumes
        const double vol = interGeo.volume();
        currVol += vol;

        const int quadOrd = discreteModel_.hasPhysical() ? faceQuadratureOrder( en ) : 0;

        // flag to trigger inflow intersection
        bool inflowIntersection = false;
        // conforming case
        if( intersection.conforming() )
        {
          FaceQuadratureType faceQuadInner(gridPart_,intersection, quadOrd, FaceQuadratureType::INSIDE);
          if( checkIntersection( intersection, faceQuadInner, inflowIntersection ) )
          {
            shockIndicator = -1;
            return true;
          }

        }
        else
        { // non-conforming case
          typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
          NonConformingQuadratureType faceQuadInner(gridPart_,intersection, quadOrd, NonConformingQuadratureType::INSIDE);
          if( checkIntersection( intersection, faceQuadInner, inflowIntersection ) )
          {
            shockIndicator = -1;
            return true;
          }
        }

        // check all neighbors
        // if we have an outflow intersection check other side too
        if (intersection.neighbor() && ! inflowIntersection )
        {
          // get neighbor entity
          const EntityType& nb = intersection.outside();

          // set neighbor to caller
          caller().setNeighbor( nb );

          if( intersection.conforming() )
          {
            FaceQuadratureType faceQuadOuter(gridPart_,intersection, quadOrd, FaceQuadratureType::OUTSIDE);
            if( checkIntersection( intersection, faceQuadOuter, inflowIntersection , false ) )
            {
              shockIndicator = -1;
              return true;
            }
          }
          else
          { // non-conforming case
            typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
            NonConformingQuadratureType faceQuadOuter(gridPart_,intersection, quadOrd,
                                                      NonConformingQuadratureType::OUTSIDE);
            if( checkIntersection( intersection, faceQuadOuter, inflowIntersection , false ) )
            {
              shockIndicator = -1;
              return true;
            }
          }

          // invert direction
          inflowIntersection = ! inflowIntersection;
        }

        // calculate indicator on inflow intersections
        if( inflowIntersection )
        {
          // add face vol to circume
          circume += vol;

          // order of quadrature
          const int jumpQuadOrd = spc_.order();

          // check all neighbors
          if (intersection.neighbor())
          {
            // get neighbor entity
            const EntityType& nb = intersection.outside();

            // conforming case
            if( ! conformingGridPart && ! intersection.conforming() )
            { // non-conforming case
              if (applyLocalNeighbor< false > (intersection, nb, shockIndicator, adaptIndicator))
              {
                shockIndicator = -1;
                return true;
              }
            }
            else
            {
              if (applyLocalNeighbor< true > (intersection, nb, shockIndicator, adaptIndicator))
              {
                shockIndicator = -1;
                return true;
              }
            }
          }

          // check all neighbors
          if ( intersection.boundary() )
          {
            FaceQuadratureType faceQuadInner(gridPart_,intersection, jumpQuadOrd, FaceQuadratureType::INSIDE);

            // initialize intersection
            caller().initializeBoundary( intersection, faceQuadInner );

            if (applyBoundary(en, intersection, faceQuadInner, enVal, shockIndicator, adaptIndicator))
            {
              shockIndicator = -1;
              return true;
            }
          }
        }

      } // end intersection iterator

      // calculate max face volume
      {
        //    min faceVol
        // faceVol = std::min( faceVol, currVol * refFaceVol );
        // omit ref vol in newest version
        faceVol = std::min( faceVol, currVol );
      }

      // volume factor is scaled with volume of reference element
      // the formula is refVol Elem / refVol face (for cubes this is 1)
      // for simplices this is 1/dim
      //const double elemVol = geo.volume() * geoInfo_.referenceVolume( geo.type() );

      // do not consider ref vol in newest version
      const double elemVol = geo.volume();

      // calculation of 1 / (h^orderPower)
      const double hVal = elemVol / faceVol;

      const double hPowPolOrder = std::pow( hVal , orderPower_ );

      // multiply h pol ord with circume
      const double circFactor = (circume > 0.0) ? (hPowPolOrder/(circume * tolFactor_ )) : 0.0;

      for(int r=0; r<dimRange; ++r)
      {
        // only scale shock indicator with tolerance
        shockIndicator[r] = std::abs(shockIndicator[r]) * circFactor * tol_1_;
        adaptIndicator[r] = std::abs(adaptIndicator[r]) * circFactor;
        if(shockIndicator[r] > 1.)
        {
          limit[r] = true;
          limiter = true;
        }
      }

      return limiter ;
    }

    template <class QuadratureType>
    bool checkIntersection(const IntersectionType& intersection,
                           const QuadratureType& quad,
                           bool& inflowIntersection,
                           const bool checkPhysical = true ) const
    {
      // check whether we have an inflow intersection or not
      const int quadNop = quad.nop();
      for(int l=0; l<quadNop; ++l)
      {
        // check physicality of value
        const bool physical = caller().checkPhysical( intersection, quad, l );

        if( checkPhysical && ! physical )
        {
          // return notPhysical
          return true ;
        }

        if( ! inflowIntersection )
        {
          // check intersection
          if( physical && caller().checkDirection(intersection, quad, l) )
          {
            inflowIntersection = true;
            // in case of physicality check is disabled
            // just break
            if ( ! checkPhysical ) break;
          }
        }
      }
      // return physical
      return false;
    }

    // make private
    LimitDGPass();
    LimitDGPass(const LimitDGPass&);
    LimitDGPass& operator=(const LimitDGPass&);

    const CornerPointSetType& cornerPointSet( const GeometryType& geomType ) const
    {
      return cornerPointSetContainer_[ geomType ];
    }

  protected:
    DiscreteModelCallerType &caller () const
    {
      assert( caller_ );
      return *caller_;
    }

  private:
    mutable DiscreteModelCallerType *caller_;
    DiscreteModelType& discreteModel_;
    mutable double currentTime_;

    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    GridPartType& gridPart_;

#if WANT_DUNE_OPTIM
    mutable LinearProgramming linProg_;
#endif

    const IndexSetType& indexSet_;
    const LocalIdSetType& localIdSet_;

    CornerPointSetContainerType cornerPointSetContainer_;

    mutable TemporaryLocalFunctionType uTmpLocal_;

    const double orderPower_;
    const DofConversionUtilityType dofConversion_;
    mutable int faceQuadOrd_;
    mutable int volumeQuadOrd_;
    mutable int argOrder_;

    mutable ComboSetMapType storedComboSets_;

    // tolerance to scale shock indicator
    const double tolFactor_;
    const double tol_1_;

    // if true scheme is TVD
    const GeometryInformationType geoInfo_;
    const FaceGeometryInformationType faceGeoInfo_;

    mutable GradientType deoMod_;
    mutable RangeType    phi0_ ;

    mutable std::vector< GradientType > deoMods_;
    mutable std::vector< CheckType >    comboVec_;

    mutable std::vector< RangeType >  aver_ ;
    mutable std::vector< DomainType > barys_;
    mutable std::vector< RangeType >  nbVals_;
    mutable std::vector< MatrixCacheType > matrixCacheVec_;

    mutable std::vector< RangeType  > factors_;
    mutable std::vector< std::vector< int > > numbers_;

    // vector for stroing the information which elements have been computed already
    mutable std::vector< bool > visited_;

    LocalMassMatrixType localMassMatrix_;
    //! true if limiter is used in adaptive scheme
    const bool adaptive_;
    //! true if grid is cartesian like
    const bool cartesianGrid_;
    mutable int limitedElements_, notPhysicalElements_;
    mutable std::vector<double> stepTime_;
    mutable size_t elementCounter_;

    //! true if indicator should be calculated
    mutable bool calcIndicator_;

    //! true if limiter is used as finite volume scheme of higher order
    mutable bool reconstruct_;

    // choice of admissible linear functions
    const AdmissibleFunctions admissibleFunctions_;
    mutable AdmissibleFunctions usedAdmissibleFunctions_ ;

    mutable std::vector< RangeType  > values_;
    mutable std::vector< GradientType > gradients_;

    mutable int counter_;
    mutable double computeTime_;

  }; // end DGLimitPass

} // namespace
} // namespace Dune

#endif // #ifndef DUNE_LIMITERPASS_HH
