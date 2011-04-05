#ifndef DUNE_LIMITERPASS_HH
#define DUNE_LIMITERPASS_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/fem/misc/utility.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/dgmodelcaller.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

#include <dune/fem/io/parameter.hh>

//*************************************************************
namespace Dune {  
/*! @addtogroup PassLimit
*/

  // base class for Limiters, mainly reading limitEps 
  struct LimiterFunctionBase
  {
    const double limitEps_; 
    LimiterFunctionBase() 
      : limitEps_( getEpsilon() )
    {}
  protected:  
    //! get tolerance for shock detector 
    double getEpsilon() const 
    {
      // default value 
      double eps = 1e-8;
      eps = Parameter::getValue("LimitEps", eps );
      return eps;
    }

    void printInfo(const std::string& name ) const 
    {
      if( Parameter::verbose() ) 
      {
        std::cout << "LimiterFunction: " << name << " with limitEps = " << limitEps_ << std::endl;
      }
    }
  };

  template <class Field>
  struct MinModLimiter : public LimiterFunctionBase 
  {
    using LimiterFunctionBase :: limitEps_ ;
    typedef Field FieldType;

    MinModLimiter() : LimiterFunctionBase() 
    {
      printInfo( "minmod" );
    }

    FieldType operator ()(const FieldType& g, const FieldType& d ) const 
    {
      const FieldType absG = std::abs( g );
      const FieldType absD = std::abs( d );
      
      // avoid rounding errors 
      if ( (absG < limitEps_) && (absD < limitEps_) ) return 1;
      
      // if product smaller than zero set localFactor to zero 
      // otherwise d/g until 1 
      const FieldType localFactor = 
            ( (g * d) < 0.0 ) ? 0.0 : 
            ( (absG > absD) ) ? (d/g) : 1.0;

      return localFactor;
    }
  };

  template <class Field>
  struct SuperBeeLimiter : public LimiterFunctionBase 
  {
    using LimiterFunctionBase :: limitEps_ ;
    typedef Field FieldType;

    SuperBeeLimiter() : LimiterFunctionBase()
    {
      printInfo( "superbee" );
    }

    FieldType operator ()(const FieldType& g, const FieldType& d ) const 
    {
      const FieldType absG = std::abs( g );
      const FieldType absD = std::abs( d );

      // avoid rounding errors 
      if ( (absG < limitEps_) && (absD < limitEps_) ) return 1;
      if ( absG < limitEps_ ) return 1;
      
      const FieldType r = d / g ;

      // if product smaller than zero set localFactor to zero 
      // otherwise d/g until 1 
      const FieldType localFactor = 
            ( (g * d) < 0.0 ) ? 0.0 :
            ( std::max( std::min( 2.0 * r , 1.0 ), std::min( r, 2.0 ) ) );

      return localFactor;
    }
  };

  template <class Field>
  struct VanLeerLimiter : public LimiterFunctionBase 
  {
    using LimiterFunctionBase :: limitEps_ ;
    typedef Field FieldType;

    VanLeerLimiter() : LimiterFunctionBase()
    {
      printInfo( "vanLeer" );
    }

    FieldType operator ()(const FieldType& g, const FieldType& d ) const 
    {
      const FieldType absG = std::abs( g );
      const FieldType absD = std::abs( d );

      // avoid rounding errors 
      if ( (absG < limitEps_) && (absD < limitEps_) ) return 1;
      if ( absG < limitEps_ ) return 1;
      
      const FieldType r = d / g ;
      const FieldType absR = std::abs( r );

      return ( absR + r ) / ( 1.0 + absR );
    }
  };

  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class LimiterDiscreteModelCaller 
    : public DGDiscreteModelCaller< DiscreteModelImp, ArgumentImp, SelectorImp > 
  {
    typedef DGDiscreteModelCaller< DiscreteModelImp, ArgumentImp, SelectorImp> BaseType ;

  public:
    LimiterDiscreteModelCaller( DiscreteModelImp& problem)
      : BaseType( problem ) 
#ifndef NDEBUG 
      , quadId_(size_t(-1)) 
      , quadPoint_(-1)
#endif
    {}
  
    // check whether we have inflow or outflow direction 
    template< class Intersection, class QuadratureType >
    bool checkPhysical( const Intersection& intersection,
                        QuadratureType& faceQuad, 
                        const  int quadPoint) 
    {
#ifndef NDEBUG 
      // store quadature info
      quadId_    = faceQuad.id();
      quadPoint_ = quadPoint;
#endif
      // evaluate data 
      BaseType :: evaluateQuad( faceQuad, quadPoint,
                                data_->localFunctionsSelf(), valuesEn_ );

      // call problem checkDirection 
      return problem_.checkPhysical( *(intersection.inside()),
                                     faceQuad.point( quadPoint ), valuesEn_ ); 
    }

    // check whether we have inflow or outflow direction 
    template <class IntersectionIterator, class QuadratureType>
    bool checkDirection(const IntersectionIterator& nit,
                        QuadratureType& faceQuad, 
                        const  int quadPoint)
    {
      // check quadature info
      assert( quadId_    == faceQuad.id() );
      assert( quadPoint_ == quadPoint );

      // call problem checkDirection 
      return problem_.checkDirection(nit, time_, 
                                     faceQuad.localPoint(quadPoint),
                                     valuesEn_ );
    }
  protected:
#ifndef NDEBUG 
    size_t quadId_ ;
    int quadPoint_ ;
#endif

    using BaseType :: time_ ;
    using BaseType :: data_ ;
    using BaseType :: problem_ ;
    using BaseType :: valuesEn_ ;
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
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef GlobalTraitsImp Traits; 
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainFieldType DomainFieldType;
    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;
    typedef FieldVector<DomainFieldType, dimDomain - 1> FaceDomainType;

    typedef LimiterDefaultDiscreteModel <GlobalTraitsImp,Model,passId> DGDiscreteModelType;

    //typedef typename ExistsLimiterFunction< Model, DomainFieldType > :: LimiterFunctionType  LimiterFunctionType;
    typedef MinModLimiter< DomainFieldType > LimiterFunctionType;
    //typedef typename Model :: Traits :: LimiterFunctionType  LimiterFunctionType;
  };
  
  
  // **********************************************
  // **********************************************
  // **********************************************
  template <class GlobalTraitsImp, class Model, int passId >
  class LimiterDefaultDiscreteModel :
    public DGDiscreteModelDefaultWithInsideOutside< LimiterDefaultTraits<GlobalTraitsImp,Model,passId >, passId >  
  {
  public:
    typedef LimiterDefaultTraits<GlobalTraitsImp,Model,passId> Traits;
    
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::FaceDomainType FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType :: Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    typedef typename DomainType :: field_type DomainFieldType;

    typedef typename Traits :: LimiterFunctionType  LimiterFunctionType;

    enum { dimRange = RangeType :: dimension };
    
  public:
    /** \brief default limiter discrete model */
    LimiterDefaultDiscreteModel(const Model& mod, 
                                const DomainFieldType veloEps = 1e-8) 
      : model_(mod) , velocity_(0) 
      , veloEps_( veloEps ) 
    {
    }

    /** \brief copy constructor */
    LimiterDefaultDiscreteModel(const LimiterDefaultDiscreteModel& org) 
      : model_(org.model_) , velocity_(org.velocity_) 
      , veloEps_(org.veloEps_) 
    {
    }
    
    //! \brief returns false 
    bool hasSource() const { return false; }
    //! \brief returns true 
    bool hasFlux() const   { return true;  }
    
    /** \brief numericalFlux of for limiter evaluateing the difference
         of the solution in the current integration point if we are a an
         inflow intersection.
         This is needed for the shock detection.
    */
    template <class ArgumentTuple>
    double numericalFlux(const IntersectionType& it,
                         const double time, 
                         const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight) const
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
 
      if( checkDirection(it,time,x, uLeft) )
      {
        gLeft  = argURight;
        gLeft -= argULeft;
        gRight = gLeft;
        return it.integrationOuterNormal( x ).two_norm();
      }
      else 
      {
        gLeft = gRight = 0.0;
        return 0.0;
      }
    }

    /** \brief boundaryFlux to evaluate the difference of the interior solution 
        with the boundary value. This is needed for the limiter. 
        The default returns 0, meaning that we use the interior value as ghost value. 
    */
    template <class ArgumentTuple>
    double boundaryFlux(const IntersectionType& it,
                        const double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft, 
                        RangeType& gLeft) const
    { 
      gLeft = 0 ;
      return 0.0;
    }

    /** \brief returns true if model provides boundary values for this
        intersection */
    inline bool hasBoundaryValue(const IntersectionType& it,
                                 const double time, 
                                 const FaceDomainType& x) const
    {
      return true; 
    }

    //! returns difference between internal value and boundary 
    //! value 
    inline void boundaryValue(const IntersectionType& it,
                              const double time,
                              const FaceDomainType& x,
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
                        const DomainType& x,
                        const ArgumentTuple& u ) const 
    { 
      // take first component 
      return physical( en, x, Element<0> :: get( u ) );
    }

    /** \brief check physical values 
     *  \param[in] xglBary Point in the global coordinates
     */ 
    bool physical( const EntityType& en, 
                   const DomainType& x,
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
                        const double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft) const 
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);

      // evaluate velocity 
      model_.velocity(this->inside(),time,it.geometryInInside().global(x),
                      argULeft,velocity_);
      return checkDirection(it, x, velocity_);
    }
    
    // g = grad L ( w_E,i - w_E ) ,  d = u_E,i - u_E  
    // default limiter function is minmod 
    DomainFieldType limiterFunction( const DomainFieldType& g, const DomainFieldType& d ) const
    {
      return limiterFunction_( g, d );
    }
    
  protected:
    //! returns true, if we have an inflow boundary
    bool checkDirection(const IntersectionType& it,
                        const FaceDomainType& x,
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
      The implemented Shock detection is described in detail in: 
        L. Krivodonova and J. Xin and J.-F. Remacle and N. Chevaugeon and J. E. Flaherty
        Shock detection and limiting with discontinuous Galerkin methods for hyperbolic conservation laws.
        Appl. Numer. Math., 48(3-4), pages 323-338, 2004.
    
      Link to paper:
        http://www.scorec.rpi.edu/REPORTS/2003-3.pdf

      Limiting is done by simply setting the polynomial order to zero.
  */
  template <class DiscreteModelImp, class PreviousPassImp, int passId = 0 >
  class LimitDGPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp , passId > 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp, passId > BaseType;

    typedef LimitDGPass<DiscreteModelImp, PreviousPassImp, passId> ThisType;
    
    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;
    
    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef const EntityType ConstEntityType;

    typedef typename BaseType::ArgumentType ArgumentType;
    typedef typename PreviousPassType::GlobalArgumentType ArgumentFunctionType;
    typedef typename ArgumentFunctionType :: LocalFunctionType LocalFunctionType;
    
    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType;
    typedef typename LocalIdSetType :: IdType IdType;
    
    // Types extracted from the underlying grids
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType :: Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;
    
    // Various other types
    typedef typename DestinationType::LocalFunctionType DestLocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;

    typedef CombinedSelector< ThisType , SelectorType > CombinedSelectorType;
    typedef LimiterDiscreteModelCaller<
      DiscreteModelType, ArgumentType, CombinedSelectorType> DiscreteModelCallerType;

    // type of Communication Manager 
    typedef CommunicationManager< DiscreteFunctionSpaceType > CommunicationManagerType;

    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::DimRange,
     dimDomain = DiscreteFunctionSpaceType::DimDomain};
    enum { dimension = GridType :: dimension };
    typedef typename GridType :: ctype ctype; 
    typedef FieldVector<ctype, dimDomain-1> FaceDomainType;

    typedef PointBasedDofConversionUtility< dimRange > DofConversionUtilityType;

    //! is true if grid is structured grid 
    enum { StructuredGrid = ! Capabilities::IsUnstructured<GridType>::v };

    typedef FieldVector< DomainType , dimRange > DeoModType; 
    typedef FieldMatrix< DomainFieldType, dimension , dimension > MatrixType;

    typedef AllGeomTypes< typename GridPartType :: IndexSetType,
                          GridType> GeometryInformationType;

    typedef GeometryInformation< GridType, 1 > FaceGeometryInformationType;

    // get LagrangePointSet of pol order 1 
    typedef LagrangeDiscreteFunctionSpace < typename
      DiscreteFunctionSpaceType :: FunctionSpaceType, GridPartType, 1 >
          LagrangeSpaceType; 
    
    typedef DGFEntityKey<int> KeyType;
    typedef std::vector<int> CheckType;
    typedef std::pair< KeyType, CheckType > VectorCompType;
    typedef std::set< VectorCompType > ComboSetType;

    typedef std::pair < MatrixType , bool > MatrixCacheEntry; 
    typedef std::map< KeyType, MatrixCacheEntry > MatrixCacheType;

    //! type of local mass matrix 
    typedef LocalDGMassMatrix< DiscreteFunctionSpaceType,
                  VolumeQuadratureType > LocalMassMatrixType;

    //! type of used adaptation method 
    typedef AdaptationMethod<GridType> AdaptationMethodType;

    //! id for choosing admissible linear functions 
    enum AdmissibleFunctions { DGFunctions = 0, ReconstructedFunctions = 1 , BothFunctions = 2 };

  public:
    //- Public methods
    /** \brief constructor
     * 
     *  \param  problem    Actual problem definition (see problem.hh)
     *  \param  pass       Previous pass
     *  \param  spc        Space belonging to the discrete function local to
     *                     this pass
     */
    LimitDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc,
                const int vQ = -1,
                const int fQ = -1) :
      BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      currentTime_(0.0),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      localIdSet_( gridPart_.grid().localIdSet()),
      lagrangeSpace_(gridPart_),
      orderPower_( -((spc_.order()+1.0) * 0.25)),
      dofConversion_(dimRange),
      faceQuadOrd_( (fQ < 0) ? (2 * spc_.order() + 1) : fQ ),
      volumeQuadOrd_( (vQ < 0) ? (2 * spc_.order()) : vQ ),
      argOrder_( spc_.order() ),
      conformingComboSet_(),
      comboSet_(),
      tolFactor_( getTolFactor() ),
      tol_1_( 1.0/getTol() ),
      geoInfo_( gridPart_.indexSet() ),
      faceGeoInfo_( geoInfo_.geomTypes(1) ),
      matrixCacheVec_( gridPart_.grid().maxLevel() + 1 ),
      localMassMatrix_( spc_ , volumeQuadOrd_ ),
      adaptive_((AdaptationMethodType(gridPart_.grid())).adaptive()),
      cartesianGrid_( CheckCartesian::check( gridPart_ ) ),
      stepTime_(3, 0.0),
      calcIndicator_(true),
      adaptMode_(false),
      reconstruct_(false),
      applyLimiter_(true),
      admissibleFunctions_( getAdmissibleFunctions() )
    {
      if( gridPart_.grid().comm().rank() == 0 )
      {
        std::cout << "LimitPass: Grid is ";
        if( cartesianGrid_ ) 
          std::cout << "cartesian";
        else 
          std::cout << "unstructured";
        //std::cout << "! LimitEps: "<< limitEps_ << ", LimitTol: "<< 1./tol_1_ << std::endl;
        std::cout << "! LimitTol: "<< 1./tol_1_ << std::endl;
      }

      // we need the flux here 
      assert(problem.hasFlux());
    }
    
    //! Destructor
    virtual ~LimitDGPass() {
    }

    void disableFirstCall() const 
    {
      applyLimiter_ = false; 
    }
    
    void enableFirstCall() const 
    {
      applyLimiter_ = true; 
    }

  protected:    
    //! get tolerance factor for shock detector 
    double getTolFactor() const 
    {
      const double dim = dimension; 
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
    double getEpsilon() const 
    {
      // default value 
      double eps = 1e-8;
      eps = Parameter::getValue("femdg.limiter.limiteps", eps );
      return eps;
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
      static bool assign(const ArgImp& arg, DestImp& dest)
      {
        // reconstruct if this combination of orders has been given
        return (arg.space().order() == 0) && (dest.space().order() == 1);
      }
    };

    template <class S1> 
    struct AssignFunction<S1,S1>
    {
      template <class ArgImp, class DestImp> 
      static bool assign(const ArgImp& arg, DestImp& dest) 
      {
        dest.assign(arg);
        return false;
      }
    };
    
    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      // get stopwatch 
      Timer timer; 
      
      // get reference to U 
      const ArgumentFunctionType& U = (*(Element<0>::get(arg)));
      
      // initialize dest as copy of U 
      // if reconstruct_ false then only reconstruct in some cases 
      reconstruct_ =
          AssignFunction<typename ArgumentFunctionType ::
          DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>::
               assign( U , dest );

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

      // if polOrder of destination is > 0 then we have to do something 
      if( spc_.order() > 0 && applyLimiter_ )
      {
        // prepare, i.e. set argument and destination 
        prepare(arg, dest);

        elementCounter_ = 0;
        // dod limitation 
        IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it) 
        {
          Timer localTime; 
          applyLocalImp(*it);
          stepTime_[2] += localTime.elapsed();
          ++elementCounter_;
        }

        // finalize
        finalize(arg, dest);
      }

      // accumulate time 
      this->computeTime_ += timer.elapsed();
    }

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
      limitedElements_ = 0;
      notPhysicalElements_ = 0;
      
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;
      
      // initialize arg in caller 
      caller_.setArgument(*arg_);

      // time initialisation
      currentTime_ = this->time();

      // set time to caller 
      caller_.setTime(currentTime_);

      // calculate maximal indicator (if necessary)
      problem_.indicatorMax();

      const int numLevels = gridPart_.grid().maxLevel() + 1;
      // check size of matrix cache vec
      if( (int) matrixCacheVec_.size() < numLevels )
      {
        matrixCacheVec_.resize( numLevels );
      }
    }
    
    //! Some management.
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      //if( notPhysicalElements_ ) 
      /*
      if( limitedElements_ > 0 )
      {
        std::cout << " Time: " << currentTime_
                  << " Elements limited: " << limitedElements_
                  << " due to side effects: " << notPhysicalElements_
                  << std::endl;
      }
      */

      // communicate dest 
      dest.communicate();

      // finalize caller 
      caller_.finalize();
    }

  protected:
    //! apply local is virtual 
    void applyLocal(ConstEntityType& en) const
    {
      applyLocalImp(en);
    }

    //! Perform the limitation on all elements.
    void applyLocalImp(ConstEntityType& en) const
    {
      Timer indiTime; 

      // true if entity has boundary intersections 
      bool boundary = false;

      // true if entity has at least one non-conforming intersections 
      bool nonConforming = false ;

      // true if current element is a cartesian element 
      bool cartesian = cartesianGrid_;

      // extract types 
      enum { dim = EntityType :: dimension };
      typedef typename EntityType :: ctype coordType;
      
      // check argument is not zero
      assert( arg_ );
      
      //- statements
      // set entity to caller 
      caller_.setEntity( en );

      // get function to limit 
      const ArgumentFunctionType& U = (*(Element<0>::get(*arg_)));

      // get U on entity
      const LocalFunctionType uEn = U.localFunction(en);

      // get geometry 
      const Geometry& geo = en.geometry();
      
      // cache geometry type 
      const GeometryType geomType = geo.type();
      // get bary center of element 
      const DomainType enBary = geo.global( geoInfo_.localCenter( geomType ) );

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
      // when we want ro reconstruct in any case then 
      // limiter is true but indicator is calculated 
      // because of adaptation 
      bool limiter = reconstruct_; 

      // check physicality of data 
      // evaluate average returns true if not physical 
      if ( evalAverage(en, uEn, enVal ) ) 
      {
        limiter = true;
        // enable adaptation because of not physical 
        adaptIndicator = -1;
      }
      else if ( calcIndicator_ )
      {
        // check shock indicator 
        limiter = calculateIndicator(en, uEn, geo, limiter, limit, shockIndicator, adaptIndicator);
      }

      // check average value
      if (! problem_.physical(en, enBary, enVal) ) 
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
      assert( (problem_.hasPhysical()) ? (problem_.physical(en, enBary, enVal)) : true );

      stepTime_[0] += indiTime.elapsed();
      indiTime.reset();

      // evaluate function  
      if( limiter ) 
      {
        // get local references 
        std::vector< DomainType >& barys = barys_;
        std::vector< RangeType >&  nbVals = nbVals_;
        barys.reserve( dim * dim );
        nbVals.reserve( dim * dim );
        // set to size zero since values are determined new
        barys.resize( 0 );
        nbVals.resize( 0 );
        
        // loop over all neighbors 
        for ( ; niter != endnit; ++niter ) 
        {
          const IntersectionType& intersection = *niter; 

          // check cartesian
          if( ! StructuredGrid && cartesian ) 
          {
            // check whether this element is really cartesian 
            if ( ! CheckCartesian::checkIntersection(intersection) ) cartesian = false;
          }
          
          // check all neighbors 
          if (intersection.neighbor()) 
          {
            EntityPointerType ep = intersection.outside();
            EntityType& nb = *ep;
            
            // get U on entity
            const LocalFunctionType uNb = U.localFunction(nb);

            // non-conforming case 
            if( ! intersection.conforming() ) nonConforming = true ;
                
            RangeType nbVal;

            // get geometry of neighbor 
            const Geometry& nbGeo = nb.geometry();

            DomainType nbBary = nbGeo.global( geoInfo_.localCenter( nbGeo.type() ) );

            // evaluate function  
            if ( evalAverage(nb, uNb, nbVal ) )
            {
              limiter = true;
            }

            // calculate difference 
            nbVal -= enVal;

            // store difference of mean values 
            nbVals.push_back(nbVal);
            
            // calculate difference 
            nbBary -= enBary;

            // store value 
            barys.push_back(nbBary);

          } // end neighbor 

          ////////////////////////////
          // --boundary
          ////////////////////////////
          // use ghost cell approach for limiting 
          if( intersection.boundary() )
          {
            // we have entity with boundary intersections 
            boundary = true ; 
            
            typedef typename IntersectionType :: Geometry LocalGeometryType;
            const LocalGeometryType& interGeo = intersection.geometry();

            /////////////////////////////////////////
            // construct bary center of ghost cell 
            /////////////////////////////////////////
            const FaceDomainType& faceMid = faceGeoInfo_.localCenter( interGeo.type() );

            // get unit normal 
            DomainType lambda ( intersection.unitOuterNormal(faceMid) );

            // get one point of intersection 
            DomainType point ( interGeo.corner( 0 ) );
            point -= enBary; 

            const double length = (point * lambda);
            const double factorLength = (cartesianGrid_) ? 2.0 * length : length ;
            //const double factorLength = length ; //(cartesianGrid_) ? 2.0 * length : length ;
            lambda *= factorLength;

            assert( lambda.two_norm () > 0 );

            // store difference between bary centers 
            barys.push_back(lambda);

            /////////////////////////////////////////////////
            /////////////////////////////////////////////////
            
            RangeType nbVal( 0 );
            
            // check for boundary Value 
            const FaceDomainType localPoint 
                  ( interGeo.local( lambda + enBary ) );
            
            if( problem_.hasBoundaryValue(intersection, currentTime_, 
                                          localPoint) )
            {
              problem_.boundaryValue(intersection, currentTime_, 
                                     localPoint, 
                                     enVal, nbVal);
              // calculate difference 
              nbVal -= enVal;
            }

            // store value 
            nbVals.push_back( nbVal );
            /////////////////////////////////

          } //end boundary

        } // end intersection iterator 
      }

      // if limit, then limit all components 
      limit = limiter;
      {
        // check whether not physical occured 
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
        problem_.adaptation( gridPart_.grid() , en, shockIndicator, adaptIndicator );
      }

      // if nothing to limit then just return here
      if ( ! limiter ) return ;

      // increase number of limited elements 
      ++limitedElements_;
      
      // get local funnction for limited values
      DestLocalFunctionType limitEn = dest_->localFunction(en);

      // create combination set 
      const ComboSetType& comboSet = setupComboSet( nbVals_.size() , geomType , nonConforming );
      
      // initialize combo vecs 
      const size_t comboSize = comboSet.size();
      deoMods_.reserve( comboSize );
      comboVec_.reserve( comboSize );

      // reset values 
      deoMods_.resize( 0 );
      comboVec_.resize( 0 );

      if( admissibleFunctions_ >= ReconstructedFunctions ) 
      {
        // calculate linear functions 
        calculateLinearFunctions(en.level(), comboSet, geomType, 
                                 boundary, nonConforming , cartesian );
      }
      
      // add DG Function
      if( (admissibleFunctions_ % 2) == DGFunctions )
      {
        addDGFunction( en, geo, uEn, enVal, enBary );
      }
      
      // Limiting 
      limitFunctions(comboVec_,barys_,nbVals_,geomType,deoMods_);

      // take maximum of limited functions 
      getMaxFunction(deoMods_, deoMod_);

      // project deoMod_ to limitEn
      L2project(en, geo, enBary, enVal, limit, deoMod_, limitEn); 
      
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
      DeoModType D;
      FieldMatrix<double,dimDomain,dimDomain> A;
      RangeType b[dimDomain];
      TemporaryLocalFunction< DiscreteFunctionSpaceType > uTmp(spc_,en);
      uTmp.clear();
      for (int r=0;r<dimRange;++r) 
      {
        for (int i=0; i<dimDomain+1; ++i) 
        { 
          const int idx = dofConversion_.combinedDof(i,r);
          uTmp[ idx ] = uEn[ idx ];
        }
      }
      typedef typename LagrangeSpaceType :: LagrangePointSetType  LagrangePointSetType;
      // use LagrangePointSet to evaluate on cornners of the 
      // geometry and also use caching 
      const LagrangePointSetType& quad = lagrangeSpace_.lagrangePointSet( geo.type() );
      for(int i=0; i<dimDomain; ++i) 
      {
        uTmp.evaluate( quad[ i ], b[ i ]);   
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
    
    // get linear function from reconstruction of the average values
    void calculateLinearFunctions(const int level, 
                                  const ComboSetType& comboSet, 
                                  const GeometryType& geomType, 
                                  const bool hasBoundaryIntersection,
                                  const bool nonConforming,
                                  const bool cartesian ) const 
    {
      assert( (StructuredGrid) ? (cartesian) : true );  
      // use matrix cache in case of structured grid 
      const bool useCache = cartesian 
                            && ! nonConforming 
                            && ! hasBoundaryIntersection;

      assert( level < (int) matrixCacheVec_.size() );
      MatrixCacheType& matrixCache = matrixCacheVec_[level]; 

      const bool initMatrixCache = (matrixCache.size () == 0);
      //const bool useCache = false;
      const bool nonCachedMatrix = ! useCache || initMatrixCache;

      enum { dim = dimension };

      // get references 
      MatrixType matrix;
      MatrixCacheEntry localEntry;
      localEntry.first  = 0;
      localEntry.second = true;

      DomainType rhs;

      typedef typename ComboSetType :: iterator iterator; 

      // create matrix cache entries 
      if ( useCache && initMatrixCache )
      {
        const iterator endit = comboSet.end();
        for(iterator it = comboSet.begin(); it != endit; ++it) 
        {
          // get tuple of numbers 
          const KeyType& v = (*it).first; 
          MatrixCacheEntry& entry = matrixCache[ v ];
          entry.first  = 0;
          entry.second = true;
        }
      }
      
      // calculate linear functions 
      // D(x) = U_i + D_i * (x - w_i)
      const iterator endit = comboSet.end();
      for(iterator it = comboSet.begin(); it != endit; ++it) 
      {
        // get tuple of numbers 
        const KeyType& v = (*it).first; 

        assert( (useCache) ? 
            (matrixCache.find( v ) != matrixCache.end()) : 1 );
        
        // get matrix from cache if ok 
        MatrixCacheEntry& entry = (useCache) ? matrixCache[v] : localEntry;
        MatrixType& inverse = entry.first; 

        RangeFieldType det = (entry.second) ? 1 : 0;
       
        const bool calculateMatrix = nonCachedMatrix || ! entry.second ;

        if( calculateMatrix )
        {
          // setup matrix 
          for(int i=0; i<dim; ++i) 
          {
            matrix[i] = barys_[ v[i] ];
          }

          // invert matrix 
          det = FMatrixHelp :: invertMatrix(matrix,inverse);
        }

        // create new instance of limiter coefficients  
        DeoModType& dM = deoMod_;

        // if matrix is regular 
        if( std::abs( det ) > 0 )
        {
          // non-singular case 
          entry.second = true;
          
          // calculate D
          for(int r=0; r<dimRange; ++r)
          {
            for(int i=0; i<dim; ++i) 
            {
              rhs[i] = nbVals_[ v[i] ][r];
            }

            // get solution 
            inverse.mv( rhs, dM[r] );
          }

          // store linear function 
          deoMods_.push_back( dM );
          comboVec_.push_back( (*it).second );
        }
        else 
        {
          // we are in singular situation 
          entry.second = false;
          
          // apply least square by adding another point 
          // this should make the linear system solvable 
           
          // creare vector with size = dim+1
          // the first dim components are equal to v 
          std::vector<int> nV( dim+1 );
          for(int i=0; i<dim; ++i) nV[i] = v[i];

          // take first point of list of points to check 
          CheckType check ( (*it).second );
          assert( check.size() > 0 );

          // get check iterator 
          typedef typename CheckType :: iterator CheckIteratorType; 
          const CheckIteratorType checkEnd = check.end();

          // take first entry as start value 
          CheckIteratorType checkIt = check.begin();
          
          // matrix should be regular now 
          if( checkIt == checkEnd )
          {
            // should work for 1 or 2 otherwise error 
            DUNE_THROW(InvalidStateException,"Check vector has no entries!");
          }

          CheckType checkNums ( check );
          checkNums.reserve( checkNums.size() + dim );
          for(int i=0; i<dim; ++i) 
          {
            checkNums.push_back( v[i] );
          }

          for( ; checkIt != checkEnd; ++ checkIt )
          {
            // assign last element 
            nV[dim] = *checkIt ;
          
            // apply least square again 
            bool matrixSingular = 
                  applyLeastSquare( barys_,
                                    nbVals_,
                                    nV,
                                    dM,
                                    matrix,
                                    inverse,
                                    rhs,
                                    calculateMatrix );

            // if matrix was valid add to functions 
            if( ! matrixSingular ) 
            {
              // store linear function
              deoMods_.push_back( dM );

              // store vector with points to check 
              comboVec_.push_back( checkNums );
            }
          }
        }
      } // end solving 
    }

    // matrix = A^T * A 
    template <class NewMatrixType, class MatrixType> 
    void multiply_AT_A(const NewMatrixType& A, MatrixType& matrix) const
    {
      assert( (int) MatrixType :: rows == (int) NewMatrixType :: cols );
      typedef typename MatrixType :: field_type value_type;

      for(int row=0; row< MatrixType :: rows; ++row)
      {
        for(int col=0; col< MatrixType :: cols; ++col)
        {
          matrix[row][col] = 0;
          for(int k=0; k<NewMatrixType :: rows;  ++k)
          {
            matrix[row][col] += A[k][row] * A[k][col];
          }
        }
      }
    }

    // solve system by applying least square method 
    template <class BaryVectorType, class NbValVectorType,
              class FunctionType, class MatrixType, class VectorType> 
    bool applyLeastSquare(const BaryVectorType& barys,
                          const NbValVectorType& nbVals,
                          const std::vector<int> &nV,
                          FunctionType& dM, 
                          MatrixType& matrix, 
                          MatrixType& inverse, 
                          VectorType& rhs,
                          const bool calculateInverse ) const 
    {
      // dimension 
      enum { dim = dimension };
      // new dimension is dim + 1 
      enum { newDim = dim + 1 };

      // apply least square by adding another point 
      // this should make the linear system solvable 
       
      // need new matrix type containing one row more  
      typedef FieldMatrix<DomainFieldType, newDim , dim> NewMatrixType;
      typedef FieldVector<DomainFieldType, newDim > NewVectorType;

      // new matrix 
      NewMatrixType A ;

      assert( (int) nV.size() == newDim );

      // create matrix 
      for(int k=0; k<newDim; ++k) 
      {
        A[k] = barys[ nV[k] ];
      }

      RangeFieldType det = 1;

      // only calculate new if not from cache 
      if( calculateInverse ) 
      {
        // matrix = A^T * A 
        multiply_AT_A(A, matrix);

        // invert matrix 
        det = FMatrixHelp :: invertMatrix(matrix,inverse);
      }

      if( std::abs( det ) > 0 )
      {
        // need new right hand side 
        NewVectorType newRhs;
        
        // calculate D
        for(int r=0; r<dimRange; ++r)
        {
          // get right hand side 
          for(int i=0; i<newDim; ++i) 
          {
            newRhs[i] = nbVals[ nV[i] ][r];
          }
          
          // convert newRhs to matrix 
          rhs = 0;
          A.umtv(newRhs, rhs);
          
          // get solution 
          inverse.mv( rhs, dM[r] );
        }
        
        // return false because matrix is invertable 
        return false;
      }

      // matrix is still singular 
      return true;
    }

    //! limit all functions 
    void limitFunctions(const std::vector< CheckType >& comboVec,
                        const std::vector< DomainType>& barys,
                        const std::vector< RangeType >& nbVals,
                        const GeometryType& geomType,
                        std::vector< DeoModType >& deoMods) const
    {
      const size_t numFunctions = deoMods.size();
      // for all functions check with all values 
      for(size_t j=0; j<numFunctions; ++j) 
      {
        const std::vector<int> & v = comboVec[j];

        // loop over dimRange 
        for(int r=0; r<dimRange; ++r) 
        {
          RangeFieldType minimalFactor = 1;
          DomainType& D = deoMods[j][r];
          
          const size_t vSize = v.size();
          for(size_t m=0; m<vSize; ++m)
          {
            // get current number of entry 
            const size_t k = v[m];

            // evaluate values for limiter function 
            const DomainFieldType g = D * barys[k];
            const DomainFieldType d = nbVals[k][r];

            // call limiter function 
            const DomainFieldType localFactor = problem_.limiterFunction( g, d );

            // take minimum 
            minimalFactor = std::min( localFactor , minimalFactor );

            // this cannot be smaller, so break if we got this value 
            if( localFactor <= 0.0 ) break; 
          }

          // scale linear function 
          D *= minimalFactor;
        }
      }
    }

    // chose function with maximal gradient 
    void getMaxFunction(const std::vector< DeoModType >& deoMods,
                        DeoModType& deoMod) const
    {
      RangeType max (0);
      const size_t numFunctions = deoMods.size();
      //const int startFunc = numFunctions-1;
      const int startFunc = 0;
      std::vector< size_t > number(dimRange, startFunc);
      for(int r=0; r<dimRange; ++r)
      {
        max[r] = deoMods[ startFunc ][r].two_norm2();
      }

      for(size_t l=1; l<numFunctions; ++l)
      //for(int l= startFunc-1 ; l>=0; --l)
      {
        for(int r=0; r<dimRange; ++r)
        {
          RangeFieldType D_abs = deoMods[l][r].two_norm2();
          if( D_abs > max[r] )
          {
            number[r] = l;
            max[r] = D_abs;
          }
        }
      }

      for(int r=0; r<dimRange; ++r)
      {
        deoMod[r] = deoMods[ number[r] ][r];
      }
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
        if ( ! problem_.physical( en, quad.point(l), u ) ) 
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
      enum { dim = dimension };
      if( problem_.hasPhysical() )
      {  
#if 1
        // use LagrangePointSet to evaluate on corners of the 
        // geometry and also use caching 
        return checkPhysicalQuad( lagrangeSpace_.lagrangePointSet(geo.type()), uEn );
#else 
        {
          VolumeQuadratureType volQuad(en, volumeQuadOrd_ );
          if( ! checkPhysicalQuad(volQuad, uEn) ) return false;
        }

        const IntersectionIteratorType endnit = gridPart_.iend(en); 
        for (IntersectionIteratorType nit = gridPart_.ibegin(en); 
             nit != endnit; ++nit) 
        {
          if( intersection.neighbor() ) 
          {
            // check whether we have an inflow intersection or not 
            typedef TwistUtility<GridType> TwistUtilityType;
        
            // conforming case 
            if( TwistUtilityType::conforming(gridPart_.grid(),intersection) )
            {
              FaceQuadratureType faceQuadInner(gridPart_,intersection, faceQuadOrd_, FaceQuadratureType::INSIDE);
              if( ! checkPhysicalQuad( faceQuadInner, uEn ) ) return false;
            }
            else 
            { // non-conforming case 
              typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
              NonConformingQuadratureType faceQuadInner(gridPart_,intersection, faceQuadOrd_, FaceQuadratureType::INSIDE);
              if( ! checkPhysicalQuad( faceQuadInner, nitBary, uEn ) ) return false;
            }
          }
          else 
          {
            FaceQuadratureType faceQuadInner(gridPart_,intersection, faceQuadOrd_, FaceQuadratureType::INSIDE);
            if( ! checkPhysicalQuad( faceQuadInner, nitBary, uEn ) ) return false;
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
        return lf.numScalarDofs();
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
        return dimension + 1;
      }
    };

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
        return dimension + 1;
      }
    };

    // L2 projection 
    template <class LocalFunctionImp>
    void L2project(const EntityType& en,
       const Geometry& geo,
       const DomainType& enBary,
       const RangeType& enVal,
       const FieldVector<bool,dimRange>& limit,
       const DeoModType& deoMod,
       LocalFunctionImp& limitEn, 
       const bool constantValue = false ) const 
    {
      enum { dim = dimension };
      
      // true if geometry mapping is affine 
      const bool affineMapping = localMassMatrix_.affine();

      // set zero dof to zero
      limitEn.clear();

      RangeType retVal, phi;
      const RangeType& ret = (constantValue) ? enVal : retVal;
      
      // get quadrature for destination space order  
      VolumeQuadratureType quad( en, spc_.order() + 1 );
      
      typedef typename LocalFunctionImp :: BaseFunctionSetType BaseFunctionSetType;

      const BaseFunctionSetType& baseset = limitEn.baseFunctionSet();
      const int quadNop = quad.nop();
      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        // get quadrature weight 
        const double intel = (affineMapping) ?
          quad.weight(qP) : // affine case 
          quad.weight(qP) * geo.integrationElement( quad.point(qP) ); // general case
  
        // if we don't have only constant value then evaluate function 
        if( ! constantValue ) 
        {
          // get global coordinate 
          DomainType point = geo.global( quad.point( qP ) );
          point -= enBary;
    
          // evaluate linear function 
          for(int r=0; r<dimRange; ++r) {
            retVal[r] = enVal[r] + (deoMod[r] * point);
          }
        }

        // assume PointBased and CombinedSpace  
        for(int i=0; 
            i< NumLinearBasis<LocalFunctionImp, typename
                              LocalFunctionImp :: DiscreteFunctionSpaceType>::numBasis( limitEn ); 
            ++i) 
        {
          const int base = dofConversion_.combinedDof(i,0);
          baseset.evaluate(base, quad[qP] , phi);
          
          // project linear function 
          for(int r=0; r<dimRange; ++r) 
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            // here evaluateScalar could be used 
            limitEn[dofIdx] += intel * (ret[r] * phi[0]) ;
          }
        }
      }

      // apply local inverse mass matrix for non-linear mappings 
      if( ! affineMapping )
      {
        localMassMatrix_.applyInverse( en, limitEn );
      }
      
      // check physicality of projected data
      if ( (! constantValue) && (! checkPhysical(en, geo, limitEn)) )
      {
        // for affine mapping we only need to set higher moments to zero 
        if( affineMapping ) 
        {
          const int numBasis = limitEn.numScalarDofs(); 
          for(int i=1; i<numBasis; ++i) 
          {
            for(int r=0; r<dimRange; ++r) 
            {
              const int dofIdx = dofConversion_.combinedDof(i,r);
              limitEn[dofIdx] = 0;
            }
          }
        }
        else 
        {  
          // if not physical project to mean value 
          L2project(en,geo,enBary,enVal,limit,deoMod_, limitEn, true);
        }
      }
    }

    // evaluate average of local function lf on entity en 
    bool evalAverage(const EntityType& en, 
                     const LocalFunctionType& lf,
                     RangeType& val) const 
    {
      if( localMassMatrix_.affine() )
      {
        // get point quadrature 
        VolumeQuadratureType quad( en, 0 );

        RangeType phi; 
        
        lf.baseFunctionSet().evaluate(0, quad[0], phi);
        for(int r=0; r<dimRange; ++r) 
        {
          const int dofIdx = dofConversion_.combinedDof(0, r);
          // here evaluateScalar could be used 
          val[r] = lf[dofIdx] * phi[0];
        }
        // return whether value is physical 
        return ( problem_.hasPhysical() && ! problem_.physical(en, quad.point(0), val) );
      }
      else 
      {
        bool notphysical = false;
        const Geometry& geo = en.geometry();
        
        // get quadrature 
        VolumeQuadratureType quad( en, volumeQuadOrd_);

        RangeType tmp; 

        // set value to zero
        val = 0;

        const int quadNop = quad.nop();
        for(int qp=0; qp<quadNop; ++qp) 
        {
          // evaluate function  
          lf.evaluate( quad[qp] , tmp );

          // check whether value is physical 
          if ( problem_.hasPhysical() && ! problem_.physical(en, quad.point(qp), tmp) )
          {
            notphysical = true;
          }

          // apply integration weight 
          tmp *= quad.weight(qp) * geo.integrationElement(quad.point(qp));
          // sum up 
          val += tmp;
        }

        // mean value, i.e. devide by volume 
        val *= 1.0/geo.volume();
        return notphysical;
      }
    }

    // fill combination vector recursive 
    template <class SetType, int dim> 
    struct FillVector
    {
      static void fill(const int neighbors, 
                       const int start,
                       SetType& comboSet,
                       std::vector<int>& v) 
      {
        assert( (int) v.size() > dim );
        for(int n = start; n<neighbors; ++n)
        {
          v[dim] = n;
          FillVector<SetType, dim+1> :: fill( neighbors, 
                                              n + 1,
                                              comboSet,
                                              v);
        }
      } 
    };
    
    // termination of fill combination vector 
    template <class SetType> 
    struct FillVector<SetType, dimension-1 >
    {
      enum { dim = dimension-1 };
      static void fill(const int neighbors, 
                       const int start,
                       SetType& comboSet,
                       std::vector<int>& v) 
      {
        assert( (int) v.size() == dimension );
        for(int n = start; n<neighbors; ++n)
        {
          v[dim] = n;
          comboSet.insert( VectorCompType( v, std::vector<int> () ) );
        }
      } 
    };
    
    // setup set storing combinations of linear functions 
    const ComboSetType& setupComboSet(const int neighbors, const GeometryType& geomType, const bool nonConforming ) const 
    {
      // check for new build (this set is constant)
      if( conformingComboSet_.size() == 0 ) 
      {
        buildComboSet(neighbors, geomType, conformingComboSet_);
      }
      
      // in case of non-conforming grid or grid with more than one
      // element type this has to be re-build on every element
      if( nonConforming || 
          spc_.multipleGeometryTypes() )
      {
        buildComboSet(neighbors, geomType, comboSet_);
        return comboSet_;
      }
      else 
      {
        return conformingComboSet_;
      }
    }

    // build combo set 
    void buildComboSet(const int neighbors, 
                       const GeometryType& geomType,
                       ComboSetType& comboSet) const 
    {
      // clear set 
      comboSet.clear();

      // maximal number of neighbors 
      std::vector<int> v(dimension,0);
      FillVector<ComboSetType, 0> :: fill(neighbors,0, comboSet, v );

      // create set containing all numbers 
      std::set<int> constNumbers;
      typedef typename std::set<int> :: iterator el_iterator;
      for(int i = 0; i<neighbors; ++i)
      {
        constNumbers.insert(i);
      }
        
      const int checkSize = neighbors - dimension ;
      typedef typename ComboSetType :: iterator iterator; 
      const iterator endit = comboSet.end();
      for(iterator it = comboSet.begin(); it != endit; ++it) 
      {
        const KeyType& v = (*it).first;
        CheckType& check = const_cast<CheckType&> ((*it).second);

        // reserve memory 
        check.reserve ( checkSize );

        // get set containing all numbers 
        std::set<int> numbers (constNumbers);
        
        // remove the key values 
        for(int i=0; i<dimension; ++i) 
        {
          el_iterator el = numbers.find( v[i] );
          numbers.erase( el );
        }
        
        // generate check vector 
        el_iterator endel = numbers.end(); 
        for(el_iterator elit = numbers.begin(); 
            elit != endel; ++ elit)
        {
          check.push_back( *elit );
        }
      }
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

      // create intersection quadrature 
      IntersectionQuadratureType interQuad( gridPart_, intersection, faceQuadOrd_ );

      // get appropriate references 
      const QuadratureImp &faceQuadInner = interQuad.inside();
      const QuadratureImp &faceQuadOuter = interQuad.outside();
                         
      // set neighbor to caller 
      caller_.setNeighbor( nb, faceQuadInner, faceQuadOuter );

      typedef typename IntersectionType :: Geometry LocalGeometryType;
      const LocalGeometryType& interGeo = intersection.geometry();   

      RangeType jump, adapt;
      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l) 
      {
        // calculate jump 
        const double val = caller_.numericalFlux( intersection, 
                                                  faceQuadInner, faceQuadOuter, l,
                                                  jump , adapt );

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
    bool applyBoundary(const IntersectionType & intersection,
                      const QuadratureImp & faceQuadInner,
                      RangeType& shockIndicator,
                      RangeType& adaptIndicator) const 
    {
      typedef typename IntersectionType :: Geometry LocalGeometryType;
      const LocalGeometryType& interGeo = intersection.geometry();   

      RangeType jump;
      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l) 
      {
        // calculate jump 
        const double val = caller_.boundaryFlux(intersection, faceQuadInner,l,jump);

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

      const IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType niter = gridPart_.ibegin(en); 
           niter != endnit; ++niter) 
      {
        const IntersectionType& intersection = *niter ;

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

        // check whether we have an inflow intersection or not 
        typedef TwistUtility<GridType> TwistUtilityType;
        
        const int quadOrd = problem_.hasPhysical() ? faceQuadOrd_ : 0;
        
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
          NonConformingQuadratureType faceQuadInner(gridPart_,intersection, quadOrd, FaceQuadratureType::INSIDE);
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
          EntityPointerType ep = intersection.outside();
          EntityType& nb = *ep; 

          // set neighbor to caller 
          caller_.setNeighbor( nb );

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
                                                      FaceQuadratureType::OUTSIDE);
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
            EntityPointerType ep = intersection.outside();
            const EntityType& nb = *ep; 

            // conforming case 
            if( ! GridPartType::conforming && ! intersection.conforming() )
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

            // set neighbor entity to inside entity 
            caller_.setBoundary(en, faceQuadInner);

            if (applyBoundary(intersection, faceQuadInner,shockIndicator,adaptIndicator))
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
        const bool physical = caller_.checkPhysical( intersection, quad, l );

        if( checkPhysical && ! physical ) 
        {
          // return notPhysical 
          return true ; 
        }

        if( ! inflowIntersection ) 
        {
          // check intersection 
          if( physical && caller_.checkDirection(intersection, quad, l) ) 
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
    
  private:
    mutable DiscreteModelCallerType caller_;
    const DiscreteModelType& problem_; 
    mutable double currentTime_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    GridPartType& gridPart_;
    const LocalIdSetType& localIdSet_;
    LagrangeSpaceType lagrangeSpace_;

    const double orderPower_;
    const DofConversionUtilityType dofConversion_; 
    mutable int faceQuadOrd_;
    mutable int volumeQuadOrd_;
    mutable int argOrder_;

    mutable ComboSetType conformingComboSet_;
    mutable ComboSetType comboSet_;

    // tolerance to scale shock indicator 
    const double tolFactor_;
    const double tol_1_;

    // if true scheme is TVD
    const GeometryInformationType geoInfo_;
    const FaceGeometryInformationType faceGeoInfo_;

    mutable DeoModType deoMod_;

    mutable std::vector< DeoModType > deoMods_;
    mutable std::vector< CheckType >  comboVec_;

    mutable std::vector< DomainType > barys_;
    mutable std::vector< RangeType >  nbVals_;
    mutable std::vector< MatrixCacheType > matrixCacheVec_;

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
    
    //! true if limiter has been called after adaptation 
    mutable bool adaptMode_;
    
    //! true if limtier is used as finite volume scheme of higher order 
    mutable bool reconstruct_;

    //! true if limiting operator is applied 
    mutable bool applyLimiter_;

    // choice of admissible linear functions 
    const AdmissibleFunctions admissibleFunctions_;
    bool dgAdded_ ;
  }; // end DGLimitPass 

} // end namespace Dune 
#endif
