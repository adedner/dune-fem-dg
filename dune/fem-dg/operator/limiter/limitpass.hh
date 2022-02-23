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

#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/basesetlocalkeystorage.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/interpolate.hh>

#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/common/capabilities.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/localfunction/bindable.hh>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/fem-dg/pass/pass.hh>
#include <dune/fem-dg/pass/context.hh>
#include <dune/fem-dg/pass/discretemodel.hh>
#include <dune/fem-dg/pass/modelcaller.hh>

#include <dune/fem-dg/operator/limiter/limiterutility.hh>
#include <dune/fem-dg/operator/limiter/limiterdiscretemodel.hh>
#include <dune/fem-dg/operator/limiter/lpreconstruction.hh>

#include <dune/fem-dg/operator/limiter/smoothness.hh>
#include <dune/fem-dg/operator/limiter/indicatorbase.hh>



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

      // call problem checkDirection
      const typename BaseType::EntityType &entity = intersection.inside();
      //return discreteModel().checkPhysical( entity, entity.geometry().local( intersection.geometry().global( quadrature.localPoint( qp ) ) ), ranges_ );
      return discreteModel().checkPhysical( entity, quadrature.point( qp ), ranges_ );
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
      return discreteModel().checkDirection( intersection, quadrature.localPoint( qp ), quadrature.point( qp ), ranges_ );
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
    typedef ConstLocalFunction< ArgumentFunctionType >                    LocalFunctionType;
    //typedef typename ArgumentFunctionType::LocalFunctionType              LocalFunctionType;

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
    typedef FieldVector< typename GridType::ctype, dimGrid-1>              FaceLocalDomainType;

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
    //! _default refers to muscl(1) on simplices and LP limiter (3) on cubes and polygons
    enum AdmissibleFunctions { dgonly = 0, muscl = 1 , muscldg = 2, lp = 3, _default = 4 };

    //! returns true if pass is currently active in the pass tree
    using BaseType :: active ;

    //! type of Cartesian grid checker
    typedef CheckCartesian< GridPartType >  CheckCartesianType;

    //! function describing an external troubled cell indicator
    typedef std::shared_ptr< TroubledCellIndicatorBase<ArgumentFunctionType> > TroubledCellIndicatorType;

  protected:
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

    struct ConstantFunction :
      public Dune::Fem::BindableGridFunction< GridPartType, Dune::Dim<dimRange> >
    {
      typedef Dune::Fem::BindableGridFunction<GridPartType, Dune::Dim<dimRange> > Base;
      typedef typename Base::RangeType RangeType;
      using Base::Base;

      RangeType value_;

      template <class Quadrature, class RangeArray>
      void evaluateQuadrature(const Quadrature& quadrature, RangeArray &values) const
      {
        const int nop = quadrature.nop();
        values.resize( nop );
        std::fill( values.begin(), values.end(), value_ );
      }

      template <class Point>
      void evaluate(const Point &x, typename Base::RangeType &ret) const
      {
        ret = value_;
      }

      bool isConstant() const { return true; }

      unsigned int order() const { return 0; }
      std::string name() const { return "LimmitPass::ConstantFunction"; } // needed for output
    };

    struct LinearFunction : public ConstantFunction
    {
      DomainType   center_;
      GradientType gradient_;

      typedef ConstantFunction Base;
      typedef typename Base::RangeType RangeType;
      using Base::Base;
      using Base::value_;

      template <class Quadrature, class RangeArray>
      void evaluateQuadrature(const Quadrature& quadrature, RangeArray &values) const
      {
        const int nop = quadrature.nop();
        values.resize( nop );
        for( int qp = 0; qp<nop; ++qp )
        {
          evaluate( quadrature[ qp ], values[ qp ] );
        }
      }

      template <class Point>
      void evaluate(const Point &x, RangeType &ret) const
      {
        Base::evaluate( x, ret );

        // get global coordinate
        auto point = Base::global( x );
        point -= center_;

        // evaluate linear function
        for( int r = 0; r < dimRange; ++r )
        {
          ret[ r ] += (gradient_[ r ] * point);
        }
      }

      bool isConstant() const { return gradient_.infinity_norm() < 1e-10; }

      unsigned int order() const { return 1; }
      std::string name() const { return "LimmitPass::LinearFunction"; } // needed for output
    };


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
                const int fQ = -1,
                const bool verbose = Dune::Fem::Parameter::verbose() )
      : LimitDGPass(problem, pass, spc, Dune::Fem::Parameter::container(), vQ, fQ, verbose )
    {}

    //- Public methods
    /** \brief constructor
     *
     *  \param  problem    Actual problem definition (see problem.hh)
     *  \param  pass       Previous pass
     *  \param  spc        Space belonging to the discrete function local to this pass
     *  \param  parameter  Parameter reader for dynamic parameters
     *  \param  vQ         order of volume quadrature
     *  \param  fQ         order of face quadrature
     */

    LimitDGPass(DiscreteModelType& problem,
                PreviousPassType& pass,
                const DiscreteFunctionSpaceType& spc,
                const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container(),
                const int vQ = -1,
                const int fQ = -1,
                const bool verbose = Dune::Fem::Parameter::verbose() ) :
      BaseType(pass, spc),
      caller_(),
      discreteModel_(problem),
      currentTime_(0.0),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      linearFunction_( gridPart_ ),
      linProg_(),
      indexSet_( gridPart_.indexSet() ),
      localIdSet_( gridPart_.grid().localIdSet()),
      uTmpLocal_( spc_ ),
      uEn_(),
      uEnAvg_(),
      limitEn_(),
      interpolEn_( spc_ ),
      orderPower_( -((spc_.order()+1.0) * 0.25)),
      dofConversion_(dimRange),
      faceQuadOrd_( fQ ),
      volumeQuadOrd_( vQ ),
      argOrder_( spc_.order() ),
      storedComboSets_(),
      tolFactor_( getTolFactor() ),
      indicator_( getIndicator(parameter) ),
      tol_1_( 2./ parameter.getValue("femdg.limiter.tolerance", double(1.0) ) ),
      geoInfo_( gridPart_.indexSet() ),
      faceGeoInfo_( geoInfo_.geomTypes(1) ),
      phi0_( 0 ),
      matrixCacheVec_( gridPart_.grid().maxLevel() + 1 ),
      factors_(),
      numbers_(),
      adaptive_((AdaptationMethodType(gridPart_.grid())).adaptive()),
      cartesianGrid_( CheckCartesianType::check( gridPart_ ) ),
      stepTime_(3, 0.0),
      calcIndicator_(discreteModel_.calculateIndicator()),
      reconstruct_(false),
      admissibleFunctions_( getAdmissibleFunctions( parameter ) ),
      usedAdmissibleFunctions_( usedAdmissibleFunctions() ),
      extTroubledCellIndicator_( indicator_ == 3
           ? new ModalSmoothnessIndicator< ArgumentFunctionType >() : nullptr ),
      counter_( 0 ),
      overallLimited_( 0 ),
      verbose_( verbose )
    {
      if( usedAdmissibleFunctions_ == lp )
      {
        linProg_.reset( new LinearProgramming( gridPart_, BoundaryValue( *this ),
                        parameter.getValue<double>("finitevolume.linearprogramming.tol", 1e-8 )) );
      }

      if( verbose_ )
      {
        std::cout << "LimitDGPass (Grid is ";
        if( cartesianGrid_ )
          std::cout << "cartesian) ";
        else
          std::cout << "unstructured) ";

        std::cout << "parameters: " << std::endl;
        std::cout << "  femdg.limiter.tolerance: " << 1./tol_1_ << std::endl;
        std::cout << "  femdg.limiter.indicator: " << indicator_ << std::endl;
        std::string adFctName( admissibleFctNames()[ admissibleFunctions_ ] );
        if( admissibleFunctions_ != usedAdmissibleFunctions_ )
          adFctName += "("+admissibleFctNames()[ usedAdmissibleFunctions_ ] + ")";
        std::cout << "  femdg.limiter.admissiblefunctions: " << adFctName << std::endl;
      }

      // we need the flux here
      assert(problem.hasFlux());
    }

    //! Destructor
    virtual ~LimitDGPass()
    {
      /*
      if( verbose_ )
      {
        std::cout << "~LimitDGPass: op calls " << counter_ << " limEl: " << overallLimited_ << std::endl;
        std::cout << "~LimitDGPass: lp: " << stepTime_[2] << "  indi: " << stepTime_[0] << "  proj: " << stepTime_[1] << std::endl;
      }
      */
    }

    void setTroubledCellIndicator(TroubledCellIndicatorType indicator)
    {
      std::abort();
      extTroubledCellIndicator_ = indicator;
    }

    //! return default face quadrature order
    static int defaultVolumeQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >::volumeOrder( space.order( entity ) );
    }

    //! return default face quadrature order
    static int defaultFaceQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >::surfaceOrder( space.order( entity ) );
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
    double getIndicator( const Dune::Fem::ParameterReader &parameter ) const
    {
      static const std::string indicators[]  = { "user", "none", "jump" ,"modal", "always" };
      std::string key( "femdg.limiter.indicator" );
      return parameter.getEnum( key, indicators, 2);
    }

    static const std::string (&admissibleFctNames())[5]
    {
      static const std::string admissiblefct[]  = { "dgonly" , "muscl", "muscl+dg", "lp", "default" };
      return admissiblefct;
    }

    //! get tolerance for shock detector
    AdmissibleFunctions getAdmissibleFunctions( const Dune::Fem::ParameterReader &parameter ) const
    {
      std::string key( "femdg.limiter.admissiblefunctions" );
      return (AdmissibleFunctions) parameter.getEnum( key, admissibleFctNames(), _default );
    }

    AdmissibleFunctions usedAdmissibleFunctions() const
    {
      if( admissibleFunctions_ == _default )
      {
        if( ! geoInfo_.multipleGeomTypes() && geoInfo_.geomTypes(0)[0].isSimplex() )
        {
          // default for simplices is muscl reconstruction because it
          // seems to work better for simplex grids
          return muscl;
        }
        else
        {
          // default for all other element types is lp reconstruction
          return lp;
        }
      }
      return admissibleFunctions_;
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

        // do limitation
        const auto endit = spc_.end();
        for( auto it = spc_.begin(); (it != endit); ++it )
        {
          // for initialization of thread passes for only a few iterations
          const auto& entity = *it;
          //Dune::Timer localTime;
          applyLocalImp( entity );
          //stepTime_[2] += localTime.elapsed();
        }

        // finalize
        finalize(arg, dest);
        this->computeTime_ += timer.elapsed();
      }
      else
      {
        // get reference to U and pass on to dest
        const ArgumentFunctionType &U = *(std::get< argumentPosition >( arg ));
        dest.assign( U );
      }

      //std::cout << std::endl;

      // accumulate time
      //this->computeTime_ += timer.elapsed();
    }

  protected:
    struct EvalAverage
    {
      const ThisType& op_;
      const GridPartType& gridPart_;
      LocalFunctionType &uLocal_;
      const DiscreteModelType& discreteModel_;
      const double volume_;

      typedef typename IntersectionType::Geometry IntersectionGeometry;

      EvalAverage( const ThisType& op, LocalFunctionType& uLocal, const DiscreteModelType& model, const double volume = -1.0 )
        : op_( op ), gridPart_( uLocal.discreteFunction().space().gridPart() ),
          uLocal_( uLocal ), discreteModel_( model ), volume_( volume )
      {}

      // return true is average value is non-physical
      bool evaluate( const EntityType& entity, RangeType& value ) const
      {
        // get U on entity
        auto guard = bindGuard( uLocal_, entity );
        return op_.evalAverage( entity, uLocal_, value );
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
      //stepTime_[0] = stepTime_[1] = stepTime_[2] = 0.0;
      return tmp;
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
      if( firstThread ) ++counter_;

      // get reference to U
      const ArgumentFunctionType &U = *(std::get< argumentPosition >( arg ));

      // initialize dest as copy of U
      // if reconstruct_ false then only reconstruct in some cases
      reconstruct_ =
          AssignFunction<typename ArgumentFunctionType ::
          DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>::
               assign( U , dest, firstThread );

      // if case of finite volume scheme set admissible functions to reconstructions
      if( reconstruct_ && ! linProg_ )
      {
        // if linProg is not available then the only other FV reconstruction is muscl
        usedAdmissibleFunctions_ = muscl;
      }

      // in case of reconstruction
      if( reconstruct_ )
      {
        // for FV scheme in case of non-adaptive scheme indicator not needed
        calcIndicator_ = adaptive_;

        // adjust quadrature orders
        argOrder_ = U.space().order();
        faceQuadOrd_ = 2 * argOrder_ + 1;
        volumeQuadOrd_ = 2 * argOrder_;
      }

      limitedElements_ = 0;
      notPhysicalElements_ = 0;
      this->numberOfElements_ = 0;

      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      // time initialisation
      currentTime_ = this->time();

      // initialize caller
      caller_.reset( new DiscreteModelCallerType( *arg_, discreteModel_ ) );
      caller_->setTime(currentTime_);

      uEn_.reset( new LocalFunctionType( U ) );
      uEnAvg_.reset( new LocalFunctionType( U ) );
      limitEn_.reset( new DestLocalFunctionType( *dest_ ) );

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

      if( linProg_ )
      {
        values_.resize( size );
        valuesComputed_.resize( size );
        std::fill( valuesComputed_.begin(), valuesComputed_.end(), false );
      }
    }

    //! Some management (interface version)
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      finalize( arg, dest, true );
    }

    //! Some management (thread parallel version)
    void finalize(const ArgumentType& arg, DestinationType& dest, const bool doCommunicate) const
    {
      overallLimited_ += limitedElements_;
      /*
      if( limitedElements_ > 0 && counter_ >0 )
      {
        std::cout << " Time: " << currentTime_
                  << " Elements limited: " << limitedElements_
                  << " due to side effects: " << notPhysicalElements_
                  << " overall: " << this->numberOfElements_
                  << std::endl;
      }
      */

      if( doCommunicate )
      {
        // communicate dest
        dest.communicate();
      }

      // finalize caller
      caller_.reset();
      uEn_.reset();
      uEnAvg_.reset();
      limitEn_.reset();

      // reset usedAdmissibleFunction value to class default
      usedAdmissibleFunctions_ = usedAdmissibleFunctions();
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
      // increase element counter
      ++this->numberOfElements_;

      // timer for shock detection
      Dune::Timer indiTime;

      // check argument is not zero
      assert( arg_ );

      //- statements
      // set entity to caller
      caller().setEntity( en );

      // initialize linear function
      linearFunction_.bind( en );

      // get function to limit
      const ArgumentFunctionType &U = *(std::get< argumentPosition >( *arg_ ));

      assert( uEn_ );
      // get U on entity
      auto guard = bindGuard( *uEn_, en );
      const LocalFunctionType& uEn = *uEn_;

      // get geometry
      const Geometry& geo = linearFunction_.geometry();

      // cache geometry type
      const GeometryType geomType = geo.type();
      // get bary center of element
      const LocalDomainType& wLocal = geoInfo_.localCenter( geomType );

      DomainType& enBary = linearFunction_.center_;

      // compute barycenter of element
      enBary = geomType.isNone() ? geo.center() : geo.global( wLocal );

      const IntersectionIteratorType endnit = gridPart_.iend( en );
      IntersectionIteratorType niter = gridPart_.ibegin(en);
      if( niter == endnit ) return ;

      // if a component is true, then this component has to be limited
      FieldVector<bool,dimRange> limit(false);
      // total jump vector
      RangeType shockIndicator(0);
      RangeType adaptIndicator(0);

      RangeType& enVal = linearFunction_.value_;

      // if limiter is true then limitation is done
      // when we want to reconstruct in any case then
      // limiter is true but indicator is calculated
      // because of adaptation
      bool limiter = reconstruct_;

      ///////////////////////////////////////////////////
      //
      //  Troubled Cell Indicator
      //
      ///////////////////////////////////////////////////

      // check physicality of data at quadrature points
      // evaluate average returns true if one quadrature point holds unphysical values
      const bool oneNotPhysical = evalAverage( en, uEn, enVal ); // returns enVal

      // make sure that the average value is physical,
      // otherwise everything may be corrupted
      {
        // get barycenter of entity
        const DomainType& enBaryLocal = (int(dimGrid) == int(dimDomain) && ! geomType.isNone() ) ?
                wLocal : geo.local( enBary ) ;

        // check that average value is physical
        if(  discreteModel_.hasPhysical() &&
            !discreteModel_.physical( en, enBaryLocal, enVal ) )
        {
          std::cerr << "Average Value "
                    << enVal
                    << " in point "
                    << enBary
                    << " is Unphysical!"
                    << std::endl << "ABORTED" << std::endl;
          assert( false );
          std::abort();
        }
      }

      if( oneNotPhysical ) // already enough to trigger limitation process
      {
        limiter = true;
        // enable adaptation because of not physical
        adaptIndicator = -1;
      }
      else if ( calcIndicator_ ) // otherwise compute shock indicator
      {
        // check shock indicator
        limiter = calculateIndicator(U, uEn, geo, limiter,                   // parameter
                                     limit, shockIndicator, adaptIndicator); // return values
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

      stepTime_[0] += indiTime.elapsed();
      indiTime.reset();

      // if limit, then limit all components
      limit = limiter;
      {
        // check whether not physical occurred
        if (limiter && shockIndicator[0] < 1.)
        {
          shockIndicator = -1;
          ++notPhysicalElements_;
        }

        // if limiter also adapt if not finite volume scheme
        if( limiter && ! reconstruct_ )
        {
          adaptIndicator = -1;
        }

        // call problem adaptation for setting refinement marker
        // discreteModel_.adaptation( gridPart_.grid() , en, shockIndicator, adaptIndicator );
      }

      ///////////////////////////////////////////////////
      //
      //  Reconstruction and Limiting
      //
      ///////////////////////////////////////////////////

      // boundary is true if boundary segment was found
      // nonconforming is true if entity has at least one non-conforming intersections
      // cartesian is true if the grid is cartesian and no nonconforming refinement present
      typename LimiterUtilityType::Flags flags( cartesianGrid_, limiter );

      // evaluate function
      if( limiter )
      {
        // helper class for evaluation of average value of discrete function
        assert( uEnAvg_ );
        EvalAverage average( *this, *uEnAvg_, discreteModel_, geo.volume() );

        // setup neighbors barycenter and mean value for all neighbors
        LimiterUtilityType::setupNeighborValues( gridPart_, en, average, enBary, enVal, uEn,
                                                 StructuredGrid, flags, barys_, nbVals_ );
      }

      const unsigned int enIndex = indexSet_.index( en );

      // mark entity as finished, even if not limited everything necessary was done
      assert( indexSet_.index( en ) < int(visited_.size()) );
      visited_[ indexSet_.index( en ) ] = true ;

      // if nothing to limit then just return here
      if ( ! limiter ) return ;

      double lpTime = 0;

      // increase number of limited elements
      ++limitedElements_;

      GradientType& gradient = linearFunction_.gradient_;

      // use linear function from LP reconstruction
      if( usedAdmissibleFunctions_ == lp )
      {
        Dune::Timer timer;
        // compute average values on neighboring elements
        fillAverageValues( en, enIndex, *uEnAvg_, enVal );
        assert( linProg_ );
        // compute optimal linear reconstruction
        linProg_->applyLocal( en, gridPart_.indexSet(), values_, gradient );
        lpTime = timer.elapsed();
        stepTime_[2] += lpTime;
      }
      else
      {
        // obtain combination set
        ComboSetType& comboSet = storedComboSets_[ nbVals_.size() ];
        if( comboSet.empty() )
        {
          // create combination set
          LimiterUtilityType::buildComboSet( nbVals_.size(), comboSet );
        }

        // reset values
        gradients_.clear();
        comboVec_.clear();

        if( usedAdmissibleFunctions_ == muscl || usedAdmissibleFunctions_ == muscldg )
        {
          // level is only needed for Cartesian grids to access the matrix caches
          const int matrixCacheLevel = ( flags.cartesian ) ? en.level() : 0 ;
          assert( matrixCacheLevel < (int) matrixCacheVec_.size() );
          MatrixCacheType& matrixCache = matrixCacheVec_[ matrixCacheLevel ];

          // calculate linear functions, stored in gradients_ and comboVec_
          LimiterUtilityType::calculateLinearFunctions( comboSet, geomType, flags,
                                    barys_, nbVals_,
                                    matrixCache,
                                    gradients_,
                                    comboVec_ );
        }

        // add DG Function
        if( usedAdmissibleFunctions_ == dgonly || usedAdmissibleFunctions_ == muscldg )
        {
          addDGFunction( en, geo, uEn, enVal, enBary );
        }

        // Limiting
        std::vector< RangeType > factors;
        LimiterUtilityType::limitFunctions(
            discreteModel_.limiterFunction(), comboVec_, barys_, nbVals_, gradients_, factors );

        // take maximum of limited functions
        LimiterUtilityType::getMaxFunction(gradients_, gradient, factors_[ enIndex ], numbers_[ enIndex ], factors );
      } // end if linProg

      ////////////////////////////////////////////////////////////////////
      //
      //  L2 projection of limited function to destination
      //
      ////////////////////////////////////////////////////////////////////

      // get local funnction for limited values
      assert( limitEn_ );
      DestLocalFunctionType& limitEn = *limitEn_;
      auto limguard = bindGuard( limitEn, en );
      auto interpolGuad = bindGuard( interpolEn_, en );

      // project linearFunction onto limitEn
      interpolate( en, geo, limit, linearFunction_, limitEn );

      assert( checkPhysical(en, geo, limitEn) );

      stepTime_[1] += (indiTime.elapsed() - lpTime);

      //end limiting process
    }

  protected:
    void fillAverageValues( const EntityType& en, const unsigned int enIndex,
                            LocalFunctionType &uLocal, const RangeType& enVal ) const
    {
      // helper class for evaluation of average value of discrete function
      EvalAverage average( *this, uLocal, discreteModel_);

      if( ! valuesComputed_[ enIndex ] )
      {
        values_[ enIndex ] = enVal;
        valuesComputed_[ enIndex ] = true;
      }

      const auto& indexSet = gridPart_.indexSet();
      for (const auto& intersection : intersections(gridPart_, en) )
      {
        if( intersection.neighbor() )
        {
          const auto neighbor = intersection.outside();
          const unsigned int nbIndex = indexSet.index( neighbor );
          if( ! valuesComputed_[ nbIndex ] )
          {
            average.evaluate( neighbor, values_[ nbIndex ] );
            valuesComputed_[ nbIndex ] = true;
          }
        }
      }
    }

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
      assert( Dune::Fem::Capabilities::isHierarchic< DiscreteFunctionSpaceType > :: v );

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

      gradients_.push_back( D );
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
      for(int l=0; l<quadNop; ++l)
      {
        uEn.evaluate( quad[l] , u );
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
      if( discreteModel_.hasPhysical() )
      {
        if( en.type().isNone() )
        {
          RangeType u;
          // for polyhedral cells check corners manually
          // since caching with CornerPointSet won't work
          const int nCorners = geo.corners();
          for( int i=0; i<nCorners; ++i )
          {
            const auto local = geo.local( geo.corner( i  ) );
            uEn.evaluate( local, u );
            if ( ! discreteModel_.physical( en, local, u ) )
            {
              // return notPhysical
              return false;
            }
          }

          // solution is physical
          return true;
        }
        else
        {
          // use LagrangePointSet to evaluate on corners of the
          // geometry and also use caching
          return checkPhysicalQuad( CornerPointSetType( en ), uEn );
        }
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

    template <class LocalFunctionImp, class FunctionSpaceImp, class GridPartImp,
              int polOrd, class StrorageImp >
    struct NumLinearBasis<LocalFunctionImp,
              DiscontinuousGalerkinSpace<FunctionSpaceImp, GridPartImp, polOrd,
                                         StrorageImp> >
    {
      inline static int numBasis(const LocalFunctionImp& lf)
      {
        return dimGrid + 1;
      }
    };


    // L2 projection
    template <class LocalFunctionImp>
    void interpolate(const EntityType& en,
                     const Geometry& geo,
                     const FieldVector<bool,dimRange>& limit,
                     const LinearFunction& linearFunction,
                     LocalFunctionImp& limitEn ) const
    {
      // true if geometry mapping is affine
      const bool affineMapping = geo.affine();

      // set zero dof to zero
      uTmpLocal_.init( en );
      uTmpLocal_.clear();

      const bool constantValue = linearFunction.isConstant();
      if( constantValue )
      {
        const ConstantFunction& constFct = linearFunction;
        interpolEn_( constFct, uTmpLocal_.localDofVector() );
      }
      else
      {
        interpolEn_( linearFunction, uTmpLocal_.localDofVector() );
      }

      // check physicality of projected data
      if ( (! constantValue) && (! checkPhysical(en, geo, uTmpLocal_)) )
      {
        // for affine mapping we only need to set higher moments to zero
        if( Dune::Fem::Capabilities::isHierarchic< DiscreteFunctionSpaceType > :: v &&
            affineMapping )
        {
          // note: the following assumes an ONB made up of moments and affine mapping
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
          // general interpolation of constant value
          const ConstantFunction& constFct = linearFunction;
          interpolEn_( constFct, uTmpLocal_.localDofVector() );
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
        for(int i=0; i<numBasis; ++i)
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
      const Geometry& geo = en.geometry();

      if( Dune::Fem::Capabilities::isHierarchic< DiscreteFunctionSpaceType > :: v && geo.affine() )
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
        // get quadrature
        VolumeQuadratureType quad( en, volumeQuadratureOrder( en ) );

        // set value to zero
        val = 0;

        RangeType aver;
        const int quadNop = quad.nop();
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

      // return true if average value is not physical
      return notphysical;
    }

    bool integrateIntersection(const IntersectionType & intersection,
                               const EntityType& nb,
                               RangeType& shockIndicator,
                               RangeType& adaptIndicator) const
    {
      if( ! conformingGridPart && ! intersection.conforming() )
        return integrateIntersectionImpl< false >(intersection, nb, shockIndicator, adaptIndicator);
      else
        return integrateIntersectionImpl< true  >(intersection, nb, shockIndicator, adaptIndicator);
    }

    template <bool conforming>
    bool integrateIntersectionImpl(const IntersectionType & intersection,
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
    bool integrateBoundary(const EntityType& entity,
                           const IntersectionType & intersection,
                           const QuadratureImp & faceQuadInner,
                           RangeType& shockIndicator,
                           RangeType& adaptIndicator) const
    {
      typedef typename IntersectionType :: Geometry LocalGeometryType;
      const LocalGeometryType& interGeo = intersection.geometry();

      RangeType jump;
      JacobianRangeType dummy ;

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
    bool calculateIndicator(const ArgumentFunctionType& U,
                            const LocalFunctionType& uEn,
                            const Geometry& geo,
                            const bool initLimiter,
                            FieldVector<bool,dimRange>& limit,
                            RangeType& shockIndicator,
                            RangeType& adaptIndicator ) const
    {
      const EntityType& en = uEn.entity();

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

          // internal shock indicator based on Krivodonova et al.
          if( ! extTroubledCellIndicator_ )
          {
            // order of quadrature
            const int jumpQuadOrd = spc_.order();

            // check all neighbors
            if (intersection.neighbor())
            {
              // get neighbor entity
              const EntityType& nb = intersection.outside();

              if (integrateIntersection(intersection, nb, shockIndicator, adaptIndicator))
              {
                shockIndicator = -1;
                return true;
              }
            }

            // check all neighbors
            if ( intersection.boundary() )
            {
              FaceQuadratureType faceQuadInner(gridPart_,intersection, jumpQuadOrd, FaceQuadratureType::INSIDE);

              // initialize intersection
              caller().initializeBoundary( intersection, faceQuadInner );

              if (integrateBoundary(en, intersection, faceQuadInner, shockIndicator, adaptIndicator))
              {
                shockIndicator = -1;
                return true;
              }
            }
          }
        }
      } // end intersection iterator

      if (indicator_ == 1)
        return false;

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

      assert(indicator_!=0 || extTroubledCellIndicator_);
      if( extTroubledCellIndicator_ ) // also true for indicator_=3
      {
        shockIndicator = ( tol_1_ * (*extTroubledCellIndicator_)( U, uEn ) );
      }
      else if( indicator_ == 2 )
      {
        for(int r=0; r<dimRange; ++r)
        {
          shockIndicator[r] = std::abs(shockIndicator[r]) * circFactor * tol_1_;
        }
      }
      else if ( indicator_ == 4 )
      {
        // always limit on each cell
        shockIndicator = 2.;
      }
      else
      {
        std::cout << "wrong indicator selection\n";
        assert(0);
        std::abort();
      }

      for(int r=0; r<dimRange; ++r)
      {
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
      for(int qp=0; qp<quadNop; ++qp)
      {
        // check physicality of value
        const bool physical = caller().checkPhysical( intersection, quad, qp );

        if( checkPhysical && ! physical )
        {
          // return notPhysical
          return true ;
        }

        if( ! inflowIntersection )
        {
          // check intersection
          if( physical && caller().checkDirection(intersection, quad, qp) )
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

  protected:
    DiscreteModelCallerType &caller () const
    {
      assert( caller_ );
      return *caller_;
    }

  private:
    mutable std::unique_ptr< DiscreteModelCallerType > caller_;
    DiscreteModelType& discreteModel_;
    mutable double currentTime_;

    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    GridPartType& gridPart_;

    mutable LinearFunction linearFunction_;

    std::unique_ptr< LinearProgramming > linProg_;

    const IndexSetType& indexSet_;
    const LocalIdSetType& localIdSet_;

    mutable TemporaryLocalFunctionType uTmpLocal_;

    mutable std::unique_ptr< LocalFunctionType > uEn_;
    mutable std::unique_ptr< LocalFunctionType > uEnAvg_;
    mutable std::unique_ptr< DestLocalFunctionType > limitEn_;

    mutable LocalInterpolation<DiscreteFunctionSpaceType> interpolEn_;

    const double orderPower_;
    const DofConversionUtilityType dofConversion_;
    mutable int faceQuadOrd_;
    mutable int volumeQuadOrd_;
    mutable int argOrder_;

    mutable ComboSetMapType storedComboSets_;

    // tolerance to scale shock indicator
    const double tolFactor_;
    std::size_t indicator_;
    const double tol_1_;

    // if true scheme is TVD
    const GeometryInformationType geoInfo_;
    const FaceGeometryInformationType faceGeoInfo_;

    mutable RangeType    phi0_ ;

    mutable std::vector< GradientType > gradients_;
    mutable std::vector< CheckType >    comboVec_;

    mutable std::vector< DomainType > barys_;
    mutable std::vector< RangeType >  nbVals_;
    mutable std::vector< MatrixCacheType > matrixCacheVec_;

    mutable std::vector< RangeType  > factors_;
    mutable std::vector< std::vector< int > > numbers_;

    // vector for stroing the information which elements have been computed already
    mutable std::vector< bool > visited_;

    //! true if limiter is used in adaptive scheme
    const bool adaptive_;
    //! true if grid is cartesian like
    const bool cartesianGrid_;
    mutable int limitedElements_, notPhysicalElements_;
    mutable std::vector<double> stepTime_;

    //! true if troubled cell indicator should be calculated
    mutable bool calcIndicator_;

    //! true if limiter is used as finite volume scheme of higher order
    mutable bool reconstruct_;

    // choice of admissible linear functions
    const AdmissibleFunctions    admissibleFunctions_;
    mutable AdmissibleFunctions  usedAdmissibleFunctions_ ;

    mutable std::vector< RangeType  > values_;
    mutable std::vector< bool >       valuesComputed_;

    // function pointer to external troubled cell indicator, can be nullptr
    TroubledCellIndicatorType extTroubledCellIndicator_;

    mutable int counter_;
    mutable int overallLimited_;
    const bool verbose_;
  }; // end DGLimitPass

} // namespace
} // namespace Dune

#endif // #ifndef DUNE_LIMITERPASS_HH
