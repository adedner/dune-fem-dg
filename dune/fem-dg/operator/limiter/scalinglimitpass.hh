#ifndef DUNE_FEMDG_SCALINGLIMITPASS_HH
#define DUNE_FEMDG_SCALINGLIMITPASS_HH

#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem-dg/operator/limiter/limitpass.hh>

//*************************************************************
namespace Dune
{
namespace Fem
{
  /** \brief Concrete implementation of Pass for Limiting.
   *
   *  \ingroup Pass
   *
   *  \note: A detailed description can be found in:
   *
   *  X. Zhang and C.-W. Shu,
   *  Maximum-principle-satisfying and positivity-preserving high order schemes
   *  for conservation laws: Survey and new developments.
   */
  template <class DiscreteModelImp, class PreviousPassImp, int passId >
  class ScalingLimitDGPass
  : public LocalPass<DiscreteModelImp, PreviousPassImp , passId >
  {
    typedef ScalingLimitDGPass< DiscreteModelImp, PreviousPassImp, passId > ThisType;
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

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType           DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType      VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType        FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::IteratorType              IteratorType;

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

    // Various other types
    typedef MutableLocalFunction< DestinationType >                       DestLocalFunctionType;

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

    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef AllGeomTypes< IndexSetType, GridType> GeometryInformationType;

    typedef GeometryInformation< GridType, 1 > FaceGeometryInformationType;

    // get LagrangePointSet of pol order 1
    typedef CornerPointSet< GridPartType >                                 CornerPointSetType;
    // get lagrange point set of order 1
    typedef std::map< Dune::GeometryType, CornerPointSetType >             CornerPointSetContainerType;

    //! type of local mass matrix
    typedef LocalMassMatrix< DiscreteFunctionSpaceType,
                  VolumeQuadratureType > LocalMassMatrixType;

    //! returns true of pass is currently active in the pass tree
    using BaseType :: active ;

    //! type of cartesian grid checker
    typedef CheckCartesian< GridPartType >  CheckCartesianType;

    //! function describing an external troubled cell indicator
    typedef std::shared_ptr< TroubledCellIndicatorBase<ArgumentFunctionType> > TroubledCellIndicatorType;

    typedef Capabilities::DefaultQuadrature<DiscreteFunctionSpaceType >   DefaultQuadratureType;

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
    ScalingLimitDGPass(DiscreteModelType& problem,
                       PreviousPassType& pass,
                       const DiscreteFunctionSpaceType& spc,
                       const int vQ = -1,
                       const int fQ = -1,
                       const bool verbose = Dune::Fem::Parameter::verbose()) :
      ScalingLimitDGPass( problem, pass, spc, Dune::Fem::Parameter::container() )
    {}

    //- Public methods
    /** \brief constructor
     *
     *  \param  problem    Actual problem definition (see problem.hh)
     *  \param  pass       Previous pass
     *  \param  spc        Space belonging to the discrete function local to this pass
     *  \param  vQ         order of volume quadrature
     *  \param  fQ         order of face quadrature
     */
    ScalingLimitDGPass(DiscreteModelType& problem,
                       PreviousPassType& pass,
                       const DiscreteFunctionSpaceType& spc,
                       const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container(),
                       const int vQ = -1,
                       const int fQ = -1,
                       const bool verbose = Dune::Fem::Parameter::verbose()) :
      BaseType(pass, spc),
      caller_(),
      discreteModel_(problem),
      currentTime_(0.0),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      scaledFunction_( gridPart_, discreteModel_.model(), spc_.order() ),
      uEn_(),
      limitEn_(),
      interpolEn_( spc_ ),
      indexSet_( gridPart_.indexSet() ),
      cornerPointSetContainer_(),
      dofConversion_(dimRange),
      faceQuadOrd_( fQ ),
      volumeQuadOrd_( vQ ),
      argOrder_( spc_.order() ),
      geoInfo_( gridPart_.indexSet() ),
      phi0_( 0 ),
      localMassMatrix_( spc_, [](const int order) { return DefaultQuadratureType::volumeOrder(order); } ),
      stepTime_(3, 0.0)
    {
      // we need the flux here
      assert(problem.hasFlux());
    }

    //! Destructor
    virtual ~ScalingLimitDGPass() {}

    //! return default face quadrature order
    static int defaultVolumeQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return DefaultQuadratureType::volumeOrder( space.order( entity ) );
    }

    //! return default face quadrature order
    static int defaultFaceQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return DefaultQuadratureType::surfaceOrder( space.order( entity ) );
    }

    //! return appropriate quadrature order, default is 2 * order(entity)
    int volumeQuadratureOrder( const EntityType& entity ) const
    {
      return ( volumeQuadOrd_ < 0 ) ? defaultVolumeQuadratureOrder( spc_, entity ) : volumeQuadOrd_ ;
    }

    //! return appropriate quadrature order, default is 2 * order( entity ) + 1
    int faceQuadratureOrder( const EntityType& entity ) const
    {
      return ( faceQuadOrd_ < 0 ) ? defaultFaceQuadratureOrder( spc_, entity ) : faceQuadOrd_ ;
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

      //std::cout << "ScalingLimitPass::compute " << std::endl;

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
          const auto& en = *it;
          Dune::Timer localTime;
          applyLocalImp(en);
          stepTime_[2] += localTime.elapsed();
        }

        // finalize
        finalize(arg, dest);
      }
      else
      {
        std::cout << "ScalingLimitPass::compute deactive " << std::endl;
        // get reference to U and pass on to dest
        const ArgumentFunctionType &U = *(std::get< argumentPosition >( arg ));
        dest.assign( U );
      }

      //std::cout << std::endl;

      // accumulate time
      this->computeTime_ += timer.elapsed();
    }

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
        const typename IntersectionGeometry::LocalCoordinate localPoint = interGeo.local( globalPoint );
        const double currentTime = op_.time();

        // check for boundary Value
        if( discreteModel_.hasBoundaryValue( intersection, currentTime, localPoint ) )
        {
          /*
          FaceQuadratureType faceQuadInner( gridPart_, intersection, 0, FaceQuadratureType::INSIDE);
          IntersectionQuadraturePointContext< IntersectionType, EntityType,
            FaceQuadratureType, RangeType, RangeType > local( intersection, entity, faceQuadInner, entityValue, entityValue,
                                                              0, currentTime, volume_);
          discreteModel_.boundaryValue( local, entityValue, neighborValue );
          */
          return true ;
        }
        return false ;
      }
    };


    struct ScaledFunction :
      public Dune::Fem::BindableGridFunction< GridPartType, Dune::Dim<dimRange> >
    {
      typedef Dune::Fem::BindableGridFunction<GridPartType, Dune::Dim<dimRange> > Base;
      typedef typename Base::RangeType   RangeType;
      typedef typename Base::EntityType  EntityType;
      using Base::Base;

      typedef typename DiscreteModelType :: ModelType  ModelType;

    protected:
      const ModelType& model_;
      int order_;

      mutable unsigned int quadId_;
      mutable std::vector< RangeType > tmpVal_ ;
      mutable RangeType theta_;
    public:
      using Base::bind;
      mutable RangeType enVal_;
      mutable int count_ ;

      ScaledFunction( const GridPartType& gridPart,
                      const ModelType& model, const int order )
        : Base( gridPart ),
          model_( model ),
          order_( order ),
          quadId_( -1 ),
          tmpVal_(), theta_(1), enVal_(0), count_(0)
      {}

      void bind( const EntityType& entity, const int order )
      {
        Base::bind( entity );
        order_ = order;
      }

      template <class Quadrature>
      std::vector< RangeType >& getValues( const Quadrature& quadrature ) const
      {
        // adjust size if necessary
        tmpVal_.resize( quadrature.nop() );

        // store quadrature id for later check during evaluate
        quadId_ = quadrature.id();
        count_  = 0;
        return tmpVal_;
      }

      // return reference to scaling factor
      RangeType& resetTheta () const
      {
        // reset theta to default value which does no scaling
        theta_ = 1;
        return theta_;
      }

      template <class Quadrature, class RangeArray>
      void evaluateQuadrature(const Quadrature& quadrature, RangeArray &values) const
      {
        assert( quadrature.id() == quadId_ );
        const int quadNop = quadrature.nop();
        assert( quadNop == int(tmpVal_.size()) );

        for(int qP = 0; qP < quadNop ; ++qP)
        {
          evaluate( qP, values[ qP ] );
        }
      }

      template <class Point>
      void evaluate(const Point &x, typename Base::RangeType &value) const
      {
        evaluate( count_, value );
        ++count_;
        return ;

        std::cout << "Wrong quadrature " << std::endl;
        assert( false );
        std::abort();
      }

      template <class Quadrature>
      void evaluate(const QuadraturePointWrapper< Quadrature > &x, typename Base::RangeType &value) const
      {
        assert( x.quadrature().id() == quadId_ );
        evaluate( x.index(), value );
      }

      unsigned int order() const { return order_; }
      std::string name() const { return "ScalingLimmitPass::ScaledFunction"; } // needed for output

    protected:
      void evaluate(const int qP, typename Base::RangeType &value) const
      {
        assert( qP < int(tmpVal_.size() ));
        // copy value
        value = tmpVal_[ qP ];

        // \tilde{p}(x) = \theta (p(x) - \bar{u}) + \bar{u}
        // limitedRange contains all components that should be modified
        // default is 0,...,dimRange-1
        for( const auto& d : model_.limitedRange() )
        {
          value[ d ] = theta_[ d ] * ( value[ d ] - enVal_[ d ]) + enVal_[ d ];
        }
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
      AssignFunction<typename ArgumentFunctionType ::
                     DiscreteFunctionSpaceType,DiscreteFunctionSpaceType>::assign( U , dest, firstThread );

      limitedElements_ = 0;
      this->numberOfElements_ = 0;

      discreteModel_.clearIndicator();

      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      // time initialisation
      currentTime_ = this->time();

      // initialize caller
      caller_.reset( new DiscreteModelCallerType( *arg_, discreteModel_ ) );
      caller_->setTime(currentTime_);

      uEn_.reset( new LocalFunctionType( U ));
      limitEn_.reset( new DestLocalFunctionType( *dest_ ));

      // calculate maximal indicator (if necessary)
      discreteModel_.indicatorMax();
      discreteModel_.obtainBounds( globalMin_, globalMax_ );

      const size_t size = indexSet_.size( 0 ) ;
      // reset visited vector
      visited_.resize( size );
      std::fill( visited_.begin(), visited_.end(), false );
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

        std::cout << "ScalingLimitPass: Elements limited = " << limitedElements_
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
      limitEn_.reset();
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

    void setTroubledCellIndicator(TroubledCellIndicatorType indicator)
    { // don't need a troubled cell indicator here
    }

  protected:
    //! Perform the limitation on all elements.
    void applyLocalImp(const EntityType& en) const
    {
      // timer for shock detection
      Dune::Timer indiTime;

      // check argument is not zero
      assert( arg_ );

      ++this->numberOfElements_;

      //- statements
      // set entity to caller
      caller().setEntity( en );

      // bind function to entity
      scaledFunction_.bind( en, spc_.order( en ) );

      // get U on entity
      assert( uEn_ );
      auto guard = bindGuard( *uEn_, en );
      const LocalFunctionType& uEn = *uEn_;
      auto interpolGuard = bindGuard( interpolEn_, en );

      // get reference to cell average
      RangeType& enVal = scaledFunction_.enVal_;

      bool limiter = false;

      // check physicality of data
      // evaluate average returns true if not physical
      if( evalAverage( en, uEn, enVal ) )
      {
        limiter = true;
      }

      RangeType& theta = scaledFunction_.resetTheta();

      CornerPointSetType cornerquad( en );
      if( ! checkPhysicalQuad( cornerquad, uEn, enVal, theta ) )
      {
        limiter = true;
      }

      // evaluate uEn on all quadrature points on the intersections
      for (const auto& intersection : intersections(gridPart_, en) )
      {
        int faceQuadOrd = faceQuadratureOrder( en );
        if( intersection.neighbor() )
        {
          faceQuadOrd = std::max( faceQuadOrd, faceQuadratureOrder( intersection.outside() ) );
        }

        FaceQuadratureType faceQuadInner(gridPart_,intersection, faceQuadOrd, FaceQuadratureType::INSIDE);
        if( !checkPhysicalQuad( faceQuadInner, uEn, enVal, theta ) )
        {
          limiter = true;
        }
      }

      // evaluate uEn on all interior quadrature points
      VolumeQuadratureType quad( en, volumeQuadratureOrder( en ) );
      if( ! checkPhysicalQuad( quad, uEn, enVal, theta ) )
      {
        limiter = true;
      }

      for (auto t : theta)
      {
        if (t<1) // if scaling would be done
          limiter = true;
      }

      // scale function
      if( limiter )
      {
        // get local funnction for limited values
        assert( limitEn_ );
        auto lguard = bindGuard( *limitEn_, en );
        DestLocalFunctionType& limitEn = *limitEn_;

        // set all dofs to zero
        limitEn.clear();

        interpolEn_( scaledFunction_, limitEn.localDofVector() );

        // set indicator 1
        discreteModel_.markIndicator();

        // increase number of limited elements
        ++limitedElements_;
      }

      stepTime_[1] += indiTime.elapsed();
      //end limiting process
    }

  protected:
    // check physicality on given quadrature
    template <class QuadratureType, class LocalFunctionImp>
    bool checkPhysicalQuad(const QuadratureType& quad,
                           const LocalFunctionImp& uEn,
                           const RangeType& enVal,
                           RangeType& theta ) const
    {
      // geometry and also use caching
      const EntityType& en = uEn.entity();

      auto& tmpVal = scaledFunction_.getValues( quad );

      bool physical = true;

      // evaluate uEn on all quadrature points
      uEn.evaluateQuadrature( quad , tmpVal );

      const int quadNop = quad.nop();
      for(int l=0; l<quadNop; ++l)
      {
        const RangeType& u = tmpVal[ l ];

        for( const auto& d : discreteModel_.model().limitedRange() )
        {
          const double denominator = std::abs( enVal[ d ] - u[ d ] );
          if( denominator < 1e-12 ) continue ;

          const double upper = std::min( std::abs( enVal[ d ] - globalMax_[ d ] ) / denominator, 1.0 );
          const double lower = std::min( std::abs( enVal[ d ] - globalMin_[ d ] ) / denominator, 1.0 );

          theta[ d ] = std::min( std::min( upper, lower ), theta[d ]);
        }

        // check data
        if ( ! discreteModel_.physical( en, quad.point(l), u ) )
        {
          // notPhysical
          physical = false ;
        }
      }

      // solution is physical
      return physical;
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
      if( Dune::Fem::Capabilities::isHierarchic< DiscreteFunctionSpaceType > :: v
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
        if( int(aver_.size()) < quadNop )
        {
          // resize value vector
          aver_.resize( quadNop );
        }

        // evaluate quadrature at once
        lf.evaluateQuadrature( quad, aver_ );

        for(int qp=0; qp<quadNop; ++qp)
        {
          // check whether value is physical
          notphysical |= (discreteModel_.hasPhysical() && !discreteModel_.physical( en, quad.point( qp ), aver_[ qp ] ) );

          // possibly adjust average value, e.g. calculate primitive vairables and so on
          discreteModel_.adjustAverageValue( en, quad.point( qp ), aver_[ qp ] );

          // apply integration weight
          aver_[ qp ] *= quad.weight(qp) * geo.integrationElement( quad.point(qp) );
          // sum up
          val += aver_[ qp ];
        }

        // mean value, i.e. devide by volume
        val *= 1.0/geo.volume();
      }

      return notphysical;
    }

    // make private
    ScalingLimitDGPass();
    ScalingLimitDGPass(const ScalingLimitDGPass&);
    ScalingLimitDGPass& operator=(const ScalingLimitDGPass&);

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
    mutable std::unique_ptr< DiscreteModelCallerType > caller_;
    DiscreteModelType& discreteModel_;
    mutable double currentTime_;

    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    GridPartType& gridPart_;

    mutable ScaledFunction scaledFunction_;

    mutable std::unique_ptr< LocalFunctionType > uEn_;
    mutable std::unique_ptr< DestLocalFunctionType > limitEn_;

    mutable LocalInterpolation< DiscreteFunctionSpaceType > interpolEn_;

    const IndexSetType& indexSet_;

    CornerPointSetContainerType cornerPointSetContainer_;

    const DofConversionUtilityType dofConversion_;
    mutable int faceQuadOrd_;
    mutable int volumeQuadOrd_;
    mutable int argOrder_;

    // if true scheme is TVD
    const GeometryInformationType geoInfo_;

    mutable RangeType    phi0_ ;

    mutable RangeType    globalMin_;
    mutable RangeType    globalMax_ ;

    mutable std::vector< RangeType >  aver_ ;

    // vector for stroing the information which elements have been computed already
    mutable std::vector< bool > visited_;

    LocalMassMatrixType localMassMatrix_;

    //! true if grid is cartesian like
    mutable int limitedElements_;
    mutable std::vector<double> stepTime_;

  }; // end DGLimitPass

} // namespace
} // namespace Dune

#endif // #ifndef DUNE_LIMITERPASS_HH
