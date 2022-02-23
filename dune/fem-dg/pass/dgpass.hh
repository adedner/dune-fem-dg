#ifndef DUNE_FEM_DG_DGPASS_HH
#define DUNE_FEM_DG_DGPASS_HH

#include <dune/common/fvector.hh>

#include <dune/grid/common/grid.hh>


#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/storage/dynamicarray.hh>

#include <dune/fem-dg/pass/pass.hh>
#include <dune/fem-dg/pass/modelcaller.hh>
#include <dune/fem-dg/pass/discretemodel.hh>

#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{
 /**
  * \brief Concrete implementation of Pass for first hyperbolic systems
  *
  * \ingroup Pass
  * Description: Solver for equations of the form
  * \f{eqnarray*}
  *   v + \nabla\cdot(F(x,u) + A(x,u)\nabla u) &=& S(x,u)  \quad\mbox{in}\quad \Omega    \\
  * \f}
  * where \f$ u \f$ is the argument and \f$ v \f$ is computed.
  * Weak formulation on a cell T:
  * \f[
  * \int_T v \varphi = -\int_{\partial T} g \varphi + \int_T F \cdot \nabla \varphi + \int_T Q \varphi
  * \f]
  * with
  * \f[ g \approx F \cdot n + \tilde{A}[u] \cdot n \f]
  * and
  * \f[ Q \approx S - A \nabla u \f]
  * where \f$ \tilde{A} \f$ denotes the arithmetic average and \f$ [u] \f$ the jump of
  * \f$ u \f$ over the cell interface.
  *
  * The discrete model provides:
  * * \f$F\f$: `analyticalFlux()`,
  * * \f$Q\f$: `source()`,
  * * \f$g\f$: `numericalFlux()`.
  */
  template< class DiscreteModelImp, class PreviousPassImp , int passIdImp = -1 >
  class LocalCDGPass :
    public LocalPass< DiscreteModelImp , PreviousPassImp , passIdImp >
  {
    typedef LocalCDGPass< DiscreteModelImp, PreviousPassImp, passIdImp > ThisType;
  public:

    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp , PreviousPassImp , passIdImp > BaseType;
    typedef typename BaseType::PassIds                                  PassIds;

    //! Repetition of template arguments
    typedef DiscreteModelImp                                            DiscreteModelType;
    typedef PreviousPassImp                                             PreviousPassType;

    // Types from the base class
    typedef typename BaseType::EntityType                               EntityType;
    typedef typename BaseType::ArgumentType                             ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType         DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType    VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType      FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Iterator over the space
    typedef typename DiscreteFunctionSpaceType::IteratorType            IteratorType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType                GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType            GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType              DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType               RangeType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType         DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType          RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType       JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType    BasisFunctionSetType;

    // Types extracted from the underlying grid part types
    typedef typename GridPartType::IntersectionIteratorType             IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType                     IntersectionType;
    typedef typename GridPartType::template Codim<0>::GeometryType      Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType                 LocalFunctionType;

    typedef CDGDiscreteModelCaller< DiscreteModelType, ArgumentType, PassIds > DiscreteModelCallerType;

    typedef typename DestinationType::DofBlockPtrType                   DofBlockPtrType;

    // type of local id set
    typedef typename GridPartType::IndexSetType                         IndexSetType;

    typedef Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType >    TemporaryLocalFunctionType;
    typedef Fem::MutableLocalFunction< DestinationType >                MutableLocalFunctionType;

#ifdef USE_CACHED_INVERSE_MASSMATRIX
    // type of communication manager object which does communication
    typedef typename DiscreteFunctionSpaceType::template ToNewDimRange< 1 >::Type ScalarDiscreteFunctionSpaceType;
    typedef Fem::DGMassInverseMassImplementation< ScalarDiscreteFunctionSpaceType, true > MassInverseMassType ;
    typedef typename MassInverseMassType::KeyType                        MassKeyType;
    typedef Fem::SingletonList< MassKeyType, MassInverseMassType >       InverseMassProviderType;
    // store a reference to the mass matrix implementation
    typedef MassInverseMassType&                                         LocalMassMatrixStorageType;
#else
    //! type of local mass matrix
    typedef Fem::LocalMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType > LocalMassMatrixType;
    typedef LocalMassMatrixType                                          LocalMassMatrixStorageType;
#endif

    // true if all intersections between element in this grid part are conforming
    static const bool conformingGridPart =
      Fem::GridPartCapabilities::isConforming< GridPartType >::v ;

    typedef Capabilities::DefaultQuadrature<DiscreteFunctionSpaceType >  DefaultQuadratureType;

  public:
    //- Public methods
    //! Constructor
    //! \param discreteModel Actual discrete model definition (see dgdiscretemodels.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order
    LocalCDGPass( DiscreteModelType& discreteModel,
                  PreviousPassType& pass,
                  const DiscreteFunctionSpaceType& spc,
                  const int volumeQuadOrd = -1,
                  const int faceQuadOrd = -1,
                  const bool verbose = Dune::Fem::Parameter::verbose() )
      : LocalCDGPass( discreteModel, pass, spc, Dune::Fem::Parameter::container(), volumeQuadOrd, faceQuadOrd, verbose )
    {}

    //- Public methods
    //! Constructor
    //! \param discreteModel Actual discrete model definition (see dgdiscretemodels.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order
    LocalCDGPass( DiscreteModelType& discreteModel,
                  PreviousPassType& pass,
                  const DiscreteFunctionSpaceType& spc,
                  const Dune::Fem::ParameterReader& /* parameter */,
                  const int volumeQuadOrd = -1,
                  const int faceQuadOrd = -1,
                  const bool verbose = Dune::Fem::Parameter::verbose() )
      : BaseType(pass, spc),
        caller_(),
        discreteModel_(discreteModel),
        arg_( nullptr ),
        dest_( nullptr ),
        spc_(spc),
        gridPart_(spc_.gridPart()),
        indexSet_(gridPart_.indexSet()),
        visited_(),
        updEn_(spc_),
        updNeigh_(spc_),
        localFunc_(),
        valEnVec_( 20 ),
        valNbVec_( 20 ),
        valJacEn_( 20 ),
        valJacNb_( 20 ),
        dtMin_(std::numeric_limits<double>::max()),
        minLimit_(2.0*std::numeric_limits<double>::min()),
        volumeQuadOrd_( volumeQuadOrd ),
        faceQuadOrd_( faceQuadOrd ),
#ifdef USE_CACHED_INVERSE_MASSMATRIX
        localMassMatrix_( InverseMassProviderType :: getObject( MassKeyType( gridPart_ ) ) ),
#else
        localMassMatrix_( spc_ , [](const int order) { return DefaultQuadratureType::volumeOrder(order); } ),
#endif
        reallyCompute_( true )
    {
      // make sure that either a ghost layer or an overlap layer is there for
      // communication of the data, otherwise the scheme will not work
      if( spc_.gridPart().comm().size() > 1 )
      {
        if( spc_.grid().ghostSize( 0 ) <= 0 && spc_.grid().overlapSize( 0 ) <= 0 )
        {
          DUNE_THROW(InvalidStateException, "No ghost or overlap layer present, which is needed for communication");
        }
      }

      valEnVec_.setMemoryFactor( 1.1 );
      valNbVec_.setMemoryFactor( 1.1 );
      valJacEn_.setMemoryFactor( 1.1 );
      valJacNb_.setMemoryFactor( 1.1 );
    }

    //! print tex info
    void printTexInfo(std::ostream& out) const {
      BaseType::printTexInfo(out);
      out << "LocalCDGPass: "
          << "\\\\ \n";
    }

    //! switch upwind if necessary
    void switchUpwind()
    {
      discreteModel_.switchUpwind();
    }

    void setTime ( double time )
    {
      discreteModel_.setTime( time );
      BaseType::setTime( time );
    }

    template <class AdaptationType>
    void setAdaptation( AdaptationType& adHandle, double weight )
    {
      discreteModel_.setAdaptation( adHandle, weight );
    }

    //! Estimate for the timestep size
    double timeStepEstimateImpl() const
    {
      // return time step size valid for this pass
      return dtMin_ ;
    }

    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      compute( arg, dest, std::numeric_limits<int>::max() );
    }

    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest, const size_t breakAfter ) const
    {
      // get stopwatch
      Dune::Timer timer;

      prepare(arg, dest);

      if( reallyCompute_ )
      {
        size_t count = 0;
        const auto endit = spc_.end();
        for( auto it = spc_.begin(); (it != endit) && (count < breakAfter) ; ++it, ++count )
        {
          applyLocal( *it );
        }

        finalize(arg, dest);
      }

      // accumulate time
      this->computeTime_ += timer.elapsed();
    }

    void applyForBlock( const EntityType& entity,
                        const int l,
                        DofBlockPtrType& block,
                        DestinationType& matrixRow ) const
    {
      // set base function dof to 1
      (*block)[ l ] = 1;

      // clear target
      matrixRow.clear();

      // apply operator locally
      applyLocal( entity );

      // apply -1
      matrixRow *= -1.0;

      // calc op( 0 )
      (*block)[ l ] = 0;

      // add right hand side
      applyLocal( entity );
    }

    template <class Entity, class Intersection, class Quadrature>
    inline void numericalFlux(const DestinationType &u,
                              const Entity &entity, const Entity &nb,
                              const Intersection &intersection,
                              const Quadrature &faceQuadInner, const Quadrature &faceQuadOuter,
                              const int l,
                              typename DestinationType::RangeType &fluxEn,
                              typename DestinationType::RangeType &fluxNb) const
    {
      // add u to argument tuple
      this->previousPass_.pass(u);
      typename PreviousPassType::NextArgumentType prevArg=this->previousPass_.localArgument();
      typename BaseType::TotalArgumentType totalArg(&u, prevArg);

      DiscreteModelCallerType caller( totalArg, discreteModel_ );
      caller.setTime( this->time() );

      caller.setEntity( entity );
      caller.initializeIntersection( nb, intersection, faceQuadInner, faceQuadOuter );

      JacobianRangeType diffFluxEn, diffFluxNb;

      caller.numericalFlux(intersection, faceQuadInner, faceQuadOuter, l,
                           fluxEn, fluxNb, diffFluxEn, diffFluxNb);

      arg_ = nullptr;
    }

    //! In the preparations, store pointers to the actual arguments and
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      prepare( arg, dest, true );
    }

    //! In the preparations, store pointers to the actual arguments and
    //! destinations. Filter out the "right" arguments for this pass.
    //! This is the version use with ThreadPass
    void prepare(const ArgumentType& arg, DestinationType& dest, const bool firstThread ) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      this->numberOfElements_ = 0;

      if( firstThread && dest_ )
      {
        // clear destination (only for first pass in thread parallel version)
        dest_->clear();
      }

      // set arguments to caller
      caller_.reset( new DiscreteModelCallerType( *arg_, discreteModel_ ) );
      caller_->setTime( this->time() );

      localFunc_.reset( new MutableLocalFunctionType( *dest_ ) );

      // resize indicator function
      visited_.resize( indexSet_.size(0) );
      // set all values to false
      const size_t indSize = visited_.size();
      for(size_t i=0; i<indSize; ++i) visited_[i] = false;

      // time initialisation to max value
      dtMin_ = std::numeric_limits<double>::max();
    }

    //! Some timestep size management.
    void doFinalize(DestinationType& dest, const bool doCommunicate) const
    {
      if( doCommunicate )
      {
        // communicate calculated function (not in thread parallel version)
        dest.communicate();
      }

      // call finalize
      caller_.reset();
      localFunc_.reset();

      arg_  = nullptr;
      dest_ = nullptr;
    }

    //! Some timestep size management.
    //! This is the version use with ThreadPass
    void finalize(const ArgumentType& arg, DestinationType& dest, const bool doCommunicate ) const
    {
      doFinalize( dest, doCommunicate );
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      doFinalize( dest, true );
    }

    //! this pass needs communication only when hasFlux on discrete model is true
    bool requireCommunication () const { return discreteModel_.hasFlux(); }

  protected:
    struct DefaultNBChecker
    {
      bool operator ()(const EntityType& , const EntityType& nb ) const
      {
#if HAVE_MPI
        // only update when neighbor is interior
        return nb.partitionType() == InteriorEntity ;
#else
        return true ;
#endif
      }

      //! returns true if the intersection with neighbor nb should be skipped
      bool skipIntersection( const EntityType& nb ) const
      {
        return false ;
      }
    };

  public:
    // compatibility
    void applyLocal(EntityType& en) const
    {
      const EntityType& entity = en;
      this->applyLocal( entity );
    }

    void applyLocal(const EntityType& entity ) const
    {
      applyLocal( entity, DefaultNBChecker() );
    }

    template <class NeighborChecker>
    void applyLocal(const EntityType& entity,
                    const NeighborChecker& nbChecker ) const
    {
      // init local function
      initLocalFunction( entity , updEn_ );

      // call real apply local
      applyLocal(entity, updEn_, nbChecker );

      // add update to real function
      updateFunctionAndApplyMass(entity, updEn_ );
    }

    void applyLocalMass( const EntityType& en) const
    {
      if( dest_ )
      {
        assert( localFunc_ );
        MutableLocalFunctionType& function = *localFunc_;

        // get local function and add update
        auto guard = bindGuard( function, en );

        // apply local inverse mass matrix
        localMassMatrix_.applyInverse( caller(), en, function );
      }
    }

    // only calculate element integral part and interior fluxes
    void elementIntegral(const EntityType& entity ) const
    {
      // init local function
      initLocalFunction( entity , updEn_ );

      // call real apply local
      elementIntegral(entity, updEn_);

      // add update to real function (do not apply local inverse mass yet)
      updateFunction(entity, updEn_ );
    }

    // only calculate element integral part and interior fluxes
    template <class NeighborChecker>
    void applyLocalInterior(const EntityType& entity, const NeighborChecker& nbChecker ) const
    {
      // init local function
      initLocalFunction( entity , updEn_ );

      // calculate element integral, but also set entity in any case
      elementIntegral(entity, updEn_, true );

      // calculate surface integral for interior and boundary intersections
      surfaceIntegral(entity, updEn_, nbChecker );

      // add update to real function (do not apply local inverse mass yet)
      updateFunction(entity, updEn_ );
    }

    // only calculate integral that need data from other processes
    template <class NeighborChecker>
    void applyLocalProcessBoundary(const EntityType& entity,
                                   const NeighborChecker& nbChecker ) const
    {
      // init local function
      initLocalFunction( entity , updEn_ );

      // calculate surface integral for intersections at the process boundary
      surfaceIntegral(entity, updEn_, nbChecker, true );

      // add update to real function
      updateFunctionAndApplyMass(entity, updEn_ );
    }

  protected:
    //! local integration
    template <class NeighborChecker>
    void applyLocal(const EntityType& entity,
                    TemporaryLocalFunctionType& updEn,
                    const NeighborChecker& nbChecker ) const
    {
      // calculate element integral
      elementIntegral( entity, updEn, true );

      // calculate surface integral
      surfaceIntegral( entity, updEn, nbChecker );
    }

    // compute element integral for given entity
    void elementIntegral(const EntityType& entity,
                         TemporaryLocalFunctionType& updEn,
                         const bool alsoSetEntity = false ) const
    {
      // increase element counter
      ++this->numberOfElements_ ;

      const bool hasSource = discreteModel_.hasSource();
      const bool hasMass   = discreteModel_.hasMass();
      // only apply volumetric integral if order > 0
      // otherwise this contribution is zero
      if( (spc_.order() > 0) || hasSource || hasMass )
      {
        assert( volumeQuadratureOrder( entity ) >=0 );
        VolumeQuadratureType volQuad(entity, volumeQuadratureOrder( entity ) );

        // only call geometry once, who know what is done in this function
        const Geometry & geo = entity.geometry();


        caller().setEntity(entity, volQuad);

        // if only flux, evaluate only flux
        if ( discreteModel_.hasFlux() && ! hasSource )
        {
          evalVolumetricPartFlux(entity, geo, volQuad, updEn);
        }
        else if( hasSource )
        {
          // evaluate flux and source
          evalVolumetricPartBoth(entity, geo, volQuad, updEn);
        }
      }
      else if( alsoSetEntity )
      {
        caller().setEntity( entity );
      }
    }

    // compute surface integral for given entity
    template <class NeighborChecker>
    void surfaceIntegral(const EntityType& entity,
                         TemporaryLocalFunctionType& updEn,
                         const NeighborChecker& nbChecker,
                         const bool setEntity = false ) const
    {
      // only compute boundary integral, when setEntity is false,
      // this is the case in applyLocal and interiorIntergra,
      // but not in processBoundaryIntegral
      const bool computeBoundary = ! setEntity ;

      if ( discreteModel_.hasFlux() )
      {
        if( setEntity )
        {
          if( caller().hasMass() )
          {
            // if mass term is present then we need to evaluate the
            // volume quadrature here for later mass term evaluation
            VolumeQuadratureType volQuad(entity, volumeQuadratureOrder( entity ) );
            caller().setEntity(entity, volQuad);
          }
          else
          {
            // set entity for caller, no quadrature needed
            caller().setEntity( entity );
          }
        }

        /////////////////////////////
        // Surface integral part
        /////////////////////////////
        // get volume of element divided by the DG polynomial factor
        const double envol = entity.geometry().volume() / ( 2.0 * spc_.order( entity ) + 1.0 ) ;

        for (const auto& intersection : intersections(gridPart_, entity) )
        {
          double nbvol = envol;
          double wspeedS = 0.0;

          if( intersection.neighbor() )
          {
            // get neighbor
            const EntityType& nb = intersection.outside();

            // check whether we have to skip this intersection
            if( nbChecker.skipIntersection( nb ) )
            {
              continue ;
            }

            // true if neighbor values can be updated (needed for thread parallel version)
            const bool canUpdateNeighbor = nbChecker( entity, nb ) && dest_ ;

            if( ! visited_[ indexSet_.index( nb ) ] )
            {
              // for conforming situations apply Quadrature given
              if( ! conformingGridPart && ! intersection.conforming() )
              {
                // occurs in a non-conforming grid
                assert( conformingGridPart == false );

                // apply neighbor part, return is volume of neighbor which is
                // needed below
                nbvol = applyLocalNeighbor< false >
                            (intersection, nb, faceQuadratureOrder( entity, nb ),
                             updEn, updNeigh_, envol,
                             wspeedS,
                             canUpdateNeighbor );
              }
              else
              {
                // apply neighbor part, return is volume of neighbor which is
                // needed below
                nbvol = applyLocalNeighbor< true >
                            (intersection, nb, faceQuadratureOrder( entity, nb ),
                             updEn, updNeigh_ , envol,
                             wspeedS,
                             canUpdateNeighbor );
              }
            } // end if do something
          } // end if neighbor
          // compute boundary only in applyLocal and interiorIntegral
          else if( computeBoundary && intersection.boundary() )
          {
            FaceQuadratureType faceQuadInner(gridPart_, intersection, faceQuadratureOrder( entity ),
                                             FaceQuadratureType::INSIDE);

            // initialize intersection
            caller().initializeBoundary( entity, intersection, faceQuadInner );

            const unsigned int faceQuadInner_nop = faceQuadInner.nop();

            if( valEnVec_.size() < faceQuadInner_nop )
              valEnVec_.resize( faceQuadInner_nop );

            if( valJacEn_.size() < faceQuadInner_nop )
              valJacEn_.resize( faceQuadInner_nop );

            for (unsigned int l = 0; l < faceQuadInner_nop; ++l)
            {
              RangeType& fluxEn = valEnVec_[ l ];
              JacobianRangeType& diffFluxEn = valJacEn_[ l ];

#ifndef NDEBUG
              fluxEn = 0;
              diffFluxEn = 0;
#endif

              // eval boundary Flux
              wspeedS += caller().boundaryFlux(intersection,
                                               faceQuadInner,
                                               l,
                                               fluxEn,
                                               diffFluxEn)
                       * faceQuadInner.weight(l);

              // apply weights
              fluxEn     *= -faceQuadInner.weight( l );
              diffFluxEn *= -faceQuadInner.weight( l );

            }

            if( DiscreteModelCallerType :: evaluateJacobian )
            {
              // update local functions at once
              updEn.axpyQuadrature( faceQuadInner, valEnVec_, valJacEn_ );
            }
            else
            {
              // update only with range values
              updEn.axpyQuadrature( faceQuadInner, valEnVec_ );
            }
          } // end if boundary

          if (wspeedS > minLimit_ )
          {
            double minvolS = std::min(envol , nbvol);
            dtMin_ = std::min(dtMin_, (minvolS/wspeedS) );
          }
        } // end intersection loop
      } // end discreteModel_.hasFlux()

      // this entity is finised by now
      visited_[ indexSet_.index( entity ) ] = reallyCompute_ ;
    }

    // initialize local update function
    template <class LocalFunctionImp>
    void initLocalFunction(const EntityType& entity, LocalFunctionImp& update) const
    {
      // init local function
      update.init( entity );
      // clear dof values
      update.clear();
    }

    //! add update to destination
    template <class LocalFunctionImp>
    void updateFunction(const EntityType& entity,
                        LocalFunctionImp& update) const
    {
      if( dest_ )
      {
        // add update to destination
        dest_->addLocalDofs( entity, update );
      }
    }

    //! add update to destination
    template <class LocalFunctionImp>
    void updateFunctionAndApplyMass(
                        const EntityType& entity,
                        LocalFunctionImp& update ) const
    {
      if( dest_ )
      {
        assert( localFunc_ );
        MutableLocalFunctionType& function = *localFunc_;

        // get local function and add update
        auto guard = bindGuard( function, entity );

        // add update
        function += update;

        // apply local inverse mass matrix
        if (DiscreteModelImp::ApplyInverseMassOperator)
        {
          localMassMatrix_.applyInverse( caller(), entity, function );
        }
      }
    }

    //////////////////////////////////////////
    // Volumetric integral part only flux
    //////////////////////////////////////////
    template <class LocalFunctionImp>
    void evalVolumetricPartFlux(const EntityType& entity,
                                const Geometry& geo,
                                const VolumeQuadratureType& volQuad,
                                LocalFunctionImp& updEn) const
    {
      const unsigned int volQuad_nop = volQuad.nop();
      if( valJacEn_.size() < volQuad_nop )
      {
        valJacEn_.resize( volQuad_nop );
      }

      for (unsigned int l = 0; l < volQuad_nop; ++l)
      {
        JacobianRangeType& flux = valJacEn_[ l ];
#ifndef NDEBUG
        flux = 0;
#endif
        // evaluate analytical flux and source
        caller().analyticalFlux(entity, volQuad, l, flux );

        const double intel = geo.integrationElement(volQuad.point(l))
                           * volQuad.weight(l);

        // apply integration weights
        flux *= intel;

      }

      // add values to local function
      updEn.axpyQuadrature( volQuad, valJacEn_ );
    }

    //////////////////////////////////////////
    // Volumetric integral part only flux
    //////////////////////////////////////////
    template <class LocalFunctionImp>
    void evalVolumetricPartBoth(const EntityType& entity,
                                const Geometry& geo,
                                const VolumeQuadratureType& volQuad,
                                LocalFunctionImp& updEn) const
    {
      const unsigned int volQuad_nop = volQuad.nop();
      if( valEnVec_.size() < volQuad_nop )
      {
        valEnVec_.resize( volQuad_nop );
      }

      if( valJacEn_.size() < volQuad_nop )
      {
        valJacEn_.resize( volQuad_nop );
      }

      for (unsigned int l = 0; l < volQuad_nop; ++l)
      {
        RangeType& source       = valEnVec_[ l ];
        JacobianRangeType& flux = valJacEn_[ l ];

#ifndef NDEBUG
        source = 0 ;
        flux = 0;
#endif

        // evaluate analytical flux and source
        const double dtEst =
          caller().analyticalFluxAndSource(entity, volQuad, l, flux, source );

        const double intel = geo.integrationElement(volQuad.point(l))
                           * volQuad.weight(l);

        // apply integration weights
        source *= intel;
        flux   *= intel;

        if( dtEst > minLimit_ )
          dtMin_ = std::min(dtMin_, dtEst);

      }
      updEn.axpyQuadrature(volQuad, valEnVec_, valJacEn_ );
    }

    template <bool conforming, class LocalFunctionImp >
    double applyLocalNeighbor(const IntersectionType & intersection,
                              const EntityType & nb,
                              const int faceQuadratureOrder,
                              LocalFunctionImp & updEn,
                              LocalFunctionImp & updNb,
                              const double enVol,
                              double & wspeedS,
                              const bool canUpdateNeighbor ) const
    {
      // make sure correct method is called
      assert( intersection.conforming() == conforming );

      // check on quadrature order
      assert( faceQuadratureOrder >= 0 );

      // use IntersectionQuadrature to create appropriate face quadratures
      typedef Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
      typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;

      // create intersection quadrature (without neighbor check)
      IntersectionQuadratureType interQuad( gridPart_, intersection, faceQuadratureOrder, true );

      // get appropriate references
      const QuadratureImp &faceQuadInner = interQuad.inside();
      const QuadratureImp &faceQuadOuter = interQuad.outside();

      // get geometry of neighbor
      const Geometry & nbGeo = nb.geometry();

      // get volume of neighbor divided by the DG polynomial factor
      const double nbVol = nbGeo.volume() / ( 2.0 * spc_.order( nb ) + 1.0 ) ;

      // set neighbor and initialize intersection
      caller().initializeIntersection( nb, intersection, faceQuadInner, faceQuadOuter );

      const unsigned int faceQuadInner_nop = faceQuadInner.nop();

      if( valNbVec_.size() < faceQuadInner_nop )
      {
        valEnVec_.resize( faceQuadInner_nop );
        valNbVec_.resize( faceQuadInner_nop );

        valJacEn_.resize( faceQuadInner_nop );
        valJacNb_.resize( faceQuadInner_nop );
      }

      for (unsigned int l = 0; l < faceQuadInner_nop; ++l)
      {
        RangeType& fluxEn = valEnVec_[ l ];
        RangeType& fluxNb = valNbVec_[ l ];

        JacobianRangeType& diffFluxEn = valJacEn_[ l ];
        JacobianRangeType& diffFluxNb = valJacNb_[ l ];

#ifndef NDEBUG
        fluxEn = 0;
        fluxNb = 0;

        diffFluxEn = 0;
        diffFluxNb = 0;
#endif

        // calculate num flux for multiplication with the basis
        // wspeedS = fastest wave speed
        wspeedS += caller().numericalFlux(intersection,
                                         faceQuadInner,
                                         faceQuadOuter,
                                         l,
                                         fluxEn, fluxNb,
                                         diffFluxEn, diffFluxNb)
                 * faceQuadInner.weight(l);

        // apply weights
        fluxEn     *= -faceQuadInner.weight(l);
        fluxNb     *=  faceQuadOuter.weight(l);

        diffFluxEn *= -faceQuadInner.weight(l);
        diffFluxNb *=  faceQuadOuter.weight(l);
      }

      // update local functions at once
      if( DiscreteModelCallerType :: evaluateJacobian )
        updEn.axpyQuadrature( faceQuadInner, valEnVec_, valJacEn_ );
      else
        updEn.axpyQuadrature( faceQuadInner, valEnVec_ );

      // if we can also update the neighbor
      // this can be different in thread parallel programs
      if( canUpdateNeighbor )
      {
        // init local function
        initLocalFunction( nb, updNb );

        // add update of the neighbor to the new discrete
        // function which represents L(u_h)
        // this update is convenient, so that we
        // don't have to visit neighbor directly
        if( DiscreteModelCallerType :: evaluateJacobian )
          updNb.axpyQuadrature( faceQuadOuter, valNbVec_, valJacNb_ );
        else
          updNb.axpyQuadrature( faceQuadOuter, valNbVec_ );

        // do the update to global values
        updateFunction( nb, updNb );
      }
      return nbVol;
    }

  private:
    LocalCDGPass();
    LocalCDGPass(const LocalCDGPass&);
    LocalCDGPass& operator=(const LocalCDGPass&);

  public:
    //! return default face quadrature order
    static int defaultVolumeQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return DefaultQuadratureType::volumeOrder( space.order( entity ) );
    }

    //! return default face quadrature
    static int defaultFaceQuadratureOrder( const DiscreteFunctionSpaceType& space, const EntityType& entity )
    {
      return DefaultQuadratureType::surfaceOrder( space.order( entity ) );
    }

  protected:
    //! return appropriate quadrature order
    int volumeQuadratureOrder( const EntityType& entity ) const
    {
      return ( volumeQuadOrd_ < 0 ) ? defaultVolumeQuadratureOrder( spc_, entity ) : volumeQuadOrd_ ;
    }

    //! return appropriate quadrature order
    int faceQuadratureOrder( const EntityType& entity ) const
    {
      return ( faceQuadOrd_ < 0 ) ? defaultFaceQuadratureOrder( spc_, entity ) : faceQuadOrd_ ;
    }

    //! return appropriate quadrature order
    int faceQuadratureOrder( const EntityType& entity, const EntityType& neighbor ) const
    {
      return ( faceQuadOrd_ < 0 ) ?
        std::max( defaultFaceQuadratureOrder( spc_, entity ), defaultFaceQuadratureOrder( spc_, neighbor ) ) :
        faceQuadOrd_ ;
    }

  protected:
    DiscreteModelCallerType &caller () const
    {
      assert( caller_ );
      return *caller_;
    }

    mutable std::unique_ptr< DiscreteModelCallerType > caller_;
    DiscreteModelType& discreteModel_;

    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType & gridPart_;
    const IndexSetType& indexSet_;

    // indicator for grid walk
    mutable Fem::DynamicArray<bool> visited_;

    mutable TemporaryLocalFunctionType updEn_;
    mutable TemporaryLocalFunctionType updNeigh_;

    mutable std::unique_ptr< MutableLocalFunctionType > localFunc_;

    //! Some helper variables
    mutable Fem::DynamicArray< RangeType > valEnVec_;
    mutable Fem::DynamicArray< RangeType > valNbVec_;

    mutable Fem::DynamicArray< JacobianRangeType > valJacEn_;
    mutable Fem::DynamicArray< JacobianRangeType > valJacNb_;

    mutable double dtMin_;
    const double minLimit_;

    const int volumeQuadOrd_, faceQuadOrd_;
    LocalMassMatrixStorageType localMassMatrix_;
    mutable bool reallyCompute_;
  };

} // end namespace
} // end namespace Dune

#endif
