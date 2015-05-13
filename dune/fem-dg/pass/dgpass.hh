#ifndef DUNE_FEM_DG_DGPASS_HH
#define DUNE_FEM_DG_DGPASS_HH

#include <dune/common/fvector.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/pass/common/local.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/arrays.hh>

#include <dune/fem-dg/pass/dgmodelcaller.hh>
#include <dune/fem/misc/compatibility.hh>

namespace Dune {
/*! @addtogroup PassHyp
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   v + div(f(x,u)) + A(x,u)\nabla u &=& S(x,u)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ u \f$ is the argument and \f$ v \f$ is computed.
** Weak formulation on a cell T:
** \f[
** \int_T v \phi = -\int_{\partial T} g \phi + \int_T f \cdot \nabla \phi + \int_T Q \phi
** \f]
** with \f$ g \approx f \cdot n + \tilde{A}[u] \cdot n \f$ and \f$ Q \approx S - A \nabla u \f$
** where \f$ \tilde{A} \f$ denotes the arithmetic average and \f$ [u] \f$ the jump of
** \f$ u \f$ over the cell interface.\\
** The discrete model provides the \b analyticalFlux f, the \b source Q and the \b numericalFlux g.
** @{
**************************************************************************/

  //! Concrete implementation of Pass for first hyperbolic systems using
  //! LDG
  template< class DiscreteModelImp, class PreviousPassImp , int passIdImp = -1 >
  class LocalCDGPass :
    public Fem::LocalPass< DiscreteModelImp , PreviousPassImp , passIdImp >
  {
    typedef LocalCDGPass< DiscreteModelImp, PreviousPassImp, passIdImp > ThisType;
  public:

    //- Typedefs and enums
    //! Base class
    typedef Fem::LocalPass< DiscreteModelImp , PreviousPassImp , passIdImp > BaseType;
    typedef typename BaseType::PassIds PassIds;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::EntityType  EntityType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! Iterator over the space
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType  RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType:: BasisFunctionSetType
      BasisFunctionSetType;

    // Types extracted from the underlying grid part types
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridPartType :: template Codim<0>:: GeometryType  Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;

    typedef CDGDiscreteModelCaller< DiscreteModelType, ArgumentType, PassIds > DiscreteModelCallerType;

    typedef typename DestinationType :: DofBlockPtrType DofBlockPtrType;

    // type of local id set
    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > TemporaryLocalFunctionType;

#ifdef USE_CACHED_INVERSE_MASSMATRIX
    // type of communication manager object which does communication
    typedef typename DiscreteFunctionSpaceType::template ToNewDimRange< 1 >::Type ScalarDiscreteFunctionSpaceType;
    typedef Fem::DGMassInverseMassImplementation< ScalarDiscreteFunctionSpaceType, true > MassInverseMassType ;
    typedef typename MassInverseMassType :: KeyType MassKeyType;
    typedef Fem::SingletonList< MassKeyType, MassInverseMassType >  InverseMassProviderType;
    // store a reference to the mass matrix implementation
    typedef MassInverseMassType&  LocalMassMatrixStorageType ;
#else
    //! type of local mass matrix
    typedef Fem::LocalMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType > LocalMassMatrixType;
    typedef LocalMassMatrixType  LocalMassMatrixStorageType;
#endif

    // true if all intersections between element in this grid part are conforming
    static const bool conformingGridPart =
      Fem::GridPartCapabilities::isConforming< GridPartType >::v ;

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
                  const int faceQuadOrd = -1 )
      : BaseType(pass, spc),
        caller_( 0 ),
        discreteModel_(discreteModel),
        arg_(0),
        dest_(0),
        spc_(spc),
        gridPart_(spc_.gridPart()),
        indexSet_(gridPart_.indexSet()),
        visited_(0),
        updEn_(spc_),
        updNeigh_(spc_),
        valEnVec_( 20 ),
        valNbVec_( 20 ),
        valJacEn_( 20 ),
        valJacNb_( 20 ),
        dtMin_(std::numeric_limits<double>::max()),
        minLimit_(2.0*std::numeric_limits<double>::min()),
        volumeQuadOrd_( volumeQuadOrd ),
        faceQuadOrd_( faceQuadOrd ),
        numberOfElements_( 0 ),
#ifdef USE_CACHED_INVERSE_MASSMATRIX
        localMassMatrix_( InverseMassProviderType :: getObject( MassKeyType( gridPart_ ) ) ),
#else
        localMassMatrix_( spc_ , volumeQuadOrd_ ),
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

    //! Destructor
    virtual ~LocalCDGPass()
    {}

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
      // get stopwatch
      Timer timer;

      prepare(arg, dest);

      if( reallyCompute_ )
      {
        const IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it)
        {
          applyLocal(*it);
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
    inline void flux(const DestinationType &u,
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

      arg_  = 0;
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

      numberOfElements_ = 0;

      if( firstThread && dest_ )
      {
        // clear destination (only for first pass in thread parallel version)
        dest_->clear();
      }

      assert( ! caller_ );
      // set arguments to caller
      caller_ = new DiscreteModelCallerType( *arg_, discreteModel_ );
      caller_->setTime( this->time() );

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
      if( doCommunicate && (&dest) )
      {
        // communicate calculated function (not in thread parallel version)
        dest.communicate();
      }

      // call finalize
      if( caller_ )
      {
        delete caller_;
        caller_ = 0;
      }

      arg_  = 0;
      dest_ = 0;
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

    size_t numberOfElements() const
    {
      return numberOfElements_;
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
        // get local function and add update
        LocalFunctionType function = dest_->localFunction(en);

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

      // calcuate element integral
      elementIntegral(entity, updEn_);

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
      // calcuate element integral
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
      ++numberOfElements_ ;

      // only call geometry once, who know what is done in this function
      const Geometry & geo = entity.geometry();

      // only apply volumetric integral if order > 0
      // otherwise this contribution is zero
      if( (spc_.order() > 0) || discreteModel_.hasSource() )
      {
        assert( volumeQuadratureOrder( entity ) >=0 );
        VolumeQuadratureType volQuad(entity, volumeQuadratureOrder( entity ) );
        caller().setEntity(entity, volQuad);

        // if only flux, evaluate only flux
        if ( discreteModel_.hasFlux() && ! discreteModel_.hasSource() )
        {
          evalVolumetricPartFlux(entity, geo, volQuad, updEn);
        }
        else
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
          // set entity for caller
          caller().setEntity( entity );
        }

        /////////////////////////////
        // Surface integral part
        /////////////////////////////
        // get volume of element divided by the DG polynomial factor
        const double envol = entity.geometry().volume() / ( 2.0 * spc_.order( entity ) + 1.0 ) ;

        const IntersectionIteratorType endnit = gridPart_.iend(entity);
        for (IntersectionIteratorType nit = gridPart_.ibegin(entity); nit != endnit; ++nit)
        {
          const IntersectionType& intersection = *nit;

          double nbvol = envol;
          double wspeedS = 0.0;

          if( intersection.neighbor() )
          {
            // get neighbor
            const EntityType nb = Fem::make_entity( intersection.outside() );

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
            caller().initializeBoundary( intersection, faceQuadInner );

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
        // get local function and add update
        LocalFunctionType function = dest_->localFunction( entity );
        function += update;
      }
    }

    //! add update to destination
    template <class LocalFunctionImp>
    void updateFunctionAndApplyMass(
                        const EntityType& entity,
                        LocalFunctionImp& update) const
    {
      if( dest_ )
      {
        // get local function and add update
        LocalFunctionType function = dest_->localFunction( entity );
        function += update;

        // apply local inverse mass matrix
        if (DiscreteModelImp::ApplyInverseMassOperator)
          localMassMatrix_.applyInverse( caller(), entity, function );
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

  protected:
    //! return appropriate quadrature order, default is 2 * order(entity)
    int volumeQuadratureOrder( const EntityType& entity ) const
    {
      return ( volumeQuadOrd_ < 0 ) ? ( spc_.order( entity ) * 2 ) : volumeQuadOrd_ ;
    }

    //! return default face quadrature order
    int defaultFaceQuadOrder( const EntityType& entity ) const
    {
      return (2 * spc_.order( entity )) + 1;
    }

    //! return appropriate quadrature order, default is 2 * order( entity ) + 1
    int faceQuadratureOrder( const EntityType& entity ) const
    {
      return ( faceQuadOrd_ < 0 ) ? defaultFaceQuadOrder( entity ) : faceQuadOrd_ ;
    }

    //! return appropriate quadrature order, default is 2 * order( entity ) + 1
    int faceQuadratureOrder( const EntityType& entity, const EntityType& neighbor ) const
    {
      return ( faceQuadOrd_ < 0 ) ?
        std::max( defaultFaceQuadOrder( entity ), defaultFaceQuadOrder( neighbor ) ) :
        faceQuadOrd_ ;
    }

  protected:
    DiscreteModelCallerType &caller () const
    {
      assert( caller_ );
      return *caller_;
    }

    mutable DiscreteModelCallerType *caller_;
    DiscreteModelType& discreteModel_;

    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType & gridPart_;
    const IndexSetType& indexSet_;

    // indicator for grid walk
    mutable Fem::MutableArray<bool> visited_;

    mutable TemporaryLocalFunctionType updEn_;
    mutable TemporaryLocalFunctionType updNeigh_;

    //! Some helper variables
    mutable Fem::MutableArray< RangeType > valEnVec_;
    mutable Fem::MutableArray< RangeType > valNbVec_;

    mutable Fem::MutableArray< JacobianRangeType > valJacEn_;
    mutable Fem::MutableArray< JacobianRangeType > valJacNb_;

    mutable double dtMin_;
    const double minLimit_;

    const int volumeQuadOrd_, faceQuadOrd_;
    mutable size_t numberOfElements_ ;
    LocalMassMatrixStorageType localMassMatrix_;
    mutable bool reallyCompute_;
  };
//! @}
} // end namespace Dune

#endif
