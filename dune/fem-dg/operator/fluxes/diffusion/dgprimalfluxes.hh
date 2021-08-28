#ifndef DUNE_FEM_DG_DGPRIMALFLUXES_HH
#define DUNE_FEM_DG_DGPRIMALFLUXES_HH

#include <dune/common/version.hh>

#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>
#include <dune/fem-dg/pass/context.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/storage/dynamicarray.hh>

#include <dune/fem-dg/pass/dgmasspass.hh>

#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{

  // DGPrimalDiffusionFluxImpl
  //----------------------

  /**
   * \brief diffusion flux
   *
   * \ingroup DiffusionFluxes
   */
  template <class DiscreteFunctionSpaceImp,
            class Model,
            class FluxParameterImp >
  class DGPrimalDiffusionFluxImpl
   : public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model, FluxParameterImp >
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model, FluxParameterImp >      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;
    static const int dimRange  = DiscreteFunctionSpaceType::dimRange;
    static const int dimGradRange = dimDomain * dimRange;
    static const int polOrd = DiscreteFunctionSpaceType::polynomialOrder;

    typedef typename DiscreteFunctionSpaceType::RangeFieldType        RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType       DomainFieldType;
    typedef FieldVector< DomainFieldType, dimDomain-1 >               FaceDomainType;
    typedef typename DiscreteFunctionSpaceType::DomainType            DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType             RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType     JacobianRangeType;


    typedef typename DiscreteFunctionSpaceType::GridPartType          GridPartType;
    typedef typename GridPartType::IntersectionIteratorType           IntersectionIterator;
    typedef typename GridPartType::IntersectionType                   IntersectionType;
    typedef typename GridPartType::GridType                           GridType;
    typedef typename DiscreteFunctionSpaceType::EntityType            EntityType;
    typedef typename GridPartType::template Codim< 0 >::IteratorType  IteratorType;
    typedef typename GridPartType::IntersectionIteratorType           IntersectionIteratorType;

    typedef typename BaseType::DiscreteGradientSpaceType              DiscreteGradientSpaceType;
    typedef typename DiscreteGradientSpaceType::RangeType             GradientType;
    typedef Fem::TemporaryLocalFunction< DiscreteGradientSpaceType >  LiftingFunctionType;

    typedef Fem::CachingQuadrature< GridPartType, 0,
         Capabilities::DefaultQuadrature<DiscreteGradientSpaceType>::template DefaultQuadratureTraits>   VolumeQuadratureType ;

    typedef Fem::LocalMassMatrix
      < DiscreteGradientSpaceType, VolumeQuadratureType >             LocalMassMatrixType;

    class Lifting
    {
    protected:
#ifdef USE_CACHED_INVERSE_MASSMATRIX
#warning "Using cached inverse local mass matrix"
      // type of communication manager object which does communication
      typedef typename DiscreteGradientSpaceType::template ToNewDimRange< 1 >::Type         ScalarDiscreteFunctionSpaceType;
      typedef Fem::DGMassInverseMassImplementation< ScalarDiscreteFunctionSpaceType, true > MassInverseMassType;
      typedef typename MassInverseMassType::KeyType                                         KeyType;
      typedef Fem::SingletonList< KeyType, MassInverseMassType >                            InverseMassProviderType;
      typedef MassInverseMassType&                                                          LocalMassMatrixStorageType;
#else
      typedef LocalMassMatrixType                                                           LocalMassMatrixStorageType;
#endif

      const DiscreteGradientSpaceType& gradSpc_;
      LiftingFunctionType r_e_;
      LocalMassMatrixStorageType localMassMatrix_;
      unsigned char isInitialized_;

      // prohibit copying
      Lifting( const Lifting& );
    public:
      explicit Lifting( const DiscreteGradientSpaceType& grdSpace )
        : gradSpc_( grdSpace )
        , r_e_( gradSpc_ )
#ifdef USE_CACHED_INVERSE_MASSMATRIX
        , localMassMatrix_( InverseMassProviderType :: getObject( KeyType( gradSpc_.gridPart() ) ) )
#else
        , localMassMatrix_( gradSpc_, [](const int order) { return Capabilities::DefaultQuadrature<DiscreteGradientSpaceType>::volumeOrder(order); } )
#endif
        , isInitialized_( 0 )
      {
      }

      bool isInitialized() const { return isInitialized_ == 2 ; }

      void initialize( const EntityType& entity )
      {
        assert( isInitialized_ != 1 );
        r_e_.init( entity );
        r_e_.clear();
        isInitialized_ = 1;
      }

      LiftingFunctionType& function()
      {
        return r_e_;
      }

      void finalize()
      {
        assert( isInitialized_ == 1 );
        localMassMatrix_.applyInverse( r_e_ );
        isInitialized_ = 2;
      }
    };

  protected:
    using BaseType :: determineDirection;
    using BaseType :: model_;
    using BaseType :: cflDiffinv_;
    using BaseType :: dimensionFactor_;
    using BaseType :: nonconformingFactor_;
    using BaseType :: numericalFlux ;
    using BaseType :: upwind_ ;

  public:
    typedef typename BaseType::ParameterType  ParameterType;
  protected:
    typedef typename BaseType::IdEnum         EnumType;
    typedef typename ParameterType::LiftingEnum    LiftingEnum;
  public:

    using BaseType :: parameter ;

    bool initAreaSwitch() const
    {
      if( method_ == EnumType::cdg2 )
      {
        // when default value is used, then use areSwitch
        if( ( upwind_ - BaseType :: upwindDefault() ).two_norm2() < 1e-10 )
          return true ;
      }
      return upwind_.two_norm2() < 1e-10 ;
    }

  public:
    static const bool evaluateJacobian = true;

    /**
     * \brief constructor reading parameters
     */
    DGPrimalDiffusionFluxImpl( GridPartType& gridPart,
                               const Model& model,
                               const ParameterType& parameters,
                               const EnumType staticMethod ) :
      BaseType( model, true, parameters ),
      gridPart_( gridPart ),
      method_( staticMethod == EnumType::primal ? parameters.getMethod() : staticMethod ),
      penalty_( parameter().penalty() ),
      nipgFactor_( (method_ == EnumType::nipg) ||
                   (method_ == EnumType::bo)
                   ? 0.5 : -0.5 ),
      liftFactor_( parameter().liftfactor() ),
      liftingMethod_( parameter().getLifting() ),
      penaltyTerm_( method_ == EnumType::ip || ((std::abs(  penalty_ ) > 0) &&
                    method_ != EnumType::br2 &&
                    method_ != EnumType::bo )),
      gradSpc_( gridPart_ ),
      LeMinusLifting_( hasLifting() ? new Lifting( gradSpc_ ) : 0 ),
      LePlusLifting_( ( method_ == EnumType::br2 ) ? new Lifting( gradSpc_ ) : 0 ),
#ifdef LOCALDEBUG
      LeMinusLifting2_( ( method_ <= EnumType::cdg ) ? new Lifting( gradSpc_ ) : 0 ),
#endif
      insideIsInflow_ ( true ),
      areaSwitch_( initAreaSwitch() ),
      useTheoryParams_( false ),
      initialized_ ( false )
    {
      // calculate maxNeighborVolumeRatio_
      maxNeighborsVolumeRatio_ = 1.;

      double theoryFactor = parameter().theoryparameters();
      useTheoryParams_ = (theoryFactor > 0.);

      double n_k = DiscreteFunctionSpaceType :: polynomialOrder ;
      ainsworthFactor_ = theoryFactor * 0.5 * n_k * ( n_k + 1.0 );

      int maxNumFaces = 0 ;
      int maxNumOutflowFaces = 0;
      if ( useTheoryParams_ )
      {
        for( const auto& entity : elements( gridPart ) )
        {
          const double insideVol = entity.geometry().volume();
          int numFaces = 0;
          int numOutflowFaces = 0;

          for (const auto& intersection : intersections(gridPart, entity) )
          {
            ++numFaces ;
            if ( intersection.neighbor() )
            {
              double outsideVol = intersection.outside().geometry().volume();
              numOutflowFaces += (determineDirection(areaSwitch_, insideVol,outsideVol,intersection) ? 1 : 0);
              if ( !areaSwitch_ || insideVol/outsideVol < 1)
                maxNeighborsVolumeRatio_ = std::max( maxNeighborsVolumeRatio_, insideVol/outsideVol );
            }
            else
              ++numOutflowFaces;
          }
          maxNumFaces = std::max( maxNumFaces, numFaces );
          maxNumOutflowFaces = std::max( maxNumOutflowFaces, numOutflowFaces );
        }

        initialized_ = true;
        liftFactor_ = 0.0;
        penalty_ = 0.;
        if (method_ == EnumType::cdg2)
        {
          liftFactor_ = theoryFactor * 0.25* ((double) maxNumFaces); // max number of faces here
          //if( ! areaSwitch_ )
          liftFactor_ *= (1.+maxNeighborsVolumeRatio_);
        }
        else if (method_ == EnumType::cdg)
        {
          liftFactor_ = theoryFactor * maxNumOutflowFaces;
        }
        else if (method_ == EnumType::br2)
        {
          liftFactor_ = theoryFactor * maxNumFaces;
        }
        else if( method_ == EnumType::nipg )
        {
          std::cerr << "ERROR: No theory parameters for NIPG" << std::endl;
          DUNE_THROW(InvalidStateException,"No theory parameters for NIPG");
        }
      }

      if( Fem::Parameter :: verbose() )
      {
        std::cout <<"Diff. flux: ";
        diffusionFluxName( std::cout );

        std::cout <<", penalty: ";
        if ( useTheoryParams_ && (method_ == EnumType::ip) )
        {
          std::cout <<"theory (";
          diffusionFluxPenalty( std::cout );
          std::cout << ")";
        }
        else
          diffusionFluxPenalty( std::cout );

        std::cout <<", max neigh. vol. ratio: " <<maxNeighborsVolumeRatio_;
        std::cout <<", liftfactor: " <<liftFactor_;
        std::cout <<", max inflow faces: " <<maxNumOutflowFaces;
        std::cout <<std::endl;
      }
    }

    //! copy constructor (needed for thread parallel programs)
    DGPrimalDiffusionFluxImpl( const DGPrimalDiffusionFluxImpl& other ) :
      BaseType( other ),
      gridPart_( other.gridPart_ ),
      method_( other.method_ ),
      penalty_( other.penalty_ ),
      nipgFactor_( other.nipgFactor_ ),
      liftFactor_( other.liftFactor_ ),
      liftingMethod_( other.liftingMethod_ ),
      penaltyTerm_( other.penaltyTerm_ ),
      gradSpc_( gridPart_ ),
      LeMinusLifting_( hasLifting() ? new Lifting( gradSpc_ ) : 0 ),
      LePlusLifting_( ( method_ == EnumType::br2 ) ? new Lifting( gradSpc_ ) : 0 ),
#ifdef LOCALDEBUG
      LeMinusLifting2_( ( method_ <= EnumType::cdg ) ? new Lifting( gradSpc_ ) : 0 ),
#endif
      maxNeighborsVolumeRatio_( other.maxNeighborsVolumeRatio_ ),
      ainsworthFactor_( other.ainsworthFactor_ ),
      insideIsInflow_ ( other.insideIsInflow_ ),
      areaSwitch_( other.areaSwitch_ ), // used area based switch
      useTheoryParams_( other.useTheoryParams_ ),
      initialized_( other.initialized_ )
    {
    }

    // return reference to gradient discrete function space
    const DiscreteGradientSpaceType& gradientSpace() const { return gradSpc_; }

    double maxNeighborsVolumeRatio() const
    {
      assert( initialized_ );
      return maxNeighborsVolumeRatio_;
    }

    void diffusionFluxName( std::ostream& out ) const
    {
      out << ParameterType::methodNames( method_ );
      if( areaSwitch_ )
        out <<"(area)";
      else
        out <<"(upwind)";
    }

    void diffusionFluxPenalty( std::ostream& out ) const
    {
      out << penalty_;
    }

    void diffusionFluxLiftFactor( std::ostream& out ) const
    {
      out <<liftFactor_;
    }

    //! returns true if lifting has to be calculated
    bool hasLifting () const { return ( method_ <= EnumType::br2 ); }

  protected:
    Lifting& LePlusLifting() const
    {
      assert( LePlusLifting_ );
      return *LePlusLifting_;
    }

    Lifting& LeMinusLifting() const
    {
      assert( LeMinusLifting_ );
      return *LeMinusLifting_;
    }
#ifdef LOCALDEBUG
    Lifting& LeMinusLifting2() const
    {
      assert( LeMinusLifting2_ );
      return *LeMinusLifting2_;
    }
#endif

  public:
    void initialize( const DiscreteFunctionSpaceType &space )
    {
    }

    template <class LocalEvaluationVec, class RangeVector >
    void initializeIntersection(const LocalEvaluationVec& left,
                                const LocalEvaluationVec& right,
                                const RangeVector& uLeftVec,
                                const RangeVector& uRightVec)
    {
      if( hasLifting() )
      {
        computeLiftings( left, right, uLeftVec, uRightVec, (method_ == EnumType::br2) );
      }
    }

    template <class LocalEvaluationVec, class RangeVector >
    void computeLiftings(const LocalEvaluationVec& left,
                         const LocalEvaluationVec& right,
                         const RangeVector& uLeftVec,
                         const RangeVector& uRightVec,
                         const bool computeBoth )
    {
      if( hasLifting() || computeBoth )
      {
        if ( ! LeMinusLifting_ )
          LeMinusLifting_.reset( new Lifting( gradSpc_ ) );

        // define for an intersection e
        //  Ke+ := { e in bnd(Ke+), s * n_Ke+ < 0 }
        //  Ke- := { e in bnd(Ke-), s * n_Ke- > 0 }
        // Notice
        //  l_e = r_e on Ke-
        //  l_e = -r_e on Ke+
        // so that
        //  L_e = 2*r_e on Ke-
        //  L_e = 0 elsewhere

        // get Ke- in entity
        insideIsInflow_ = determineDirection( areaSwitch_, left.entity().geometry().volume(),
                                              right.entity().geometry().volume(),
                                              left.intersection() );

        // get right local evaluation
        const LocalEvaluationVec& inflow = insideIsInflow_ ? right : left;

        const RangeVector& uInflow = insideIsInflow_ ? uRightVec : uLeftVec;

        const size_t quadNoInp = left.quadrature().nop();
        liftingEvalLeMinus_.resize( quadNoInp );

        // back to the real computation, entity=Ke-
        LeMinusLifting().initialize( inflow.entity() );

        // calculate real lifting
        for(size_t qp = 0; qp < quadNoInp; ++qp )
        {
          addLifting( inflow[qp], uInflow[qp], uLeftVec[qp], uRightVec[qp], liftingEvalLeMinus_[qp] );
        }
        // add to local function
        LeMinusLifting().function().axpyQuadrature( inflow.quadrature(), liftingEvalLeMinus_ );

        // LeMinusLifting_ has L_e=2*r_e on Ke-
        LeMinusLifting().finalize( );

        // already evaluate for all quadrature points
        LeMinusLifting().function().evaluateQuadrature( inflow.quadrature(), liftingEvalLeMinus_ );

        if ( computeBoth )
        {
          if ( ! LePlusLifting_ )
            LePlusLifting_.reset( new Lifting( gradSpc_ ) );

          // get right local evaluation
          const LocalEvaluationVec& outflow = ( insideIsInflow_ ) ? left : right;
          const RangeVector& uOutflow = insideIsInflow_ ? uLeftVec : uRightVec;

          // calculate 2*r_e on Ke+
          LePlusLifting().initialize( outflow.entity() );

          const size_t quadNoOutp = right.quadrature().nop();
          liftingEvalLePlus_.resize( quadNoOutp );
          for(size_t qp = 0; qp < quadNoOutp; ++qp )
          {
            // get value of 2*r_e in quadrature point
            addLifting( outflow[qp], uOutflow[qp], uLeftVec[qp], uRightVec[qp], liftingEvalLePlus_[qp] );
          }

          // add to local function
          LePlusLifting().function().axpyQuadrature( outflow.quadrature(), liftingEvalLePlus_ );

          // LePlusLifting_ carries 2*r_e on Ke+
          LePlusLifting().finalize( );

          // already evaluate for all quadrature points, reuse liftTmp here
          LeMinusLifting().function().evaluateQuadrature( outflow.quadrature(), liftingEvalLePlus_ );
        }
      }
    }

#ifdef LOCALDEBUG
    template <class LiftingFunction , class Geometry >
    double integrateLifting( const LiftingFunction& lifting, const Geometry& geometry ) const
    {
      typedef typename LiftingFunction::RangeType RangeType;
      VolumeQuadratureType quad( lifting.entity(), 2 * lifting.order() + 2 );
      const int quadNop = quad.nop();
      RangeType val;
      double sum = 0.0;
      for( int qp = 0; qp < quadNop; ++qp )
      {
        const double weight = quad.weight( qp ) *
          geometry.integrationElement( quad.point( qp ) );
        lifting.evaluate( quad[ qp ], val );
        sum += weight * (val * val);
      }
      return sum;
    }
#endif

    template <class LocalEvaluationVec, class RangeVector, class RangeVector2 >
    void initializeBoundary(const LocalEvaluationVec& local,
                            const RangeVector& uLeftVec,
                            const RangeVector2& uRightVec)
    {
      if( hasLifting() )
      {
        insideIsInflow_ = true;

        LeMinusLifting().initialize( local.entity() );

        const size_t quadNop = local.quadrature().nop();
        liftingEvalLeMinus_.resize( quadNop );
        for(size_t qp = 0; qp < quadNop; ++qp )
        {
          addLifting( local[qp], uLeftVec[qp], uLeftVec[qp], uRightVec[qp], liftingEvalLeMinus_[qp] );
        }
        // add to local function
        LeMinusLifting().function().axpyQuadrature( local.quadrature(), liftingEvalLeMinus_ );
        // finalize function
        LeMinusLifting().finalize();
        // already evaluate all liftings
        LeMinusLifting().function().evaluateQuadrature( local.quadrature(), liftingEvalLeMinus_ );
      }
    }

  protected:
    template <class LocalEvaluation, class LiftingFunction >
    void addLifting(const LocalEvaluation& inside,
                    const RangeType& uInflow,
                    const RangeType& uLeft,
                    const RangeType& uRight,
                    LiftingFunction& func ) const
    {
      DomainType normal = inside.intersection().integrationOuterNormal( inside.localPosition() );

      JacobianRangeType jumpUNormal;
      for(int r = 0; r < dimRange; ++r)
      {
        for(int j=0; j<dimDomain; ++j)
          jumpUNormal[r][j] = normal[j] * (uLeft[r] - uRight[r]);
      }

      Fem::FieldMatrixConverter< GradientType, JacobianRangeType> func1( func );
      if (liftingMethod_ == LiftingEnum::id_id)
        func1 = jumpUNormal;
      else
      {
        JacobianRangeType AJumpUNormal;
        model_.diffusion( inside, uInflow, jumpUNormal, AJumpUNormal );
        func1 = AJumpUNormal;
      }
      func *= -inside.quadrature().weight( inside.index() );
    }

    /** \return A(u)L_e*n */
    template <class LocalEvaluation>
    void applyLifting(const LocalEvaluation& inside,
                      const DomainType& normal,
                      const RangeType& u,
                      const GradientType& sigma,
                      RangeType& lift) const
    {
      // convert sigma into JacobianRangeType
      Fem::FieldMatrixConverter< GradientType, JacobianRangeType> gradient( sigma );

      if (liftingMethod_ != LiftingEnum::id_A)
      {
        JacobianRangeType mat;
        // set mat = A(u)L_e
        model_.diffusion( inside, u, gradient, mat );
        // set lift = A(u)L_e*n
        mat.mv( normal, lift );
      }
      else
      {
        // just apply gradient
        gradient.mv( normal, lift );
      }
    }

    /** \brief calculate \f$\sum_{e\in\partial K} \Lambda_e |e|^2\f$
     *
     *  \note \f$\Lambda_e = 1\f$ for Dirichlet face @a e, \f$\Lambda_e = 0.5\f$
     *        for interface, and \f$\Lambda_e = 0\f$ for Neumann face
     */
    double sumFaceVolumeSqr( const EntityType& entity ) const
    {
      double sumFaceVolSqr  = 0.0;

      for (const auto& intersection : intersections(gridPart_, entity) )
      {
        const double faceVol = intersection.geometry().volume();

        // !!!!! forget about Neumann for now
        // 1/2 for interior intersections
        if ( intersection.neighbor() )
          sumFaceVolSqr += 0.5 * faceVol * faceVol;
        else
          sumFaceVolSqr += faceVol * faceVol;
      }

      return sumFaceVolSqr;
    }

  public:
    /**
     * \brief flux function on interfaces between cells
     *
     * \param intersection intersection
     * \param time current time given by TimeProvider
     * \param x coordinate of required evaluation local to \c intersection
     * \param uLeft DOF evaluation on this side of \c intersection
     * \param uRight DOF evaluation on the other side of \c intersection
     * \param gLeft result for this side of \c intersection
     * \param gRight result for the other side of \c intersection
     *
     * \note The total numerical (diffusive) flux for multiplication with phi
     *       is given with
     *        BO:
     *          gLeft = - {A(u)grad(u)}*n
     *        NIPG
     *          gLeft = - {A(u)grad(u)}*n
     *                  + C_nipg/h {A(u)}[u]*n
     *        CDG2:
     *          gLeft = - {A(u)grad(u)}*n
     *                  + C_cdg2/h {A(u)}[u]*n
     *                  + liftFactor*(A(u)L_e)|Ke-
     *        CDG:
     *          gLeft = - (A(u)grad(u)*n)|Ke-
     *                  + C_cdg/h {A(u)}[u]*n - beta*n[A(u)grad(u)]
     *                  + liftFactor*(A(u)L_e)|Ke-
     *        BR2:
     *          gLeft = - {A(u)grad(u)}*n + C_br2 {A(u)r_e([u])}*n
     *                = - {A(u)grad(u)}*n + C_br2 0.5*{A(u)L_e}
     *        IP:
     *          gLeft = - {A(u)grad(u)}*n + C_ip/h {A(u)}[u]*n
     *
     *       The total numerical flux for multiplication with grad(phi) is
     *        NIPG, BO:
     *          gDiffLeft = 0.5*A(u-)[u]
     *        IP, BR2, CDG2:
     *          gDiffLeft = -0.5*A(u-)[u]
     *        CDG:
     *          gDiffLeft = -0.5*A(u-)[u] - beta*[A(u)][u]
     *       where h = min(|entity+|,|entity-|) / |intersection|.
     *       In this method we need to return
     *          in gLeft:      = gRight
     *          in gDiffLeft:  = gDiffRight
     *                           0.5*A(u+)[u]

     *
     * \return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const JacobianRangeType& jacLeft,
                         const JacobianRangeType& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight )
    {
      gLeft  = 0;
      gRight = 0;

      const FaceDomainType& x = left.localPosition();
      const DomainType normal = left.intersection().integrationOuterNormal( x );

      const double faceLengthSqr = normal.two_norm2();

      JacobianRangeType diffmatrix;
      RangeType diffflux;

      // for all methods except CDG we need to evaluate {A(u)grad(u)}
      if (method_ != EnumType::cdg)
      {
        // A(u-)grad(u-) for multiplication with phi
        // call on inside
        model_.diffusion( left, uLeft, jacLeft, diffmatrix );

        // diffflux=A(u-)grad(u-)*n
        diffmatrix.mv( normal, diffflux );

        // A(u+)grad(u+) for multiplication with phi
        // call on outside
        model_.diffusion( right, uRight, jacRight, diffmatrix );

        // diffflux = 2{A(u)grad(u)}*n
        diffmatrix.umv( normal, diffflux );

        // ===== CDG2, BR2, IP, NIPG, BO =====
        // gLeft = gRight = -{A(u)grad(u)}*n
        gLeft.axpy ( -0.5, diffflux );
        gRight.axpy( -0.5, diffflux );
      }
      else
      // for CDG we need A(u)grad(u) on Ke-
      {
        ////////////////////////////////
        // [ phi ] * [ A(u^out)grad(u^out) ] ...
        ///////////////////////////////
        if ( insideIsInflow_ )
          model_.diffusion( right, uRight, jacRight, diffmatrix );
        else
          model_.diffusion( left, uLeft, jacLeft, diffmatrix );

        // diffflux=(A(u)grad(u)*n)|Ke-
        diffmatrix.mv( normal, diffflux );

        // ===== CDG =====
        // gLeft = gRight = -(A(u)grad(u)*n)|Ke-
        gLeft.axpy ( -1., diffflux );
        gRight.axpy( -1., diffflux );
      }

      // jumpUNormal = [u]*n
      RangeType jumpUNormal( uLeft );
      jumpUNormal -= uRight ;

      // set jumpU = [u]
      JacobianRangeType jumpU;
      for( int r = 0; r < dimRange; ++r )
      {
        jumpU[ r ]  = normal;
        jumpU[ r ] *= jumpUNormal[ r ];
      }

      // get A(u-)[u] in gDiffLeft
      // get A(u+)[u] in gDiffRight
      // this is not the final value for gDiffLeft, gDiffRight
      // A(u-)[u], A(u+)[u] are needed in the penalty term
      model_.diffusion( left,  uLeft,  jumpU, gDiffLeft );
      model_.diffusion( right, uRight, jumpU, gDiffRight );

      ////////////////////////////////////////////////
      //  start penalty term
      ///////////////////////////////////////////////
      // BR2 has its own special and BO doesn't have penalty term
      // every other method has this penalty term
      if ( penaltyTerm_ )
      {
        RangeType penaltyTerm ;

        if( (method_ == EnumType::ip) && useTheoryParams_ )
        {
          // penaltyTerm
          // = ainsworthFactor * maxEigenValue(A(u)) * FaceEntityVolumeRatio * [u] * n
          RangeType maxLeft, maxRight;
          model_.eigenValues( left,  uLeft, maxLeft );
          model_.eigenValues( right, uRight, maxRight );
          RangeFieldType maxEigenValue = 0;
          for( int r = 0; r<dimRange; ++r )
          {
            maxEigenValue = std::max( maxEigenValue, maxLeft[r] );
            maxEigenValue = std::max( maxEigenValue, maxRight[r] );
          }

          // calculate penalty factor
          const double penaltyFactorInside = sumFaceVolumeSqr( left.entity() ) / left.volume();
          const double penaltyFactorOutside = sumFaceVolumeSqr( right.entity() ) / right.volume();
          penalty_ = std::max( penaltyFactorInside, penaltyFactorOutside );
          penalty_ *= ainsworthFactor_ * maxEigenValue ;

          jumpU.mv( normal, penaltyTerm );
        }
        else
        {
          // \int C_IP {A(u)}[u][phi] dx
          //    = \int C_IP (A(u-)[u]*n + A(u+)[u]*n)/2 phi
          // so penaltyTerm = C_IP (A(u-)[u]*n + A(u+)[u]*n)/2

          // apply with normal
          gDiffLeft.mv( normal, penaltyTerm );
          gDiffRight.umv( normal, penaltyTerm );

          penaltyTerm *= 0.5 ;
        }

        double minvol =
          std::min( left.volume(), right.volume() );
        penaltyTerm /= minvol;

        // ===== IP, CDG2+, CDG+, NIPG+ =====
        // add to fluxes
        gLeft.axpy( penalty_, penaltyTerm );
        gRight.axpy( penalty_, penaltyTerm );
      }
      ////////////////////////////////////////////////
      //  end penalty term
      ///////////////////////////////////////////////


      // gDiffLeft = -0.5 A(u-)[u] (IP,BR2)
      // gDiffLeft = 0.5 A(u-)[u] (NIPG,BO)
      gDiffLeft *= nipgFactor_;

      // current entity gets (i.e. for IP)
      // {A(u)grad(u)}* n  -C11[u]*n)phi    + 0.5A(u-)[u]*grad(phi)
      // and neighbor gets
      // {A(u)grad(u)}*(-n)-C11[u]*(-n))phi + 0.5A(u-)[u]*grad(phi)
      // so term 0.5A(u-)[u]*grad(phi) stays the same therefore:

      // gDiffRight = 0.5 A(u-)[u] (IP,BR2)
      // gDiffRight = -0.5 A(u-)[u] (NIPG,BO)
      gDiffRight *= (-nipgFactor_);

      ////////////////////////////////////////////////
      //  begin lifting terms
      ///////////////////////////////////////////////
      if( hasLifting() )
      {
        // ... and the values of u from Ke-
        const RangeType& u = insideIsInflow_ ? uRight : uLeft;
        const LocalEvaluation& inflow = insideIsInflow_ ? right : left ;

        RangeType lift;
        // LiftingFunctionType& LeMinus = LeMinusLifting().function();
        // get value of A(u^out)L_e^out*n into liftTotal
        applyLifting( inflow, normal, u, liftingEvalLeMinus_[inflow.index()], lift );

        // only for CDG-type methods
        if (method_ != EnumType::br2)
        {
          lift   *= liftFactor_ ;

          // ===== CDG2, CDG =====
          gLeft  -= lift;
          gRight -= lift;
        }

        if( method_ == EnumType::cdg )
        {
          const RangeFieldType C_12 = insideIsInflow_ ? 0.5 : -0.5;
          JacobianRangeType resU;

          ////////////////////////////////
          // [ u ] * [ A(u)grad(phi) ]
          ///////////////////////////////


          resU = gDiffLeft;
          resU *= C_12/nipgFactor_;

          // save gDiffLeft
          gDiffLeft += resU;

          // MINUS MINUS MINUS MINUS
          resU  = gDiffRight;
          resU *= C_12/(-nipgFactor_);

          // save gDiffLeft
          gDiffRight += resU;
        }

        if (method_ == EnumType::br2)
        {
          // BR2 hasn't had penalty term until now
          // so we add it at this place.
          //             \int C_BR2 {A(u)r_e([u])}[phi] dx
          //           = \int C_BR2 (A(u_out)r_e([u])*n + A(u_in)r_e([u])*n)/2 phi
          // so penaltyTerm = C_BR2 (A(u_out)r_e([u])*n + A(u_in)r_e([u])*n)/2

          // get correct quadrature, the one from Ke+
          // ... and the values of u from Ke+
          const RangeType& uPlus = insideIsInflow_ ? uLeft : uRight;
          const LocalEvaluation& inflow = insideIsInflow_ ? left : right;

          RangeType liftTotal;

          // get value of A(u^in)L_e^in*n into liftTotal
          applyLifting( inflow, normal, uPlus, liftingEvalLePlus_[inflow.index()], liftTotal );

          // set liftTotal = {A(u)r_e}*n = 0.25*(A(u^in)L_e^in*n + A(u^out)L_e^out*n)
          liftTotal += lift;
          // add penalty coefficient
          liftTotal *= 0.25 * liftFactor_;

          // ===== BR2 =====
          gLeft -= liftTotal;
          gRight -= liftTotal;
        }
      }
      ////////////////////////////////////////////////
      //  end lifting terms
      ///////////////////////////////////////////////

      //////////////////////////////////////////////////////////
      //
      //  --Time step calculation
      //
      //////////////////////////////////////////////////////////
      const double faceVolumeEstimate = dimensionFactor_ *
        (left.intersection().conforming() ? faceLengthSqr
          : (nonconformingFactor_ * faceLengthSqr));

      const double diffTimeLeft = diffusionTimeStep( left, faceVolumeEstimate, uLeft );

      const double diffTimeRight = diffusionTimeStep( right, faceVolumeEstimate, uRight );

      // take minimum to proceed
      const double diffTimeStep = std::max( diffTimeLeft, diffTimeRight );

      // timestep restict to diffusion timestep
      // WARNING: reconsider this
      return diffTimeStep * cflDiffinv_;
    }


    /**
     * \brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        const RangeType& uLeft,
                        const RangeType& uRight,
                        const JacobianRangeType& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft )
    {
      // get local point
      const FaceDomainType& x = left.localPosition();
      const DomainType normal = left.intersection().integrationOuterNormal( x );

      const double faceLengthSqr = normal.two_norm2();

      JacobianRangeType diffmatrix;

      // diffmatrix = A(u-)grad(u-)
      model_.diffusion( left, uLeft, jacLeft, diffmatrix);

      // gLeft = -A(u-)grad(u-)*n
      diffmatrix.mv( normal, gLeft );
      gLeft *= -1.0;

      if( hasLifting() )
      {
        RangeType lift;
        // get value of A(u)L_e*n into lift
        applyLifting( left, normal, uRight, liftingEvalLeMinus_[left.index()], lift );

        if( method_ == EnumType::br2 )
        {
          // set liftTotal = A(u)r_e*n = 0.5*A(u_in)L_e_in*n
          lift *= (0.5*liftFactor_);
        }
        else
        {
          // only for CDG-type methods
          lift *= liftFactor_ ;
        }

        gLeft -= lift;
      }

      /****************************/
      /* Diffusion                 *
       ****************************/
      const double bndNipgFactor = 2.0 * nipgFactor_ ;

      // A(u)[u] = A(u)(u-g_D)*n (SIPG)
      // -A(u)[u] = -A(u)(u-g_D)*n (NIPG)
      // for multiplication with grad(phi)
      JacobianRangeType bndJumpU ;
      for( int r = 0; r < dimRange; ++r )
      {
        bndJumpU[r]  = normal;
        bndJumpU[r] *= (uLeft[r] - uRight[r]);
      }


      // get A(u-)[u] in gDiffLeft
      // this is not hte final value for gDiffLeft
      // but it's used in the penalty term
      model_.diffusion( left, uLeft, bndJumpU, gDiffLeft );

      // add penalty term
      if ( penaltyTerm_ )
      {
        // penalty term for IP
        RangeType penaltyTerm;
        const double enVolInv = 1./left.volume();

        if( (method_ == EnumType::ip) && useTheoryParams_ )
        {
          // penaltyTerm
          // = ainsworthFactor * maxEigenValue(A(u)) * FaceEntityVolumeRatio * [u] * n
          RangeType maxInside;
          model_.eigenValues( left, uLeft, maxInside );
          RangeFieldType maxEigenValue = 0;
          for( int r = 0; r<dimRange; ++r )
            maxEigenValue = std::max( maxEigenValue, maxInside[r] );

          // calculate penalty factor
          penalty_  = sumFaceVolumeSqr( left.entity() ) * enVolInv;
          penalty_ *= ainsworthFactor_ * maxEigenValue ;

          bndJumpU.mv( normal, penaltyTerm );
        }
        else
        {
          /*
          double penFac = model_.penaltyBoundary( left, uLeft );
          bndJumpU.mv( normal, penaltyTerm );
          penaltyTerm *= penFac ;
          */
          // \int C_IP {A(u)}[u][phi] dx
          //    = \int C_IP A(u_L)[u]*n phi
          // so penaltyTerm = C_IP A(u_L)[u]*n

          // apply with normal
          gDiffLeft.mv( normal, penaltyTerm );
        }

        // scale with 1/|e|, because normal is not unit normal
        penaltyTerm *= enVolInv;

        gLeft.axpy( penalty_, penaltyTerm );
      }

      // gDiffLeft = A(u-)[u]  (SIPG)
      // gDiffLeft = -A(u-)[u] (NIPG)
      gDiffLeft *= bndNipgFactor;

      ////////////////////////////////////////////////////
      //
      //  --Time step boundary
      //
      ////////////////////////////////////////////////////
      return diffusionTimeStep( left, faceLengthSqr, uLeft ) * cflDiffinv_;
    }

  protected:
    template <class LocalEvaluation>
    double diffusionTimeStep(const LocalEvaluation& local,
                             const double faceVolumeEstimate,
                             const RangeType& u) const
    {
      // apply proper scaling to obtain h^2 scaling of time step estimate
      return (faceVolumeEstimate / local.volume()) *
        model_.diffusionTimeStep( local, u ); // model only provides scaling with diffusion coefficient
    }

    GridPartType&                 gridPart_;
    const EnumType                method_;
    double                        penalty_;
    const double                  nipgFactor_;
    double                        liftFactor_;
    LiftingEnum                   liftingMethod_;
    const bool                    penaltyTerm_;
    DiscreteGradientSpaceType  gradSpc_;
    std::unique_ptr< Lifting > LeMinusLifting_;
    std::unique_ptr< Lifting > LePlusLifting_;
#ifdef LOCALDEBUG
    std::unique_ptr< Lifting > LeMinusLifting2_;
#endif
    mutable Fem::DynamicArray< GradientType > liftingEvalLeMinus_ ;
    mutable Fem::DynamicArray< GradientType > liftingEvalLePlus_ ;

    double            maxNeighborsVolumeRatio_; // for CDG2 only
    double            ainsworthFactor_;
    bool              insideIsInflow_;
    const bool        areaSwitch_;
    bool              useTheoryParams_;
    bool              initialized_;
  }; // end DGPrimalDiffusionFluxImpl



  //////////////////////////////////////////////////////////
  //
  //  extended flux for matrix assembly
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp,
            class Model,
            class FluxParametersImp = DGPrimalDiffusionFluxParameters >
  class ExtendedDGPrimalDiffusionFlux
   : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, FluxParametersImp >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model, FluxParametersImp >      BaseType;

  public:
    typedef typename BaseType::GridPartType        GridPartType;
    typedef typename BaseType::IntersectionType    IntersectionType;
    typedef typename BaseType::EntityType          EntityType;
    typedef typename BaseType::RangeType           RangeType;
    typedef typename BaseType::JacobianRangeType   JacobianRangeType;
    typedef typename BaseType::GradientType        GradientType;
    typedef typename BaseType::DomainType          DomainType;

    typedef typename BaseType::ParameterType       ParameterType;

    ExtendedDGPrimalDiffusionFlux( GridPartType& gridPart,
                                   const Model& model,
                                   const ParameterType& parameters )
      : BaseType( gridPart, model, parameters, BaseType::EnumType::primal )
    { }

    ExtendedDGPrimalDiffusionFlux( GridPartType& gridPart,
                                   const Model& model )
      : ExtendedDGPrimalDiffusionFlux( gridPart, model, ParameterType() )
    { }

    //! copy constructor (needed for thread parallel programs)
    ExtendedDGPrimalDiffusionFlux( const ExtendedDGPrimalDiffusionFlux& other ) :
      BaseType( other )
    {
    }

    using BaseType::initializeIntersection;

    template <class LocalEvaluationVec, class RangeVector >
    void initializeIntersection(const LocalEvaluationVec& left,
                                const LocalEvaluationVec& right,
                                const RangeVector& uLeftVec,
                                const RangeVector& uRightVec,
                                bool computeBoth)
    {
      this->computeLiftings(left, right, uLeftVec, uRightVec, computeBoth );
    }

    // return AL_e.n on element and neighbor
    template <class LocalEvaluation>
    void evaluateLifting(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         const RangeType& uEn,
                         const RangeType& uNb,
                         JacobianRangeType& liftEn,
                         JacobianRangeType& liftNb) const
    {
      assert( this->LePlusLifting().isInitialized() );
      assert( this->LePlusLifting().entity() == left.entity() );
      assert( this->LeMinusLifting().isInitialized() );
      assert( this->LeMinusLifting().entity() == right.entity() );

      const int qp = left.index();
      if ( this->insideIsInflow_)
      {
        applyLifting( left, left.intersection().normal(), uEn, liftingEvalLePlus_[qp], liftEn );
        applyLifting( right, right.intersection().normal(), uNb, liftingEvalLeMinus_[qp], liftNb );
      }
      else
      {
        applyLifting( right, left.intersection().normal(), uEn, liftingEvalLeMinus_[qp], liftNb );
        applyLifting( left, right.intersection().normal(), uNb, liftingEvalLePlus_[qp], liftEn );

      }
      assert( liftEn == liftEn );
      assert( liftNb == liftNb );
    }

    // return AL_e.n on element and neighbor
    const typename BaseType::LiftingFunctionType &getInsideLifting() const
    {
      assert( this->LePlusLifting().isInitialized() );
      assert( this->LeMinusLifting().isInitialized() );
      if ( this->insideIsInflow_ )
      {
        return this->LePlusLifting().function();
      }
      else
      {
        return this->LeMinusLifting().function();
      }
    }

  protected:
    using BaseType :: liftingEvalLeMinus_;
    using BaseType :: liftingEvalLePlus_;

    /*
    template <class QuadratureImp>
    void applyLifting(const EntityType& entity,
                      const double time,
                      const QuadratureImp& quad,
                      const int qp,
                      const DomainType& normal,
                      const RangeType& u,
                      const GradientType& sigma,
                      JacobianRangeType& mat) const
    {
      Fem::FieldMatrixConverter< GradientType, JacobianRangeType> gradient( sigma );
      // set mat = A(u)L_e
      this->model_.diffusion( entity, time, x, u, gradient, mat );
    }
    */
  }; // end ExtendedDGPrimalDiffusionFlux

} // end namespace
} // end namespace
#endif
