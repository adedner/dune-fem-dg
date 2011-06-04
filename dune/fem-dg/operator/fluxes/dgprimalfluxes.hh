#ifndef DUNE_FEM_DG_DGPRIMALFLUXES_HH
#define DUNE_FEM_DG_DGPRIMALFLUXES_HH

// Dune-Fem includes
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

// local includes
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

namespace Dune {


  // DGPrimalDiffusionFluxImpl
  //----------------------

  template <class DiscreteFunctionSpaceImp, 
            class Model>
  class DGPrimalDiffusionFluxImpl 
   : public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model > 
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model >      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    enum { dimRange  = DiscreteFunctionSpaceType :: dimRange };
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType :: polynomialOrder };

    typedef typename DiscreteFunctionSpaceType :: RangeFieldType        RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType       DomainFieldType;
    typedef FieldVector< DomainFieldType, dimDomain-1 >                 FaceDomainType;
    typedef typename DiscreteFunctionSpaceType :: DomainType            DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType             RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType     JacobianRangeType;


    typedef typename DiscreteFunctionSpaceType :: GridPartType          GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType           IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection               Intersection;
    typedef typename GridPartType :: GridType                           GridType;
    typedef typename GridType::template Codim< 0 >::Entity              EntityType;
    typedef typename GridPartType::template Codim< 0 >::IteratorType    IteratorType;
    typedef typename GridPartType::IntersectionIteratorType
                                                                        IntersectionIteratorType;

    typedef typename BaseType :: DiscreteGradientSpaceType              DiscreteGradientSpaceType;
    typedef typename DiscreteGradientSpaceType :: RangeType             GradientType;
    typedef TemporaryLocalFunction< DiscreteGradientSpaceType >         LiftingFunctionType;

    typedef CachingQuadrature< GridPartType, 0>                         VolumeQuadratureType ;
    
    typedef LocalDGMassMatrix
      < DiscreteGradientSpaceType, VolumeQuadratureType >               LocalMassMatrixType;

    class Lifting 
    {
      DiscreteGradientSpaceType gradSpc_;
      LiftingFunctionType r_e_;
      LocalMassMatrixType localMassMatrix_;

    public:  
      Lifting(GridPartType &gridPart) 
        : gradSpc_( gridPart )
        , r_e_( gradSpc_ )
        , localMassMatrix_( gradSpc_, 2*gradSpc_.order() )
      {}

      void initialize(const EntityType& entity ) 
      { 
        r_e_.init( entity ); 
        r_e_.clear();
      }

      LiftingFunctionType& function() 
      {
        return r_e_;
      }

      void finalize()
      {
        localMassMatrix_.applyInverse( r_e_ );
      }

    private:
      Lifting( const Lifting& );
    };

  private:
  protected:
    using BaseType :: determineDirection;
    using BaseType :: model_;
    using BaseType :: cflDiffinv_;
    using BaseType :: dimensionFactor_;
    using BaseType :: nonconformingFactor_;
    using BaseType :: numericalFlux ;
    using BaseType :: upwind_ ;
    using BaseType :: methodNames ;

    typedef typename BaseType :: MethodType  MethodType;

  public:

    bool initAreaSwitch() const 
    {
      if( method_ == method_cdg2 )
      {
        // when default value is used, then use areSwitch 
        if( ( upwind_ - BaseType :: upwindDefault() ).two_norm2() < 1e-10 ) 
          return true ;
      }
      return upwind_.two_norm2() < 1e-10 ;
    }
  public:
    enum { evaluateJacobian = true };

    /**
     * @brief constructor reading parameters 
     */
    DGPrimalDiffusionFluxImpl( GridPartType& gridPart,
                               const Model& model,
                               const MethodType method ) :
      BaseType( model, true ),
      gridPart_( gridPart ),
      method_( method ),
      penalty_( Parameter::getValue<double>("dgdiffusionflux.penalty") ),
      nipgFactor_( (method_ == method_nipg) || 
                   (method_ == method_bo) 
                   ? 0.5 : -0.5 ),
      liftFactor_( Parameter::getValue<double>("dgdiffusionflux.liftfactor") ),
      penaltyTerm_( method_ip || ((std::abs(  penalty_ ) > 0) && 
                    method_ != method_br2 && 
                    method_ != method_bo )),
      LeMinusLifting_( hasLifting() ? new Lifting( gridPart ) : 0 ),
      LePlusLifting_( ( method_ == method_br2 ) ? new Lifting( gridPart ) : 0 ),
#ifdef LOCALDEBUG
      LeMinusLifting2_( ( method_ <= method_cdg ) ? new Lifting( gridPart ) : 0 ),
#endif
      insideIsInflow_ ( true ),
      areaSwitch_( initAreaSwitch() ),
      useTheoryParams_( false ),
      initialized_ ( false )
    {
      
      // calculate maxNeighborVolumeRatio_
      maxNeighborsVolumeRatio_ = 1.;

      double theoryFactor = Parameter::getValue< double >( "dgdiffusionflux.theoryparameters", 0. );
      useTheoryParams_ = (theoryFactor > 0.);

      double n_k = DiscreteFunctionSpaceType :: polynomialOrder ; 
      ainsworthFactor_ = theoryFactor * 0.5 * n_k * ( n_k + 1.0 );

      int maxNumFaces = 0 ;
      int maxNumOutflowFaces = 0;
      if ( useTheoryParams_ )
      {
        const IteratorType itend = gridPart.template end<0>();
        for( IteratorType it = gridPart.template begin<0>(); it != itend; ++it ) 
        {
          const EntityType& entity = * it ;
          const double insideVol = entity.geometry().volume();
          int numFaces = 0;
          int numOutflowFaces = 0;
          const IntersectionIteratorType intitend = gridPart.iend( entity );
          for(IntersectionIteratorType intit = gridPart.ibegin( entity ); 
              intit != intitend; ++intit )
          {
            const Intersection& intersection = * intit ;

            ++numFaces ;
            if ( intersection.neighbor() )
            {
              numOutflowFaces += (determineDirection(intersection) ? 1 : 0);
              double outsideVol = intersection.outside()->geometry().volume();
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
        if (method_ == method_cdg2)
        {
          liftFactor_ = theoryFactor * 0.25* ((double) maxNumFaces); // max number of faces here
          //if( ! areaSwitch_ ) 
          liftFactor_ *= (1.+maxNeighborsVolumeRatio_);
        }
        else if (method_ == method_cdg) 
        { 
          liftFactor_ = theoryFactor * maxNumOutflowFaces; 
        }
        else if (method_ == method_br2) 
        { 
          liftFactor_ = theoryFactor * maxNumFaces;
        }
        else 
        {
          std::cerr << "ERROR: No theory parameters for NIPG" << std::endl;
          DUNE_THROW(InvalidStateException,"No theory parameters for NIPG");
        }
      }

      if (Parameter :: verbose()) 
      {
        std::cout <<"Diff. flux: ";
        diffusionFluxName( std::cout );

        std::cout <<", penalty: ";
        if ( useTheoryParams_ && (method_ == method_ip) )
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
      penaltyTerm_( other.penaltyTerm_ ),
      LeMinusLifting_( hasLifting() ? new Lifting( gridPart_ ) : 0 ),
      LePlusLifting_( ( method_ == method_br2 ) ? new Lifting( gridPart_ ) : 0 ),
#ifdef LOCALDEBUG
      LeMinusLifting2_( ( method_ <= method_cdg ) ? new Lifting( gridPart_ ) : 0 ),
#endif
      maxNeighborsVolumeRatio_( other.maxNeighborsVolumeRatio_ ),
      ainsworthFactor_( other.ainsworthFactor_ ),
      insideIsInflow_ ( other.insideIsInflow_ ),
      areaSwitch_( other.areaSwitch_ ), // used area based switch
      useTheoryParams_( other.useTheoryParams_ ),
      initialized_( other.initialized_ )
    {
    }

    //! destructor deleting liftings 
    ~DGPrimalDiffusionFluxImpl() 
    {
      delete LeMinusLifting_; 
      LeMinusLifting_ = 0;
      delete LePlusLifting_; 
      LePlusLifting_ = 0;
#ifdef LOCALDEBUG
      delete LeMinusLifting2_; 
      LeMinusLifting2_ = 0;
#endif
    }

    double maxNeighborsVolumeRatio() const
    {
      assert( initialized_ );
      return maxNeighborsVolumeRatio_;
    }

    void diffusionFluxName( std::ostream& out ) const
    {
      out << methodNames( method_ );
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
    const bool hasLifting () const { return ( method_ <= method_br2 ); }

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
    template <class QuadratureImp, class ArgumentTupleVector > 
    void initializeIntersection(const Intersection& intersection,
                                const EntityType& inside,
                                const EntityType& outside,
                                const double time,
                                const QuadratureImp& quadInner, 
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec) 
    {
      if( hasLifting() ) 
      {
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
        if (areaSwitch_)
          insideIsInflow_ = determineDirection( inside.geometry().volume(),
                                             outside.geometry().volume() ); 
        else
          insideIsInflow_ = determineDirection( intersection ); 

        const EntityType& entity = ( insideIsInflow_ ) ? outside : inside;

        const size_t quadNoInp = quadInner.nop();
        liftTmp_.resize( quadNoInp );

        // get the right quadrature for the lifting entity
        const QuadratureImp& faceQuad = ( insideIsInflow_ ) ? quadOuter : quadInner; 

#ifdef LOCALDEBUG
        const size_t quadNoOutp = quadOuter.nop();
        double sum = 0;
        double sum2 = 0;

        // get Ke+ in entity2
        // calculate L_e=2*r_e on Ke+
        const EntityType& entity2 = ( insideIsInflow_ ) ? inside : outside;
        LeMinusLifting().initialize( entity2 );

        {
          // get the right quadrature for the lifting entity
          const QuadratureImp& faceQuad2 = ( insideIsInflow_ ) ? quadInner : quadOuter; 

          for(size_t qp = 0; qp < quadNoInp; ++qp )
          {
            // get value of 2*r_e in quadrature point
            addLifting(intersection, time, faceQuad2, qp, 
                       uLeftVec[ qp ], uRightVec[ qp ], liftTmp_[ qp ] );
          }
          // add to local function 
          LeMinusLifting().function().axpyQuadrature( faceQuad, liftTmp_ );

          LeMinusLifting().finalize();
        }

        // calculate 4*\int_Ke+(r_e*r_e)
        double term1;
        term1 = integrateLifting( LeMinusLifting().function(), 
                                  LeMinusLifting().function().entity().geometry() );
        const double interiorFactor = 
          ( insideIsInflow_ ? -1. : 1. );

        // set sum += -4*\int_Ke+(r_e*r_e)
        sum += interiorFactor*term1;

        // calculate L_e=2*r_e on Ke-
        LeMinusLifting2().initialize( entity );

        for(size_t qp = 0; qp < quadNoOutp; ++qp )
        {
          // get value of 2*r_e in quadrature point
          addLifting(intersection, time, faceQuad, qp, 
                     uLeftVec[ qp ], uRightVec[ qp ], liftTmp_ );
        }
        LeMinusLifting2().function().axpyQuadrature( faceQuad, liftTmp_ );
        LeMinusLifting2().finalize();

        // calculate 4*\int_Ke-(r_e*r_e)
        double term2;
        term2 = integrateLifting( LeMinusLifting2().function(),
                                  LeMinusLifting2().function().entity().geometry() );

        sum2 = sum;

        // put 4*(\int_Ke-(r_e*r_e) - \int_Ke+(r_e*r_e))
        //  = 4*\int_e(r_e*l_e) in sum
        // put 4*(\int_Ke-(r_e*r_e) + \int_Ke+(r_e*r_e))
        //  = 4*\int_e(r_e*r_e) in sum2
        sum -= interiorFactor*term2;
        sum2 += std::abs( term2 );

        if (sim2 > 1e-20)
        {
          localMaxRatio_ = std::max( localMaxRatio_, std::abs(sum/sum2) );
          localMinRatio_ = std::min( localMinRatio_, std::abs(sum/sum2) );
        } else
          localMaxRatio_ = std::numeric_limits<double>();

        if ( (std::abs(sum) > sum2) )
        {
          std::cout <<"sum: " <<sum <<" sum2: " <<sum2 <<std::endl;
          std::cout <<"term1: " <<term1 <<std::endl;
          std::cout <<"term2: " <<term2 <<std::endl <<std::endl;
          abort();
        }

        // add to final sum 
        //    \sum_e (\int_Ke+ r_e*r_e - \int_Ke- r_e*r_e)
        //    \sum_e (\int_Ke+ r_e*r_e + \int_Ke- r_e*r_e)
        sum_ += sum;
        sum2_ += sum2;
#endif

        // back to the real computation, entity=Ke-
        LeMinusLifting().initialize( entity );

        // calculate real lifting
        for(size_t qp = 0; qp < quadNoInp; ++qp )
        {
          addLifting(intersection, time, faceQuad,  qp, 
                     uLeftVec[ qp ], uRightVec[ qp ], liftTmp_[ qp ] );
        }
        // add to local function 
        LeMinusLifting().function().axpyQuadrature( faceQuad, liftTmp_ );

        // LeMinusLifting_ has L_e=2*r_e on Ke-
        LeMinusLifting().finalize();

        if (method_ == method_br2)
        {
          // get Ke+ in entity2
          // calculate 2*r_e on Ke+
          const EntityType& entity2 = ( insideIsInflow_ ) ? inside : outside;
          LePlusLifting().initialize( entity2 );

          // get the right quadrature for the lifting entity
          const QuadratureImp& faceQuad2 = ( insideIsInflow_ ) ? quadInner : quadOuter; 

          const size_t quadNoOutp = quadOuter.nop();
          for(size_t qp = 0; qp < quadNoOutp; ++qp )
          {
            // get value of 2*r_e in quadrature point
            // use correct order on interface quadratures!
            addLifting(intersection, time, faceQuad2,  qp, 
                       uLeftVec[ qp ], uRightVec[ qp ], liftTmp_[ qp ] );
          }

          // add to local function 
          LePlusLifting().function().axpyQuadrature( faceQuad2, liftTmp_ );

          // LePlusLifting_ carries 2*r_e on Ke+
          LePlusLifting().finalize();
        }
      }
    }


#ifdef LOCALDEBUG
    template <class LiftingFunction , class Geometry >
    double integrateLifting( const LiftingFunction& lifting, const Geometry& geometry ) const 
    {
      typedef typename LiftingFunction :: RangeType RangeType;
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


    template <class QuadratureImp, class ArgumentTupleVector>
    void initializeBoundary(const Intersection& intersection,
                            const EntityType& entity,
                            const double time,
                            const QuadratureImp& quadInner, 
                            const ArgumentTupleVector& uLeftVec,
                            const std::vector< RangeType >& uRight) 
    {
      if( hasLifting() ) 
      {
        insideIsInflow_ = true; 

        LeMinusLifting().initialize( entity );

        const size_t quadNop = quadInner.nop();
        liftTmp_.resize( quadNop );
        for(size_t qp = 0; qp < quadNop; ++qp )
        {
          addLifting(intersection, time, quadInner, qp, 
                     uLeftVec[ qp ], uRight[ qp ] , liftTmp_[ qp ] );
        }
        // add to local function 
        LeMinusLifting().function().axpyQuadrature( quadInner, liftTmp_ );
        LeMinusLifting().finalize();
      }
    }

  protected:  
    template <class QuadratureImp, class LiftingFunction >
    void addLifting(const Intersection& intersection,
                    const double time,
                    const QuadratureImp& faceQuad,
                    const int quadPoint, 
                    const RangeType& uLeft,
                    const RangeType& uRight,
                    LiftingFunction& func ) const
    {
      const FaceDomainType& x = faceQuad.localPoint( quadPoint );
      DomainType normal = intersection.integrationOuterNormal( x );

      // calculate jump 
      RangeType jumpUNormal ( uLeft );
      jumpUNormal -= uRight;

      addLifting( normal, jumpUNormal, faceQuad, 
                  quadPoint, func );
    }

    template <class QuadratureImp>
    void addLifting( const DomainType& normal,
                     const RangeType& jumpUNormal,
                     const QuadratureImp& faceQuad,
                     const int quadPoint,
                     GradientType& factor) const 
    {
      FieldMatrixConverter< GradientType, JacobianRangeType> factor1( factor );
      for(int r = 0; r < dimRange; ++r)
      {
        for(int j=0; j<dimDomain; ++j) 
          factor1[ r ][ j ] = normal[ j ] * jumpUNormal[ r ];
      }
      factor *= -faceQuad.weight( quadPoint );
    }

    /** @return A(u)L_e*n */
    template <class QuadratureImp>
    void applyLifting(const QuadratureImp& faceQuad,
                      const int quadPoint,
                      const double time,
                      const DomainType& normal,
                      const RangeType& u,
                      LiftingFunctionType& r_e,
                      RangeType& lift) const
    { 
      GradientType sigma;

      // evaluate lifting at quadrature point 
      r_e.evaluate( faceQuad[ quadPoint ], sigma );

      FieldMatrixConverter< GradientType, JacobianRangeType> gradient( sigma );

      JacobianRangeType mat; 

      // set mat = G(u)L_e
      model_.diffusion( r_e.entity(),
                        time, faceQuad.point( quadPoint ),
                        u, gradient, mat );

      // set lift = G(u)L_e*n
      mat.mv( normal, lift );
    }

    /** \brief calculate \f$\sum_{e\in\partial K} \Lambda_e |e|^2\f$
     *  
     *  \note \f$\Lambda_e = 1\f$ for Dirichlet face @a e, \f$\Lambda_e = 0.5\f$
     *        for interface, and \f$\Lambda_e = 0\f$ for Neumann face
     */
    double sumFaceVolumeSqr( const EntityType& entity ) const
    {
      double sumFaceVolSqr  = 0.0;

      const IntersectionIteratorType intitend = gridPart_.iend( entity );
      for(IntersectionIteratorType intit = gridPart_.ibegin( entity ); 
          intit != intitend; ++intit )
      {
        const Intersection& intersection = * intit ;
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
     * @brief flux function on interfaces between cells
     *
     * @param intersection intersection
     * @param time current time given by TimeProvider
     * @param x coordinate of required evaluation local to \c intersection
     * @param uLeft DOF evaluation on this side of \c intersection
     * @param uRight DOF evaluation on the other side of \c intersection
     * @param gLeft result for this side of \c intersection
     * @param gRight result for the other side of \c intersection
     *
     * @note The total numerical flux for multiplication with phi
     *       is given with
     *        CDG2:
     *          gLeft = numflux(f(u)) - {G(u)grad(u)}*n 
     *            + C_cdg2/h {G(u)}[u]*n
     *            + liftFactor*(G(u)L_e)|Ke-
     *        CDG:
     *          gLeft = numflux(f(u)) - {G(u)grad(u)}*n 
     *            - beta*n[G(u)grad(u)] + C_cdg/h {G(u)}[u]*n
     *            + liftFactor*(G(u)L_e)|Ke-
     *        BR2:
     *          gLeft = numflux(f(u)) - {G(u)grad(u)}*n + C_br2 {G(u)r_e([u])}*n
     *        IP:
     *          gLeft = numflux(f(u)) - {G(u)grad(u)}*n + C_ip/h {G(u)}[u]*n
     *
     *       The total numerical flux for multiplication with grad(phi) is
     *        NIPG, BO:
     *          gDiffLeft = 0.5*G(u-)[u]
     *        IP, BR2, CDG2:
     *          gDiffLeft = -0.5*G(u-)[u]
     *        CDG:
     *          gDiffLeft = -0.5*G(u-)[u] - beta*G(u)[u]
     *       where h = min(|entity+|,|entity-|) / |intersection|.
     *       In this method we need to return
     *          in gLeft:      same like above, but without numflux(f(u))
     *          in gDiffLeft:  same like above
     *
     * @return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class DiscreteModelImp, class QuadratureImp>
    double numericalFlux(const Intersection& intersection,
                         const DiscreteModelImp& discreteModel,
                         const double time,
                         const QuadratureImp& faceQuadInner,
                         const QuadratureImp& faceQuadOuter,
                         const int quadPoint, 
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

      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
      const DomainType normal = intersection.integrationOuterNormal( x );

      const DomainType& xglInside  = faceQuadInner.point( quadPoint );
      const DomainType& xglOutside = faceQuadOuter.point( quadPoint );

      const EntityType& inside  = discreteModel.inside();
      const EntityType& outside = discreteModel.outside();

      const double faceLengthSqr = normal.two_norm2();

      JacobianRangeType diffmatrix ;
      RangeType diffflux ;

      // for all methods except CDG we need to evaluate {G(u)grad(u)}
      if (method_ != method_cdg) 
      {
        // G(u-)grad(u-) for multiplication with phi
        // call on inside 
        model_.diffusion( inside, time, xglInside, 
                          uLeft, jacLeft, diffmatrix );

        // diffflux=G(u-)grad(u-)*n
        diffmatrix.mv( normal, diffflux );

        // G(u+)grad(u+) for multiplication with phi
        // call on outside 
        model_.diffusion( outside, time, xglOutside,
                          uRight, jacRight, diffmatrix );

        // diffflux = 2{G(u)grad(u)}*n
        diffmatrix.umv( normal, diffflux );

        // gLeft = gRight = -{G(u)grad(u)}*n
        gLeft.axpy ( -0.5, diffflux );
        gRight.axpy( -0.5, diffflux );
      } 
      else
      // for CDG we need G(u)grad(u) on Ke-
      {
        ////////////////////////////////
        // [ phi ] * [ G(u)grad(u) ] ...
        ///////////////////////////////
        if ( ! insideIsInflow_ )
          model_.diffusion( inside, time, xglInside,
                            uLeft, jacLeft, diffmatrix );
        else
          model_.diffusion( outside, time, xglOutside,
                            uRight, jacRight, diffmatrix );
                            
        // diffflux=(G(u)grad(u)*n)|Ke-
        diffmatrix.mv( normal, diffflux );

        // gLeft = gRight = -(G(u)grad(u)*n)|Ke-
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

      // get G(u-)[u] in gDiffLeft
      // get G(u+)[u] in gDiffRight
      // this is not the final value for gDiffLeft, gDiffRight
      // G(u-)[u], G(u+)[u] are needed in the penalty term
      model_.diffusion( inside, time, xglInside,
                        uLeft, jumpU, gDiffLeft );
      model_.diffusion( outside, time, xglOutside,
                        uRight, jumpU, gDiffRight );

      ////////////////////////////////////////////////
      //  start penalty term 
      ///////////////////////////////////////////////
      // BR2 has its own special and BO doesn't have penalty term
      // every other method has this penalty term
      if ( penaltyTerm_ ) 
      {
        RangeType penaltyTerm ;

        if( (method_ == method_ip) && useTheoryParams_ ) 
        {
          // penaltyTerm 
          // = ainsworthFactor * maxEigenValue(A(u)) * FaceEntityVolumeRatio * [u] * n
          RangeType maxLeft, maxRight; 
          model_.eigenValues( inside, time, xglInside, uLeft, maxLeft );   
          model_.eigenValues( outside, time, xglOutside, uRight, maxRight );   
          double maxEigenValue = 0.;
          for( int r = 0; r<dimRange; ++r )
          {
            maxEigenValue = std::max( maxEigenValue, maxLeft[r] );
            maxEigenValue = std::max( maxEigenValue, maxRight[r] );
          }

          // calculate penalty factor
          const double penaltyFactorInside = sumFaceVolumeSqr( inside ) / discreteModel.enVolume();
          const double penaltyFactorOutside = sumFaceVolumeSqr( outside ) / discreteModel.nbVolume();
          penalty_ = std::max( penaltyFactorInside, penaltyFactorOutside ); 
          penalty_ *= ainsworthFactor_ * maxEigenValue ;

          jumpU.mv( normal, penaltyTerm );
        }
        else
        {
          // \int C_IP {G(u)}[u][phi] dx
          //    = \int C_IP (G(u_L)[u]*n + G(u_R)[u]*n)/2 phi
          // so penaltyTerm = C_IP (G(u_L)[u]*n + G(u_R)[u]*n)/2
          
          // apply with normal
          gDiffLeft.mv( normal, penaltyTerm );
          gDiffRight.umv( normal, penaltyTerm );
    
          penaltyTerm *= 0.5 ;

          /*
          double penFac = model_.penaltyFactor( inside,    outside, time, 
                                                xglInside, uLeft, uRight );
          penaltyTerm = jumpUNormal ;
          penaltyTerm *= penFac ;
          */
        }

        double minvol = 
          std::min( discreteModel.enVolume(), discreteModel.nbVolume() );
        penaltyTerm /= minvol;

        // add to fluxes 
        gLeft.axpy( penalty_, penaltyTerm );
        gRight.axpy( penalty_, penaltyTerm );
      }
      ////////////////////////////////////////////////
      //  end penalty term 
      ///////////////////////////////////////////////


      // gDiffLeft = -0.5 G(u-)[u] (IP,BR2)
      // gDiffLeft = 0.5 G(u-)[u] (NIPG,BO)
      gDiffLeft *= nipgFactor_;

      // current entity gets (i.e. for IP)
      //    (numflux(f(u))*n+{G(u)grad(u)}*n-C11[u]*n)phi
      //    + 0.5G(u-)[u]*grad(phi)
      // and neighbor gets
      //    (numflux(f(u))*(-n)+{G(u)grad(u)}*(-n)-C11[u]*(-n))phi
      //    + 0.5G(u-)[u]*grad(phi)
      // so term 0.5G(u-)[u]*grad(phi) stays the same therefore:
      gDiffRight *= (-nipgFactor_);

      ////////////////////////////////////////////////
      //  begin lifting terms 
      ///////////////////////////////////////////////
      if( hasLifting() )
      {
        // get correct quadrature, the one from Ke-
        const QuadratureImp& faceQuad = ( insideIsInflow_ ) ? faceQuadOuter : faceQuadInner; 
        // ... and the values of u from Ke-
        const RangeType& u = ( insideIsInflow_ ) ? uRight : uLeft;

        RangeType lift; 
        LiftingFunctionType& LeLeft = LeMinusLifting().function();
        // get value of G(u-)L_e-*n into liftTotal
        applyLifting( faceQuad, quadPoint, time, normal, u, LeLeft, lift );

        // only for CDG-type methods
        if (method_ != method_br2)
        {
          lift   *= liftFactor_ ;
          gLeft  -= lift;
          gRight -= lift;
        }

        if( method_ == method_cdg ) 
        {
          const RangeFieldType C_12 = ( insideIsInflow_ ) ? 0.5 : -0.5;
          JacobianRangeType resU;

          ////////////////////////////////
          // [ u ] * [ G(u)grad(phi) ] ...
          ///////////////////////////////

          resU = gDiffLeft; 
          resU *= C_12/nipgFactor_;

          // save gDiffLeft 
          gDiffLeft += resU;
          
          resU  = gDiffRight;
          resU *= C_12/(-nipgFactor_);

          // save gDiffLeft 
          gDiffRight += resU;
        }

        if (method_ == method_br2)
        {
          // BR2 hasn't had penalty term until now
          // so we add it at this place.
          // \int C_BR2 {G(u)r_e([u])}[phi] dx
          //    = \int C_BR2 (G(u_L)r_e([u])*n + G(u_R)r_e([u])*n)/2 phi
          // so penaltyTerm = C_BR2 (G(u_L)r_e([u])*n + G(u_R)r_e([u])*n)/2

          // get correct quadrature, the one from Ke+
          const QuadratureImp& faceQuadPlus = 
            ( ! insideIsInflow_ ) ? faceQuadOuter : faceQuadInner; 
          // ... and the values of u from Ke+
          const RangeType& uPlus = ( ! insideIsInflow_ ) ? uRight : uLeft;

          RangeType liftTotal;
          LiftingFunctionType& LeRight = LePlusLifting().function();

          // get value of G(u+)L_e+*n into liftTotal
          applyLifting( faceQuadPlus, quadPoint, time, 
                        normal, uPlus, LeRight, liftTotal );
          
          // set liftTotal = {G(u)r_e}*n = 0.25*(G(u+)L_e+*n + G(u-)L_e-*n)
          liftTotal += lift;
          // add penalty coefficient
          liftTotal *= 0.25 * liftFactor_; 

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
        (intersection.conforming() ? faceLengthSqr 
          : (nonconformingFactor_ * faceLengthSqr));

      const double diffTimeLeft = 
        model_.diffusionTimeStep( intersection,
            discreteModel.enVolume(),
            faceVolumeEstimate, 
            time, x, uLeft );
                                  
      const double diffTimeRight = 
        model_.diffusionTimeStep( intersection,
            discreteModel.nbVolume(),
            faceVolumeEstimate, 
            time, x, uRight );
                                  
      // take minimum to proceed 
      const double diffTimeStep = std::max( diffTimeLeft, diffTimeRight );

      // timestep restict to diffusion timestep
      // WARNING: reconsider this
      return diffTimeStep * cflDiffinv_;
    }


    /**
     * @brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class DiscreteModelImp, class QuadratureImp> 
    double boundaryFlux(const Intersection& intersection,
                        const DiscreteModelImp& discreteModel,
                        const double time, 
                        const QuadratureImp& faceQuadInner,
                        const int quadPoint,
                        const RangeType& uLeft,
                        const RangeType& uRight,
                        const JacobianRangeType& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft )   /*@LST0E@*/
    {
      // get local point 
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
      const DomainType normal = intersection.integrationOuterNormal(x);
      const DomainType& xglInside = faceQuadInner.point( quadPoint );

      const EntityType& inside  = discreteModel.inside();

      const double faceLengthSqr = normal.two_norm2();

      JacobianRangeType diffmatrix;

      // diffmatrix = G(u-)grad(u-)
      model_.diffusion( inside, time,
                        xglInside, 
                        uLeft, jacLeft, diffmatrix);

      // gLeft = -G(u-)grad(u-)*n
      diffmatrix.mv( normal, gLeft );
      gLeft *= -1.0;

      if( hasLifting() )
      {
        LiftingFunctionType& LeLeft = LeMinusLifting().function();

        RangeType lift; 
        // get value of G(u)L_e*n into lift
        applyLifting( faceQuadInner, quadPoint, time, normal, uRight, LeLeft, lift );

        if( method_ == method_br2 )
        {
          // set liftTotal = G(u)r_e*n = 0.5*G(u_in)L_e_in*n
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

      // G(u)[u] = G(u)(u-g_D)*n (SIPG)
      // -G(u)[u] = -G(u)(u-g_D)*n (NIPG)
      // for multiplication with grad(phi)
      JacobianRangeType bndJumpU ;
      for( int r = 0; r < dimRange; ++r )
      {
        bndJumpU[r]  = normal;
        bndJumpU[r] *= (uLeft[r] - uRight[r]);
      }


      // get G(u-)[u] in gDiffLeft
      // this is not hte final value for gDiffLeft
      // but it's used in the penalty term
      model_.diffusion( inside, 
                        time,
                        xglInside, 
                        uLeft, 
                        bndJumpU, gDiffLeft );

      // add penalty term
      if ( penaltyTerm_ ) 
      {
        // penalty term for IP
        RangeType penaltyTerm;
        const double enVolInv = 1./discreteModel.enVolume();

        if( (method_ == method_ip) && useTheoryParams_ ) 
        {
          // penaltyTerm 
          // = ainsworthFactor * maxEigenValue(A(u)) * FaceEntityVolumeRatio * [u] * n
          RangeType maxInside;
          model_.eigenValues( inside, time, xglInside, uLeft, maxInside );   
          double maxEigenValue = 0.;
          for( int r = 0; r<dimRange; ++r )
            maxEigenValue = std::max( maxEigenValue, maxInside[r] );

          // calculate penalty factor
          penalty_  = sumFaceVolumeSqr( inside ) * enVolInv;
          penalty_ *= ainsworthFactor_ * maxEigenValue ;

          bndJumpU.mv( normal, penaltyTerm );
        }
        else
        {
          /*
          double penFac = model_.penaltyBoundary( inside, time, 
                                                  xglInside, uLeft );
          bndJumpU.mv( normal, penaltyTerm );
          penaltyTerm *= penFac ;
          */
          // \int C_IP {G(u)}[u][phi] dx
          //    = \int C_IP G(u_L)[u]*n phi
          // so penaltyTerm = C_IP G(u_L)[u]*n
          
          // apply with normal
          gDiffLeft.mv( normal, penaltyTerm );
        }

        // scale with 1/|e|, because normal is not unit normal
        penaltyTerm *= enVolInv;

        gLeft.axpy( penalty_, penaltyTerm );
      }

      // gDiffLeft = G(u-)[u]  (SIPG)
      // gDiffLeft = -G(u-)[u] (NIPG)
      gDiffLeft *= bndNipgFactor;

      ////////////////////////////////////////////////////
      //
      //  --Time step boundary 
      //
      ////////////////////////////////////////////////////
      const double diffTimeStep = 
        model_.diffusionTimeStep( intersection,
            discreteModel.enVolume(),
            faceLengthSqr, 
            time, x, uLeft );
                                  
      return diffTimeStep * cflDiffinv_; 
    }

  protected:
    GridPartType&     gridPart_;
    const MethodType  method_; 
    double            penalty_;
    const double      nipgFactor_;
    double            liftFactor_;
    const bool        penaltyTerm_; 
    Lifting*          LeMinusLifting_;
    Lifting*          LePlusLifting_;
#ifdef LOCALDEBUG
    Lifting*          LeMinusLifting2_;
#endif
    mutable MutableArray< GradientType > liftTmp_ ;
    double            maxNeighborsVolumeRatio_; // for CDG2 only
    double            ainsworthFactor_;
    bool              insideIsInflow_; 
    const bool        areaSwitch_;
    bool              useTheoryParams_;
    bool              initialized_;
  }; // end DGPrimalDiffusionFluxImpl                        

  
  //! DG primal diffusion flux 
  template <class DiscreteFunctionSpaceImp, 
            class Model, 
            DGDiffusionFluxIdentifier >
  class DGPrimalDiffusionFlux;

  //////////////////////////////////////////////////////////
  //
  //  general diffusion flux allows choice of method via Parameter
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp, 
            class Model> 
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, method_general >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  protected:  
    typedef typename BaseType :: MethodType  MethodType;
    using BaseType :: getMethod ;

  public:  
    /**
      * @brief constructor reading parameters 
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model)
      : BaseType( gridPart, model, getMethod() )
    {
    }
  };

  //////////////////////////////////////////////////////////
  //
  //  specialization for CDG2 
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp, 
            class Model> 
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, method_cdg2 >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  public:  
    /**
      * @brief constructor reading parameters 
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model)
      : BaseType( gridPart, model, method_cdg2 )
    {
    }
  };


  //////////////////////////////////////////////////////////
  //
  //  specialization for CDG 
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp, 
            class Model> 
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, method_cdg >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  public:  
    /**
      * @brief constructor reading parameters 
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model)
      : BaseType( gridPart, model, method_cdg )
    {
    }
  };


  //////////////////////////////////////////////////////////
  //
  //  specialization for BR2 
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp, 
            class Model> 
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, method_br2 >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  public:  
    /**
      * @brief constructor reading parameters 
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model)
      : BaseType( gridPart, model, method_br2 )
    {
    }
  };


  //////////////////////////////////////////////////////////
  //
  //  specialization for IP 
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp, 
            class Model> 
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, method_ip >
    : public DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGPrimalDiffusionFluxImpl< DiscreteFunctionSpaceImp, Model >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  public:  
    /**
      * @brief constructor reading parameters 
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model)
      : BaseType( gridPart, model, method_ip )
    {
    }
  };

  //////////////////////////////////////////////////////////
  //
  //  specialization for no-diffusion
  //
  //////////////////////////////////////////////////////////
  template <class DiscreteFunctionSpaceImp, 
            class Model> 
  class DGPrimalDiffusionFlux<  DiscreteFunctionSpaceImp, Model, method_none >
    : public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model >
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model >
      BaseType;

  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  public:  
    /**
      * @brief constructor reading parameters 
      */
    DGPrimalDiffusionFlux( GridPartType& gridPart,
                           const Model& model)
      : BaseType( model, false )
    {
    }
  };

  //////////////////////////////////////////////////////////
  //
  //  specialization for NIPG and BO are missing since these methods are not so
  //  interesting, use method_general for this 
  //
  //////////////////////////////////////////////////////////

} // end namespace 
#endif
