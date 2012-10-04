#ifndef PRIMALMATRIXASSEMBLY_HH
#define PRIMALMATRIXASSEMBLY_HH

#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

#include <dune/fem-dg/operator/fluxes/dgprimalfluxes.hh>
#include <dune/fem-dg/operator/dg/primaloperator.hh>

namespace Dune {

template <class ModelType>
class UpwindFlux {
public:
  typedef ModelType Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model :: DomainType DomainType;
  typedef typename Model :: RangeType RangeType;
  typedef typename Model :: FluxRangeType FluxRangeType;
  // typedef typename Model :: DiffusionRangeType DiffusionRangeType;
  typedef typename Model :: FaceDomainType  FaceDomainType;
  typedef typename Model :: EntityType  EntityType;
  typedef typename Model :: IntersectionType  IntersectionType;
protected:
  template <class Model, bool constVelo>
    struct Velocity
    {
      /**
       * @brief computes and returns the wind direction
       */
      static inline double upwind(const Model& model,
                                  const IntersectionType& it,
                                  const double time,
                                  const FaceDomainType& x,
                                  const RangeType& uLeft)
      {
        const DomainType normal = it.integrationOuterNormal(x);
        DomainType velocity;
        model.velocity(*it.inside(),time,
                       it.geometryInInside().global(x),
                       velocity);
        return normal*velocity;
      }
    };

  template <class Model>
    struct Velocity<Model,true>
    {
      /**
       * @brief computes and returns the wind direction for models with
       * constant velocity
       */
      static inline double upwind( const Model& model,
                                   const IntersectionType& it,
                                   const double time,
                                   const FaceDomainType& x,
                                   const RangeType& uLeft )
      {
        const DomainType normal = it.integrationOuterNormal(x);
        return normal * model.velocity_;
      }
    };

public:
  /**
   * @brief constructor
   */
  UpwindFlux(const Model& mod) : model_(mod) {}

  static std::string name () { return "UpwindFlux"; }

  const Model& model() const {return model_;}

  /**
   * @brief evaluates the flux \f$g(u,v)\f$
   *
   * @return maximum wavespeed * normal
   */
  template <class QuadratureImp>
  inline double numericalFlux( const IntersectionType& it,
                               const EntityType& inside,
                               const EntityType& outside,
                               const double time,
                               const QuadratureImp& faceQuadInner,
                               const QuadratureImp& faceQuadOuter,
                               const int quadPoint,
                               const RangeType& uLeft,
                               const RangeType& uRight,
                               RangeType& gLeft,
                               RangeType& gRight ) const
  {
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
    const double upwind = Velocity<Model,Model::ConstantVelocity>::
      upwind(model_,it,time,x,uLeft);

    RangeType uFlux;
    if (upwind>0)
      uFlux = uLeft;
    else
      uFlux = uRight;
    uFlux *= upwind;

    /*
    RangeType cFlux;
    cFlux = uLeft;
    cFlux += uRight;
    cFlux *= 0.5;
    cFlux *= upwind;
    */

    gLeft  = uFlux;
    gRight = gLeft;
    return std::abs(upwind);
  }
protected:
  const Model& model_;
};

template <class ModelType>
class LLFAdvFlux
{
public:
  typedef ModelType Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model :: DomainType DomainType;
  typedef typename Model :: RangeType RangeType;
  typedef typename Model :: FluxRangeType FluxRangeType;
  typedef typename Model :: FaceDomainType  FaceDomainType;
  typedef typename Model :: EntityType  EntityType;
  typedef typename Model :: IntersectionType  IntersectionType;
  /**
   * @brief Constructor
   */
  LLFAdvFlux(const Model& mod) : model_(mod) {}

  static std::string name () { return "LaxFriedrichsFlux"; }

  const Model& model() const {return model_;}

  template <class QuadratureImp>
  inline double numericalFlux( const IntersectionType& intersection,
                               const EntityType& inside,
                               const EntityType& outside,
                               const double time,
                               const QuadratureImp& faceQuadInner,
                               const QuadratureImp& faceQuadOuter,
                               const int quadPoint,
                               const RangeType& uLeft,
                               const RangeType& uRight,
                               RangeType& gLeft,
                               RangeType& gRight ) const
  {
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
    DomainType normal = intersection.integrationOuterNormal(x);  
    const double len = normal.two_norm();
    normal *= 1./len;
    
    RangeType visc;
    FluxRangeType anaflux;
    model_.advection( inside, time, faceQuadInner.point( quadPoint ),
                      uLeft, anaflux );
    anaflux.mv( normal, gLeft );
    model_.advection( outside, time, faceQuadOuter.point( quadPoint ),
                      uRight, anaflux );
    anaflux.umv( normal, gLeft );

    double maxspeedl, maxspeedr, maxspeed;
    double viscparal, viscparar, viscpara;
    
    const DomainType xGlobal = intersection.geometry().global(x);

    model_.maxSpeed( normal, time, xGlobal, 
                     uLeft, viscparal, maxspeedl );
    model_.maxSpeed( normal, time, xGlobal,
                     uRight, viscparar, maxspeedr );

    maxspeed = (maxspeedl > maxspeedr) ? maxspeedl : maxspeedr;
    viscpara = (viscparal > viscparar) ? viscparal : viscparar;
    visc = uRight;
    visc -= uLeft;
    visc *= viscpara;
    gLeft -= visc;
    
    gLeft *= 0.5*len;
    gRight = gLeft;

    return 0.;
  }
protected:
  const Model& model_;
};


/*************************
 * Assemble system matrix using DGPrimalDiffusionFlux implementation
 *
 * Remark: at the moment this class sets up the transposed of the system
 *         matrix similar to DgPass::operator2Matrix
 ************************/

template <class DestinationImp,
          class M, 
          class AdvFlux,
          DGDiffusionFluxIdentifier = method_general >
class DGPrimalMatrixAssembly
{
  public:
  typedef DestinationImp DestinationType;
  typedef M Model;
  typedef typename DestinationType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType; 
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType; 
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType; 
  typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename DestinationType::LocalFunctionType LocalFunctionType;

  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry IntersectionGeometryType;

  typedef Fem :: Parameter ParameterType ;

  // need treatment of non conforming grids
  // CACHING
  // typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType; 
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType; 
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef typename GridPartType :: GridType :: template Codim< 0 > :: EntityPointer EntityPointerType;

  typedef ExtendedDGPrimalDiffusionFlux<DiscreteFunctionSpaceType,Model> FluxType;
  typedef AdvFlux AdvFluxType;

  struct ZeroFunction
  {
    typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;

    template <class PointType> 
    void evaluate( const PointType &x, RangeType &phi ) const
    {
      phi = 0;
    }
  };

  typedef Dune::Fem::GridFunctionAdapter< ZeroFunction, GridPartType > ZeroFunctionType;

  struct RangeValues
  {
    typedef std::vector< std::vector< RangeType > > VectorType;
    RangeValues(int col, const VectorType &vec) : col_(col), vec_(vec), zero_(0) 
    {}
    const RangeType &operator[](int row) const
    {
      if (col_==-1) 
        return zero_;
      else 
        return vec_[row][col_];
    }
  private:
    const int col_;
    const VectorType& vec_;
    const RangeType zero_;
  };
  struct JacobianRangeValues
  {
    typedef std::vector< std::vector< JacobianRangeType > > VectorType;
    JacobianRangeValues(int col, const VectorType &vec) : col_(col), vec_(vec), zero_(0) 
    {}
    const JacobianRangeType &operator[](int row) const
    {
      if (col_==-1) 
        return zero_;
      else 
        return vec_[row][col_];
    }
  private:
    const int col_;
    const VectorType& vec_;
    const JacobianRangeType zero_;
  };
  struct IntersectionStorage
  {
    IntersectionStorage( 
        const IntersectionType& intersection,
        const EntityType &entity, const EntityType &neighbor, const double areaEn, const double areaNb )
      : intersection_( intersection ),
        en_(entity), 
        nb_(neighbor), 
        areaEn_(areaEn), areaNb_(areaNb)
    {}

    const IntersectionType& intersection() const { return intersection_; }
    const EntityType &inside() const { return en_; }
    const EntityType &outside() const { return nb_; }
    const double enVolume() const { return areaEn_; }
    const double nbVolume() const { return areaNb_; }

    private:
    const IntersectionType& intersection_;
    const EntityType &en_,&nb_;
    const double areaEn_, areaNb_;
  };

  public:

  DGPrimalMatrixAssembly( GridPartType& gridPart,
                          const Model& model )
    : model_(model),
      space_(gridPart),
      zero_(), 
      advFlux_(model_),
      flux_(gridPart, model),
      calculateFluxes_( ParameterType::getValue<bool>( "use_dgstabilization",true ) ),
      useStrongBoundaryCondition_( ParameterType::getValue<bool>( "use_strongbnd",false ) )
  {
  }

  const DiscreteFunctionSpaceType &space() const
  {
    return space_;
  }

  const typename FluxType::DiscreteGradientSpaceType &gradientSpace() const
  {
    return flux_.gradientSpace();
  }

  const void operator()(const DestinationType &arg, DestinationType &dest) const
  {
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId = Dune :: method_general ;
    typedef DGAdvectionDiffusionOperator< Model, AdvFluxType, PrimalDiffusionFluxId, DiscreteFunctionSpaceType::polynomialOrder >   
            DgOperatorType; 
    static const AdvFluxType flux( model_ );
    static const DgOperatorType dgOperator( space().gridPart(), flux );
    dgOperator(arg,dest);
  }
  template <class Matrix> 
  void operator2Matrix( Matrix& matrix, DestinationType& rhs ) const 
  {
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId = Dune :: method_general ;
    typedef DGAdvectionDiffusionOperator< Model, AdvFluxType, PrimalDiffusionFluxId, DiscreteFunctionSpaceType::polynomialOrder >   
            DgOperatorType; 
    const AdvFluxType flux( model_ );
    const DgOperatorType dgOperator( space().gridPart(), flux );

    DestinationType vector("Op2Mat::arg", space() ); 
    DestinationType matrixRow("Op2Mat::matixRow", space() ); 
    DestinationType result("Op2Mat::result", space() ); 
    vector.clear();  

    matrix.clear();
      
    // get right hand side, = op( 0 )
    this->operator()( vector, rhs );

    const int size = space().size();
    for(int i=0; i<size; ++i) 
    {
      vector.leakPointer()[ i ] = 1;
      this->operator()( vector, matrixRow );
      vector.leakPointer()[ i ] = 0;
      matrixRow *= -1.0;
      matrixRow += rhs;
      for(int j=0; j<size; ++j) 
      {
        const double value = matrixRow.leakPointer()[ j ];
        if( std::abs( value ) > 0 )
        {
          matrix.add( i, j, value );  
        }
      }
    }
  }
  /* 
   * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
   * implementation.
   */
  template <class Matrix> 
  void assemble( const double time, 
                 Matrix& matrix, DestinationType& rhs ) const 
  {
    typedef typename Matrix::LocalMatrixType LocalMatrixType;
    // matrix.reserve();
    matrix.clear();
    rhs.clear();

    const DiscreteFunctionSpaceType &dfSpace = rhs.space();
    const size_t maxNumBaseFunctions = dfSpace.mapper().maxNumDofs();


    std::vector< RangeType > phi( maxNumBaseFunctions );
    std::vector< JacobianRangeType > dphi( maxNumBaseFunctions );

    flux_.initialize( dfSpace );

    const RangeType uZero(0);
    const JacobianRangeType uJacZero(0);


    const IteratorType end = dfSpace.end();
    for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      LocalMatrixType localOpEn = matrix.localMatrix( entity, entity );
      LocalFunctionType rhsLocal = rhs.localFunction( entity );

      const BaseFunctionSetType &baseSet = localOpEn.domainBaseFunctionSet();
      const unsigned int numBaseFunctionsEn = baseSet.size();
            
      QuadratureType quadrature( entity, 
                                 elementQuadOrder( dfSpace.order( entity ) ) );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        const typename GeometryType::Jacobian &gjit = geometry.jacobianInverseTransposed( x );

        // resize of phi and dphi is done in evaluate methods 
        baseSet.evaluateAll( quadrature[ pt ], phi );
        baseSet.jacobianAll( quadrature[ pt ], gjit, dphi );

        RangeType arhs(0);
        // first assemble rhs (evaluate source with u=0)
        if ( model_.hasStiffSource() )
          model_.stiffSource( entity, time, x, uZero, uJacZero, arhs );
        if ( model_.hasNonStiffSource() )
        {
          RangeType sNonStiff (0);
          model_.nonStiffSource( entity, time, x, uZero, uJacZero, sNonStiff );
          arhs += sNonStiff;
        }

        JacobianRangeType arhsdphi;
        model_.diffusion(entity, time, x, uZero, uJacZero, arhsdphi);
        JacobianRangeType brhsphi;
        model_.advection(entity, time, x, uZero, brhsphi);
        arhsdphi -= brhsphi;

        for( unsigned int localCol = 0; localCol < numBaseFunctionsEn; ++localCol )
        {
          // now assemble source part depending on u (mass part)
          RangeType aphi(0);
          if ( model_.hasStiffSource() )
            model_.stiffSource( entity, time, x, phi[localCol], dphi[localCol], aphi );
          if ( model_.hasNonStiffSource() )
          {
            RangeType sNonStiff (0);
            model_.nonStiffSource( entity, time, x, phi[localCol], dphi[localCol], sNonStiff );
            aphi += sNonStiff;
          }
          // subtract affine part and move to left hand side
          aphi -= arhs;
          aphi *= -1;

          JacobianRangeType adphi;
          model_.diffusion(entity, time, x, phi[localCol], dphi[localCol], adphi);

          JacobianRangeType bphi;
          model_.advection(entity, time, x, phi[localCol], bphi);

          adphi -= bphi;

          adphi -= arhsdphi;

          // get column object and call axpy method 
          localOpEn.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
        }
        // assemble rhs
        arhs *= weight;
        arhsdphi *= -weight;
        rhsLocal.axpy( quadrature[ pt ], arhs, arhsdphi );
      }

      const IntersectionIteratorType endiit = dfSpace.gridPart().iend( entity );
      for ( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity );
            iit != endiit ; ++ iit )
      {
        const IntersectionType& intersection = *iit ;
        
        if( intersection.neighbor() && calculateFluxes_ )
        {
          if ( dfSpace.continuous(intersection) ) continue;
          if( intersection.conforming() ) 
          {
            assembleIntersection< true > ( time, entity, geometry, intersection, dfSpace, baseSet, matrix, localOpEn, rhsLocal );
          }
          else 
          {
            assembleIntersection< false > ( time, entity, geometry, intersection, dfSpace, baseSet, matrix, localOpEn, rhsLocal );
          }
        }
        else if ( intersection.boundary() && ! useStrongBoundaryCondition_ )
        {
          IntersectionStorage intersectionStorage( intersection, entity, entity, geometry.volume(), geometry.volume() );
        
          FaceQuadratureType faceQuadInside(dfSpace.gridPart(), intersection,
                                            faceQuadOrder( dfSpace.order( entity ) ),
                                            FaceQuadratureType::INSIDE);

          const size_t numFaceQuadPoints = faceQuadInside.nop();
          resize( numFaceQuadPoints );

          // evalute base functions
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
            baseSet.jacobianAll( faceQuadInside[ pt ], 
                                 geometry.jacobianInverseTransposed( faceQuadInside.point( pt ) ), dphiFaceEn[pt] );
          }

          // storage for all flux values
          //
          boundaryValues(dfSpace.gridPart(),intersection,entity,time,faceQuadInside,bndValues);

          // first compute affine part of boundary flux
          boundaryFlux(dfSpace.gridPart(),intersection,entity,time,faceQuadInside,
                       RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
                       bndValues,
                       valueNb,dvalueNb);

          // for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          //   bndValues[pt] = RangeType(0);
          // compute boundary fluxes depending on u
          for( unsigned int localCol = 0; localCol < numBaseFunctionsEn; ++localCol )
          {
            boundaryFlux(dfSpace.gridPart(),intersection,entity,time,faceQuadInside,
                         RangeValues(localCol,phiFaceEn), JacobianRangeValues(localCol,dphiFaceEn),
                         bndValues,
                         valueEn,dvalueEn);
            for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            {
              const double weight = faceQuadInside.weight( pt );
              valueEn[pt] -= valueNb[pt];
              dvalueEn[pt] -= dvalueNb[pt];
              localOpEn.column( localCol ).axpy( phiFaceEn[pt], dphiFaceEn[pt], 
                                                 valueEn[pt], dvalueEn[pt], weight );
            }
          }
          // now move affine part to right hand side
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            RangeType& rhsFlux          = valueNb[ pt ];
            JacobianRangeType& drhsFlux = dvalueNb[ pt ];

            const double weight = faceQuadInside.weight( pt );
            rhsFlux *= -weight;
            drhsFlux *= -weight;
          }
          rhsLocal.axpyQuadrature( faceQuadInside, valueNb, dvalueNb );
        }
      }
    }
  }
  // assemble vector containing boundary fluxes for right hand side
  void assemble( const double time, 
                 DestinationType& rhs ) const 
  {
    rhs.clear();

    const DiscreteFunctionSpaceType &dfSpace = rhs.space();
    const size_t maxNumBaseFunctions = dfSpace.mapper().maxNumDofs();

    flux_.initialize(dfSpace);

    const RangeType uZero(0);
    const JacobianRangeType uJacZero(0);

    const IteratorType end = dfSpace.end();
    for( IteratorType it = dfSpace.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      LocalFunctionType rhsLocal = rhs.localFunction( entity );

      const BaseFunctionSetType &baseSet = rhsLocal.baseFunctionSet();
      const unsigned int numBaseFunctionsEn = baseSet.size();
            
      const IntersectionIteratorType endiit = dfSpace.gridPart().iend( entity );
      for ( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity );
            iit != endiit ; ++ iit )
      {
        const IntersectionType& intersection = *iit ;
        
        if( intersection.neighbor() && calculateFluxes_ )
        {
        }
        else if ( intersection.boundary() && ! useStrongBoundaryCondition_ )
        {
          IntersectionStorage intersectionStorage( intersection, entity, entity, geometry.volume(), geometry.volume() );
        
          FaceQuadratureType faceQuadInside(dfSpace.gridPart(), intersection,
                                            faceQuadOrder( dfSpace.order( entity ) ),
                                            FaceQuadratureType::INSIDE);

          const size_t numFaceQuadPoints = faceQuadInside.nop();
          resize( numFaceQuadPoints );

          // store all basis functions

          // evalute base functions
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
            baseSet.jacobianAll( faceQuadInside[ pt ], 
                                 geometry.jacobianInverseTransposed( faceQuadInside.point( pt ) ), dphiFaceEn[pt] );
          }

          boundaryValues(dfSpace.gridPart(),intersection,entity,time,faceQuadInside,bndValues);

          // first compute affine part of boundary flux
          boundaryFlux(dfSpace.gridPart(),intersection,entity,time,faceQuadInside,
                       RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
                       bndValues,
                       valueNb,dvalueNb);

          // now move affine part to right hand side
          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            RangeType& rhsFlux          = valueNb[ pt ];
            JacobianRangeType& drhsFlux = dvalueNb[ pt ];

            const double weight = faceQuadInside.weight( pt );
            rhsFlux *= -weight;
            drhsFlux *= -weight;
          }
          rhsLocal.axpyQuadrature( faceQuadInside, valueNb, dvalueNb );
        }
      }
    }
  }

  template <bool conforming, class Matrix>
  void assembleIntersection( const double time, 
                             const EntityType& entity,
                             const GeometryType& geometry,
                             const IntersectionType& intersection,      
                             const DiscreteFunctionSpaceType& dfSpace,
                             const BaseFunctionSetType& baseSet,
                             Matrix& matrix,
                             typename Matrix::LocalMatrixType& localOpEn,
                             LocalFunctionType &rhsLocal) const
  {
    typedef typename Matrix::LocalMatrixType LocalMatrixType;
    // make sure we got the right conforming statement
    assert( intersection.conforming() == conforming );

    EntityPointerType ep = intersection.outside();
    const EntityType& neighbor = *ep ;

    const int entityOrder   = dfSpace.order( entity );
    const int neighborOrder = dfSpace.order( neighbor );

    // get local matrix for face entries 
    LocalMatrixType localOpNb = matrix.localMatrix( entity, neighbor );
    const BaseFunctionSetType &baseSetNb = localOpNb.rangeBaseFunctionSet();

    // get neighbor's base function set 
    const unsigned int numBaseFunctionsEn = baseSet.size();
    const unsigned int numBaseFunctionsNb = baseSetNb.size();

    // only do one sided evaluation if the polynomial orders match 
    const bool oneSidedEvaluation = ( numBaseFunctionsEn == numBaseFunctionsNb );

    const bool updateOnNeighbor   = 
      dfSpace.gridPart().indexSet().index(entity) >
      dfSpace.gridPart().indexSet().index(neighbor) ;

    // only do one sided evaluation if the polynomial orders match 
    if( updateOnNeighbor && oneSidedEvaluation )
      return;

    const int polOrdOnFace = std::max( entityOrder, neighborOrder );

    // use IntersectionQuadrature to create appropriate face quadratures 
    typedef Fem :: IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
    typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;

    // create intersection quadrature 
    IntersectionQuadratureType interQuad( dfSpace.gridPart(), intersection, faceQuadOrder( polOrdOnFace ));

    // get appropriate references 
    const QuadratureImp &faceQuadInside  = interQuad.inside();
    const QuadratureImp &faceQuadOutside = interQuad.outside();

    const GeometryType& nbGeometry = neighbor.geometry();
    const size_t numFaceQuadPoints = faceQuadInside.nop();
    resize( numFaceQuadPoints );

    // evalute base functions
    for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
    {
      // resize is done in evaluateAll and jacobianAll 
      baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
      baseSet.jacobianAll( faceQuadInside[ pt ], 
                           geometry.jacobianInverseTransposed( faceQuadInside.point( pt ) ), dphiFaceEn[pt] );
      baseSetNb.evaluateAll( faceQuadOutside[ pt ], phiFaceNb[pt] );
      baseSetNb.jacobianAll( faceQuadOutside[ pt ], 
                             nbGeometry.jacobianInverseTransposed( faceQuadOutside.point( pt ) ), dphiFaceNb[pt] );
    }

    IntersectionStorage intersectionStorage( intersection,
                  entity, neighbor, 
                  entity.geometry().volume(), neighbor.geometry().volume() );

    flux(dfSpace.gridPart(),
         intersectionStorage,
         time, faceQuadInside, faceQuadOutside,
         RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
         RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
         rhsValueEn, rhsDValueEn, rhsValueNb, rhsDValueNb );

    // compute fluxes and assemble matrix
    for( unsigned int localCol = 0; localCol < numBaseFunctionsEn; ++localCol )
    {
      // compute flux for one base function, i.e.,
      // - uLeft=phiFaceEn[.][localCol]
      // - uRight=0
      flux(dfSpace.gridPart(),
           intersectionStorage,
           time, faceQuadInside, faceQuadOutside,
           RangeValues(localCol,phiFaceEn), JacobianRangeValues(localCol,dphiFaceEn),
           RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
           valueEn, dvalueEn, valueNb, dvalueNb );

      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        const double weight = faceQuadInside.weight( pt );
        valueEn[pt] -= rhsValueEn[pt];
        dvalueEn[pt] -= rhsDValueEn[pt];
        valueNb[pt] -= rhsValueNb[pt];
        dvalueNb[pt] -= rhsDValueNb[pt];
        localOpEn.column( localCol ).axpy( phiFaceEn[pt], dphiFaceEn[pt], 
                                           valueEn[pt], dvalueEn[pt], weight );
        localOpNb.column( localCol ).axpy( phiFaceNb[pt], dphiFaceNb[pt], 
                                           valueNb[pt], dvalueNb[pt], -weight );
      }
    }

    // assemble part from neighboring row
    if( oneSidedEvaluation )
    {
      LocalMatrixType localOpNbNb = matrix.localMatrix( neighbor, neighbor );
      LocalMatrixType localOpNbEn = matrix.localMatrix( neighbor, entity );
      for( unsigned int localCol = 0; localCol < numBaseFunctionsEn; ++localCol )
      {
        // compute flux for one base function, i.e.,
        // - uLeft=phiFaceEn[.][localCol]
        // - uRight=0
        flux(dfSpace.gridPart(),
             intersection, entity, neighbor, time, faceQuadInside, faceQuadOutside,
             RangeValues(-1,phiFaceNb), JacobianRangeValues(-1,dphiFaceNb),
             RangeValues(localCol,phiFaceNb), JacobianRangeValues(localCol,dphiFaceNb),
             valueEn, dvalueEn, valueNb, dvalueNb );
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const double weight = faceQuadInside.weight( pt );
          valueEn[pt] -= rhsValueEn[pt];
          dvalueEn[pt] -= rhsDValueEn[pt];
          valueNb[pt] -= rhsValueNb[pt];
          dvalueNb[pt] -= rhsDValueNb[pt];
          localOpNbNb.column( localCol ).axpy( phiFaceNb[pt], dphiFaceNb[pt],  // +
                                               valueNb[pt], dvalueNb[pt], -weight );
          localOpNbEn.column( localCol ).axpy( phiFaceEn[pt], dphiFaceEn[pt],  // -
                                               valueEn[pt], dvalueEn[pt], weight );
        }
      }
    }

    // now move affine part to right hand side
    for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
    {
      const double weight = faceQuadInside.weight( pt );
      rhsValueEn[pt] *= -weight;
      rhsDValueEn[pt] *= -weight;
    }
    rhsLocal.axpyQuadrature( faceQuadInside, rhsValueEn, rhsDValueEn );
  }

  template <class Matrix> 
  void testSymmetrie(const Matrix &matrix) const
  {
    // NOTE: this is bad coding style 
    // not check at the moment 
    return;


    typedef typename Matrix::LocalMatrixType LocalMatrixType;
    double maxSymError_diagonal = 0;
    double maxSymError_offdiagonal = 0;
    const IteratorType end = space().end();
    for( IteratorType it = space().begin(); it != end; ++it )
    {
      const EntityType& entity = *it;
      LocalMatrixType localOpEn = matrix.localMatrix( entity, entity  );
      const BaseFunctionSetType &baseSet = localOpEn.domainBaseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.size();

      for( unsigned int r = 0; r < numBaseFunctions; ++r )
      {
        for( unsigned int c = 0; c < numBaseFunctions; ++c )
        {
          double v1 = localOpEn.get(r,c);
          double v2 = localOpEn.get(c,r);
          // if (std::abs(v1-v2)>1e-8)
          //  std::cout << "symmetry error on diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
          maxSymError_diagonal = std::max(maxSymError_diagonal,std::abs(v1-v2));
        }
      }
      const IntersectionIteratorType endiit = space().gridPart().iend( entity );
      for ( IntersectionIteratorType iit = space().gridPart().ibegin( entity );
            iit != endiit ; ++ iit )
      {
        const IntersectionType& intersection = *iit ;
        if ( !intersection.neighbor() )
          continue;
        EntityPointerType ep = intersection.outside();
        const EntityType& neighbor = *ep ;

        LocalMatrixType localOpNb1 = matrix.localMatrix( entity, neighbor );
        LocalMatrixType localOpNb2 = matrix.localMatrix( neighbor, entity );
        for( unsigned int r = 0; r < numBaseFunctions; ++r )
        {
          for( unsigned int c = 0; c < numBaseFunctions; ++c )
          {
            double v1 = localOpNb1.get(r,c);
            double v2 = localOpNb2.get(c,r);
            // if (std::abs(v1-v2)>1e-8)
            //    std::cout << "symmetry error on off diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
            maxSymError_offdiagonal = std::max(maxSymError_offdiagonal,std::abs(v1-v2));
          }
        }
      }
    }
    if (maxSymError_diagonal > 1e-8)
      std::cout << "non symmetric diagonal: " << maxSymError_diagonal << std::endl;
    if (maxSymError_offdiagonal > 1e-8)
      std::cout << "non symmetric off diagonal: " << maxSymError_offdiagonal << std::endl;
  }

  template <class Quadrature,class Value,class DValue,class RetType, class DRetType>
  void flux(const GridPartType &gridPart,
            const IntersectionStorage& intersectionStorage,
            const double time,
            const Quadrature &faceQuadInside, const Quadrature &faceQuadOutside,
            const Value &valueEn, const DValue &dvalueEn,
            const Value &valueNb, const DValue &dvalueNb,
            RetType &retEn, DRetType &dretEn, 
            RetType &retNb, DRetType &dretNb) const
  {
    RangeType gLeft,gRight;
    flux_.initializeIntersection( intersectionStorage.intersection(), 
                                  intersectionStorage.inside(),
                                  intersectionStorage.outside(), time,
                                  //zero_, zero_,
                                  faceQuadInside, faceQuadOutside,
                                  valueEn, valueNb);

    const size_t numFaceQuadPoints = faceQuadInside.nop();
    for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
    {
      flux_.numericalFlux( intersectionStorage.intersection(), 
                           intersectionStorage,
                           time, faceQuadInside, faceQuadOutside, pt, 
                           valueEn[ pt ], valueNb[ pt ], dvalueEn[ pt ], dvalueNb[ pt ],
                           retEn[ pt ], retNb[ pt ],
                           dretEn[ pt ], dretNb[ pt ]);
      advFlux_.numericalFlux(intersectionStorage.intersection(),
                             intersectionStorage.inside(),
                             intersectionStorage.outside(),
                             time,faceQuadInside,faceQuadOutside,pt,
                             valueEn[ pt ],valueNb[ pt ], 
                             gLeft, gRight);
      retEn[pt] += gLeft;
      retNb[pt] += gRight;
    }
  }

  template <class Quadrature,class Value,class DValue,class RetType, class DRetType>
  void flux(const GridPartType &gridPart,
            const IntersectionType& intersection,
            const EntityType &entity, const EntityType &neighbor,
            const double time,
            const Quadrature &faceQuadInside, const Quadrature &faceQuadOutside,
            const Value &valueEn, const DValue &dvalueEn,
            const Value &valueNb, const DValue &dvalueNb,
            RetType &retEn, DRetType &dretEn, 
            RetType &retNb, DRetType &dretNb) const
  {
    IntersectionStorage intersectionStorage( intersection,
                  entity, neighbor, 
                  entity.geometry().volume(), neighbor.geometry().volume() );

    flux(gridPart,intersectionStorage,
         time,faceQuadInside,faceQuadOutside,
         valueEn,dvalueEn,valueNb,dvalueNb,
         retEn,dretEn,retNb,dretNb);
  }

  template <class FaceQuadrature,class Value,class DValue,class RetType, class DRetType>
  void fluxAndLift(const GridPartType &gridPart,
                   const IntersectionType &intersection,
                   const EntityType &entity, const EntityType &neighbor,
                   const double time,
                   const FaceQuadrature &faceQuadInside, const FaceQuadrature &faceQuadOutside,
                   const Value &valueEn, const DValue &dvalueEn,
                   const Value &valueNb, const DValue &dvalueNb,
                   RetType &retEn, DRetType &dretEn, 
                   RetType &retNb, DRetType &dretNb,
                   DRetType &liftEn, DRetType &liftNb) const
  {
    IntersectionStorage intersectionStorage( intersection,
                  entity, neighbor, 
                  entity.geometry().volume(), neighbor.geometry().volume() );

    flux(gridPart,intersectionStorage,
         time,faceQuadInside,faceQuadOutside,
         valueEn,dvalueEn,valueNb,dvalueNb,
         retEn,dretEn,retNb,dretNb);

    const size_t numFaceQuadPoints = faceQuadInside.nop();
    for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
    {
      flux_.evaluateLifting(faceQuadInside,faceQuadOutside,pt,time,
                            valueEn[pt],valueNb[pt],
                            liftEn[pt],liftNb[pt]);
    }
  }
  template <class FaceQuadrature,class Quadrature,class Value,class DValue,class RetType, class DRetType>
  void fluxAndLift(const GridPartType &gridPart,
                   const IntersectionType &intersection,
                   const EntityType &entity, const EntityType &neighbor,
                   const double time,
                   const FaceQuadrature &faceQuadInside, const FaceQuadrature &faceQuadOutside,
                   const Value &faceValueEn, const DValue &dfaceValueEn,
                   const Value &faceValueNb, const DValue &dfaceValueNb,
                   RetType &retEn, DRetType &dretEn, 
                   RetType &retNb, DRetType &dretNb,
                   const Quadrature &quadInside, const Quadrature &quadOutside, 
                   const Value &valueEn, const Value &valueNb,
                   DRetType &liftEn, DRetType &liftNb) const
  {
    IntersectionStorage intersectionStorage( intersection,
                  entity, neighbor, 
                  entity.geometry().volume(), neighbor.geometry().volume() );

    flux(gridPart,intersectionStorage,time,faceQuadInside,faceQuadOutside,
         faceValueEn,dfaceValueEn,faceValueNb,dfaceValueNb,
         retEn,dretEn,retNb,dretNb);
    const size_t numQuadPoints = quadInside.nop();
    for( size_t pt = 0; pt < numQuadPoints; ++pt )
    {
      flux_.evaluateLifting(quadInside,quadOutside,pt,time,
                            valueEn[pt],valueNb[pt],
                            liftEn[pt],liftNb[pt]);
    }
  }
  template <class FaceQuadrature,class Value,class LiftingFunction>
  void lifting(const GridPartType &gridPart,
               const IntersectionType &intersection,
               const EntityType &entity, const EntityType &neighbor,
               const double time,
               const FaceQuadrature &faceQuadInside, const FaceQuadrature &faceQuadOutside,
               const Value &valueEn, const Value &valueNb,
               LiftingFunction &lifting) const
  {
    flux_.initializeIntersection( intersection, entity, neighbor, time,
                                  // zero_, zero_, 
                                  faceQuadInside, faceQuadOutside,
                                  valueEn, valueNb, true );
    lifting += flux_.getInsideLifting();
  }

  template <class Quadrature,class RetType>
  void boundaryValues(const GridPartType &gridPart,
            const IntersectionType &intersection,
            const EntityType &entity, 
            const double time,
            const Quadrature &faceQuadInside, 
            RetType &bndValues) const
  {
    const RangeType uZero(0);
    const size_t numFaceQuadPoints = faceQuadInside.nop();
    for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      model_.boundaryValue(intersection, time,  faceQuadInside.localPoint( pt ), 
                           uZero, bndValues[ pt ]);
  }
  template <class Quadrature,class Value,class DValue,class GValue,class RetType, class DRetType>
  void boundaryFlux(const GridPartType &gridPart,
            const IntersectionType &intersection,
            const EntityType &entity, 
            const double time,
            const Quadrature &faceQuadInside, 
            const Value &valueEn, const DValue &dvalueEn,
            const GValue &valueNb,
            RetType &retEn, DRetType &dretEn) const
  {
    RangeType gLeft,gRight;
    flux_.initializeBoundary( intersection, entity, time, 
                              // zero_,
                              faceQuadInside, 
                              valueEn, valueNb );
    IntersectionStorage intersectionStorage( intersection, entity, entity, entity.geometry().volume(), entity.geometry().volume() );
    const size_t numFaceQuadPoints = faceQuadInside.nop();
    for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
    {
      if ( model_.hasBoundaryValue(intersection,time,faceQuadInside.localPoint(pt)) )
      {
        flux_.boundaryFlux( intersection, intersectionStorage, 
                            time, faceQuadInside, pt, 
                            valueEn[ pt ], valueNb[ pt ],
                            dvalueEn[ pt ], 
                            retEn[ pt ], dretEn[ pt ]);
        advFlux_.numericalFlux(intersection,entity,entity,
                               time,faceQuadInside,faceQuadInside,pt,
                               valueEn[ pt ],valueNb[ pt ],
                               gLeft,gRight);
        retEn[pt] += gLeft;
      }
      else
      {
        model_.boundaryFlux(intersection,time,faceQuadInside.localPoint(pt),valueEn[pt],retEn[pt]);
        dretEn[pt] = 0;
      }
    }

  }

  const Model &model() const
  {
    return model_;
  }

private:
  int elementQuadOrder( int polOrder ) const 
  {
    return 2*polOrder;
  }
  int faceQuadOrder( int polOrder ) const
  {
    return 2*polOrder;
  }

  void resize( unsigned int numFaceQuadPoints ) const
  {
    if (phiFaceEn.size() >= numFaceQuadPoints) 
      return;

    phiFaceEn.resize( numFaceQuadPoints );   
    dphiFaceEn.resize( numFaceQuadPoints );
    phiFaceNb.resize( numFaceQuadPoints );
    dphiFaceNb.resize( numFaceQuadPoints );
    valueEn.resize( numFaceQuadPoints );  
    dvalueEn.resize( numFaceQuadPoints );
    valueNb.resize( numFaceQuadPoints );
    dvalueNb.resize( numFaceQuadPoints );
    rhsValueEn.resize( numFaceQuadPoints );
    rhsDValueEn.resize( numFaceQuadPoints );
    rhsValueNb.resize( numFaceQuadPoints );
    rhsDValueNb.resize( numFaceQuadPoints );
    bndValues.resize( numFaceQuadPoints );
  }

  const Model &model_;
  DiscreteFunctionSpaceType space_;
  ZeroFunction zero_;

  AdvFluxType advFlux_;
  mutable FluxType flux_;
  const bool calculateFluxes_ ;
  const bool useStrongBoundaryCondition_ ;

  // storage for all flux values
  mutable std::vector< RangeType > valueEn;
  mutable std::vector< JacobianRangeType > dvalueEn;
  mutable std::vector< RangeType > valueNb;
  mutable std::vector< JacobianRangeType > dvalueNb;
  mutable std::vector< RangeType > rhsValueEn;
  mutable std::vector< JacobianRangeType > rhsDValueEn;
  mutable std::vector< RangeType > rhsValueNb;
  mutable std::vector< JacobianRangeType > rhsDValueNb;

  // store all basis functions
  mutable std::vector< std::vector< RangeType > > phiFaceEn;
  mutable std::vector< std::vector< JacobianRangeType > > dphiFaceEn;
  mutable std::vector< std::vector< RangeType > > phiFaceNb;
  mutable std::vector< std::vector< JacobianRangeType > > dphiFaceNb;
  mutable std::vector< RangeType > bndValues;
};

};

#endif
