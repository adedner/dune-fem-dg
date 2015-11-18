#ifndef PRIMALMATRIXASSEMBLY_HH
#define PRIMALMATRIXASSEMBLY_HH

#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/compatibility.hh>

#include <dune/fem-dg/operator/fluxes/diffusion/dgprimalfluxes.hh>
#include <dune/fem-dg/operator/dg/primaloperator.hh>

namespace Dune
{
namespace Fem
{


  /*************************
   * Assemble system matrix using DGPrimalDiffusionFlux implementation
   *
   * Remark: at the moment this class sets up the transposed of the system
   *         matrix similar to DgPass::operator2Matrix
   ************************/

  template <class OperatorImp>
  class DGPrimalMatrixAssembly
  {
    public:
    typedef OperatorImp OperatorType;
    typedef typename OperatorType::Traits::ModelType ModelType;
    typedef typename OperatorType::Traits::DestinationType DestinationType;
    //typedef typename OperatorType::Traits::DiscreteModelType DiscreteModelType;
    //typedef typename DiscreteModelType::Selector   SelectorType;
    //static const Dune::DGDiffusionFluxIdentifier DGDiffusionFluxIdentifier = OperatorType::Traits::PrimalDiffusionFluxId;
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

    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename IntersectionType::Geometry IntersectionGeometryType;

    //typedef DiscreteModelType::ExtraParameterType   ExtraParameterType;
    //typedef DiscreteModelType::ModelParameter   MoelParameterType;
    //typedef CDGDiscreteModelCaller< DiscreteModelType, ExtraParameterType, PassIds > DiscreteModelCallerType;

    //typedef typename OperatorType::DiscreteModelCallerType DiscreteModelCallerType;

    typedef Fem :: Parameter ParameterType ;

    // need treatment of non conforming grids
    // CACHING
    // typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;
    typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > VolumeQuadratureType;

    typedef ExtendedDGPrimalDiffusionFlux<DiscreteFunctionSpaceType,ModelType> DiffusionFluxType;
    // deprecated FluxType
    typedef DiffusionFluxType  FluxType;

    typedef typename OperatorType :: AdvectionFluxType AdvectionFluxType;

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

      const RangeType& tuple( int i ) const { return this->operator[] ( i ); }
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

    struct VectorToTupleVector
    {
      typedef std::vector< RangeType > VectorType;
      VectorToTupleVector(const VectorType &vec) : vec_(vec)
      {}

      const RangeType& tuple( int i ) const { return this->operator[] ( i ); }
      const RangeType &operator[](int i) const
      {
        return vec_[ i ];
      }
    private:
      const VectorType& vec_;
    };
    struct JacobianRangeValues
    {
      typedef std::vector< std::vector< JacobianRangeType > > VectorType;
      JacobianRangeValues(int col, const VectorType &vec) : col_(col), vec_(vec), zero_(0)
      {}
      const JacobianRangeType& tuple( int i ) const { return this->operator[] ( i ); }
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

    public:

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



    //! constructor for DG matrix assembly
    DGPrimalMatrixAssembly( GridPartType& gridPart,
                            const OperatorType& op,
                            const bool calculateFluxes = true,
                            const bool strongBC = false )
      : op_(op),
        model_(op_.model()),
        space_(gridPart),
        zero_(),
        advFlux_(model_),
        diffusionFlux_(gridPart, model_, DGPrimalFormulationParameters( ParameterKey::generate( model_.problem().name(), "dgdiffusionflux." ) ) ),
        calculateFluxes_( calculateFluxes ),
        useStrongBoundaryCondition_( strongBC )
    {
    }

    const DiscreteFunctionSpaceType &space() const
    {
      return space_;
    }

    const typename DiffusionFluxType::DiscreteGradientSpaceType &gradientSpace() const
    {
      return diffusionFlux_.gradientSpace();
    }

    const void operator()(const DestinationType &arg, DestinationType &dest) const
    {
      op_( arg, dest );
    }

    size_t maxNumScalarBasisFunctions( const DiscreteFunctionSpaceType& space ) const
    {
      return space.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType :: localBlockSize ;
    }

    /*
     * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
     * implementation.
     */
    template <class Matrix>
    void assemble( const double time,
                   Matrix& matrix, DestinationType& rhs ) const
    {
      doAssemble( time, matrix, &rhs );
    }

    template <class Matrix>
    void assemble( const double time,
                   Matrix& matrix ) const
    {
      doAssemble( time, matrix, (DestinationType * ) 0 );
    }


    /*
     * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
     * implementation.
     */
    template <class Matrix>
    void doAssemble( const double time,
                     Matrix& matrix, DestinationType* rhs ) const
    {
      Dune::Timer timer ;

      typedef RangeType           RangeTuple;
      typedef JacobianRangeType   JacobianTuple;
      typedef ElementQuadraturePointContext< EntityType, VolumeQuadratureType, RangeTuple, JacobianTuple >      LocalEvaluationType;

      typedef typename Matrix::LocalMatrixType LocalMatrixType;
      matrix.clear();
      if( rhs )
      {
        rhs->clear();
      }

      const DiscreteFunctionSpaceType &dfSpace = op_.space();
      const size_t maxNumBasisFunctions = maxNumScalarBasisFunctions( dfSpace );

      std::vector< RangeType > phi( maxNumBasisFunctions );
      std::vector< JacobianRangeType > dphi( maxNumBasisFunctions );

      diffusionFlux_.initialize( dfSpace );

      const RangeType uZero(0);
      const JacobianRangeType uJacZero(0);

      Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > rhsLocal( dfSpace );

      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        const GeometryType &geometry = entity.geometry();
        const double volume = geometry.volume();

        LocalMatrixType localOpEn = matrix.localMatrix( entity, entity );

        if( rhs )
        {
          rhsLocal.init( entity );
          rhsLocal.clear();
        }

        const BasisFunctionSetType &baseSet = localOpEn.domainBasisFunctionSet();
        const unsigned int numBasisFunctionsEn = baseSet.size();

        VolumeQuadratureType quadrature( entity, elementQuadOrder( dfSpace.order( entity ) ) );
        const size_t numQuadraturePoints = quadrature.nop();
        for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
        {
          LocalEvaluationType local( entity, quadrature, uZero, uJacZero, pt, time, volume );

          const double weight = quadrature.weight( pt ) * geometry.integrationElement( local.point() );

          // resize of phi and dphi is done in evaluate methods
          baseSet.evaluateAll( quadrature[ pt ], phi );
          baseSet.jacobianAll( quadrature[ pt ], dphi );

          RangeType arhs(0);
          // first assemble rhs (evaluate source with u=0)
          if ( model_.hasStiffSource() )
          model_.stiffSource( local, uZero, uJacZero, arhs );
          if ( model_.hasNonStiffSource() )
          {
            RangeType sNonStiff (0);
            model_.nonStiffSource( local, uZero, uJacZero, sNonStiff );
            arhs += sNonStiff;
          }

          //if model_.hasFlux()
          //else()
          //endif
          //...

          JacobianRangeType arhsdphi;
          model_.diffusion( local, uZero, uJacZero, arhsdphi);
          JacobianRangeType brhsphi;
          model_.advection( local, uZero, uJacZero, brhsphi);
          arhsdphi -= brhsphi;

          for( unsigned int localCol = 0; localCol < numBasisFunctionsEn; ++localCol )
          {
            // now assemble source part depending on u (mass part)
            RangeType aphi(0);
            if ( model_.hasStiffSource() )
            {
              model_.stiffSource( local, phi[localCol], dphi[localCol], aphi );
            }
            if ( model_.hasNonStiffSource() )
            {
              RangeType sNonStiff (0);
              model_.nonStiffSource( local, phi[localCol], dphi[localCol], sNonStiff );
              aphi += sNonStiff;
            }
            // subtract affine part and move to left hand side
            aphi -= arhs;
            aphi *= -1;

            JacobianRangeType adphi;
            model_.diffusion( local, phi[localCol], dphi[localCol], adphi);

            JacobianRangeType bphi;
            model_.advection( local, phi[localCol], dphi[localCol], bphi);

            adphi -= bphi;

            adphi -= arhsdphi;

            // get column object and call axpy method
            localOpEn.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
          }

          if( rhs )
          {
            // assemble rhs
            arhs     *=  weight;
            arhsdphi *= -weight;
            rhsLocal.axpy( quadrature[ pt ], arhs, arhsdphi );
          }
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
              assembleIntersection< true > ( time, entity, geometry, intersection, dfSpace, baseSet, matrix, localOpEn, rhsLocal, rhs != 0 );
            }
            else
            {
              assembleIntersection< false > ( time, entity, geometry, intersection, dfSpace, baseSet, matrix, localOpEn, rhsLocal, rhs != 0 );
            }
          }
          else if ( intersection.boundary() && ! useStrongBoundaryCondition_ )
          {
            FaceQuadratureType faceQuadInside(dfSpace.gridPart(), intersection,
                                              faceQuadOrder( dfSpace.order( entity ) ),
                                              FaceQuadratureType::INSIDE);

            const size_t numFaceQuadPoints = faceQuadInside.nop();
            resize( numFaceQuadPoints, maxNumBasisFunctions );

            // evalute base functions
            for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            {
              baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
              baseSet.jacobianAll( faceQuadInside[ pt ], dphiFaceEn[pt] );
            }

            // storage for all flux values
            boundaryValues(intersection, entity, time, faceQuadInside, volume, bndValues);

            // first compute affine part of boundary flux
            boundaryFlux(intersection, entity, time, volume, faceQuadInside,
                         RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
                         bndValues, valueNb,dvalueNb);

            // for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            //   bndValues[pt] = RangeType(0);
            // compute boundary fluxes depending on u
            for( unsigned int localCol = 0; localCol < numBasisFunctionsEn; ++localCol )
            {
              boundaryFlux(intersection, entity, time, volume, faceQuadInside,
                           RangeValues(localCol,phiFaceEn), JacobianRangeValues(localCol,dphiFaceEn),bndValues,
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
              rhsFlux  *= -weight;
              drhsFlux *= -weight;
            }

            if( rhs )
            {
              rhsLocal.axpyQuadrature( faceQuadInside, valueNb, dvalueNb );
            }
          }
        }

        // accumulate right hand side
        if( rhs )
        {
          rhs->localFunction( entity ) += rhsLocal ;
        }

      } // end grid iteration

      // finish matrix build process
      matrix.communicate();
      //matrix.systemMatrix().matrix().print( std::cout );
      //rhs.print( std::cout );
      //abort();
      //
      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout << "DG( " << dfSpace.grid().size( 0 ) << " ) matrix assemble took " << timer.elapsed() << " sec." << std::endl;
      }
    }
    // assemble vector containing boundary fluxes for right hand side
    void assemble( const double time,
                   DestinationType& rhs ) const
    {
      rhs.clear();

      const DiscreteFunctionSpaceType &dfSpace = rhs.space();
      const size_t maxNumBasisFunctions = maxNumScalarBasisFunctions( dfSpace );

#ifndef EULER
      diffusionFlux_.initialize(dfSpace);
#endif

      const RangeType uZero(0);
      const JacobianRangeType uJacZero(0);

      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        const GeometryType &geometry = entity.geometry();
        const double volume = geometry.volume();

        LocalFunctionType rhsLocal = rhs.localFunction( entity );

        const BasisFunctionSetType &baseSet = rhsLocal.baseFunctionSet();
        const unsigned int numBasisFunctionsEn = baseSet.size();

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
            FaceQuadratureType faceQuadInside(dfSpace.gridPart(), intersection,
                                              faceQuadOrder( dfSpace.order( entity ) ),
                                              FaceQuadratureType::INSIDE);

            const size_t numFaceQuadPoints = faceQuadInside.nop();
            resize( numFaceQuadPoints, maxNumBasisFunctions );

            // store all basis functions

            // evalute base functions
            for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            {
              baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
              baseSet.jacobianAll( faceQuadInside[ pt ], dphiFaceEn[pt] );
            }

            // store for all flux values
            boundaryValues(intersection, entity, time, faceQuadInside, volume, bndValues);

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
              rhsFlux  *= -weight;
              drhsFlux *= -weight;
            }
            rhsLocal.axpyQuadrature( faceQuadInside, valueNb, dvalueNb );
          }
        }
      }
    }

    template <bool conforming, class Matrix, class LocalFunction>
    void assembleIntersection( const double time,
                               const EntityType& entity,
                               const GeometryType& geometry,
                               const IntersectionType& intersection,
                               const DiscreteFunctionSpaceType& dfSpace,
                               const BasisFunctionSetType& baseSet,
                               Matrix& matrix,
                               typename Matrix::LocalMatrixType& localOpEn,
                               LocalFunction& rhsLocal,
                               const bool assembleRHS ) const
    {
      const size_t maxNumBasisFunctions = maxNumScalarBasisFunctions( dfSpace );
      typedef typename Matrix::LocalMatrixType LocalMatrixType;
      // make sure we got the right conforming statement
      assert( intersection.conforming() == conforming );

      const EntityType& neighbor = Fem::make_entity( intersection.outside() );

      const int entityOrder   = dfSpace.order( entity );
      const int neighborOrder = dfSpace.order( neighbor );

      // get local matrix for face entries
      LocalMatrixType localOpNb = matrix.localMatrix( entity, neighbor );
      const BasisFunctionSetType &baseSetNb = localOpNb.rangeBasisFunctionSet();

      // get neighbor's base function set
      const unsigned int numBasisFunctionsEn = baseSet.size();
      const unsigned int numBasisFunctionsNb = baseSetNb.size();

      // only do one sided evaluation if the polynomial orders match
      const bool oneSidedEvaluation = ( numBasisFunctionsEn == numBasisFunctionsNb );

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

      // create intersection quadrature (without neighbor check)
      IntersectionQuadratureType interQuad( dfSpace.gridPart(), intersection, faceQuadOrder( polOrdOnFace ), true);

      // get appropriate references
      const QuadratureImp &faceQuadInside  = interQuad.inside();
      const QuadratureImp &faceQuadOutside = interQuad.outside();


      IntersectionStorage intersectionStorage( intersection,
                                               entity, neighbor,
                                               entity.geometry().volume(),
                                               neighbor.geometry().volume() );

      // const GeometryType& nbGeometry = neighbor.geometry();
      const size_t numFaceQuadPoints = faceQuadInside.nop();
      resize( numFaceQuadPoints, maxNumBasisFunctions );

      // evalute base functions
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        // resize is done in evaluateAll and jacobianAll
        baseSet.evaluateAll( faceQuadInside[ pt ], phiFaceEn[pt] );
        baseSet.jacobianAll( faceQuadInside[ pt ], dphiFaceEn[pt] );
        baseSetNb.evaluateAll( faceQuadOutside[ pt ], phiFaceNb[pt] );
        baseSetNb.jacobianAll( faceQuadOutside[ pt ], dphiFaceNb[pt] );
      }

      flux(dfSpace.gridPart(), intersectionStorage, time,
           faceQuadInside, faceQuadOutside,
           RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
           RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
           rhsValueEn, rhsDValueEn, rhsValueNb, rhsDValueNb, true );

      // compute fluxes and assemble matrix
      for( unsigned int localCol = 0; localCol < numBasisFunctionsEn; ++localCol )
      {
        // compute flux for one base function, i.e.,
        // - uLeft=phiFaceEn[.][localCol]
        // - uRight=0
        flux(dfSpace.gridPart(), intersectionStorage, time,
             faceQuadInside, faceQuadOutside,
             RangeValues(localCol,phiFaceEn), JacobianRangeValues(localCol,dphiFaceEn),
             RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
             valueEn, dvalueEn, valueNb, dvalueNb, true );

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
        for( unsigned int localCol = 0; localCol < numBasisFunctionsEn; ++localCol )
        {
          // compute flux for one base function, i.e.,
          // - uLeft=phiFaceEn[.][localCol]
          // - uRight=0
          flux(dfSpace.gridPart(), intersectionStorage, time,
               faceQuadInside, faceQuadOutside,
               RangeValues(-1,phiFaceNb), JacobianRangeValues(-1,dphiFaceNb),
               RangeValues(localCol,phiFaceNb), JacobianRangeValues(localCol,dphiFaceNb),
               valueEn, dvalueEn, valueNb, dvalueNb, true );
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

      if( assembleRHS )
      {
        rhsLocal.axpyQuadrature( faceQuadInside, rhsValueEn, rhsDValueEn );
      }
    }

    template <class Matrix>
    void testSymmetrie(const Matrix &matrix) const
    {
      // NOTE: this is bad coding style
      // not check at the moment
      return;
  //
  //
  //    typedef typename Matrix::LocalMatrixType LocalMatrixType;
  //    double maxSymError_diagonal = 0;
  //    double maxSymError_offdiagonal = 0;
  //    const IteratorType end = space().end();
  //    for( IteratorType it = space().begin(); it != end; ++it )
  //    {
  //      const EntityType& entity = *it;
  //      LocalMatrixType localOpEn = matrix.localMatrix( entity, entity  );
  //      const BasisFunctionSetType &baseSet = localOpEn.domainBasisFunctionSet();
  //      const unsigned int numBasisFunctions = baseSet.size();
  //
  //      for( unsigned int r = 0; r < numBasisFunctions; ++r )
  //      {
  //        for( unsigned int c = 0; c < numBasisFunctions; ++c )
  //        {
  //          double v1 = localOpEn.get(r,c);
  //          double v2 = localOpEn.get(c,r);
  //          // if (std::abs(v1-v2)>1e-8)
  //          //  std::cout << "symmetry error on diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
  //          maxSymError_diagonal = std::max(maxSymError_diagonal,std::abs(v1-v2));
  //        }
  //      }
  //      const IntersectionIteratorType endiit = space().gridPart().iend( entity );
  //      for ( IntersectionIteratorType iit = space().gridPart().ibegin( entity );
  //            iit != endiit ; ++ iit )
  //      {
  //        const IntersectionType& intersection = *iit ;
  //        if ( !intersection.neighbor() )
  //          continue;
  //        EntityPointerType ep = intersection.outside();
  //        const EntityType& neighbor = *ep ;
  //
  //        LocalMatrixType localOpNb1 = matrix.localMatrix( entity, neighbor );
  //        LocalMatrixType localOpNb2 = matrix.localMatrix( neighbor, entity );
  //        for( unsigned int r = 0; r < numBasisFunctions; ++r )
  //        {
  //          for( unsigned int c = 0; c < numBasisFunctions; ++c )
  //          {
  //            double v1 = localOpNb1.get(r,c);
  //            double v2 = localOpNb2.get(c,r);
  //            // if (std::abs(v1-v2)>1e-8)
  //            //    std::cout << "symmetry error on off diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
  //            maxSymError_offdiagonal = std::max(maxSymError_offdiagonal,std::abs(v1-v2));
  //          }
  //        }
  //      }
  //    }
  //    if (maxSymError_diagonal > 1e-8)
  //      std::cout << "non symmetric diagonal: " << maxSymError_diagonal << std::endl;
  //    if (maxSymError_offdiagonal > 1e-8)
  //      std::cout << "non symmetric off diagonal: " << maxSymError_offdiagonal << std::endl;
    }

    template <class Quadrature,class Value,class DValue,class RetType, class DRetType>
    void flux(const GridPartType &gridPart,
              const IntersectionStorage& intersectionStorage,
              const double time,
              const Quadrature &faceQuadInside, const Quadrature &faceQuadOutside,
              const Value &valueEn, const DValue &dvalueEn,
              const Value &valueNb, const DValue &dvalueNb,
              RetType &retEn, DRetType &dretEn,
              RetType &retNb, DRetType &dretNb,
              const bool initializeIntersection = true ) const
    {
      RangeType gLeft,gRight;
#ifndef EULER
      if( initializeIntersection )
      {
        diffusionFlux_.initializeIntersection( intersectionStorage.intersection(),
                                               intersectionStorage.inside(),
                                               intersectionStorage.outside(), time,
                                               faceQuadInside, faceQuadOutside,
                                               valueEn, valueNb);
      }
#endif
      const size_t numFaceQuadPoints = faceQuadInside.nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        typedef RangeType           RangeTuple;
        typedef JacobianRangeType   JacobianTuple;
        typedef IntersectionQuadraturePointContext< IntersectionType, EntityType, Quadrature, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;
        IntersectionLocalEvaluationType left( intersectionStorage.intersection(), intersectionStorage.inside(),
                                              faceQuadInside, valueEn[ pt ], dvalueEn[ pt ], pt, time, intersectionStorage.enVolume() );
        IntersectionLocalEvaluationType right( intersectionStorage.intersection(), intersectionStorage.outside(),
                                              faceQuadInside, valueNb[ pt ], dvalueNb[ pt ], pt, time, intersectionStorage.nbVolume() );

#ifndef EULER
        diffusionFlux_.numericalFlux( left, right,
                                      valueEn[ pt ], valueNb[ pt ],
                                      dvalueEn[ pt ], dvalueNb[ pt ],
                                      retEn[ pt ], retNb[ pt ],
                                      dretEn[ pt ], dretNb[ pt ]);
#endif
        advFlux_.numericalFlux(left, right,
                               valueEn[ pt ],valueNb[ pt ],
                               dvalueEn[ pt ], dvalueNb[ pt ],
                               gLeft, gRight);
        retEn[pt] += gLeft;
        retNb[pt] += gRight;
      }
    }

    template <class LocalEvaluation,class Value,class DValue,class RetType, class DRetType>
    void fluxAndLift(const LocalEvaluation& left,
                     const LocalEvaluation& right,
                     const Value &valueEn, const DValue &dvalueEn,
                     const Value &valueNb, const DValue &dvalueNb,
                     RetType &retEn, DRetType &dretEn,
                     RetType &retNb, DRetType &dretNb,
                     DRetType &liftEn, DRetType &liftNb) const
    {
      flux(left, right, valueEn,dvalueEn,valueNb,dvalueNb,
           retEn,dretEn,retNb,dretNb);

#ifndef EULER
     const size_t numFaceQuadPoints = left.quadrature().nop();
     for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
     {
        diffusionFlux_.evaluateLifting(left, right, valueEn[pt],valueNb[pt],
                                       liftEn[pt],liftNb[pt]);
     }
#endif
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
      VectorToTupleVector valEn( valueEn );
      VectorToTupleVector valNb( valueNb );
#ifndef EULER
      diffusionFlux_.initializeIntersection( intersection, entity, neighbor, time, faceQuadInside, faceQuadOutside, valEn, valNb, true );
      lifting += diffusionFlux_.getInsideLifting();
#endif
    }

    template <class Quadrature,class RetType>
    void boundaryValues(const IntersectionType &intersection,
                        const EntityType &entity,
                        const double time,
                        const Quadrature &faceQuadInside,
                        const double volume,
                        RetType &bndValues) const
    {
      const RangeType uZero(0);
      const JacobianRangeType uJacZero( 0 );

      typedef RangeType           RangeTuple;
      typedef JacobianRangeType   JacobianTuple;
      typedef IntersectionQuadraturePointContext< IntersectionType, EntityType, Quadrature, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;

      const size_t numFaceQuadPoints = faceQuadInside.nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        IntersectionLocalEvaluationType local( intersection, entity, faceQuadInside, uZero, uJacZero, pt, time, volume );
        model_.boundaryValue(local, uZero, bndValues[ pt ]);
      }
    }

    template <class QuadratureImp, class Value,class DValue,class GValue,class RetType, class DRetType>
    void boundaryFlux( const IntersectionType& intersection,
                       const EntityType& entity,
                       const double time,
                       const double volume,
                       const QuadratureImp& faceQuadInside,
                       const Value &valueEn,
                       const DValue &dvalueEn,
                       const GValue &valueNb,
                       RetType &retEn,
                       DRetType &dretEn) const
    {
      RangeType gLeft,gRight;
#ifndef EULER
      diffusionFlux_.initializeBoundary( intersection, entity, time, faceQuadInside, valueEn, valueNb );
#endif


      const size_t numFaceQuadPoints = faceQuadInside.nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        typedef RangeType           RangeTuple;
        typedef JacobianRangeType   JacobianTuple;
        typedef IntersectionQuadraturePointContext< IntersectionType, EntityType, QuadratureImp, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;
        IntersectionLocalEvaluationType local( intersection, entity, faceQuadInside, valueEn[ pt ], dvalueEn[ pt ], pt, time, volume );

        if ( model_.hasBoundaryValue( local ) )
        {
#ifndef EULER
          diffusionFlux_.boundaryFlux( local, valueEn[ pt ], valueNb[ pt ],  dvalueEn[ pt ],
                                       retEn[ pt ], dretEn[ pt ]);
#endif
          advFlux_.numericalFlux(local, local,
                                 valueEn[ pt ],valueNb[ pt ],
                                 dvalueEn[ pt ], dvalueEn[ pt ],
                                 gLeft,gRight);
          retEn[pt] += gLeft;
        }
        else
        {
          model_.boundaryFlux(local,valueEn[pt], dvalueEn[ pt ], retEn[pt]);
          dretEn[pt] = 0;
        }
      }

    }

    const ModelType &model() const
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

    void resize( unsigned int numFaceQuadPoints, unsigned int maxNumBasisFunctions ) const
    {
      if (phiFaceEn.size() >= numFaceQuadPoints)
        return;

      phiFaceEn.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) phiFaceEn[i].resize(maxNumBasisFunctions);
      dphiFaceEn.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) dphiFaceEn[i].resize(maxNumBasisFunctions);
      phiFaceNb.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) phiFaceNb[i].resize(maxNumBasisFunctions);
      dphiFaceNb.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) dphiFaceNb[i].resize(maxNumBasisFunctions);
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

    const OperatorType& op_;
    const ModelType &model_;
    DiscreteFunctionSpaceType space_;
    ZeroFunction zero_;

    AdvectionFluxType advFlux_;
    mutable DiffusionFluxType diffusionFlux_;
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
    mutable std::vector< RangeType > bndValues;

    // store all basis functions
    mutable std::vector< std::vector< RangeType > > phiFaceEn;
    mutable std::vector< std::vector< JacobianRangeType > > dphiFaceEn;
    mutable std::vector< std::vector< RangeType > > phiFaceNb;
    mutable std::vector< std::vector< JacobianRangeType > > dphiFaceNb;
  };

}
}

#endif
