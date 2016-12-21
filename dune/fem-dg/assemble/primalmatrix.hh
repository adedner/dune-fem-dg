#ifndef PRIMALMATRIXASSEMBLY_HH
#define PRIMALMATRIXASSEMBLY_HH

#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/compatibility.hh>

#include <dune/fem-dg/algorithm/sub/steadystate.hh>
#include <dune/fem-dg/algorithm/sub/elliptic.hh>
#include <dune/fem-dg/misc/uniquefunctionname.hh>
#include <dune/fem-dg/misc/matrixutils.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/dgprimalfluxes.hh>
#include "assemblertraits.hh"

namespace Dune
{
namespace Fem
{

  template< class, class >
  struct InsertFunctionTuple;

  template <class Tuple, unsigned long int... i >
  struct InsertFunctionTuple< Tuple, std::integer_sequence< unsigned long int, i...> >
  {
    typedef std::tuple< std::shared_ptr< typename std::tuple_element< i, Tuple >::type >... > TupleType;

    template< class ExtraArgImp >
    static decltype(auto) createTuple( const std::shared_ptr< ExtraArgImp >& tuple )
    {
      return std::make_shared<TupleType>( (*tuple)( _index<i>() )... );
    }
  };


  /**
   * \brief Assembles the primal DG matrix.
   *
   * \ingroup AssemblyOperator
   */
  template <class Traits>
  class DGPrimalMatrixAssembly
  {
    public:
    typedef typename Traits::DomainDiscreteFunctionType           DomainDiscreteFunctionType;
    typedef typename Traits::RangeDiscreteFunctionType            RangeDiscreteFunctionType;
    typedef typename Traits::MatrixContainerType                  LinearOperatorType;
    typedef LinearOperatorType                                    MatrixType;
    typedef DomainDiscreteFunctionType                            DestinationType;

    template< class Row, class Col >
    using FakeMatrixAdapter = MatrixType;

    typedef Dune::Fem::SubEllipticContainer< FakeMatrixAdapter, DestinationType, DestinationType > ContainerType;

    typedef typename Traits::ModelType                            ModelType;
    static const bool hasDiffusion = ModelType::hasDiffusion;

    typedef typename DestinationType::DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::IteratorType      IteratorType;
    typedef typename IteratorType::Entity                         EntityType;
    typedef typename EntityType::Geometry                         GeometryType;

    typedef typename DiscreteFunctionSpaceType::GridPartType      GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType        DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType         RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType   DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType    RangeFieldType;

    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    typedef typename DestinationType::LocalFunctionType           LocalFunctionType;

    typedef typename GridPartType::IntersectionIteratorType       IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType               IntersectionType;
    typedef typename IntersectionType::Geometry                   IntersectionGeometryType;

    typedef Fem::Parameter                                        ParameterType;

    // need treatment of non conforming grids
    typedef typename Traits::FaceQuadratureType                   FaceQuadratureType;
    typedef typename Traits::VolumeQuadratureType                 VolumeQuadratureType;


    typedef typename Traits::AdvectionFluxType                    AdvectionFluxType;
    typedef ExtendedDGPrimalDiffusionFlux<DiscreteFunctionSpaceType, ModelType,
                                          typename Traits::DiffusionFluxType::ParameterType >
                                                                  DiffusionFluxType;


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

      const RangeType& at( int i ) const { return this->operator[] ( i ); }
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

      const RangeType& at( int i ) const { return this->operator[] ( i ); }
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
      const JacobianRangeType& at( int i ) const { return this->operator[] ( i ); }
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
    //! constructor for DG matrix assembly
    template< class ContainerImp >
    DGPrimalMatrixAssembly( std::shared_ptr< ContainerImp > cont,
                            const ModelType& model )
      : model_( model ),
        space_( (*cont)(_0)->solution()->space() ),
        rhs_( (*cont)(_0)->rhs() ),
        matrix_( (*cont)(_0,_0)->matrix() ),
        zero_(),
        time_( 0 ),
        advFlux_( model_ ),
        diffusionFlux_( space_.gridPart(), model_, typename Traits::DiffusionFluxType::ParameterType( ParameterKey::generate( "", "dgdiffusionflux." ) ) ),
        calculateFluxes_( Dune::Fem::Parameter::getValue<bool>( "poissonassembler.calculateFluxes", true ) ),
        useStrongBoundaryCondition_( Dune::Fem::Parameter::getValue<bool>( "poissonassembler.strongBC", false ) )
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

    size_t maxNumScalarBasisFunctions( const DiscreteFunctionSpaceType& space ) const
    {
      return space.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType :: localBlockSize ;
    }

    void setTime ( double time )
    {
      time_ = time;
    }


    /*
     * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
     * implementation.
     */
    void assemble() const
    {
      Dune::Timer timer ;

      typedef RangeType           RangeTuple;
      typedef JacobianRangeType   JacobianTuple;
      typedef ExtraQuadraturePointContext< EntityType, VolumeQuadratureType, RangeTuple, JacobianTuple > LocalEvaluationType;

      typedef typename MatrixType::LocalMatrixType LocalMatrixType;
      matrix_->clear();
      if( rhs_ )
      {
        rhs_->clear();
      }

      const size_t maxNumBasisFunctions = maxNumScalarBasisFunctions( space_ );

      std::vector< RangeType > phi( maxNumBasisFunctions );
      std::vector< JacobianRangeType > dphi( maxNumBasisFunctions );

      diffusionFlux_.initialize( space_ );

      const RangeType uZero(0);
      const JacobianRangeType uJacZero(0);

      Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > rhsLocal( space_ );

      for( const auto& entity : elements( space_.gridPart() ) )
      {
        const GeometryType &geometry = entity.geometry();
        const double volume = geometry.volume();

        LocalMatrixType localOpEn = matrix_->localMatrix( entity, entity );

        if( rhs_ )
        {
          rhsLocal.init( entity );
          rhsLocal.clear();
        }

        const BasisFunctionSetType &baseSet = localOpEn.domainBasisFunctionSet();
        const unsigned int numBasisFunctionsEn = baseSet.size();

        VolumeQuadratureType quadrature( entity, elementQuadOrder( space_.order( entity ) ) );

        for( const auto qp : quadrature )
        {
          LocalEvaluationType local( entity, quadrature, qp.index(), time_, volume, uZero, uJacZero );

          const double weight = qp.weight() * geometry.integrationElement( local.position() );

          // resize of phi and dphi is done in evaluate methods
          baseSet.evaluateAll( qp, phi );
          baseSet.jacobianAll( qp, dphi );

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

          if( rhs_ )
          {
            // assemble rhs
            arhs     *=  weight;
            arhsdphi *= -weight;
            rhsLocal.axpy( qp, arhs, arhsdphi );
          }
        }

        for (const auto& intersection : intersections(space_.gridPart(), entity) )
        {
          if( intersection.neighbor() && calculateFluxes_ )
          {
            if ( space_.continuous(intersection) ) continue;
            if( intersection.conforming() )
            {
              assembleIntersection< true > ( entity, intersection, space_, baseSet, localOpEn, rhsLocal, rhs_ != 0 );
            }
            else
            {
              assembleIntersection< false > ( entity, intersection, space_, baseSet, localOpEn, rhsLocal, rhs_ != 0 );
            }
          }
          else if ( intersection.boundary() && ! useStrongBoundaryCondition_ )
          {
            FaceQuadratureType faceQuadInside(space_.gridPart(), intersection,
                                              faceQuadOrder( space_.order( entity ) ),
                                              FaceQuadratureType::INSIDE);

            const size_t numFaceQuadPoints = faceQuadInside.nop();
            resize( numFaceQuadPoints, maxNumBasisFunctions );

            // evaluate base functions
            for( const auto qp : faceQuadInside )
            {
              baseSet.evaluateAll( qp, phiFaceEn[qp.index()] );
              baseSet.jacobianAll( qp, dphiFaceEn[qp.index()] );
            }

            typedef QuadratureContext< EntityType, IntersectionType, FaceQuadratureType > ContextType;
            ContextType local( entity, intersection, faceQuadInside, time_, volume );

            // storage for all flux values
            boundaryValues(local, bndValues);

            // first compute affine part of boundary flux
            boundaryFlux(local,
                         RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
                         bndValues, valueNb,dvalueNb);

            // for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            //   bndValues[pt] = RangeType(0);
            // compute boundary fluxes depending on u
            for( unsigned int localCol = 0; localCol < numBasisFunctionsEn; ++localCol )
            {
              boundaryFlux(local,
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

            if( rhs_ )
            {
              rhsLocal.axpyQuadrature( faceQuadInside, valueNb, dvalueNb );
            }
          }
        }

        // accumulate right hand side
        if( rhs_ )
        {
          rhs_->localFunction( entity ) += rhsLocal ;
        }

      } // end grid iteration

      // finish matrix build process
      matrix_->communicate();

     // matrix.systemMatrix().matrix().print();
     // //rhs->print( std::cout );
     // //abort();
     // //

     // const int sep = 8;
     // int i;
     // int j;

     // std::cout << "###########" << std::endl;
     // for(i=0; i<space_.size(); ++i)
     // {
     //   if( i % sep == 0 )
     //   {
     //     for( int j=0; j<space_.size(); j++)
     //       std::cout << "------";
     //     std::cout << std::endl;
     //   }
     //   for(j=0;j<space_.size(); ++j)
     //   {
     //     if( j % sep == 0 )
     //       std::cout << "|";

     //     std::cout.width(5);
     //     if( std::abs(  matrix.matrix()(j,i) ) < 1e-14 )
     //       std::cout << std::setprecision(2) << "0" << " ";
     //     else
     //       std::cout << std::setprecision(2) << matrix.matrix()(j,i) << " ";
     //   }
     //   if( j % sep == 0 )
     //     std::cout << "|";

     //   std::cout << std::endl;
     // }
     // if( i % sep == 0 )
     // {
     //   for( int j=0; j<space_.size(); j++)
     //     std::cout << "------";
     //   std::cout << std::endl;
     // }
     // std::cout << std::endl;
     // std::cout << "###########" << std::endl;


      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout << "DG( " << space_.grid().size( 0 ) << " ) matrix assemble took " << timer.elapsed() << " sec." << std::endl;
      }
    }
    // assemble vector containing boundary fluxes for right hand side
    void assembleRhs() const
    {
      rhs_->clear();

      const size_t maxNumBasisFunctions = maxNumScalarBasisFunctions( space_ );

      if( hasDiffusion )
      {
        diffusionFlux_.initialize(space_);
      }

      const RangeType uZero(0);
      const JacobianRangeType uJacZero(0);

      for( const auto& entity : elements( space_.gridPart() ) )
      {
        const GeometryType &geometry = entity.geometry();
        const double volume = geometry.volume();

        LocalFunctionType rhsLocal = rhs_->localFunction( entity );

        const BasisFunctionSetType &baseSet = rhsLocal.baseFunctionSet();
        const unsigned int numBasisFunctionsEn = baseSet.size();

        for (const auto& intersection : intersections(space_.gridPart(), entity) )
        {

          if( intersection.neighbor() && calculateFluxes_ )
          {
          }
          else if ( intersection.boundary() && ! useStrongBoundaryCondition_ )
          {
            FaceQuadratureType faceQuadInside(space_.gridPart(), intersection,
                                              faceQuadOrder( space_.order( entity ) ),
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

            typedef QuadratureContext< EntityType, IntersectionType, FaceQuadratureType > ContextType;
            ContextType local( entity, intersection, faceQuadInside, time_, volume );

            // store for all flux values
            boundaryValues(local, bndValues);

            // first compute affine part of boundary flux
            boundaryFlux(local,
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

    template <bool conforming, class LocalFunction>
    void assembleIntersection( const EntityType& entity,
                               const IntersectionType& intersection,
                               const DiscreteFunctionSpaceType& dfSpace,
                               const BasisFunctionSetType& baseSet,
                               typename MatrixType::LocalMatrixType& localOpEn,
                               LocalFunction& rhsLocal,
                               const bool assembleRHS ) const
    {
      const size_t maxNumBasisFunctions = maxNumScalarBasisFunctions( dfSpace );
      typedef typename MatrixType::LocalMatrixType LocalMatrixType;
      // make sure we got the right conforming statement
      assert( intersection.conforming() == conforming );

      const EntityType& neighbor = intersection.outside();

      const int entityOrder   = dfSpace.order( entity );
      const int neighborOrder = dfSpace.order( neighbor );

      // get local matrix for face entries
      LocalMatrixType localOpNb = matrix_->localMatrix( entity, neighbor );
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

      typedef QuadratureContext< EntityType, IntersectionType, QuadratureImp > ContextType;
      ContextType localInside ( entity, intersection, faceQuadInside, time_, entity.geometry().volume() );
      ContextType localOutside( neighbor, intersection, faceQuadOutside, time_, neighbor.geometry().volume() );

      const size_t numFaceQuadPoints = faceQuadInside.nop();
      resize( numFaceQuadPoints, maxNumBasisFunctions );

      // evalute base functions
      for( const auto qp : faceQuadInside )
      {
        // resize is done in evaluateAll and jacobianAll
        const auto idx = qp.index();
        baseSet.evaluateAll( qp, phiFaceEn[idx] );
        baseSet.jacobianAll( qp, dphiFaceEn[idx] );
        baseSetNb.evaluateAll( faceQuadOutside[idx], phiFaceNb[idx] );
        baseSetNb.jacobianAll( faceQuadOutside[idx], dphiFaceNb[idx] );
      }

      numericalFlux(localInside, localOutside,
                    RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
                    RangeValues(-1,phiFaceEn), JacobianRangeValues(-1,dphiFaceEn),
                    rhsValueEn, rhsDValueEn, rhsValueNb, rhsDValueNb, true );

      // compute fluxes and assemble matrix
      for( unsigned int localCol = 0; localCol < numBasisFunctionsEn; ++localCol )
      {
        // compute flux for one base function, i.e.,
        // - uLeft=phiFaceEn[.][localCol]
        // - uRight=0
        numericalFlux(localInside, localOutside,
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
        LocalMatrixType localOpNbNb = matrix_->localMatrix( neighbor, neighbor );
        LocalMatrixType localOpNbEn = matrix_->localMatrix( neighbor, entity );
        for( unsigned int localCol = 0; localCol < numBasisFunctionsEn; ++localCol )
        {
          // compute flux for one base function, i.e.,
          // - uLeft=phiFaceEn[.][localCol]
          // - uRight=0
          numericalFlux(localInside, localOutside,
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

    void testSymmetry() const
    {
      //if( matrix_ )
      //  testMatrixSymmetry( *matrix_ );
    }

    template <class LocalEvaluation, class Value,class DValue,class RetType, class DRetType>
    void numericalFlux(const LocalEvaluation& inside, const LocalEvaluation& outside,
                       const Value &valueEn, const DValue &dvalueEn,
                       const Value &valueNb, const DValue &dvalueNb,
                       RetType &retEn, DRetType &dretEn,
                       RetType &retNb, DRetType &dretNb,
                       const bool initializeIntersection = true ) const
    {
      RangeType gLeft,gRight;
      if( hasDiffusion & initializeIntersection )
      {
        diffusionFlux_.initializeIntersection( inside, outside, valueEn, valueNb);
      }

      const size_t numFaceQuadPoints = inside.quadrature().nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        typedef RangeType           RangeTuple;
        typedef JacobianRangeType   JacobianTuple;
        typedef ExtraQuadraturePointContext< EntityType, IntersectionType, typename LocalEvaluation::QuadratureType, RangeTuple, JacobianTuple > LocalEvaluationType;

        LocalEvaluationType left( inside, pt, valueEn[ pt ], dvalueEn[ pt ] );
        LocalEvaluationType right( outside, pt, valueNb[ pt ], dvalueNb[ pt ] );

        if( hasDiffusion )
        {
          diffusionFlux_.numericalFlux( left, right,
                                        valueEn[ pt ], valueNb[ pt ],
                                        dvalueEn[ pt ], dvalueNb[ pt ],
                                        retEn[ pt ], retNb[ pt ],
                                        dretEn[ pt ], dretNb[ pt ]);
        }
        else
        {
          retEn[pt]  = RangeType(0);
          retNb[pt]  = RangeType(0);
          dretEn[pt] = JacobianRangeType(0);
          dretNb[pt] = JacobianRangeType(0);
        }

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
      numericalFlux(left, right, valueEn,dvalueEn,valueNb,dvalueNb,
                    retEn,dretEn,retNb,dretNb);

     if( hasDiffusion )
     {
       const size_t numFaceQuadPoints = left.quadrature().nop();
       for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
       {
          diffusionFlux_.evaluateLifting(left, right, valueEn[pt],valueNb[pt],
                                         liftEn[pt],liftNb[pt]);
       }
      }
    }

    template <class LocalEvaluation, class Value,class LiftingFunction>
    void lifting(const LocalEvaluation& left, const LocalEvaluation& right,
                 const Value &valueEn, const Value &valueNb,
                 LiftingFunction &lifting) const
    {
      VectorToTupleVector valEn( valueEn );
      VectorToTupleVector valNb( valueNb );
      if( hasDiffusion )
      {
        diffusionFlux_.initializeIntersection( left, right, valEn, valNb, true );
        lifting += diffusionFlux_.getInsideLifting();
      }
    }

    template <class LocalEvaluation,class RetType>
    void boundaryValues(const LocalEvaluation& local,
                        RetType &bndValues) const
    {
      const RangeType uZero(0);
      const JacobianRangeType uJacZero( 0 );

      typedef RangeType           RangeTuple;
      typedef JacobianRangeType   JacobianTuple;
      typedef ExtraQuadraturePointContext< EntityType, IntersectionType, typename LocalEvaluation::QuadratureType, RangeTuple, JacobianTuple > LocalEvaluationType;

      const size_t numFaceQuadPoints = local.quadrature().nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        LocalEvaluationType local2( local, pt, uZero, uJacZero );
        model_.boundaryValue( local2, uZero, bndValues[ pt ]);
      }
    }

    template <class LocalEvaluation, class Value,class DValue,class GValue,class RetType, class DRetType>
    void boundaryFlux( const LocalEvaluation& local,
                       const Value &valueEn,
                       const DValue &dvalueEn,
                       const GValue &valueNb,
                       RetType &retEn,
                       DRetType &dretEn) const
    {
      RangeType gLeft,gRight;
      if( hasDiffusion )
      {
        diffusionFlux_.initializeBoundary( local, valueEn, valueNb );
      }

      const size_t numFaceQuadPoints = local.quadrature().nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        typedef RangeType           RangeTuple;
        typedef JacobianRangeType   JacobianTuple;
        typedef ExtraQuadraturePointContext< EntityType, IntersectionType, typename LocalEvaluation::QuadratureType, RangeTuple, JacobianTuple > LocalEvaluationType;
        LocalEvaluationType local2( local, pt, valueEn[ pt ], dvalueEn[ pt ] );

        if ( model_.hasBoundaryValue( local2 ) )
        {
          if( hasDiffusion )
          {
            diffusionFlux_.boundaryFlux( local2, valueEn[ pt ], valueNb[ pt ],  dvalueEn[ pt ],
                                         retEn[ pt ], dretEn[ pt ]);
          }
          else
          {
            retEn[pt]  = RangeType(0);
            dretEn[pt] = JacobianRangeType(0);
          }
          advFlux_.numericalFlux(local2, local2,
                                 valueEn[ pt ],valueNb[ pt ],
                                 dvalueEn[ pt ], dvalueEn[ pt ],
                                 gLeft,gRight);
          retEn[pt] += gLeft;
        }
        else
        {
          model_.boundaryFlux(local2,valueEn[pt], dvalueEn[ pt ], retEn[pt]);
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

    const ModelType&                   model_;
    const DiscreteFunctionSpaceType&   space_;
    std::shared_ptr< DestinationType > rhs_;
    std::shared_ptr< MatrixType >      matrix_;
    ZeroFunction                       zero_;
    double                             time_;

    AdvectionFluxType                  advFlux_;
    mutable DiffusionFluxType          diffusionFlux_;
    const bool                         calculateFluxes_;
    const bool                         useStrongBoundaryCondition_ ;

    // storage for all flux values
    mutable std::vector< RangeType >         valueEn;
    mutable std::vector< JacobianRangeType > dvalueEn;
    mutable std::vector< RangeType >         valueNb;
    mutable std::vector< JacobianRangeType > dvalueNb;
    mutable std::vector< RangeType >         rhsValueEn;
    mutable std::vector< JacobianRangeType > rhsDValueEn;
    mutable std::vector< RangeType >         rhsValueNb;
    mutable std::vector< JacobianRangeType > rhsDValueNb;
    mutable std::vector< RangeType >         bndValues;

    // store all basis functions
    mutable std::vector< std::vector< RangeType > >         phiFaceEn;
    mutable std::vector< std::vector< JacobianRangeType > > dphiFaceEn;
    mutable std::vector< std::vector< RangeType > >         phiFaceNb;
    mutable std::vector< std::vector< JacobianRangeType > > dphiFaceNb;
  };

}
}

#endif
