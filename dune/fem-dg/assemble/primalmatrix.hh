#ifndef PRIMALMATRIXASSEMBLY_HH
#define PRIMALMATRIXASSEMBLY_HH

#include <dune/common/hybridutilities.hh>

#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/operator/common/localcontribution.hh>

#include <dune/fem-dg/algorithm/sub/steadystate.hh>
#include <dune/fem-dg/algorithm/sub/elliptic.hh>
#include <dune/fem-dg/misc/uniquefunctionname.hh>
#include <dune/fem-dg/misc/matrixutils.hh>
#include <dune/fem-dg/operator/fluxes/advection/parameters.hh>
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
    typedef std::tuple< std::shared_ptr< std::tuple_element_t< i, Tuple > >... > type;

    template< class... AddRanges >
    using RangeTypes = tuple_concat_t< std::tuple< typename std::tuple_element_t< i, Tuple >::RangeType... >, std::tuple< std::vector< AddRanges >... > >;

    template< class... AddRanges >
    using JacobianRangeTypes = tuple_concat_t< std::tuple< typename std::tuple_element_t< i, Tuple >::JacobianRangeType... >, std::tuple< std::vector< AddRanges >... > >;

    template< class ExtraArgImp >
    static decltype(auto) create( const std::shared_ptr< ExtraArgImp >& tuple )
    {
      return std::make_shared<type>( (*tuple)( _index<i>() )... );
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

    typedef Dune::Fem::AddLocalContribution< LinearOperatorType > LocalMatrixContributionType;
    //typedef Dune::Fem::TemporaryLocalMatrix<
    //  typename DomainDiscreteFunctionType::DiscreteFunctionSpaceType,
    //  typename RangeDiscreteFunctionType::DiscreteFunctionSpaceType > TemporaryLocalMatrixType;
    //typedef TemporaryLocalMatrixType  LocalMatrixContributionType;

    template< class Row, class Col >
    using FakeMatrixAdapter = MatrixType;

    typedef Dune::Fem::SubEllipticContainer< FakeMatrixAdapter, DestinationType, DestinationType > ContainerType;

    typedef typename Traits::ModelType                            ModelType;
    typedef typename Traits::ExtraParameterTupleType              ExtraParameterTupleType;
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
    static const int size = ModelType::modelParameterSize;

    typedef InsertFunctionTuple< ExtraParameterTupleType, std::make_integer_sequence<unsigned long int, size > >
                                                                  InsertFunctionTupleType;

    typedef std::vector< typename InsertFunctionTupleType::template RangeTypes< RangeType > > RangeEvalType;
    typedef std::vector< typename InsertFunctionTupleType::template JacobianRangeTypes< JacobianRangeType > > JacobianEvalType;

    std::integral_constant<int,size> Id;

  public:
    //! constructor for DG matrix assembly
    template< class ContainerImp, class ExtraParameterTupleImp >
    DGPrimalMatrixAssembly( std::shared_ptr< ContainerImp > cont,
                            ExtraParameterTupleImp tuple,
                            const ModelType& model )
      : model_( model ),
        space_( (*cont)(_0)->solution()->space() ),
        rhs_( (*cont)(_0)->rhs() ),
        matrix_( (*cont)(_0,_0)->matrix() ),
        extra_( InsertFunctionTupleType::create( tuple ) ),
        advFlux_( model_, AdvectionFluxParameters() ),
        diffFlux_( space_.gridPart(), model_, typename Traits::DiffusionFluxType::ParameterType( ParameterKey::generate( "", "dgdiffusionflux." ) ) ),
        calculateFluxes_( Dune::Fem::Parameter::getValue<bool>( "poissonassembler.calculateFluxes", true ) ),
        useStrongBoundaryCondition_( Dune::Fem::Parameter::getValue<bool>( "poissonassembler.strongBC", false ) ),
        maxNumBasisFunctions_( maxNumScalarBasisFunctions( space_ ) )
    {
    }

    const DiscreteFunctionSpaceType &space() const
    {
      return space_;
    }

    const typename DiffusionFluxType::DiscreteGradientSpaceType &gradientSpace() const
    {
      return diffFlux_.gradientSpace();
    }

    size_t maxNumScalarBasisFunctions( const DiscreteFunctionSpaceType& space ) const
    {
      return space.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ;
    }

    /*
     * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
     * implementation.
     */
    void assemble() const
    {
      Dune::Timer timer ;

      //typedef RangeType           RangeTuple;
      //typedef JacobianRangeType   JacobianTuple;

      matrix_->clear();
      if( rhs_ )
      {
        rhs_->clear();
      }

      std::vector< RangeType > phi( maxNumBasisFunctions_ );
      std::vector< JacobianRangeType > dphi( maxNumBasisFunctions_ );

      diffFlux_.initialize( space_ );

      const RangeType uZero(0);
      const JacobianRangeType uJacZero(0);

      Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > rhsLocal( space_ );

      LocalMatrixContributionType localOpEn( *matrix_ );
      LocalMatrixContributionType localOpNb( *matrix_ );
      LocalMatrixContributionType localOpNbNb( *matrix_ );
      LocalMatrixContributionType localOpNbEn( *matrix_ );
      //LocalMatrixContributionType localOpEn( space_, space_ );
      //LocalMatrixContributionType localOpNb( space_, space_ );
      //LocalMatrixContributionType localOpNbNb( space_, space_ );
      //LocalMatrixContributionType localOpNbEn( space_, space_ );

      for( const auto& entity : elements( space_.gridPart() ) )
      {
        const GeometryType &geometry = entity.geometry();
        const double volume = geometry.volume();

        auto guard = Dune::Fem::bindGuard( localOpEn, entity, entity );

        if( rhs_ )
        {
          rhsLocal.init( entity );
          rhsLocal.clear();
        }

        const BasisFunctionSetType &baseSet = localOpEn.domainBasisFunctionSet();
        const int numBasisFunctionsEn = baseSet.size();

        VolumeQuadratureType quadrature( entity, elementQuadOrder( space_.order( entity ) ) );

        //set entity values
        setEntity( baseSet, quadrature );

        typedef QuadratureContext< EntityType, VolumeQuadratureType > ContextType;
        typedef LocalEvaluation< ContextType, RangeEvalType, JacobianEvalType > LocalEvaluationType;
        ContextType cLocal( entity, quadrature, volume );

        LocalEvaluationType localF( cLocal, phi_, dphi_ );

        for( const auto qp : quadrature )
        {
          const auto& local = localF[qp.index()];
          const auto& local0 = local(0);

          const double weight = qp.weight() * geometry.integrationElement( local.position() );

          RangeType arhs(0);
          // first assemble rhs (evaluate source with u=0)
          if ( model_.hasStiffSource() )
            model_.stiffSource( local0, uZero, uJacZero, arhs );
          if ( model_.hasNonStiffSource() )
          {
            RangeType sNonStiff (0);
            model_.nonStiffSource( local0, uZero, uJacZero, sNonStiff );
            arhs += sNonStiff;
          }

          JacobianRangeType arhsdphi;
          model_.diffusion( local0, uZero, uJacZero, arhsdphi);
          JacobianRangeType brhsphi;
          model_.advection( local0, uZero, uJacZero, brhsphi);
          arhsdphi -= brhsphi;

          for( int i = 0; i < numBasisFunctionsEn; ++i )
          {
            const auto& locali = local[i];
            const auto& phii = locali.values()[Id];
            const auto& dphii = locali.jacobians()[Id];

            // now assemble source part depending on u (mass part)
            RangeType aphi(0);
            if ( model_.hasStiffSource() )
            {
              model_.stiffSource( locali, phii, dphii, aphi );
            }
            if ( model_.hasNonStiffSource() )
            {
              RangeType sNonStiff (0);
              model_.nonStiffSource( locali, phii, dphii, sNonStiff );
              aphi += sNonStiff;
            }
            // subtract affine part and move to left hand side
            aphi -= arhs;
            aphi *= -1;

            JacobianRangeType adphi;
            model_.diffusion( locali, phii, dphii, adphi);

            JacobianRangeType bphi;
            model_.advection( locali, phii, dphii, bphi);

            adphi -= bphi;

            adphi -= arhsdphi;

            // get column object and call axpy method
            localOpEn.matrixColumn( i ).axpy( std::get<size>(phi_[qp.index()]), std::get<size>(dphi_[qp.index()]), aphi, adphi, weight );
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
              assembleIntersection< true > ( entity, intersection, space_, baseSet,
                                             localOpEn, localOpNb, localOpNbNb, localOpNbEn,
                                             rhsLocal, rhs_ != 0 );
            }
            else
            {
              assembleIntersection< false > ( entity, intersection, space_, baseSet,
                                              localOpEn, localOpNb, localOpNbNb, localOpNbEn,
                                              rhsLocal, rhs_ != 0 );
            }
          }
          else if ( intersection.boundary() && ! useStrongBoundaryCondition_ )
          {
            FaceQuadratureType faceQuadInside(space_.gridPart(), intersection,
                                              faceQuadOrder( space_.order( entity ) ),
                                              FaceQuadratureType::INSIDE);

            // initialize
            initializeBoundary( baseSet, faceQuadInside );

            typedef QuadratureContext< EntityType, IntersectionType, FaceQuadratureType > ContextType;
            typedef LocalEvaluation< ContextType, RangeEvalType, JacobianEvalType > LocalEvaluationType;
            ContextType cLocal( entity, intersection, faceQuadInside, volume );
            LocalEvaluationType local( cLocal, phiFaceEn_, dphiFaceEn_ );

            // first compute affine part of boundary flux
            const auto& local0 = local(0);
            boundaryValues(local0, bndValues_);
            boundaryFlux(local0, local0.values()[Id], local0.jacobians()[Id],
                         bndValues_, valueNb_, dvalueNb_);

            const size_t numFaceQuadPoints = faceQuadInside.nop();
            // compute boundary fluxes depending on u
            for( int i = 0; i < numBasisFunctionsEn; ++i )
            {
              const auto& locali = local[i];
              boundaryFlux(locali, locali.values()[Id], locali.jacobians()[Id],
                           bndValues_, valueEn_, dvalueEn_);

              for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
              {
                const double weight = faceQuadInside.weight( pt );
                valueEn_[pt] -= valueNb_[pt];
                dvalueEn_[pt] -= dvalueNb_[pt];
                const auto& phiFaceEn  = std::get<size>(phiFaceEn_[pt]);
                const auto& dphiFaceEn = std::get<size>(dphiFaceEn_[pt]);
                localOpEn.matrixColumn( i ).axpy(phiFaceEn, dphiFaceEn,
                                           valueEn_[pt], dvalueEn_[pt], weight );
              }
            }
            // now move affine part to right hand side
            for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            {
              RangeType& rhsFlux          = valueNb_[ pt ];
              JacobianRangeType& drhsFlux = dvalueNb_[ pt ];

              const double weight = faceQuadInside.weight( pt );
              rhsFlux  *= -weight;
              drhsFlux *= -weight;
            }

            if( rhs_ )
            {
              rhsLocal.axpyQuadrature( faceQuadInside, valueNb_, dvalueNb_ );
            }
          }
        }

        // accumulate right hand side
        if( rhs_ )
        {
          rhs_->addLocalDofs( entity, rhsLocal.localDofVector() );
        }

      } // end grid iteration

      // finish matrix build process
      matrix_->finalize();

      //matrix.systemMatrix().matrix().print();
      //rhs->print( std::cout );
      //abort();
      //

      //const int sep = 8;
      //int i;
      //int j;

      //std::cout << "###########" << std::endl;
      //for(i=0; i<space_.size(); ++i)
      //{
      //  if( i % sep == 0 )
      //  {
      //    for( int j=0; j<space_.size(); j++)
      //      std::cout << "------";
      //    std::cout << std::endl;
      //  }
      //  for(j=0;j<space_.size(); ++j)
      //  {
      //    if( j % sep == 0 )
      //      std::cout << "|";

      //    std::cout.width(5);
      //    if( std::abs(  matrix_->matrix()(j,i) ) < 1e-14 )
      //      std::cout << std::setprecision(2) << "0" << " ";
      //    else
      //      std::cout << std::setprecision(2) << matrix_->matrix()(j,i) << " ";
      //  }
      //  if( j % sep == 0 )
      //    std::cout << "|";

      //  std::cout << std::endl;
      //}
      //if( i % sep == 0 )
      //{
      //  for( int j=0; j<space_.size(); j++)
      //    std::cout << "------";
      //  std::cout << std::endl;
      //}
      //std::cout << std::endl;
      //std::cout << "###########" << std::endl;


      if( Dune::Fem::Parameter::verbose() )
      {
        std::cout << "DG( " << space_.grid().size( 0 ) << " ) matrix assemble took " << timer.elapsed() << " sec." << std::endl;
      }
    }
    // assemble vector containing boundary fluxes for right hand side
    void assembleRhs() const
    {
      rhs_->clear();

      if( hasDiffusion )
      {
        diffFlux_.initialize(space_);
      }

      const RangeType uZero(0);
      const JacobianRangeType uJacZero(0);

      MutableLocalFunction< DestinationType > rhsLocal( *rhs_ );

      for( const auto& entity : elements( space_.gridPart() ) )
      {
        const GeometryType &geometry = entity.geometry();
        const double volume = geometry.volume();

        auto uGuard = Dune::Fem::bindGuard( rhsLocal, entity );

        const BasisFunctionSetType &baseSet = rhsLocal.baseFunctionSet();

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

            // initialize
            initializeBoundary( baseSet, faceQuadInside );

            typedef QuadratureContext< EntityType, IntersectionType, FaceQuadratureType > ContextType;
            typedef LocalEvaluation< ContextType, RangeEvalType, JacobianEvalType > LocalEvaluationType;
            ContextType cLocal( entity, intersection, faceQuadInside, volume );
            LocalEvaluationType local( cLocal, phiFaceEn_, dphiFaceEn_ );

            // store for all flux values
            boundaryValues(local, bndValues_);

            // first compute affine part of boundary flux
            const auto& local0 = local(0);
            boundaryValues(local0, bndValues_);
            boundaryFlux(local0, local0.values()[Id], local0.jacobians()[Id],
                         bndValues_, valueNb_, dvalueNb_);

            const size_t numFaceQuadPoints = faceQuadInside.nop();
            // now move affine part to right hand side
            for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
            {
              RangeType& rhsFlux          = valueNb_[ pt ];
              JacobianRangeType& drhsFlux = dvalueNb_[ pt ];

              const double weight = faceQuadInside.weight( pt );
              rhsFlux  *= -weight;
              drhsFlux *= -weight;
            }
            rhsLocal.axpyQuadrature( faceQuadInside, valueNb_, dvalueNb_ );
          }
        }
      }
    }

    template <bool conforming, class LocalFunction>
    void assembleIntersection( const EntityType& entity,
                               const IntersectionType& intersection,
                               const DiscreteFunctionSpaceType& dfSpace,
                               const BasisFunctionSetType& baseSet,
                               LocalMatrixContributionType& localOpEn,
                               LocalMatrixContributionType& localOpNb,
                               LocalMatrixContributionType& localOpNbNb,
                               LocalMatrixContributionType& localOpNbEn,
                               LocalFunction& rhsLocal,
                               const bool assembleRHS ) const
    {
      // make sure we got the right conforming statement
      assert( intersection.conforming() == conforming );

      const EntityType& neighbor = intersection.outside();

      const int entityOrder   = dfSpace.order( entity );
      const int neighborOrder = dfSpace.order( neighbor );

      // get local matrix for face entries
      auto nGuard = Dune::Fem::bindGuard( localOpNb, entity, neighbor );
      const BasisFunctionSetType &baseSetNb = localOpNb.rangeBasisFunctionSet();

      // get neighbor's base function set
      const int numBasisFunctionsEn = baseSet.size();
      const int numBasisFunctionsNb = baseSetNb.size();

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

      //initialize
      initializeIntersection( baseSet,  baseSetNb, faceQuadInside, faceQuadOutside );

      typedef QuadratureContext< EntityType, IntersectionType, QuadratureImp > ContextType;
      typedef LocalEvaluation< ContextType, RangeEvalType, JacobianEvalType > LocalEvaluationType;

      ContextType cLeft( entity, intersection, faceQuadInside, entity.geometry().volume() );
      ContextType cRight( neighbor, intersection, faceQuadInside, neighbor.geometry().volume() );
      LocalEvaluationType localInside( cLeft, phiFaceEn_, dphiFaceEn_ );
      LocalEvaluationType localOutside( cRight, phiFaceNb_, dphiFaceNb_ );

      const auto& localInside0 = localInside(0);
      const auto& localOutside0 = localOutside(0);
      numericalFlux(localInside0, localOutside0, localInside0.values()[Id], localInside0.jacobians()[Id], localOutside0.values()[Id], localOutside0.jacobians()[Id],
                    rhsValueEn_, rhsDValueEn_, rhsValueNb_, rhsDValueNb_, true );

      const size_t numFaceQuadPoints = faceQuadInside.nop();

      // compute fluxes and assemble matrix
      for( int i = 0; i < numBasisFunctionsEn; ++i )
      {
        // compute flux for one base function, i.e.,
        // - uLeft=phiFaceEn_[.][i]
        // - uRight=0
        const auto& localInsidei = localInside[i];
        const auto& localOutsidei = localOutside[i];
        numericalFlux(localInsidei, localOutsidei, localInsidei.values()[Id], localInsidei.jacobians()[Id], localInside0.values()[Id], localInside0.jacobians()[Id],
                      valueEn_, dvalueEn_, valueNb_, dvalueNb_, true );


        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const double weight = faceQuadInside.weight( pt );
          valueEn_[pt] -= rhsValueEn_[pt];
          dvalueEn_[pt] -= rhsDValueEn_[pt];
          valueNb_[pt] -= rhsValueNb_[pt];
          dvalueNb_[pt] -= rhsDValueNb_[pt];
          localOpEn.matrixColumn( i ).axpy( std::get<size>(phiFaceEn_[pt]), std::get<size>(dphiFaceEn_[pt]),
                                      valueEn_[pt], dvalueEn_[pt], weight );
          localOpNb.matrixColumn( i ).axpy( std::get<size>(phiFaceNb_[pt]), std::get<size>(dphiFaceNb_[pt]),
                                      valueNb_[pt], dvalueNb_[pt], -weight );
        }
      }

      // assemble part from neighboring row
      if( oneSidedEvaluation )
      {
        auto nnGuard = Dune::Fem::bindGuard( localOpNbNb, neighbor, neighbor );
        auto neGuard = Dune::Fem::bindGuard( localOpNbEn, neighbor, entity   );
        for( int i = 0; i < numBasisFunctionsEn; ++i )
        {
          // compute flux for one base function, i.e.,
          // - uLeft=phiFaceEn_[.][i]
          // - uRight=0
          const auto& localInsidei = localInside[i];
          const auto& localOutsidei = localOutside[i];
          numericalFlux(localInsidei, localOutsidei, localOutside0.values()[Id], localOutside0.jacobians()[Id], localOutsidei.values()[Id], localOutsidei.jacobians()[Id],
                        valueEn_, dvalueEn_, valueNb_, dvalueNb_, true );


          for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
          {
            const double weight = faceQuadInside.weight( pt );
            valueEn_[pt] -= rhsValueEn_[pt];
            dvalueEn_[pt] -= rhsDValueEn_[pt];
            valueNb_[pt] -= rhsValueNb_[pt];
            dvalueNb_[pt] -= rhsDValueNb_[pt];
            localOpNbNb.matrixColumn( i ).axpy( std::get<size>(phiFaceNb_[pt]), std::get<size>(dphiFaceNb_[pt]),  // +
                                          valueNb_[pt], dvalueNb_[pt], -weight );
            localOpNbEn.matrixColumn( i ).axpy( std::get<size>(phiFaceEn_[pt]), std::get<size>(dphiFaceEn_[pt]),  // -
                                          valueEn_[pt], dvalueEn_[pt], weight );
          }
        }
      }

      // now move affine part to right hand side
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        const double weight = faceQuadInside.weight( pt );
        rhsValueEn_[pt] *= -weight;
        rhsDValueEn_[pt] *= -weight;
      }

      if( assembleRHS )
      {
        rhsLocal.axpyQuadrature( faceQuadInside, rhsValueEn_, rhsDValueEn_ );
      }
    }

    void testSymmetry() const
    {
      //if( matrix_ )
      //  testMatrixSymmetry( *matrix_ );
    }

    template <class LocalEvaluationVec, class Value,class DValue,class Value2,class DValue2,class RetType, class DRetType>
    void numericalFlux(const LocalEvaluationVec& left, const LocalEvaluationVec& right,
                       const Value &valueEn, const DValue &dvalueEn,
                       const Value2 &valueNb, const DValue2 &dvalueNb,
                       RetType &retEn, DRetType &dretEn,
                       RetType &retNb, DRetType &dretNb,
                       const bool initializeIntersection = true ) const
    {
      RangeType gLeft,gRight;
      if( hasDiffusion & initializeIntersection )
      {
        diffFlux_.initializeIntersection( left, right, valueEn, valueNb);
      }

      const size_t numFaceQuadPoints = left.quadrature().nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        if( hasDiffusion )
        {
          diffFlux_.numericalFlux( left[pt], right[pt],
                                   valueEn[pt], valueNb[pt], dvalueEn[pt], dvalueNb[pt],
                                   retEn[pt], retNb[pt], dretEn[pt], dretNb[pt]);
        }
        else
        {
          retEn[pt]  = RangeType(0);
          retNb[pt]  = RangeType(0);
          dretEn[pt] = JacobianRangeType(0);
          dretNb[pt] = JacobianRangeType(0);
        }

        advFlux_.numericalFlux(left[pt], right[pt],
                               valueEn[pt], valueNb[pt], dvalueEn[pt], dvalueNb[pt],
                               gLeft, gRight);
        retEn[pt] += gLeft;
        retNb[pt] += gRight;
      }
    }

    template <class LocalEvaluationVec,class Value,class DValue,class RetType, class DRetType>
    void fluxAndLift(const LocalEvaluationVec& left,
                     const LocalEvaluationVec& right,
                     const Value &valueEn, const DValue &dvalueEn,
                     const Value &valueNb, const DValue &dvalueNb,
                     RetType &retEn, DRetType &dretEn,
                     RetType &retNb, DRetType &dretNb,
                     DRetType &liftEn, DRetType &liftNb) const
    {
      numericalFlux(left, right,
                    valueEn, dvalueEn, valueNb, dvalueNb,
                    retEn, dretEn, retNb, dretNb);

      if( hasDiffusion )
      {
        const size_t numFaceQuadPoints = left.quadrature().nop();
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          diffFlux_.evaluateLifting(left[pt], right[pt], valueEn[pt], valueNb[pt],
                                    liftEn[pt], liftNb[pt]);
        }
      }
    }

    template <class LocalEvaluationVec, class Value,class LiftingFunction>
    void lifting(const LocalEvaluationVec& left, const LocalEvaluationVec& right,
                 const Value &valueEn, const Value &valueNb,
                 LiftingFunction &lifting) const
    {
      if( hasDiffusion )
      {
        diffFlux_.initializeIntersection( left, right, valueEn, valueNb, true );
        lifting += diffFlux_.getInsideLifting();
      }
    }

    template <class LocalEvaluationVec,class RetType>
    void boundaryValues(const LocalEvaluationVec& local,
                        RetType &bndValues) const
    {
      const RangeType uZero(0);
      const JacobianRangeType uJacZero( 0 );

      const size_t numFaceQuadPoints = local.quadrature().nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        model_.boundaryValue( local[pt], uZero, bndValues[pt]);
      }
    }

    template <class LocalEvaluationVec, class Value,class DValue,class GValue,class RetType, class DRetType>
    void boundaryFlux( const LocalEvaluationVec& local,
                       const Value &valueEn,
                       const DValue &dvalueEn,
                       const GValue &valueNb,
                       RetType &retEn,
                       DRetType &dretEn) const
    {
      const size_t numFaceQuadPoints = local.quadrature().nop();

      RangeType gLeft,gRight;
      if( hasDiffusion )
      {
        diffFlux_.initializeBoundary( local, valueEn, valueNb );
      }

      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        if ( model_.hasBoundaryValue( local[pt] ) )
        {
          if( hasDiffusion )
          {
            diffFlux_.boundaryFlux( local[pt],
                                    valueEn[pt], valueNb[pt], dvalueEn[pt],
                                    retEn[pt], dretEn[pt]);
          }
          else
          {
            retEn[pt]  = RangeType(0);
            dretEn[pt] = JacobianRangeType(0);
          }
          advFlux_.numericalFlux(local[pt], local[pt],
                                 valueEn[pt], valueNb[pt], dvalueEn[pt], dvalueEn[pt],
                                 gLeft,gRight);
          retEn[pt] += gLeft;
        }
        else
        {
          model_.boundaryFlux(local[pt], valueEn[pt], dvalueEn[pt], retEn[pt]);
          dretEn[pt] = 0;
        }
      }
    }

    const ModelType &model() const
    {
      return model_;
    }

  private:
    template< class QuadratureImp, class BasisFunctionSetImp >
    void setEntity ( const BasisFunctionSetImp& basisSet, const QuadratureImp &quad ) const
    {
      resizeEntity( quad.nop() );

      //evaluate all base functions
      for( const auto qp : quad )
      {
        basisSet.evaluateAll( qp, std::get<size>(phi_[qp.index()]) );
        basisSet.jacobianAll( qp, std::get<size>(dphi_[qp.index()]) );
      }
      extra( basisSet, quad );
    }

    //template< class QuadratureImp, class BasisFunctionSetImp >
    //void setNeighbor ( const BasisFunctionSetImp& basisSet, const QuadratureImp &quad )
    //{
    //  resize( quad.nop() );

    //  //evaluate all base functions
    //  for( const auto qp : quad )
    //  {
    //    basisSet.evaluateAll( qp, phi_[qp.index()] );
    //    basisSet.jacobianAll( qp, dphi_[qp.index()] );
    //  }
    //  extra( basisSet, quad );
    //}

    template <class QuadratureImp, class BasisFunctionSetImp >
    void initializeIntersection( const BasisFunctionSetImp basisSetInner,
                                 const BasisFunctionSetImp basisSetOuter,
                                 const QuadratureImp &quadInner,
                                 const QuadratureImp &quadOuter ) const
    {
      resizeIntersection( quadInner.nop() );

      //evaluate all base functions
      for( const auto qp : quadInner )
      {
        basisSetInner.evaluateAll( qp, std::get<size>(phiFaceEn_[qp.index()]) );
        basisSetInner.jacobianAll( qp, std::get<size>(dphiFaceEn_[qp.index()]) );

      }
      for( const auto qp : quadOuter )
      {
        basisSetOuter.evaluateAll( qp, std::get<size>(phiFaceNb_[qp.index()]) );
        basisSetOuter.jacobianAll( qp, std::get<size>(dphiFaceNb_[qp.index()]) );
      }
      extra( basisSetInner, quadInner );
      extra( basisSetOuter, quadOuter );
    }

    template <class QuadratureImp, class BasisFunctionSetImp >
    void initializeBoundary( const BasisFunctionSetImp basisSet,
                             const QuadratureImp &quadInner ) const
    {
      resizeIntersection( quadInner.nop() );

      //evaluate all base functions
      for( const auto qp : quadInner )
      {
        basisSet.evaluateAll( qp, std::get<size>(phiFaceEn_[qp.index()]) );
        basisSet.jacobianAll( qp, std::get<size>(dphiFaceEn_[qp.index()]) );
      }
      extra( basisSet, quadInner );
    }

    template< class QuadratureImp, class BasisFunctionSetImp >
    void extra( const BasisFunctionSetImp basisSet, const QuadratureImp &quad ) const
    {
      Dune::Hybrid::forEach(*extra_,
      [&](auto i)
      {
        std::cout << "reading extra: " << i << std::endl;
        for( const auto qp : quad )
        {
          std::get<i>( *extra_ ).localFunction( basisSet.entity() ).evaluateQuadrature( quad, std::get<i>(phiFaceEn_[qp.index()]) );
          std::get<i>( *extra_ ).localFunction( basisSet.entity() ).evaluateQuadrature( quad, std::get<i>(dphiFaceEn_[qp.index()]) );
        }
      } );
    }

    int elementQuadOrder( int polOrder ) const
    {
      return 2*polOrder;
    }
    int faceQuadOrder( int polOrder ) const
    {
      return 2*polOrder;
    }

    void resizeEntity( unsigned int numFaceQuadPoints ) const
    {
      if (phi_.size() >= numFaceQuadPoints)
        return;

      phi_.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) std::get<size>(phi_[i]).resize(maxNumBasisFunctions_);
      dphi_.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) std::get<size>(dphi_[i]).resize(maxNumBasisFunctions_);
    }

    void resizeIntersection( unsigned int numFaceQuadPoints ) const
    {
      if (phiFaceEn_.size() >= numFaceQuadPoints)
        return;

      for (unsigned int i=0;i<numFaceQuadPoints;++i) std::get<size>(dphi_[i]).resize(maxNumBasisFunctions_);
      phiFaceEn_.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) std::get<size>(phiFaceEn_[i]).resize(maxNumBasisFunctions_);
      dphiFaceEn_.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) std::get<size>(dphiFaceEn_[i]).resize(maxNumBasisFunctions_);
      phiFaceNb_.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) std::get<size>(phiFaceNb_[i]).resize(maxNumBasisFunctions_);
      dphiFaceNb_.resize( numFaceQuadPoints );
      for (unsigned int i=0;i<numFaceQuadPoints;++i) std::get<size>(dphiFaceNb_[i]).resize(maxNumBasisFunctions_);
      valueEn_.resize( numFaceQuadPoints );
      dvalueEn_.resize( numFaceQuadPoints );
      valueNb_.resize( numFaceQuadPoints );
      dvalueNb_.resize( numFaceQuadPoints );
      rhsValueEn_.resize( numFaceQuadPoints );
      rhsDValueEn_.resize( numFaceQuadPoints );
      rhsValueNb_.resize( numFaceQuadPoints );
      rhsDValueNb_.resize( numFaceQuadPoints );
      bndValues_.resize( numFaceQuadPoints );
    }

    const ModelType&                   model_;
    const DiscreteFunctionSpaceType&   space_;
    std::shared_ptr< DestinationType > rhs_;
    std::shared_ptr< MatrixType >      matrix_;
    std::shared_ptr< typename InsertFunctionTupleType::type > extra_;

    AdvectionFluxType                  advFlux_;
    mutable DiffusionFluxType          diffFlux_;
    const bool                         calculateFluxes_;
    const bool                         useStrongBoundaryCondition_;
    const size_t                       maxNumBasisFunctions_;

    // storage for all flux values
    mutable Dune::DynamicVector< RangeType >         valueEn_;
    mutable Dune::DynamicVector< JacobianRangeType > dvalueEn_;
    mutable Dune::DynamicVector< RangeType >         valueNb_;
    mutable Dune::DynamicVector< JacobianRangeType > dvalueNb_;
    mutable Dune::DynamicVector< RangeType >         rhsValueEn_;
    mutable Dune::DynamicVector< JacobianRangeType > rhsDValueEn_;
    mutable Dune::DynamicVector< RangeType >         rhsValueNb_;
    mutable Dune::DynamicVector< JacobianRangeType > rhsDValueNb_;
    mutable Dune::DynamicVector< RangeType >         bndValues_;

    // store all basis functions
    mutable RangeEvalType         phi_;
    mutable JacobianEvalType      dphi_;
    mutable RangeEvalType         phiFaceEn_;
    mutable JacobianEvalType      dphiFaceEn_;
    mutable RangeEvalType         phiFaceNb_;
    mutable JacobianEvalType      dphiFaceNb_;
  };

}
}

#endif
