#ifndef DUNE_FEM_DG_DISCRETEMODELCALLER_HH
#define DUNE_FEM_DG_DISCRETEMODELCALLER_HH

#include <cassert>
#include <vector>

#include <dune/fem/pass/localdg/modelcaller.hh>
#include <dune/fem/version.hh>
#include <dune/fem/misc/compatibility.hh>
#include <dune/fem-dg/pass/context.hh>

namespace Dune
{

  // CDGDiscreteModelCaller
  // ----------------------

  /**
   * \brief Model caller for CDG pass.
   */
  template< class DiscreteModel, class Argument, class PassIds >
  class CDGDiscreteModelCaller
  : public Dune::Fem::DGDiscreteModelCaller< DiscreteModel, Argument, PassIds >
  {
    typedef CDGDiscreteModelCaller< DiscreteModel, Argument, PassIds > ThisType;
    typedef Dune::Fem::DGDiscreteModelCaller< DiscreteModel, Argument, PassIds > BaseType;

  public:
    typedef typename BaseType::DiscreteModelType DiscreteModelType;
    typedef typename BaseType::ArgumentType ArgumentType;

    typedef typename BaseType::Selector Selector;

    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    typedef typename BaseType::EntityType EntityType;
    typedef typename BaseType::IntersectionType IntersectionType;

    typedef typename BaseType::VolumeQuadratureType VolumeQuadratureType;
    typedef typename BaseType::FaceQuadratureType FaceQuadratureType;

    typedef typename BaseType::MassFactorType MassFactorType;

  protected:
    typedef typename BaseType::RangeTupleType RangeTupleType;
    typedef typename BaseType::JacobianRangeTupleType JacobianRangeTupleType;

      typedef  ElementQuadraturePointContext< EntityType, VolumeQuadratureType,
                     Dune::TypeIndexedTuple< RangeTupleType, Selector >,
                     Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > > ElementQuadratureContextType ;

  public:
    static const bool evaluateJacobian = DiscreteModelType::evaluateJacobian;

    using BaseType::setEntity;
    using BaseType::setNeighbor;
    using BaseType::time;

    CDGDiscreteModelCaller ( ArgumentType &argument, DiscreteModelType &discreteModel )
    : BaseType( argument, discreteModel )
#ifndef NDEBUG
      , quadInnerId_( 0 )
      , quadOuterId_( 0 )
      , quadId_( 0 )
#endif
    {}

    void setEntity ( const EntityType &entity, const VolumeQuadratureType &quadrature )
    {
      BaseType::setEntity( entity, quadrature );

      if( discreteModel().hasSource () || evaluateJacobian )
      {
        jacobians_.resize( quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, jacobians_ );
      }

#ifndef NDEBUG
      quadId_ = quadrature.id();
#endif
    }

    // Ensure: entities set correctly before call
    template <class QuadratureImp>
    void initializeIntersection( const EntityType &neighbor,
                                 const IntersectionType &intersection,
                                 const QuadratureImp &inside,
                                 const QuadratureImp &outside )
    {
      assert( intersection.neighbor() );

      BaseType::setNeighbor( neighbor, inside, outside );

      if( evaluateJacobian )
      {
        jacobiansInside_.resize( inside.nop() );
        localFunctionsInside_.evaluateQuadrature( inside, jacobiansInside_ );
        jacobiansOutside_.resize( outside.nop() );
        localFunctionsOutside_.evaluateQuadrature( outside, jacobiansOutside_ );
      }

#ifndef NDEBUG
      quadInnerId_ = inside.id();
      quadOuterId_ = outside.id();
#endif

      discreteModel().initializeIntersection( intersection, time(), inside, outside, valuesInside_, valuesOutside_ );
    }

    template <class QuadratureImp>
    void initializeBoundary ( const IntersectionType &intersection,
                              const QuadratureImp &quadrature )
    {
      assert( intersection.boundary() );

      BaseType::setBoundary( Fem::make_entity( intersection.inside() ), quadrature );
      if( evaluateJacobian )
      {
        jacobiansInside_.resize( quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, jacobiansInside_ );
      }

#ifndef NDEBUG
      quadInnerId_ = quadrature.id();
#endif

      discreteModel().initializeBoundary( intersection, time(), quadrature, valuesInside_ );
    }

    void analyticalFlux ( const EntityType &entity,
                          const VolumeQuadratureType &quadrature,
                          const int qp,
                          JacobianRangeType &flux )
    {
      assert( quadId_ == quadrature.id() );
      assert( (int) values_.size() > qp );

      ElementQuadratureContextType context( entity, quadrature, values_[ qp ], jacobianValue( jacobians_, qp ), qp, time(), discreteModel().enVolume() );
      discreteModel().analyticalFlux( context, flux );
    }

    double source ( const EntityType &entity,
                    const VolumeQuadratureType &quadrature,
                    const int qp,
                    RangeType &source )
    {
      assert( quadId_ == quadrature.id() );
      assert( (int) values_.size() > qp );

      ElementQuadratureContextType context( entity, quadrature, values_[ qp ], jacobianValue( jacobians_, qp ), qp, time(), discreteModel().enVolume() );
      return discreteModel().source( context, source );
    }

    double analyticalFluxAndSource( const EntityType &entity,
                                    const VolumeQuadratureType &quadrature,
                                    const int qp,
                                    JacobianRangeType &flux,
                                    RangeType &source )
    {
      ElementQuadratureContextType context( entity, quadrature, values_[ qp ], jacobianValue( jacobians_, qp ), qp, time(), discreteModel().enVolume() );
      discreteModel().analyticalFlux( context, flux );
      return discreteModel().source( context, source );
    }


    template <class QuadratureType>
    double numericalFlux ( const IntersectionType &intersection,
                           const QuadratureType &inside,
                           const QuadratureType &outside,
                           const int qp,
                           RangeType &gLeft,
                           RangeType &gRight,
                           JacobianRangeType &hLeft,
                           JacobianRangeType &hRight)
    {
      assert( valuesInside_.size() >= inside.nop() );
      assert( quadInnerId_ == inside.id() );
      assert( valuesOutside_.size() >= inside.nop() );
      assert( quadOuterId_ == outside.id() );

      typedef IntersectionQuadraturePointContext< IntersectionType, EntityType, QuadratureType,
                     Dune::TypeIndexedTuple< RangeTupleType, Selector >,
                     Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > > QuadratureContextType ;

      QuadratureContextType ctxLeft ( intersection, localFunctionsInside_.entity(),  inside,  valuesInside_[ qp ],  jacobianValue( jacobiansInside_, qp ),  qp, time(), discreteModel().enVolume() );
      QuadratureContextType ctxRight( intersection, localFunctionsOutside_.entity(), outside, valuesOutside_[ qp ], jacobianValue( jacobiansOutside_, qp ), qp, time(), discreteModel().nbVolume() );

      return discreteModel().numericalFlux( ctxLeft, ctxRight,
                                            gLeft, gRight, hLeft, hRight );
    }

    double boundaryFlux ( const IntersectionType &intersection,
                          const FaceQuadratureType &quadrature,
                          const int qp,
                          RangeType &gLeft,
                          JacobianRangeType &hLeft )
    {
      assert( valuesInside_.size() >= quadrature.nop() );
      assert( quadInnerId_ == quadrature.id() );

      typedef  IntersectionQuadraturePointContext< IntersectionType, EntityType, FaceQuadratureType,
                     Dune::TypeIndexedTuple< RangeTupleType, Selector >,
                     Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > > QuadratureContextType ;

      QuadratureContextType ctxLeft ( intersection, localFunctionsInside_.entity(), quadrature, valuesInside_[ qp ],  jacobianValue( jacobiansInside_, qp ),  qp, time(), discreteModel().enVolume() );

      return discreteModel().boundaryFlux( ctxLeft, gLeft, hLeft );
    }

  protected:
    template< class JacobianRangeTupleVectorType >
    const typename JacobianRangeTupleVectorType::value_type &jacobianValue ( const JacobianRangeTupleVectorType &jacobians, const int qp ) const
    {
      assert( ( evaluateJacobian ) ? (int) jacobians.size() > qp : true );
      return ( evaluateJacobian ) ? jacobians[ qp ] : jacobians[ 0 ];
    }

    using BaseType::discreteModel;
    using BaseType::jacobians_;
    using BaseType::localFunctionsInside_;
    using BaseType::localFunctionsOutside_;
    using BaseType::values_;
    using BaseType::valuesInside_;
    using BaseType::valuesOutside_;

    std::vector< Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > > jacobiansInside_, jacobiansOutside_;

  private:
#ifndef NDEBUG
    size_t quadInnerId_;
    size_t quadOuterId_;
    size_t quadId_;
#endif
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_DG_DISCRETEMODELCALLER_HH
