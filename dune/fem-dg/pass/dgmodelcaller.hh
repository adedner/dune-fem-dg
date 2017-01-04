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
namespace Fem
{

  // CDGDiscreteModelCaller
  // ----------------------

  /**
   * \brief Model caller for CDG pass.
   *
   * \ingroup PassBased
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

    typedef typename BaseType::MassFactorType MatrixMassFactorType;
    struct DiagonalMassFactor : public RangeType
    {
      typedef RangeType  BaseType;
      DiagonalMassFactor() : BaseType() {}
      template <class T>
      DiagonalMassFactor( const T& other ) : BaseType( other ) {}
      DiagonalMassFactor( const DiagonalMassFactor& other ) : BaseType( other ) {}

      //! multiply method needed in LocalMassMatrix
      void mv( const RangeType& arg, RangeType& dest ) const
      {
        for( int i=0; i<RangeType::dimension; ++i )
        {
          dest[ i ] = arg[ i ] * this->operator[]( i );
        }
      }
    };

    //typedef std::conditional< DiscreteModelType::scalarMassFactor,
    //          DiagonalMassFactor,
    //          MatrixMassFactorType > :: type MassFactorType;

    typedef DiagonalMassFactor MassFactorType;

  protected:
    typedef typename BaseType::RangeTupleType RangeTupleType;
    typedef typename BaseType::JacobianRangeTupleType JacobianRangeTupleType;

    template< class ContextImp >
    using LocalEval = LocalEvaluation< ContextImp,
                                       Dune::TypeIndexedTuple< RangeTupleType, Selector >,
                                       Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > >;

    template< class ContextImp >
    using LocalEvalVec = LocalEvaluation< ContextImp,
                                       std::vector< Dune::TypeIndexedTuple< RangeTupleType, Selector > >,
                                       std::vector< Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > > >;
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
                                 const QuadratureImp &quadInner,
                                 const QuadratureImp &quadOuter )
    {
      assert( intersection.neighbor() );

      BaseType::setNeighbor( neighbor, quadInner, quadOuter );

      if( evaluateJacobian )
      {
        jacobiansOutside_.resize( quadOuter.nop() );
        localFunctionsOutside_.evaluateQuadrature( quadOuter, jacobiansOutside_ );
        jacobiansInside_.resize( quadInner.nop() );
        localFunctionsInside_.evaluateQuadrature( quadInner, jacobiansInside_ );
      }

#ifndef NDEBUG
      quadInnerId_ = quadInner.id();
      quadOuterId_ = quadOuter.id();
#endif
      typedef QuadratureContext< EntityType, IntersectionType, QuadratureImp > ContextType;
      typedef LocalEvalVec< ContextType > EvalType;

      ContextType cLeft( discreteModel().inside(), intersection, quadInner, discreteModel().enVolume() );
      ContextType cRight( discreteModel().outside(), intersection, quadOuter, discreteModel().nbVolume() );
      discreteModel().initializeIntersection( EvalType( cLeft,  valuesInside_, jacobiansInside_ ),
                                              EvalType( cRight, valuesOutside_, jacobiansOutside_ ) );
    }

    template <class QuadratureImp>
    void initializeBoundary ( const EntityType& inside,
                              const IntersectionType &intersection,
                              const QuadratureImp &quadrature )
    {
      assert( intersection.boundary() );

      BaseType::setBoundary( inside, quadrature );
      if( evaluateJacobian )
      {
        jacobiansInside_.resize( quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, jacobiansInside_ );
      }

#ifndef NDEBUG
      quadInnerId_ = quadrature.id();
#endif
      typedef QuadratureContext< EntityType, IntersectionType, QuadratureImp > ContextType;
      typedef LocalEvalVec< ContextType > EvalType;


      ContextType cLocal( inside, intersection, quadrature, discreteModel().enVolume() );
      discreteModel().initializeBoundary( EvalType( cLocal,  valuesInside_, jacobiansInside_ ) );
    }

    template <class QuadratureImp>
    void initializeBoundary ( const IntersectionType &intersection,
                              const QuadratureImp &quadrature )
    {
      initializeBoundary( discreteModel().inside(), intersection, quadrature );
    }

    void analyticalFlux ( const EntityType &entity,
                          const VolumeQuadratureType &quadrature,
                          const int qp,
                          JacobianRangeType &flux )
    {
      assert( quadId_ == quadrature.id() );
      assert( (int) values_.size() > qp );

      typedef QuadratureContext< EntityType, VolumeQuadratureType > ContextType;
      typedef LocalEval< ContextType > EvalType;

      ContextType cLocal( entity, quadrature, discreteModel().enVolume() );
      discreteModel().analyticalFlux( EvalType( cLocal[qp], values_[qp], jacobianValue( jacobians_, qp ) ), flux );
    }

    double source ( const EntityType &entity,
                    const VolumeQuadratureType &quadrature,
                    const int qp,
                    RangeType &src )
    {
      assert( quadId_ == quadrature.id() );
      assert( (int) values_.size() > qp );

      typedef QuadratureContext< EntityType, VolumeQuadratureType > ContextType;
      typedef LocalEval< ContextType > EvalType;

      ContextType cLocal( entity, quadrature, discreteModel().enVolume() );
      return discreteModel().source( EvalType( cLocal[qp], values_[qp], jacobianValue( jacobians_, qp ) ), src );
    }

    double analyticalFluxAndSource( const EntityType &entity,
                                    const VolumeQuadratureType &quadrature,
                                    const int qp,
                                    JacobianRangeType &flux,
                                    RangeType &src )
    {
      analyticalFlux( entity, quadrature, qp, flux );
      return source( entity, quadrature, qp, src );
    }


    template <class QuadratureType>
    double numericalFlux ( const IntersectionType &intersection,
                           const QuadratureType &quadInner,
                           const QuadratureType &quadOuter,
                           const int qp,
                           RangeType &gLeft,
                           RangeType &gRight,
                           JacobianRangeType &hLeft,
                           JacobianRangeType &hRight)
    {
      assert( valuesInside_.size() >= quadInner.nop() );
      assert( quadInnerId_ == quadInner.id() );
      assert( valuesOutside_.size() >= quadOuter.nop() );
      assert( quadOuterId_ == quadOuter.id() );

      typedef QuadratureContext< EntityType, IntersectionType, QuadratureType > ContextType;
      typedef LocalEval< ContextType > EvalType;

      ContextType cLeft( localFunctionsInside_.entity(), intersection, quadInner, discreteModel().enVolume() );
      ContextType cRight( localFunctionsOutside_.entity(), intersection, quadOuter, discreteModel().nbVolume() );
      return discreteModel().numericalFlux(
                  EvalType( cLeft[qp], valuesInside_[qp],  jacobianValue( jacobiansInside_, qp ) ),
                  EvalType( cRight[qp], valuesOutside_[qp], jacobianValue( jacobiansOutside_, qp ) ),
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

      typedef QuadratureContext< EntityType, IntersectionType, FaceQuadratureType > ContextType;
      typedef LocalEval< ContextType > EvalType;

      ContextType cLocal( localFunctionsInside_.entity(), intersection, quadrature, discreteModel().enVolume(), qp );
      return discreteModel().boundaryFlux( EvalType( cLocal, valuesInside_[qp], jacobianValue( jacobiansInside_, qp ) ), gLeft, hLeft );
    }

    void mass ( const EntityType &entity,
                const VolumeQuadratureType &quadrature,
                const int qp,
                MassFactorType &m )
    {
      typedef QuadratureContext< EntityType, VolumeQuadratureType > ContextType;
      typedef LocalEval< ContextType > EvalType;

      ContextType cLocal( entity, quadrature, discreteModel().enVolume(), qp );
      discreteModel().mass( EvalType( cLocal, values_[ qp ], jacobianValue( jacobians_, qp ) ), m );
    }
  protected:
    template< class JacobianRangeTupleVectorType >
    const typename JacobianRangeTupleVectorType::value_type &jacobianValue ( const JacobianRangeTupleVectorType &jacobians, const int qp ) const
    {
      assert( ( evaluateJacobian ) ? (int) jacobians.size() > qp : true );
      return ( evaluateJacobian ) ? jacobians[qp] : jacobians[0];
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

} // namespace
} // namespace Dune

#endif // #ifndef DUNE_FEM_DG_DISCRETEMODELCALLER_HH
