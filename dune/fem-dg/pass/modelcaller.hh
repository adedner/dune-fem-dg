#ifndef DUNE_FEM_DG_PASS_LOCALDG_MODELCALLER_HH
#define DUNE_FEM_DG_PASS_LOCALDG_MODELCALLER_HH

#ifdef DUNE_FEM_PASS_COMMON_SELECTOR_HH
#error "selector.hh from dune/fem/pass/common is outdated!"
#endif

#include <cstddef>
#include <vector>

#include <dune/fem/common/tupletypetraits.hh>
#include <dune/fem/common/tupleutility.hh>
#include <dune/fem/common/typeindexedtuple.hh>

#include <dune/fem-dg/pass/selector.hh>
#include <dune/fem-dg/pass/localfunctiontuple.hh>
#include <dune/fem-dg/pass/pass.hh>
#include <dune/fem-dg/pass/context.hh>

namespace Dune
{

  namespace Fem
  {

    // DGDiscreteModelCaller
    // ---------------------

    /** \brief   model caller for local DG pass
     *  \class   DGDiscreteModelCaller
     *  \ingroup PassHyp
     *
     *  \todo please doc me
     */
    template< class DiscreteModel, class Argument, class PassIds >
    class DGDiscreteModelCaller
    {
      typedef DGDiscreteModelCaller< DiscreteModel, Argument, PassIds > ThisType;

    public:
      // discrete model type
      typedef DiscreteModel DiscreteModelType;
      // argument type
      typedef Argument ArgumentType;

      /** \brief selector (tuple of integral constants) */
      typedef typename DiscreteModelType::Selector Selector;

      // entity type
      typedef typename DiscreteModelType::EntityType EntityType;
      // intersection type
      typedef typename DiscreteModelType::IntersectionType IntersectionType;

      // type of volume quadrature
      typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
      // type of face quadrature
      typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;

      // type of mass factor
      typedef typename DiscreteModelType::MassFactorType MassFactorType;

      // function space type
      typedef typename DiscreteModelType::FunctionSpaceType FunctionSpaceType;
      // range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      // jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    private:
      typedef Dune::MakeSubTuple< ArgumentType, typename Dune::FirstTypeIndexTuple< PassIds, Selector >::type > FilterType;
      typedef typename FilterType::type DiscreteFunctionPointerTupleType;

      typedef typename TupleTypeTraits< DiscreteFunctionPointerTupleType >::PointeeTupleType DiscreteFunctionTupleType;

    protected:
      typedef LocalFunctionTuple< DiscreteFunctionTupleType, EntityType > LocalFunctionTupleType;

      typedef typename LocalFunctionTupleType::RangeTupleType RangeTupleType;
      typedef typename LocalFunctionTupleType::JacobianRangeTupleType JacobianRangeTupleType;

    public:
      DGDiscreteModelCaller ( ArgumentType &argument, DiscreteModelType &discreteModel )
      : discreteModel_( discreteModel ),
        time_( 0. ),
        discreteFunctions_( FilterType::apply( argument ) ),
        localFunctionsInside_( Dune::DereferenceTuple< DiscreteFunctionPointerTupleType >::apply( discreteFunctions_ ) ),
        localFunctionsOutside_( Dune::DereferenceTuple< DiscreteFunctionPointerTupleType >::apply( discreteFunctions_ ) )
      {}

      // return true, if discrete model has flux
      bool hasFlux () const { return discreteModel().hasFlux(); }
      // return true, if discrete model has mass
      bool hasMass () const { return discreteModel().hasMass(); }
      // return true, if discrete model has source
      bool hasSource () const { return discreteModel().hasSource(); }

      // set time
      void setTime ( double time ) { time_ = time; }
      // return time
      double time () const { return time_; }

      // set inside local functions to entity
      void setEntity ( const EntityType &entity )
      {
        localFunctionsInside_.init( entity );
        discreteModel().setEntity( entity );
      }

      // set outside local functions to entity
      void setNeighbor ( const EntityType &entity )
      {
        localFunctionsOutside_.init( entity );
        discreteModel().setNeighbor( entity );
      }

      // evaluate inside local functions in all quadrature points
      void setEntity ( const EntityType &entity, const VolumeQuadratureType &quadrature )
      {
        setEntity( entity );
        values_.resize( quadrature.nop() );
        jacobians_.resize( quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, values_ );
        localFunctionsInside_.evaluateQuadrature( quadrature, jacobians_ );
      }

      // evaluate outside local functions in all quadrature points
      template< class QuadratureType >
      void setNeighbor ( const EntityType &neighbor,
                         const QuadratureType &inside,
                         const QuadratureType &outside )
      {
        // we assume that setEntity() was called in advance!
        valuesInside_.resize( inside.nop() );
        localFunctionsInside_.evaluateQuadrature( inside, valuesInside_ );

        setNeighbor( neighbor );

        valuesOutside_.resize( outside.nop() );
        localFunctionsOutside_.evaluateQuadrature( outside, valuesOutside_ );
      }

      // please doc me
      template< class QuadratureType >
      void setBoundary ( const EntityType &entity, const QuadratureType &quadrature )
      {
        valuesInside_.resize( quadrature.nop() );
        localFunctionsInside_.evaluateQuadrature( quadrature, valuesInside_ );
      }

      // evaluate analytical flux
      void analyticalFlux ( const EntityType &entity,
                            const VolumeQuadratureType &quadrature,
                            const int qp,
                            JacobianRangeType &flux )
      {
        assert( hasFlux() );
        discreteModel().analyticalFlux( entity, time(), quadrature.point( qp ), values_[ qp ], flux );
      }

      // evaluate analytical flux
      double source ( const EntityType &entity,
                      const VolumeQuadratureType &quadrature,
                      const int qp,
                      RangeType &source )
      {
        assert( hasSource() );
        return discreteModel().source( entity, time(), quadrature.point( qp ), values_[ qp ], jacobians_[ qp ], source );
      }


      // evaluate both analytical flux and source
      double analyticalFluxAndSource ( const EntityType &entity,
                                       const VolumeQuadratureType &quadrature,
                                       const int qp,
                                       JacobianRangeType &flux,
                                       RangeType &source )
      {
        // we may only assume that hasSource() == true, cf. pass.hh
        if( hasFlux() )
          analyticalFlux( entity, quadrature, qp, flux );
        return ThisType::source( entity, quadrature, qp, source );
      }

      // evaluate numerical flux
      template< class QuadratureType >
      double numericalFlux ( const IntersectionType &intersection,
                             const QuadratureType &inside,
                             const QuadratureType &outside,
                             const int qp,
                             RangeType &gLeft, RangeType &gRight )
      {
        return discreteModel().numericalFlux( intersection, time(), inside.localPoint( qp ), valuesInside_[ qp ], valuesOutside_[ qp ], gLeft, gRight );
      }

      // evaluate boundary flux
      double boundaryFlux ( const IntersectionType &intersection,
                            const FaceQuadratureType &quadrature,
                            const int qp,
                            RangeType &gLeft )
      {
        return discreteModel().boundaryFlux( intersection, time(), quadrature.localPoint( qp ), valuesInside_[ qp ], gLeft );
      }

      // evaluate mass
      void mass ( const EntityType &entity,
                  const VolumeQuadratureType &quadrature,
                  const int qp,
                  MassFactorType &massFactor )
      {
        discreteModel().mass( entity, time(), quadrature.point( qp ), values_[ qp ], massFactor );
      }

      DGDiscreteModelCaller ( const ThisType & ) = delete;
      ThisType operator= ( const ThisType & ) = delete;

    protected:
      DiscreteModelType &discreteModel () { return discreteModel_; }
      const DiscreteModelType &discreteModel () const { return discreteModel_; }

    private:
      DiscreteModelType &discreteModel_;
      double time_;
      DiscreteFunctionPointerTupleType discreteFunctions_;

    protected:
      LocalFunctionTupleType localFunctionsInside_, localFunctionsOutside_;
      std::vector< Dune::TypeIndexedTuple< RangeTupleType, Selector > > values_, valuesInside_, valuesOutside_;
      std::vector< Dune::TypeIndexedTuple< JacobianRangeTupleType, Selector > > jacobians_;
    };




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
        discreteModel().analyticalFlux( EvalType( cLocal, values_[qp], jacobianValue( jacobians_, qp ), qp ), flux );
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
        return discreteModel().source( EvalType( cLocal, values_[qp], jacobianValue( jacobians_, qp ), qp ), src );
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
                    EvalType( cLeft, valuesInside_[qp],  jacobianValue( jacobiansInside_, qp ), qp ),
                    EvalType( cRight, valuesOutside_[qp], jacobianValue( jacobiansOutside_, qp ), qp ),
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

        ContextType cLocal( localFunctionsInside_.entity(), intersection, quadrature, discreteModel().enVolume() );
        return discreteModel().boundaryFlux( EvalType( cLocal, valuesInside_[qp], jacobianValue( jacobiansInside_, qp ), qp ), gLeft, hLeft );
      }

      void mass ( const EntityType &entity,
                  const VolumeQuadratureType &quadrature,
                  const int qp,
                  MassFactorType &m )
      {
        typedef QuadratureContext< EntityType, VolumeQuadratureType > ContextType;
        typedef LocalEval< ContextType > EvalType;

        ContextType cLocal( entity, quadrature, discreteModel().enVolume() );
        discreteModel().mass( EvalType( cLocal, values_[qp], jacobianValue( jacobians_, qp ), qp ), m );
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

  } // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_DG_FEM_PASS_LOCALDG_MODELCALLER_HH
