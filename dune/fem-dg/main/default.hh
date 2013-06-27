#ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
#define DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH

#ifdef USE_BASEFUNCTIONSET_CODEGEN 
#define USE_BASEFUNCTIONSET_OPTIMIZED
#endif

// C++ includes
#include <cassert>
#include <cstddef>

// dune-common includes
#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/space/common/arrays.hh>

//#define NEWBASEFCT_CACHING

#ifdef NEWBASEFCT_CACHING
//#include <dune/fem/space/shapefunctionset/caching.hh>
#include "caching2.hh"
#else 
#include "caching.hh"
#endif

// dune-fem includes
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/basisfunctionset/transformation.hh>
#include <dune/fem/space/shapefunctionset/caching.hh>
#include <dune/fem/quadrature/cachingpointlist.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
#ifdef NEWBASEFCT_CACHING
#include "codegen2.hh"
#else
#include "codegen.hh"
#endif
#endif


#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
#include <dune/fem/space/basefunctions/evaluatecaller.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    // DefaultBasisFunctionSet
    // -----------------------

    /**
     * \brief implementation of a basis function set for given entity
     *
     * \tparam  Entity            entity type
     * \tparam  ShapeFunctionSet  shape function set
     *
     * \note ShapeFunctionSet must be a copyable object. For most 
     *       non-trivial implementations, you may want to use a 
     *       proxy, see file
\code
    <dune/fem/space/shapefunctionset/proxy.hh>
\endcode
     */
    template< class Entity, class ShapeFunctionSet >
    class DefaultBasisFunctionSet
    {
      typedef DefaultBasisFunctionSet< Entity, ShapeFunctionSet > ThisType;

    public:
      //! \brief entity type
      typedef Entity EntityType;
      //! \brief shape function set type
      typedef ShapeFunctionSet ShapeFunctionSetType;

    protected:
      typedef typename ShapeFunctionSetType::FunctionSpaceType   LocalFunctionSpaceType;
      typedef typename LocalFunctionSpaceType::JacobianRangeType LocalJacobianRangeType;
      typedef typename LocalFunctionSpaceType::HessianRangeType  LocalHessianRangeType;

      typedef typename LocalFunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename EntityType::Geometry GeometryType;
      typedef typename EntityType::Geometry Geometry ;

      typedef typename GeometryType::ctype ctype;

    public:
      //  slight misuse of struct ToLocalFunctionSpace!!!
      //! \brief type of function space
      typedef typename ToLocalFunctionSpace< LocalFunctionSpaceType, GeometryType::coorddimension >::Type FunctionSpaceType;

      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

      typedef typename ScalarFunctionSpaceType::RangeType ScalarRangeType;
      typedef typename ScalarFunctionSpaceType::JacobianRangeType ScalarJacobianRangeType;

      //! \brief type of reference element
      typedef Dune::ReferenceElement< ctype, GeometryType::coorddimension > ReferenceElementType;

      enum { dimDomain = FunctionSpaceType::dimDomain };
      enum { dimRange  = FunctionSpaceType::dimRange  };

#ifdef NEWBASEFCT_CACHING
      typedef std::vector< ScalarRangeType >          RangeVectorType;
      typedef std::vector< ScalarJacobianRangeType >  JacobianRangeVectorType;
#else
      typedef MutableArray< MutableArray< ScalarRangeType > >         RangeVectorType;
      typedef MutableArray< MutableArray< ScalarJacobianRangeType > > JacobianRangeVectorType;
#endif

      //! \brief constructor
      DefaultBasisFunctionSet ()
      : entity_( nullptr )
      {}

      //! \brief constructor
      DefaultBasisFunctionSet ( const EntityType &entity, const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : entity_( &entity ),
        shapeFunctionSet_( shapeFunctionSet )
      {}

      void registerEntry() const
      {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
        // add my dimrange 
        Fem::CodegenInfo::instance().addDimRange( this, dimRange );
#endif
      }

#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
#include "vectorialbasefunctionset.hh"
#endif

      // Basis Function Set Interface Methods
      // ------------------------------------

      //! \brief return order of basis function set
      int order () const { return shapeFunctionSet().order(); }

      //! \brief return size of basis function set
      std::size_t size () const { return shapeFunctionSet().size(); }

      //! \brief return size of basis function set
      std::size_t numDifferentBaseFunctions () const { return size()/dimRange; }

      //! \brief return reference element
      const ReferenceElementType &referenceElement () const
      {
        return Dune::ReferenceElements< ctype, GeometryType::coorddimension >::general( type() );
      }

      //! \brief evaluate all basis function and multiply with given values and add to dofs 
      template< class QuadratureType, class Vector, class DofVector >
      void axpy ( const QuadratureType &quad, const Vector &values, DofVector &dofs ) const
      {
        axpyImpl( quad, values, dofs, values[ 0 ] );
      }

      /** \brief evaluate all basis function and multiply with given values and add to dofs 
          \note valuesA and valuesB can be vectors of RangeType or JacobianRangeType 
      */
      template< class QuadratureType, class VectorA, class VectorB, class DofVector >
      void axpy ( const QuadratureType &quad, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const
      {
        assert( valuesA.size() > 0 );
        assert( valuesB.size() > 0 );

        axpyImpl( quad, valuesA, dofs, valuesA[ 0 ] );
        axpyImpl( quad, valuesB, dofs, valuesB[ 0 ] );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        FunctionalAxpyFunctor< RangeType, DofVector > f( valueFactor, dofs );
        shapeFunctionSet().evaluateEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        typedef typename GeometryType::JacobianInverseTransposed GeometryJacobianInverseTransposedType;
        const GeometryJacobianInverseTransposedType &gjit = geometry().jacobianInverseTransposed( coordinate( x ) );
        LocalJacobianRangeType tmpJacobianFactor;
        for( int r = 0; r < FunctionSpaceType::dimRange; ++r )
          gjit.mtv( jacobianFactor[ r ], tmpJacobianFactor[ r ] );

        FunctionalAxpyFunctor< LocalJacobianRangeType, DofVector > f( tmpJacobianFactor, dofs );
        shapeFunctionSet().jacobianEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      /** \copydoc BasisFunctionSet::evaluateAll( quad, dofs, ranges ) */
      template< class QuadratureType, class DofVector, class RangeArray >
      void evaluateAll ( const QuadratureType &quad, const DofVector &dofs, RangeArray &ranges ) const
      {
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
        typedef Fem :: EvaluateCallerInterfaceTraits< 
                QuadratureType, RangeArray, DofVector > Traits;
        typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

        // get base function evaluate caller (calls evaluateRanges) 
        const BaseEvaluationType& baseEval = 
            BaseEvaluationType::storage( *this, rangeCache( quad ), quad );

        baseEval.evaluateRanges( quad, dofs, ranges );
#else 

        registerEntry();
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
        Fem::CodegenInfo::instance().addEntry( "evalranges", 
            Fem :: CodeGeneratorType :: evaluateCodegen, dimDomain, dimRange, quad.nop(), size()/dimRange );
#endif
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          evaluateAll( quad[ qp ], dofs, ranges[ qp ] );
        }
#endif // #ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        value = RangeType( 0 );
        AxpyFunctor< DofVector, RangeType > f( dofs, value );
        shapeFunctionSet().evaluateEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        AssignFunctor< RangeArray > f( values );
        shapeFunctionSet().evaluateEach( x, f );
      }

      /** \copydoc BasisFunctionSet::jacobianAll( quad, dofs, jacobians ) */
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quad, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        assert( jacobians.size() > 0 );
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
        typedef Fem :: EvaluateCallerInterfaceTraits< QuadratureType, 
                JacobianArray, DofVector, Geometry >  Traits;
        typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

        // get base function evaluate caller (calls axpyRanges) 
        const BaseEvaluationType& baseEval = 
          BaseEvaluationType::storage( *this, jacobianCache( quad ), quad );

        // call appropriate axpyJacobian method 
        baseEval.evaluateJacobians( quad, geometry(), dofs, jacobians );
#else 
        registerEntry();

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
        Fem::CodegenInfo::instance().addEntry( "evaljacobians", 
              Fem :: CodeGeneratorType :: evaluateJacobiansCodegen, dimDomain, dimRange, quad.nop(), size()/dimRange );
#endif
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          jacobianAll( quad[ qp ], dofs, jacobians[ qp ] );
        }
#endif  // #ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        LocalJacobianRangeType localJacobian( RangeFieldType( 0 ) );
        AxpyFunctor< DofVector, LocalJacobianRangeType > f( dofs, localJacobian );
        shapeFunctionSet().jacobianEach( x, f );

        typedef JacobianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        transformation( localJacobian, jacobian );
      }

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        typedef JacobianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        AssignFunctor< JacobianRangeArray, Transformation > f( jacobians, transformation );
        shapeFunctionSet().jacobianEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        LocalHessianRangeType localHessian( typename LocalHessianRangeType::value_type( RangeFieldType( 0 ) ) );
        AxpyFunctor< DofVector, LocalHessianRangeType > f( dofs, localHessian );
        shapeFunctionSet().hessianEach( x, f );

        typedef HessianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        transformation( localHessian, hessian );
      }

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        typedef HessianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        AssignFunctor< HessianRangeArray, Transformation > f( hessians, transformation );
        shapeFunctionSet().hessianEach( x, f );
      }

      //! \brief return entity
      const Entity &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      //! \brief return geometry type
      Dune::GeometryType type () const { return entity().type(); }

      // Non-interface methods
      // ---------------------

      //! \brief return shape function set
      const ShapeFunctionSetType &shapeFunctionSet () const { return shapeFunctionSet_; }

    protected:
      //! \brief evaluate all basis function and multiply with given values and add to dofs 
      template< class QuadratureType, class RangeArray, class DofVector >
      void axpyImpl ( const QuadratureType &quad, const RangeArray &rangeFactors, DofVector &dofs, const RangeType& ) const
      {
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
        typedef Fem :: EvaluateCallerInterfaceTraits< 
            QuadratureType, RangeArray, DofVector > Traits;
        typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

        // get base function evaluate caller (calls axpyRanges) 
        const BaseEvaluationType& baseEval = 
          BaseEvaluationType::storage( *this, rangeCache( quad ), quad );

        // call appropriate axpyRanges method 
        baseEval.axpyRanges( quad, rangeFactors, dofs );
#else 

        registerEntry();
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
        Fem::CodegenInfo::instance().addEntry( "axpyranges", 
              Fem :: CodeGeneratorType :: axpyCodegen, dimDomain, dimRange, quad.nop(), size()/dimRange );
#endif
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          axpy( quad[ qp ], rangeFactors[ qp ], dofs );
        }
#endif  // #ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      }

      //! \brief evaluate all basis function and multiply with given values and add to dofs 
      template< class QuadratureType, class JacobianArray, class DofVector >
      void axpyImpl ( const QuadratureType &quad, const JacobianArray &jacobianFactors, DofVector &dofs, const JacobianRangeType& ) const
      {
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
        typedef Fem :: EvaluateCallerInterfaceTraits< QuadratureType,
                JacobianArray, DofVector, Geometry >  Traits;
        typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

        // get base function evaluate caller (calls axpyRanges) 
        const BaseEvaluationType& baseEval = 
          BaseEvaluationType::storage( *this, jacobianCache( quad ), quad );

        // call appropriate axpyRanges method 
        baseEval.axpyJacobians( quad, geometry(), jacobianFactors, dofs );
#else 

        registerEntry();
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
        Fem::CodegenInfo::instance().addEntry( "axpyjacobians", 
                Fem :: CodeGeneratorType :: axpyJacobianCodegen, dimDomain, dimRange, quad.nop(), size()/dimRange ); 
#endif
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          axpy( quad[ qp ], jacobianFactors[ qp ], dofs );
        }
#endif  // #ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      }

      GeometryType geometry () const { return entity().geometry(); }

      template <class QuadratureType>
      const RangeVectorType& rangeCache( const QuadratureType& quad ) const 
      { 
        return shapeFunctionSet().scalarShapeFunctionSet().impl().rangeCache( quad );
      }

      template <class QuadratureType>
      const JacobianRangeVectorType& jacobianCache( const QuadratureType& quad ) const 
      { 
        return shapeFunctionSet().scalarShapeFunctionSet().impl().jacobianCache( quad );
      }
    private:
      const EntityType *entity_;
      ShapeFunctionSetType shapeFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
