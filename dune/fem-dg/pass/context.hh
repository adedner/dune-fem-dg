#ifndef DUNE_FEM_DG_QUADRATURECONTEXT_HH
#define DUNE_FEM_DG_QUADRATURECONTEXT_HH

#include <cassert>
#include <vector>
#include <dune/fem/common/typeindexedtuple.hh>


namespace Dune
{
namespace Fem
{

  /**
   * \brief Container class describing (many) information which are related to
   * the local evaluation of a discrete function.
   *
   * \note This class was introduced to allow more flexible interfaces for classes
   * needing a local evaluation.
   */
  template <class Entity,
            class Quadrature,
            class RangeTuple,
            class JacobianTuple>
  class ElementQuadraturePointContext
  {
  protected:
    const Entity& entity_;
    const Quadrature& quad_;
    const RangeTuple& values_;
    const JacobianTuple& jacobians_;

    const double time_;
    const double volume_;
    const int qp_;
  public:
    typedef Entity          EntityType;
    typedef Quadrature      QuadratureType;
    typedef RangeTuple      RangeTupleType;
    typedef JacobianTuple   JacobianTupleType;
    typedef typename QuadratureType :: QuadraturePointWrapperType   QuadraturePointWrapperType;
    typedef typename QuadratureType :: CoordinateType               CoordinateType;
    typedef typename QuadratureType :: LocalCoordinateType          LocalCoordinateType;

    ElementQuadraturePointContext( const Entity& entity,
                                   const Quadrature& quadrature,
                                   const RangeTuple& values,
                                   const JacobianTuple& jacobians,
                                   const int qp,
                                   const double time,
                                   const double volume )
     : entity_( entity ),
       quad_( quadrature ),
       values_( values ),
       jacobians_( jacobians ),
       time_( time ),
       volume_( volume ),
       qp_( qp )
    {}

    const Entity& entity() const { return entity_; }
    const Quadrature& quadrature() const { return quad_; }
    const RangeTuple& values() const { return values_; }
    const JacobianTuple& jacobians() const { return jacobians_; }
    const double time () const { return time_; }
    const double volume() const { return volume_; }
    const QuadraturePointWrapperType quadraturePoint() const { return quadrature()[ index() ]; }
    const CoordinateType& point() const { return quadrature().point( index() ); }
    const LocalCoordinateType& localPoint() const { return quadrature().localPoint( index() ); }
    const int index() const { return qp_; }

    template <class Tuple, class VarId >
    struct Contains
    {
      static const bool value = false;
    };

    template <class Tuple, class Types, class VarId >
    struct Contains< TypeIndexedTuple< Tuple, Types >, VarId >
    {
      static const bool value = TypeIndexedTuple< Tuple, Types >::template Contains< VarId > :: value;
    };

    template <class Functor, bool containedInTuple >
    struct Evaluate;

    template <class Functor>
    struct Evaluate<Functor, true>
    {
      typedef typename Functor :: VarId   VarId;
      typedef typename RangeTuple :: template Value< VarId > :: Type ReturnType;

      template< class ... Args >
      static const ReturnType& eval( const RangeTuple& tuple, const Functor& functor, const Args& ... args )
      {
        return tuple.template at< VarId > ();
      }

    };

    template <class Functor>
    struct Evaluate<Functor, false>
    {
      typedef typename Functor::ReturnType ReturnType;

      template< class ... Args >
      //static auto
      static ReturnType
      eval( const RangeTuple& tuple, const Functor& functor, const Args& ... args )
       // -> decltype(functor( args ... ))
      {
        return functor( args ... );
      }
    };

    /**
     * \brief Returns the values for a local function regarding a static id or the evaluation of a functor
     *
     * The evaluation of the local function is dependent on the VarId which is
     * given by the functor which is passed as an argument.
     *
     * There are two cases
     * - Case Functor::VarId is found RangeTuple: return values()[Functor::VarId]
     * - Case Functor::VarId is not found RangeTuple: return Functor( args )
     *
     * The idea of this class is to allow the combination of either
     * - discrete data (Function::VarId is found) or
     * - analytical data (Function::VarId is not found).
     *
     * The simplest functor should be a constant functor:
     * \code{.cpp}
     * struct ConstantData
     * {
     *   typedef varId VarId;
     *   typedef double ReturnType;
     *
     *   const RangeType& operator() () const
     *   {
     *     return 42.0;
     *   }
     * };
     * \endcode
     *
     * This structure is then called by
     *
     * \code{.cpp}
     * this->evaluate( ConstantData() );
     * \endcode
     *
     * An interesting usage of this functionality could be the Model classes.
     *
     * \param functor A functor which should be applied
     * \param args arguments which should be applied to the functor
     */
    template <class Functor, class ... Args>
    typename Evaluate< Functor, Contains< RangeTuple, typename Functor::VarId >::value >::ReturnType
    evaluate( const Functor& functor, const Args& ... args ) const
    {
      return Evaluate< Functor, Contains< RangeTuple, typename Functor::VarId >::value>::eval( values(), functor, args ... );
    }
  };

  /**
   * \brief Container class describing (many) information which are related to
   * the local evaluation of a discrete function.
   *
   * This class just adds an intersection() method to the class ElementQuadraturePointContext
   */
  template <class Intersection,
            class Entity,
            class Quadrature,
            class RangeTuple,
            class JacobianTuple>
  class IntersectionQuadraturePointContext
    : public ElementQuadraturePointContext< Entity, Quadrature, RangeTuple, JacobianTuple >
  {
    typedef ElementQuadraturePointContext< Entity, Quadrature, RangeTuple, JacobianTuple
      > BaseType;
  protected:
    const Intersection& intersection_;

  public:
      typedef Intersection    IntersectionType;

      IntersectionQuadraturePointContext( const IntersectionType& intersection,
                                          const Entity& entity,
                                          const Quadrature& quadrature,
                                          const RangeTuple& values,
                                          const JacobianTuple& jacobians,
                                          const int qp,
                                          const double time,
                                          const double volume )
       : BaseType( entity, quadrature, values, jacobians, qp, time, volume ),
         intersection_( intersection )
    {}

    const IntersectionType& intersection() const { return intersection_; }
  };

} // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_FEM_DG_QUADRATURECONTEXT_HH
