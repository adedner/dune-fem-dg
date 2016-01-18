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
   * This class collects several information which are relevant for the approximation
   * \f$ \int_E U(x,t) + \nabla U(x,t) \mathrm{d}x\f$.
   * It is common practice to approximate the last term with a quadrature rule of
   * a given order \f$ k \f$. Furthermore, we suppose that each discrete function
   * \f$ U|_E \f$ is localizable, i.e. can be evaluated by evaluation of shape functions
   * with support on the corresponding reference element \f$ \hat{E} \f$.
   *
   * \f[ U|_E(x,t) = \sum_{i\in I_E} U_i^E(t) (\varphi_i^E \circ F_E )(\hat{x}) U_i^E(t) \hat{varphi}_i^E(\hat{x}), \f]
   * where
   * - \f$ F_E:\hat{E}\rightarrow E \f$ is the reference mapping.
   * - \f$ \varphi_i^E \circ F_E \f$ are basis functions with support on the entity \f$ E \f$,
   * - \f$ \hat{varphi}_i^{\hat{E}} \f$ are shape functions with support on the reference element \f$ \hat{E} \f$,
   * - \f$ I_E \f$ is an index set, enumerating the set of basis and shape functions with support on \f$ E \f$ and \f$ \hat{E} \f$, respectively,
   * - \f$ \hat{x} \f$ is the local coordinate regarding the reference element \f$ \hat{E} \f$
   *
   * In detail, we want
   * \f{eqnarray*}{
   * \f$ \int_E U(x,t) + \nabla U(x,t) \mathrm{d}x &=& \int_E U|_E(x,t) + \nabla U|_E(x,t) \mathrm{d}x\\
   *      &\approx & \sum_{i\in I_E} \hat{w}_i U_i(t) \varphi(x_i) + \nabla U_i(t) \varphi(x_i) \\
   *      & = &      \sum_{i\in I_E} \hat{w}_i ( U(\hat{x}_i) + \nabla U(\hat{x}_i) )
   * }
   *
   *
   * \note This class was introduced to allow more flexible interfaces for classes
   * needing a local evaluation of a discrete function.
   */
  template <class Entity,
            class Quadrature,
            class RangeTuple,
            class JacobianTuple>
  class ElementQuadraturePointContext
  {
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

  public:
    typedef Entity          EntityType;
    typedef Quadrature      QuadratureType;
    typedef RangeTuple      RangeTupleType;
    typedef JacobianTuple   JacobianTupleType;
    typedef typename QuadratureType::QuadraturePointWrapperType   QuadraturePointWrapperType;
    typedef typename QuadratureType::CoordinateType               CoordinateType;
    typedef typename QuadratureType::LocalCoordinateType          LocalCoordinateType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ E \f$
     *             (to be more precise: the quadrature rule for the corresponding reference element \f$ \hat{E} \f$)
     *  \param[in] values \f$ U_i \hat{\varphi}(\hat{x}_i) \f$
     *  \param[in] jacobians \f$ U_i \nabla \hat{\varphi}(\hat{x}_i) \f$
     *  \param[in] qp the number of the quadrature point, i.e. \f$ i \f$
     *  \param[in] time the current time \f$ t \f$
     *  \param[in] volume the volume of the entity \f$ vol(E) \f$
     */
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

    /**
     *  \brief returns the entity \f$ E \f$
     */
    const Entity& entity() const { return entity_; }

    /**
     *  \brief returns the quadrature
     */
    const Quadrature& quadrature() const { return quad_; }

    /**
     *  \brief returns the value \f$ U_i \hat{\varphi}(\hat{x}_i) \f$
     */
    const RangeTuple& values() const { return values_; }

    /**
     *  \brief returns the value \f$ U_i \nabla \hat{\varphi}(\hat{x}_i) \f$
     */
    const JacobianTuple& jacobians() const { return jacobians_; }

    /**
     *  \brief return the current time \f$ t \f$, which is needed for instationary problems
     */
    const double time () const { return time_; }

    /**
     *  \brief return the volume of the entity \f$ E \f$
     */
    const double volume() const { return volume_; }

    /**
     *  \brief returns a quadrature point object containing the following information \f$ \hat{x}_i, w_i, i, x_i \f$
     */
    const QuadraturePointWrapperType quadraturePoint() const { return quadrature()[ index() ]; }

    /**
     *  \brief returns the global quadrature point \f$ x_i \f$
     */
    const CoordinateType& position() const { return quadrature().point( index() ); }

    /**
     *  \brief returns the local point \f$ \hat{x}_i \f$
     */
    const LocalCoordinateType& localPosition() const { return quadrature().localPoint( index() ); }

    /**
     *  \brief returns the number of the quadrature point \f$ i \f$
     */
    const int index() const { return qp_; }


  public:
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

  protected:
    const Entity& entity_;
    const Quadrature& quad_;
    const RangeTuple& values_;
    const JacobianTuple& jacobians_;

    const double time_;
    const double volume_;
    const int qp_;
  };

  /**
   * \brief Container class describing (many) information which are related to
   * the local evaluation of a discrete function.
   *
   * This class just adds an intersection() method to the class ElementQuadraturePointContext.
   * This additional information is needed for evaluations on intersections, i.e. numerical fluxes etc.
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

      /**
       * \brief constructor
       */
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
