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
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   */
  template <class Entity,
            class Quadrature >
  class ElementQuadratureContext
  {
  public:
    typedef Entity          EntityType;
    typedef Quadrature      QuadratureType;
    typedef typename QuadratureType::QuadraturePointWrapperType   QuadraturePointWrapperType;
    typedef typename QuadratureType::CoordinateType               CoordinateType;
    typedef typename QuadratureType::LocalCoordinateType          LocalCoordinateType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] time the current time \f$ t \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    ElementQuadratureContext( const Entity& entity,
                              const Quadrature& quadrature,
                              const double time,
                              const double volume )
     : entity_( entity ),
       quad_( quadrature ),
       time_( time ),
       volume_( volume )
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
     *  \brief return the current time \f$ t \f$, which is needed for instationary problems
     */
    const double time () const { return time_; }

    /**
     *  \brief return the volume of the entity \f$ E \f$
     */
    const double volume() const { return volume_; }

  protected:
    const Entity& entity_;
    const Quadrature& quad_;

    const double time_;
    const double volume_;
  };





  /**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   */
  template <class Intersection,
            class Entity,
            class Quadrature >
  class IntersectionQuadratureContext
    : public ElementQuadratureContext<Entity,Quadrature>
  {
    typedef ElementQuadratureContext<Entity,Quadrature> BaseType;
  public:
    typedef Intersection    IntersectionType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] time the current time \f$ t \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    IntersectionQuadratureContext( const Intersection& intersection,
                                   const Entity& entity,
                                   const Quadrature& quadrature,
                                   const double time,
                                   const double volume )
     : BaseType( entity, quadrature, time, volume ),
       intersection_( intersection )
    {}

    const IntersectionType& intersection() const { return intersection_; }

  protected:
    const Intersection& intersection_;
  };



/**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   *
   * \ingroup PassBased
   *
   * We are following the \ref Notation "general notation":
   *
   * Let \f$ u^t,v^t:\mathcal{\Omega}_\mathcal{G} \rightarrow \mathbb{R}^d\f$ be discrete functions for a fixed time
   * \f$ t\in [t_{\text{start}},t_{\text{end}}] \f$. Now, we want to be able to approximate the integral over an
   * entity \f$ E \in \mathcal{G} \f$
   *
   * \f[ \int_E u^t(x) + \nabla v^t(x) \mathrm{d}x \f]
   *
   * via a quadrature rule. This can be done in the following way
   *
   * \f{eqnarray*}{
   *    \int_E u^t(x) + \nabla v^t(x) \mathrm{d}x &=& \int_E u^t|_E(x) + \nabla v^t|_E(x) \mathrm{d}x\\
   *      &=&        \int_{E} u^t_E(F_E^{-1}(x)) + \nabla v^t_E(F_E^{-1}(x)) \mathrm{d}x\\
   *      &=&        \int_{\hat{E}} \left|\det DF_E(\hat{x})\right|\left( u^t_E(\hat{x}) + \nabla v^t_E(\hat{x}) \right)\mathrm{d}\hat{x} \\
   *      &\approx & \sum_{p\in X_p} \left|\det DF_E(\hat{x})\right| \omega_i \left(u^t_E(\hat{x}_p) + \nabla v^t_E(\hat{x}_p) \right)
   * \f}
   *
   * To archieve this goal, we collect some of the information above in this class.
   *
   * This class was introduced to allow more flexible interfaces for classes
   * needing a local evaluation of a discrete function.
   *
   * \note This class is following the pass concept and the above description is simplified in
   * the sense that we are not only interested in the approximation of one discrete function \f$ u^t \f$
   * but in a tuple of functions \f$ (u^{a,t}, {u^b,t} \ldots\f$ yielding from a pass.
   */
  template <class Entity,
            class Quadrature >
  class ElementQuadraturePointContext
    : public ElementQuadratureContext< Entity, Quadrature >
  {
    typedef ElementQuadratureContext< Entity, Quadrature > BaseType;
  public:
    using BaseType::quadrature;
    typedef typename BaseType::QuadraturePointWrapperType   QuadraturePointWrapperType;
    typedef typename BaseType::CoordinateType               CoordinateType;
    typedef typename BaseType::LocalCoordinateType          LocalCoordinateType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] qp the number of the quadrature point, i.e. \f$ p \f$
     *  \param[in] time the current time \f$ t \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    ElementQuadraturePointContext( const Entity& entity,
                                   const Quadrature& quadrature,
                                   const int qp,
                                   const double time,
                                   const double volume )
     : BaseType( entity, quadrature, time, volume ),
       qp_( qp )
    {}

    /**
     *  \brief returns a quadrature point object containing the following information \f$ \hat{x}_p, x_p, \omega_p, p \f$
     */
    const QuadraturePointWrapperType quadraturePoint() const { return quadrature()[ index() ]; }

    /**
     *  \brief returns the global quadrature point \f$ x_p \f$
     */
    const CoordinateType& position() const { return quadrature().point( index() ); }

    /**
     *  \brief returns the local point \f$ \hat{x}_p \f$
     */
    const LocalCoordinateType& localPosition() const { return quadrature().localPoint( index() ); }

    /**
     *  \brief returns the number of the quadrature point \f$ p \f$
     */
    const int index() const { return qp_; }

  protected:
    const int qp_;
  };



  /**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   */
  template <class Intersection,
            class Entity,
            class Quadrature >
  class IntersectionQuadraturePointContext
    : public ElementQuadraturePointContext<Entity,Quadrature>
  {
    typedef ElementQuadraturePointContext<Entity,Quadrature> BaseType;
  public:
    typedef Intersection    IntersectionType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] time the current time \f$ t \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    IntersectionQuadraturePointContext( const Intersection& intersection,
                                        const Entity& entity,
                                        const Quadrature& quadrature,
                                        const int qp,
                                        const double time,
                                        const double volume )
     : BaseType( entity, quadrature, qp, time, volume ),
       intersection_( intersection )
    {}

    const IntersectionType& intersection() const { return intersection_; }

  protected:
    const Intersection& intersection_;
  };



  /**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   *
   * \ingroup PassBased
   *
   * We are following the \ref Notation "general notation":
   *
   * Let \f$ u^t,v^t:\mathcal{\Omega}_\mathcal{G} \rightarrow \mathbb{R}^d\f$ be discrete functions for a fixed time
   * \f$ t\in [t_{\text{start}},t_{\text{end}}] \f$. Now, we want to be able to approximate the integral over an
   * entity \f$ E \in \mathcal{G} \f$
   *
   * \f[ \int_E u^t(x) + \nabla v^t(x) \mathrm{d}x \f]
   *
   * via a quadrature rule. This can be done in the following way
   *
   * \f{eqnarray*}{
   *    \int_E u^t(x) + \nabla v^t(x) \mathrm{d}x &=& \int_E u^t|_E(x) + \nabla v^t|_E(x) \mathrm{d}x\\
   *      &=&        \int_{E} u^t_E(F_E^{-1}(x)) + \nabla v^t_E(F_E^{-1}(x)) \mathrm{d}x\\
   *      &=&        \int_{\hat{E}} \left|\det DF_E(\hat{x})\right|\left( u^t_E(\hat{x}) + \nabla v^t_E(\hat{x}) \right)\mathrm{d}\hat{x} \\
   *      &\approx & \sum_{p\in X_p} \left|\det DF_E(\hat{x})\right| \omega_i \left(u^t_E(\hat{x}_p) + \nabla v^t_E(\hat{x}_p) \right)
   * \f}
   *
   * To archieve this goal, we collect some of the information above in this class.
   *
   * This class was introduced to allow more flexible interfaces for classes
   * needing a local evaluation of a discrete function.
   *
   * \note This class is following the pass concept and the above description is simplified in
   * the sense that we are not only interested in the approximation of one discrete function \f$ u^t \f$
   * but in a tuple of functions \f$ (u^{a,t}, {u^b,t} \ldots\f$ yielding from a pass.
   */
  template <class Entity,
            class Quadrature,
            class RangeTuple,
            class JacobianTuple>
  class ExtraElementQuadraturePointContext
    : public ElementQuadraturePointContext< Entity, Quadrature >
  {
    typedef ElementQuadraturePointContext< Entity, Quadrature > BaseType;

    template <class Tuple, class VarId >
    struct Contains
    {
      static const bool value = false;
    };

    //pass version
    template <class Tuple, class Types, class VarId >
    struct Contains< TypeIndexedTuple< Tuple, Types >, VarId >
    {
      static const bool value = TypeIndexedTuple< Tuple, Types >::template Contains< VarId > :: value;
    };

    //general version, without passes, just a simple std::tuple of arguments
    template <class... Args, class VarId >
    struct Contains< std::tuple< Args...>, VarId >
    {
      static const bool value = (VarId::value >= 0 && std::tuple_size<std::tuple<Args...> >::value < VarId::value);
    };

    template <class Functor, bool containedInTuple >
    struct Evaluate;

    template <class Functor>
    struct Evaluate<Functor, false>
    {
      template< class ... Args >
      static decltype(auto) eval( const RangeTuple& tuple, const Functor& functor, const Args& ... args )
      {
        return functor( args ... );
      }
    };

    template <class Functor>
    struct Evaluate<Functor, true>
    {
      //pass version
      template< class Tuple, class Types, class ... Args >
      static decltype(auto) eval( const TypeIndexedTuple< Tuple, Types >& tuple, const Functor& functor, const Args& ... args )
      {
        return tuple.template at< typename Functor::VarId >();
      }

      //no passes: assume a simple std::tuple of elements
      template< class... ExtraArgs, class ... Args >
      static decltype(auto) eval( const std::tuple< ExtraArgs... >& tuple, const Functor& functor, const Args& ... args )
      {
        return std::get< typename Functor::VarId::value >( tuple );
      }
    };

  public:
    typedef RangeTuple      RangeTupleType;
    typedef JacobianTuple   JacobianTupleType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] values tuple of a evaluation of local functions,
     *             i.e. \f$ ( u^{a,t}_E( \hat{x}_p ), u^{b,t}_E( \hat{x}_p )\ldots ) \f$
     *  \param[in] jacobians tuple of a evaluation of the gradient of local functions,
     *             i.e. \f$ ( \nabla u^{a,t}_E( \hat{x}_p ), \nabla u^{b,t}_E( \hat{x}_p )\ldots ) \f$
     *  \param[in] qp the number of the quadrature point, i.e. \f$ p \f$
     *  \param[in] time the current time \f$ t \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    ExtraElementQuadraturePointContext( const Entity& entity,
                                        const Quadrature& quadrature,
                                        const int qp,
                                        const double time,
                                        const double volume,
                                        const RangeTuple& values,
                                        const JacobianTuple& jacobians )
     : BaseType( entity, quadrature, qp, time, volume ),
       values_( values ),
       jacobians_( jacobians )
    {}

    /**
     *  \brief returns tuple of a evaluation of local functions,
     *         i.e. \f$ ( u^{a,t}_E( \hat{x}_p ), u^{b,t}_E( \hat{x}_p )\ldots ) \f$
     */
    const RangeTuple& values() const { return values_; }

    /**
     *  \brief returns tuple of a evaluation of the gradient of local functions,
     *         i.e. \f$ ( \nabla u^{a,t}_E( \hat{x}_p ), \nabla u^{b,t}_E( \hat{x}_p )\ldots ) \f$
     */
    const JacobianTuple& jacobians() const { return jacobians_; }

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
    decltype(auto) evaluate( const Functor& functor, const Args& ... args ) const
    {
      return Evaluate< Functor, Contains< RangeTuple, typename Functor::VarId >::value>::eval( values(), functor, args ... );
    }

  protected:
    const RangeTuple& values_;
    const JacobianTuple& jacobians_;
  };

  /**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   *
   * This class just adds an intersection() method to the class ExtraElementQuadraturePointContext.
   * This additional information is needed for evaluations on intersections, i.e. numerical fluxes etc.
   */
  template <class Intersection,
            class Entity,
            class Quadrature,
            class RangeTuple,
            class JacobianTuple>
  class ExtraIntersectionQuadraturePointContext
    : public ExtraElementQuadraturePointContext< Entity, Quadrature, RangeTuple, JacobianTuple >
  {
    typedef ExtraElementQuadraturePointContext< Entity, Quadrature, RangeTuple, JacobianTuple
      > BaseType;
  protected:
    const Intersection& intersection_;

  public:
      typedef Intersection    IntersectionType;

      /**
       * \brief constructor
       */
      ExtraIntersectionQuadraturePointContext( const IntersectionType& intersection,
                                               const Entity& entity,
                                               const Quadrature& quadrature,
                                               const int qp,
                                               const double time,
                                               const double volume,
                                               const RangeTuple& values,
                                               const JacobianTuple& jacobians )
       : BaseType( entity, quadrature, qp, time, volume, values, jacobians ),
         intersection_( intersection )
    {}

    const IntersectionType& intersection() const { return intersection_; }
  };

} // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_FEM_DG_QUADRATURECONTEXT_HH
