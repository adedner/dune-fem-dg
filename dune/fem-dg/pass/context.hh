#ifndef DUNE_FEM_DG_QUADRATURECONTEXT_HH
#define DUNE_FEM_DG_QUADRATURECONTEXT_HH

#include <cassert>
#include <vector>
#include <dune/fem/common/typeindexedtuple.hh>


namespace Dune
{
namespace Fem
{


  template < class T, typename std::enable_if < !is_vector<T>::value, int >::type = 0 >
  decltype(auto) indexIfVector(const T& value, int i)
  {
    return value;
  }

  template<class T, typename std::enable_if<is_vector<T>::value, int>::type = 0>
  decltype(auto) indexIfVector(const T& value, int i)
  {
    typedef std::remove_reference_t<decltype(std::declval<T>()[0] )> ElemType;
    return i < 0 ? ElemType(0) : value[i];
  }

  // (qp, id, base)   std::get<id>(arg[qp])[base]

  // (qp, id      )   std::get<id>(arg[qp])
  // (qp,     base)   operator[*]: std::get<*>(arg[qp])[base]
  // (    id, base)   operator[*]: std::get<id>(arg[*])[base]

  // (qp          )   operator[*]: std::get<*>(arg[qp])[*]
  // (    id      )   operator[*]: std::get<id>(arg[*])[*]
  // (        base)   operator[*]: std::get<*>(arg[*])[base]

  template< class E, bool qpAct, bool idAct, bool baseAct, int id >
  struct Eval;

  template< class E, int id >
  struct FinalEval
  {
    FinalEval( const E& e, int qp, int b ) : e_(e) {}

    static_assert( static_fail< E >::value, "alert: This should not happen: Use std::vector, std::tuple etc..." );
    decltype(auto) operator()()
    { return e_; }
  private:
    const E& e_;
  };


  //////////// Last Element
  template< class Arg, int id >
  struct FinalEval<std::vector< Arg >,id>
  {
    typedef std::vector< Arg > E;
    FinalEval( const E& e, int qp, int b ) : e_(e), b_(b) {}

    decltype(auto) operator()()
    { return indexIfVector( e_, b_ ); }
  private:
    const E& e_; int b_;
  };

  template< class Arg, int id >
  struct FinalEval<std::vector<std::vector<Arg> >,id>
  {
    typedef std::vector<std::vector<Arg> > E;
    FinalEval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    decltype(auto) operator()()
    { return indexIfVector( e_[qp_], b_ ); }
  private:
    const E& e_; int qp_; int b_;
  };

  template< class... Args, int id >
  struct FinalEval<std::vector<std::tuple<Args...> >,id>
  {
    typedef std::vector<std::tuple<Args...> > E;
    FinalEval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    decltype(auto) operator()()
    { return indexIfVector( std::get<id>(e_[qp_]), b_ ); }
  private:
    const E& e_; int qp_; int b_;
  };

  template< class... Args, int id >
  struct FinalEval<std::vector<Dune::TypeIndexedTuple<Args...> >,id>
  {
    typedef std::vector<Dune::TypeIndexedTuple<Args...> > E;
    FinalEval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    decltype(auto) operator()()
    { return indexIfVector( e_[qp_][std::integral_constant<int,id>()], b_ ); }
  private:
    const E& e_; int qp_; int b_;
  };

  template< class... Args, int id >
  struct FinalEval<Dune::TypeIndexedTuple<Args...>,id>
  {
    typedef Dune::TypeIndexedTuple<Args...> E;
    FinalEval( const E& e, int qp, int b ) : e_(e), b_(b) {}

    decltype(auto) operator()()
    { return indexIfVector( e_[std::integral_constant<int,id>()], b_ ); }
  private:
    const E& e_; int b_;
  };


  // this should not happen
  template< class E, bool b1, bool b2, bool b3, bool id, int outerId >
  struct FinalEval<Eval<E,b1,b2,b3,id>,outerId >
  {
    static_assert( static_fail< E >::value, "alert: This should not happen: Use std::vector, std::tuple etc..." );
  };

  //////////// End Last Element


  template< class E, int id >
  struct Eval<E,false,false,false,id>
    : public std::remove_cv_t<std::remove_reference_t<decltype(std::declval<FinalEval<E,id> >()())> >
  {
  private:
    typedef FinalEval<E,id> FinalEvalType;
    typedef std::remove_cv_t<std::remove_reference_t<decltype(std::declval<FinalEvalType>()())> > BaseType;
  public:
    using BaseType::size;
    Eval( const E& e, int qp, int b )
      : BaseType( FinalEvalType( e, qp, b )() )
    {}
  };


  template< class E, int id >
  struct Eval<E,true,false,false,id>
  {
    Eval( const E& e, int qp, int b ) : e_(e), b_(b) {}

    decltype(auto) operator[]( int idx ) const
    { return FinalEval<E,id>(e_,idx,b_)(); }
  private:
    const E& e_; int b_;
  };

  template< class E, int id >
  struct Eval<E,false,true,false,id>
  {
    Eval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    template< class Int, Int idx >
    decltype(auto) operator[]( const std::integral_constant<Int,idx>& ) const
    { return FinalEval<E,idx>(e_,qp_,b_)(); }
  private:
    const E& e_; int qp_; int b_;
  };

  template< class E, int id >
  struct Eval<E,false,false,true,id>
  {
    Eval( const E& e, int qp, int b ) : e_(e), qp_(qp) {}

    decltype(auto) operator[]( int idx ) const
    { return FinalEval<E,id>(e_,qp_,idx)(); }

  private:
    const E& e_; int qp_;
  };

  template< class E, int id >
  struct Eval<E,true,true,false,id>
  {
    Eval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    //access: quadrature point
    decltype(auto) operator[]( unsigned int idx ) const
    { return Eval<E,false,false,true,id>( e_, idx, b_ ); }

    //access: quadrature point II
    decltype(auto) operator[]( unsigned long int idx ) const
    { return Eval<E,false,false,true,id>( e_, idx, b_ ); }

    //access: multi values
    template< class Int, Int idx >
    decltype(auto) operator[]( const std::integral_constant<Int,idx>& ) const
    { return Eval<E,true,false,false,idx>( e_, qp_, b_ ); }

  private:
    template< class NoImplicitTypeConversion >
    void operator[]( NoImplicitTypeConversion ) const {}

    const E& e_; int qp_; int b_;
  };

  template< class E, int id >
  struct Eval<E,true,false,true,id>
  {
    Eval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    //access: quadrature point
    decltype(auto) operator[]( unsigned int idx ) const
    { return Eval<E,false,false,true,id>( e_, idx, b_ ); }

    //access: quadrature point II
    decltype(auto) operator[]( unsigned long int idx ) const
    { return Eval<E,false,false,true,id>( e_, idx, b_ ); }

    //access: basis functions
    decltype(auto) operator[]( int idx ) const
    { return Eval<E,true,false,false,id>( e_, qp_, idx ); }
  private:
    template< class NoImplicitTypeConversion >
    void operator[]( NoImplicitTypeConversion ) const {}

    const E& e_; int qp_; int b_;
  };

  template< class E, int id >
  struct Eval<E,false,true,true,id>
  {
    Eval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    //access: multi values
    template< class Int, Int idx >
    decltype(auto) operator[]( const std::integral_constant<Int,idx>& ) const
    { return Eval<E,false,false,true,idx>( e_, qp_, b_ ); }

    //acess: basis functions
    decltype(auto) operator[]( int idx ) const
    { return Eval<E,false,true,false,id>( e_, qp_, idx ); }
  private:
    template< class NoImplicitTypeConversion >
    void operator[]( NoImplicitTypeConversion ) const {}

    const E& e_; int qp_; int b_;
  };

  template< class E, int id >
  struct Eval<E,true,true,true,id>
  {
    Eval( const E& e, int qp, int b ) : e_(e), qp_(qp), b_(b) {}

    //access: quadrature point
    decltype(auto) operator[]( unsigned int idx ) const
    { return Eval<E,false,true,true,id>( e_, idx, b_ ); }

    //access: quadrature point II
    decltype(auto) operator[]( unsigned long int idx ) const
    { return Eval<E,false,true,true,id>( e_, idx, b_ ); }

    //access: multi values
    template< class Int, Int idx >
    decltype(auto) operator[]( const std::integral_constant<Int,idx>& ) const
    { return Eval<E,true,false,true,idx>( e_, qp_, b_ ); }

    //access: basis functions
    decltype(auto) operator[]( int idx ) const
    { return Eval<E,true,true,false,id>( e_, qp_, idx ); }
  private:
    template< class NoImplicitTypeConversion >
    void operator[]( NoImplicitTypeConversion ) const {}

    const E& e_; int qp_; int b_;
  };


  template< bool hasQPs = true, bool hasIds = true, bool hasBasis = true, bool isNonZero = true, int id = -1 >
  struct Access
  {
    static constexpr bool hasQuadPoints = hasQPs;
    static constexpr bool hasBasisFunctions = hasBasis;
    static constexpr bool hasMultiValues = hasIds;

    typedef Access<false,hasIds,hasBasis,isNonZero,id> QuadraturePointType;
    template< int i >
    using MultiValueType = Access<hasQPs,false,hasBasis,isNonZero,i>;
    typedef Access<hasQPs,hasIds,false,false,id>     ZeroBasisFunctionType;
    typedef Access<hasQPs,hasIds,false,isNonZero,id> BasisFunctionType;

    static const int multiValue = id;
  };


  //maximal access length based on types
  template< class EvalImp >
  struct EvalAccess
  {
    //fallback: no, we cannot select anything
    typedef Access< false, false, false > type;
  };

  template< class Arg >
  struct EvalAccess<std::vector< Arg > >
  {
    typedef Access< true, false, false > type;
  };

  template< class Arg >
  struct EvalAccess<std::vector<std::vector<Arg> > >
  {
    typedef Access< true, false, true > type;
  };

  template< class... Args >
  struct EvalAccess<std::vector<std::tuple<Args...> > >
  {
    typedef Access< true, true, true > type;
  };

  template< class... Args >
  struct EvalAccess<std::tuple<Args...> >
  {
    typedef Access< false, true, true > type;
  };


  template< class... Args >
  struct EvalAccess<std::vector<Dune::TypeIndexedTuple<Args...> > >
  {
    typedef Access< true, true, false > type;
  };

  template< class... Args >
  struct EvalAccess<Dune::TypeIndexedTuple<Args...> >
  {
    typedef Access< false, true, false > type;
  };


  struct FunctorContext
  {
    template <class Tuple, class VarId >
    struct Contains
    {
      static const bool value = false;
    };

    //pass version
    template <class Tuple, class Types, class VarId >
    struct Contains< TypeIndexedTuple< Tuple, Types >, VarId >
    {
      static const bool value = TypeIndexedTuple< Tuple, Types >::template Contains<VarId>::value;
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
      template< class RangeTuple, class ... Args >
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
  };


  //forward declarations
  template< class A, class... >
  class Context
  {
    static_assert( static_fail<A>::value, "Context expects either one or two template arguments" );
  };

  //forward declarations
  template< class A, class, class... >
  class QuadratureContext
  {
    static_assert( static_fail<A>::value, "QuadraturePointContext expects either two or three template arguments" );
  };

  //forward declarations
  template< class A, class... >
  class PointContext
  {
    static_assert( static_fail<A>::value, "PointContext expects either one or two template arguments" );
  };


  template< class Intersection >
  class IntersectionStorage
  {
  public:
    typedef Intersection IntersectionType;

    explicit IntersectionStorage( const Intersection& intersection )
      : intersection_( intersection )
    {}
    const IntersectionType& intersection() const { return intersection_; }
  protected:
    const Intersection& intersection_;
  };

  /**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   */
  template< class Entity >
  class Context< Entity >
  {
  public:
    typedef Entity          EntityType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    Context( const Entity& entity,
             const double volume )
     : entity_( entity ),
       volume_( volume )
    {}

    /**
     *  \brief returns the entity \f$ E \f$
     */
    const Entity& entity() const { return entity_; }

    /**
     *  \brief return the volume of the entity \f$ E \f$
     */
    const double volume() const { return volume_; }

    //nasty dummy, just to have an intersection() method.
    decltype(auto) intersection() const { return{}; }

  protected:
    const Entity& entity_;
    const double volume_;
  };

  /**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   */
  template< class Entity, class Intersection >
  class Context< Entity, Intersection >
    : public Context< Entity >,
      public IntersectionStorage< Intersection >
  {
    typedef Context< Entity >                   BaseType;
    typedef IntersectionStorage< Intersection > InterBaseType;
  public:
    typedef Entity          EntityType;
    typedef Intersection    IntersectionType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] intersection the intersection where the local evaluation should be done
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    Context( const Entity& entity,
             const Intersection& intersection,
             const double volume )
     : BaseType( entity, volume ),
       InterBaseType( intersection )
    {}
  };


  /**
   * \brief This class collects several information which are relevant for the approximation of
   * integrals of discrete functions via quadrature schemes.
   */
  template <class Entity >
  class PointContext< Entity >
    : public Context< Entity >
  {
    typedef Context< Entity > BaseType;
  public:
    typedef typename Entity::Geometry::LocalCoordinate  LocalCoordinateType;
    typedef typename Entity::Geometry::GlobalCoordinate CoordinateType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    PointContext( const Entity& entity,
                  const CoordinateType& position,
                  const LocalCoordinateType& localPos,
                  const double volume )
     : BaseType( entity, volume ),
       position_( position ),
       localPos_( localPos )
    {}

    /**
     *  \brief returns the global quadrature point \f$ x_p \f$
     */
    const CoordinateType& position() const { return position_; }

    /**
     *  \brief returns the local point \f$ \hat{x}_p \f$
     */
    const LocalCoordinateType& localPosition() const { return localPos_; }

  protected:
    const CoordinateType& position_;
    const LocalCoordinateType& localPos_;
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
  template <class Entity, class Quadrature >
  class QuadratureContext< Entity, Quadrature >
    : public Context< Entity >
  {
    typedef Context< Entity >                                     BaseType;
  public:
    typedef Entity                                                EntityType;
    typedef Quadrature                                            QuadratureType;
    typedef typename QuadratureType::QuadraturePointWrapperType   QuadraturePointWrapperType;
    typedef typename QuadratureType::CoordinateType               CoordinateType;
    typedef typename QuadratureType::LocalCoordinateType          LocalCoordinateType;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] qp the number of the quadrature point, i.e. \f$ p \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    QuadratureContext( const Entity& entity,
                       const Quadrature& quadrature,
                       const double volume,
                       const int qp = -1 )
     : BaseType( entity, volume ),
       quad_( quadrature ),
       qp_( qp )
    {}

    using BaseType::entity;
    using BaseType::volume;

    //short cut
    template< class LocalEvaluation >
    QuadratureContext( const LocalEvaluation& local,
                       const int qp )
     : BaseType( local.entity(), local.quadrature(), qp, local.volume() )
    {}

    /**
     *  \brief returns a quadrature point object containing the following information \f$ \hat{x}_p, x_p, \omega_p, p \f$
     */
    const QuadraturePointWrapperType quadraturePoint() const { assert( index() >= 0 ); return quadrature()[ index() ]; }

    /**
     *  \brief returns the global quadrature point \f$ x_p \f$
     */
    const CoordinateType& position() const { assert( index() >= 0 ); return quadrature().point( index() ); }

    /**
     *  \brief returns the local point \f$ \hat{x}_p \f$
     */
    const LocalCoordinateType& localPosition() const { assert( index() >= 0 ); return quadrature().localPoint( index() ); }

    /**
     *  \brief returns the number of the quadrature point \f$ p \f$
     */
    const int index() const { return qp_; }

    /**
     *  \brief returns the quadrature
     */
    const Quadrature& quadrature() const { return quad_; }

    decltype(auto) operator[]( int idx ) const
    {
      return QuadratureContext<Entity,Quadrature>( entity(), quad_, volume(), idx );
    }

    void setIndex( int qp ) const { qp_ = qp; }
  protected:
    const Quadrature& quad_;
    mutable int qp_;
  };

  template <class Entity, class Intersection, class Quadrature >
  class QuadratureContext< Entity, Intersection, Quadrature >
    : public QuadratureContext< Entity, Quadrature >,
      public IntersectionStorage< Intersection >

  {
    typedef QuadratureContext< Entity, Quadrature >               BaseType;
    typedef IntersectionStorage< Intersection >                   BaseType2;
  public:
    typedef Entity                                                EntityType;
    typedef Intersection                                          IntersectionType;
    typedef Quadrature                                            QuadratureType;
    typedef typename QuadratureType::QuadraturePointWrapperType   QuadraturePointWrapperType;
    typedef typename QuadratureType::CoordinateType               CoordinateType;
    typedef typename QuadratureType::LocalCoordinateType          LocalCoordinateType;

    using BaseType2::intersection;

    /**
     *  \brief constructor
     *
     *  \param[in] entity the entity \f$ E \f$ where the local evaluation should be done
     *  \param[in] quadrature the quadrature rule for the entity \f$ \hat{E} \f$
     *  \param[in] qp the number of the quadrature point, i.e. \f$ p \f$
     *  \param[in] volume the volume of the entity \f$ \mathrm{vol}(E) \f$
     */
    QuadratureContext( const Entity& entity,
                       const Intersection& intersection,
                       const Quadrature& quadrature,
                       const double volume,
                       const int qp = -1 )
     : BaseType( entity, quadrature, volume, qp ),
       BaseType2( intersection )
    {}

    decltype(auto) operator[]( int idx ) const
    {
      return QuadratureContext<Entity,Intersection,Quadrature>( BaseType::entity(), intersection(), BaseType::quad_, BaseType::volume(), idx );
    }
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
   */

  template< class QuadratureContextImp, class RangeType, class JacobianType = RangeType, class AccessImp = typename EvalAccess<RangeType>::type >
  class LocalEvaluation
  {
    template <class EvalImp >
    using Evaluator = Eval< EvalImp, AccessImp::hasQuadPoints, AccessImp::hasMultiValues, AccessImp::hasBasisFunctions, AccessImp::multiValue >;

  public:
    LocalEvaluation( const QuadratureContextImp& quadImp, const RangeType& values, const JacobianType& jacobians, int qp, int basis = -1 )
    : quadImp_( quadImp ),
      values_( values ),
      jacobians_( jacobians ),
      qp_( qp ),
      basis_( basis )
    {
      //small hack
      quadImp_.setIndex(qp);
    }

    //default range = jacobian
    LocalEvaluation( const QuadratureContextImp& quadImp, const RangeType& values, int qp, int basis = -1 )
    : quadImp_( quadImp ),
      values_( values ),
      jacobians_( values ),
      qp_( qp ),
      basis_( basis )
    {
      //small hack
      quadImp_.setIndex(qp);
    }

    LocalEvaluation( const QuadratureContextImp& quadImp, const RangeType& values, const JacobianType& jacobians )
    : quadImp_( quadImp ),
      values_( values ),
      jacobians_( jacobians ),
      qp_( -1 ),
      basis_( -1 )
    {}

    //default range = jacobian
    LocalEvaluation( const QuadratureContextImp& quadImp, const RangeType& values )
    : quadImp_( quadImp ),
      values_( values ),
      jacobians_( values ),
      qp_( -1 ),
      basis_( -1 )
    {}

    typedef typename QuadratureContextImp::EntityType                 EntityType;
    typedef typename QuadratureContextImp::QuadratureType             QuadratureType;
    typedef typename QuadratureContextImp::QuadraturePointWrapperType QuadraturePointWrapperType;
    typedef typename QuadratureContextImp::CoordinateType             CoordinateType;
    typedef typename QuadratureContextImp::LocalCoordinateType        LocalCoordinateType;

    /**
     *  \brief returns a quadrature point object containing the following information \f$ \hat{x}_p, x_p, \omega_p, p \f$
     */
    const QuadraturePointWrapperType quadraturePoint() const { return quadImp_.quadraturePoint(); }

    /**
     *  \brief returns the global quadrature point \f$ x_p \f$
     */
    const CoordinateType& position() const { return quadImp_.position(); }

    /**
     *  \brief returns the local point \f$ \hat{x}_p \f$
     */
    const LocalCoordinateType& localPosition() const { return quadImp_.localPosition(); }

    /**
     *  \brief returns the number of the quadrature point \f$ p \f$
     */
    const int index() const { return quadImp_.index(); }

    /**
     *  \brief returns the quadrature
     */
    const QuadratureType& quadrature() const { return quadImp_.quadrature(); }

    /**
     *  \brief returns the entity \f$ E \f$
     */
    const EntityType& entity() const { return quadImp_.entity(); }

    /**
     *  \brief returns the entity \f$ E \f$
     */
    decltype(auto) intersection() const { return quadImp_.intersection(); }

    /**
     *  \brief return the volume of the entity \f$ E \f$
     */
    const double volume() const { return quadImp_.volume(); }

    //access: quadpoints
    decltype(auto) operator[]( unsigned int idx ) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::QuadraturePointType > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, idx, basis_ );
    }

    //access: quadpoints II
    decltype(auto) operator[]( unsigned long int idx ) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::QuadraturePointType > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, idx, basis_ );
    }

    //access: multi values, i.e. tuples
    template< class Int, Int id >
    decltype(auto) operator[]( const std::integral_constant<Int,id>& idx ) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::template MultiValueType<id> > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, qp_, basis_ );
    }

    //access: basis functions, zero
    decltype(auto) operator()(int i=0) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::ZeroBasisFunctionType > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, qp_, -1 );
    }

    //access: basis functions
    decltype(auto) operator[]( int idx ) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::BasisFunctionType> NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, qp_, idx );
    }

    //access: failure
    template< class Fail >
    decltype(auto) operator[]( Fail ) const
    {
      static_assert( static_fail<Fail>::value, "This should not happen. Please check the exact type of your argument." );
      return {};
    }

    const decltype(auto) values() const
    {
      return Evaluator<RangeType>( values_, qp_, basis_ );
    }

    const decltype(auto) jacobians() const
    {
      return Evaluator<JacobianType>( jacobians_, qp_, basis_ );
    }

    template <class Functor, class ... Args>
    decltype(auto) values( const Functor& functor, const Args& ... args ) const
    {
      return FunctorContext::Evaluate< Functor, FunctorContext::Contains< RangeType, typename Functor::VarId >::value>::eval( values(), functor, args ... );
    }

    template <class Functor, class ... Args>
    decltype(auto) jacobians( const Functor& functor, const Args& ... args ) const
    {
      return FunctorContext::Evaluate< Functor, FunctorContext::Contains< JacobianType, typename Functor::VarId >::value>::eval( jacobians(), functor, args ... );
    }

  private:
    const QuadratureContextImp& quadImp_;
    const RangeType& values_;
    const JacobianType& jacobians_;
    int qp_;
    int basis_;
  };





} // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_FEM_DG_QUADRATURECONTEXT_HH
