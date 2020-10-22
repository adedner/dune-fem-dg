#ifndef DUNE_FEM_DG_QUADRATURECONTEXT_HH
#define DUNE_FEM_DG_QUADRATURECONTEXT_HH

#include <cassert>
#include <vector>
#include <dune/fem/common/typeindexedtuple.hh>
#include <dune/fem-dg/misc/tupleutility.hh>


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
    //negative index i means: set to zero
    //otherwise: return value
    return i < 0 ? ElemType(0) : value[i];
  }

  // (qp, id, base)   std::get<id>(arg[qp])[base]

  // (qp, id      )   std::get<id>(arg[qp])
  // (qp,     base)   operator[*]: std::get<*>(arg[qp])[base]
  // (    id, base)   operator[*]: std::get<id>(arg[*])[base]

  // (qp          )   operator[*]: std::get<*>(arg[qp])[*]
  // (    id      )   operator[*]: std::get<id>(arg[*])[*]
  // (        base)   operator[*]: std::get<*>(arg[*])[base]




  // Pass (instationary):         std::vector< TypeIndexedTuple<> >
  //
  // Non-Pass (stationary):       std::vector< std::tuple< std::vector< > > >
  // Non-Pass (stationary):       std::vector< std::vector< > >


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

  ////////
  //template< class... Args, int id >
  //struct FinalEval<std::tuple<Args...>,id>
  //{
  //  typedef std::tuple<Args...> E;
  //  FinalEval( const E& e, int qp, int b ) : e_(e), b_(b) {}

  //  decltype(auto) operator()()
  //  { return indexIfVector( std::get<id>(e_), b_ ); }
  //private:
  //  const E& e_; int b_;
  //};

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


  /**
   * \brief class for storing the access pattern of a local evaluation in a static way.
   *
   * \tparam hasQps true, if quadrature point has been selected
   * \tparam hasIds true, if id has been selected
   * \tparam hasBasis true, if basis function has been selected
   * \tparam isNonZero true, if basis function should not be set to zero, i.e. the evaluation will be done.
   * \tparam id the id.
   */
  template< bool hasQPs = true, bool hasIds = true, bool hasBasis = true, bool isNonZero = true, int id = -1 >
  struct Access
  {
    static constexpr bool hasQuadPoints = hasQPs;
    static constexpr bool hasBasisFunctions = hasBasis;
    static constexpr bool hasMultiValues = hasIds;
    static const int multiValue = id;

    /**
     * \brief Use this struct and typedefs therein to retrieve an Access class
     * where a basis function, a multivalue or a quadrature point is set.
     */
    struct Set
    {
      typedef Access<false,hasIds,hasBasis,isNonZero,id> QuadPoint;
      template< int i >
      using MultiValue = Access<hasQPs,false,hasBasis,isNonZero,i>;
      typedef Access<hasQPs,hasIds,false,false,id>     ZeroBasisFunction;
      typedef Access<hasQPs,hasIds,false,isNonZero,id> BasisFunction;
    };
  };


  /**
   * \brief Maps types of local evaluations arising from pass-based or non pass-based
   * frameworks to a default Access class.
   *
   * The type of a local evaluation is exchangable, but not completely arbritray.
   * Using Dune-Fem-DG, a certain structure is expected.
   * In order to map this certain structure to an Access class,
   * this class is defined via partial specialization.
   *
   * Use the `type` typedef to extract the correct Access type.
   *
   * \note This class realizises also cases of not already implemented
   * structures.
   *
   * \note This class tries to fill the gap between the Pass-based
   * approach used by instationary, matrix-free problems
   * and non Pass-based stationary problem.
   *
   *
   * Currently, there are three diffent types of local evaluations in
   * Dune-Fem-DG implemented:
   *
   * - Pass-based approach for (matrix-free) instationary problems
   * - A simple non Pass-based approach for stationary problems
   * - A non Pass-based approach for stationary problems
   *   (containing multiple values, accessible via an id,
   *   similar to the Pass-based approach)
   *
   * This non specialized class is the fallback implementation.
   * We can neither access basis functions, nor quadrature points,
   * nor the id of a multiple value.
   *
   * \note This class also implements the case: simple non Pass-based
   * for instationary (matrix-free) implementations.
   * Nevertheless, this case is currently not used by Dune-Fem-DG!
   */
  template< class EvalImp >
  struct DefaultAccess
  {
    //fallback: no, we cannot select anything
    typedef Access< false, false, false > type;
  };



  /**
   * \brief Specialization for a non Pass-based approach.
   *
   * This class knows the Access type for a simple non Pass-based approach.
   *
   * ![Simple non Pass-based access](local_eval_structure_q.png)
   *
   * \note This class also implements the case: simple non Pass-based
   * for instationary (matrix-free) implementations.
   * Nevertheless, this case is currently not used by Dune-Fem-DG!
   */
  template< class Arg >
  struct DefaultAccess<std::vector< Arg > >
  {
    //no pass (simple), 'stationary'
    typedef Access< true, false, false > type;
  };
  /**
   * \brief Specialization for a non Pass-based approach.
   *
   * This class knows the Access type for a simple non Pass-based approach.
   *
   * ![Simple non Pass-based access](local_eval_structure_qb.png)
   */
  template< class Arg >
  struct DefaultAccess<std::vector<std::vector<Arg> > >
  {
    //no pass (simple), 'stationary'
    typedef Access< true, false, true > type;
  };

  //no pass (simple), 'instationary'  [NOT USED!]
  //-> already specialized, s.a.
  //template< class Arg > struct DefaultAccess< Arg >;
  //no pass (simple), 'instationary'  [NOT USED!]
  //-> already specialized, s.a.
  //template< class Arg > struct DefaultAccess<std::vector<Arg> >;

  /**
   * \brief Specialization for a non Pass-based approach.
   *
   * This class knows the Access type for a non Pass-based approach
   * containing multiple values.
   *
   * ![Non Pass-based access](local_eval_structure_qib.png)
   */
  template< class Arg >
  struct DefaultAccess<std::vector<std::tuple<std::vector<Arg> > > >
  {
    //no pass, 'stationary'
    typedef Access< true, true, true > type;
  };
  /**
   * \brief Specialization for a non Pass-based approach.
   *
   * This class knows the Access type for a non Pass-based approach
   * containing multiple values.
   *
   * ![Non Pass-based access](local_eval_structure_ib.png)
   */
  template< class Arg >
  struct DefaultAccess<std::tuple<std::vector<Arg> > >
  {
    //no pass, 'stationary'
    typedef Access< false, true, true > type;
  };

  /**
   * \brief Specialization for a non Pass-based approach.
   *
   * This class knows the Access type for a non Pass-based approach
   * containing multiple values.
   *
   * ![Non Pass-based access](local_eval_structure_qi.png)
   *
   * \note This version is currently not used in Dune-Fem-DG.
   */
  template< class... Args >
  struct DefaultAccess<std::vector<std::tuple<Args...> > >
  {
    //no pass, 'instationary' [NOT USED!]
    typedef Access< true, true, false > type;
  };
  /**
   * \brief Specialization for a non Pass-based approach.
   *
   * This class knows the Access type for a non Pass-based approach
   * containing multiple values.
   *
   * ![Non Pass-based access](local_eval_structure_i.png)
   *
   * \note This version is currently not used in Dune-Fem-DG.
   */
  template< class... Args >
  struct DefaultAccess<std::tuple<Args...> >
  {
    //no pass, 'instationary' [NOT USED!]
    typedef Access< false, true, false > type;
  };


  /**
   * \brief Specialization for a Pass-based approach.
   *
   * This class knows the Access type for a Pass-based approach.
   *
   * ![Pass-based access](local_eval_structure_qi.png)
   */
  template< class... Args >
  struct DefaultAccess<std::vector<Dune::TypeIndexedTuple<Args...> > >
  {
    //pass, 'instationary'
    typedef Access< true, true, false > type;
  };
  /**
   * \brief Specialization for a Pass-based approach.
   *
   * This class knows the Access type for a Pass-based approach.
   *
   * ![Pass-based access](local_eval_structure_i.png)
   */
  template< class... Args >
  struct DefaultAccess<Dune::TypeIndexedTuple<Args...> >
  {
    //pass, 'instationary'
    typedef Access< false, true, false > type;
  };
  /**
   * \brief Specialization for a Pass-based approach.
   *
   * This class knows the Access type for a Pass-based approach.
   *
   * ![Pass-based access](local_eval_structure_qib.png)
   *
   * \note This version is currently not used in Dune-Fem-DG.
   */
  template< class... Args, class Arg >
  struct DefaultAccess<std::vector<Dune::TypeIndexedTuple<std::tuple<std::vector<Args>... >, Arg > > >
  {
    //pass, 'stationary' [NOT USED!]
    typedef Access< true, true, true > type;
  };
  /**
   * \brief Specialization for a Pass-based approach.
   *
   * This class knows the Access type for a Pass-based approach.
   *
   * ![Pass-based access](local_eval_structure_ib.png)
   *
   * \note This version is currently not used in Dune-Fem-DG.
   */
  template< class... Args, class Arg >
  struct DefaultAccess<Dune::TypeIndexedTuple<std::tuple<std::vector<Args>... >, Arg > >
  {
    //pass, 'stationary' [NOT USED!]
    typedef Access< false, true, true > type;
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
      //true, if VarId \in [0,#tuple-1]
      static const bool value = (VarId::value >= 0 && VarId::value < std::tuple_size<std::tuple<Args...> >::value );
    };

    /**
     * \brief This class evalu
     */
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
    int intersection() const { return 0; }

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
   * In many discretization schemes for Partial Differential Equations local data has to be evaluated.
   * This local data can be more complex than just the simple evaluation of a discrete function \f$ u\f$ in
   * some point \f$ x\in\Omega \f$. In order to provide all needed local data, we have developed this class
   * which should provide a flexible structure.
   *
   * A local evaluation simply knows, _how_ a complex local data structure has to be treated and evaluated
   * and adds a certain structure to this evaluation.
   *
   * In a mathematical way, the structure of a local evaluation \f$ \mathsf{E}\f$ can be described
   * as a function, defined by
   *
   * \f{eqnarray*}{
   *   \mathsf{E}:\mathcal{Q} \times \mathcal{I} \times \mathcal{B} \rightarrow \mathbb{K}^n
   * \f}
   *
   * with \f$ \mathcal{Q}=[0,\ldots,\text{#quadPoints}] \f$, \f$ \mathcal{I}=\mathbb{N}\f$ and
   * \f$ \mathcal{B}=[0,\ldots,\text{#basisfkt}-1] \f$.
   *
   *
   *
   * This class was introduced to allow more flexible interfaces for classes
   * needing a local evaluation of a discrete function.
   */

  template< class QuadratureContextImp, class RangeType, class JacobianType = RangeType, class AccessImp = typename DefaultAccess<RangeType>::type >
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
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::Set::QuadPoint > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, idx, basis_ );
    }

    //access: quadpoints II
    decltype(auto) operator[]( unsigned long int idx ) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::Set::QuadPoint > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, idx, basis_ );
    }

    //access: multi values, i.e. tuples
    template< class Int, Int id >
    decltype(auto) operator[]( const std::integral_constant<Int,id>& idx ) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::Set::template MultiValue<id> > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, qp_, basis_ );
    }

    //access: basis functions, zero
    decltype(auto) operator()(int i=0) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::Set::ZeroBasisFunction > NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, qp_, -1 );
    }

    //access: basis functions
    decltype(auto) operator[]( int idx ) const
    {
      typedef LocalEvaluation<QuadratureContextImp, RangeType, JacobianType, typename AccessImp::Set::BasisFunction> NewContextType;
      return NewContextType( quadImp_, values_, jacobians_, qp_, idx );
    }

    //access: failure
    template< class Fail >
    int operator[]( Fail ) const
    {
      static_assert( static_fail<Fail>::value, "This should not happen. Please check the exact type of your argument." );
      return 0;
    }

    const Evaluator<RangeType> values() const
    {
      return Evaluator<RangeType>( values_, qp_, basis_ );
    }

    const Evaluator<JacobianType> jacobians() const
    {
      return Evaluator<JacobianType>( jacobians_, qp_, basis_ );
    }

    template <class Functor, class ... Args>
    decltype(auto) values( const Functor& functor, const Args& ... args ) const
    {
      return FunctorContext::Evaluate< Functor, FunctorContext::Contains< RangeType, typename Functor::VarId >::value>::eval( values_, functor, args ... );
    }

    template <class Functor, class ... Args>
    decltype(auto) jacobians( const Functor& functor, const Args& ... args ) const
    {
      return FunctorContext::Evaluate< Functor, FunctorContext::Contains< JacobianType, typename Functor::VarId >::value>::eval( jacobians_, functor, args ... );
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
