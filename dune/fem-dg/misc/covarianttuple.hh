#ifndef DUNE_COVARIANT_TUPLE_HH
#define DUNE_COVARIANT_TUPLE_HH

#include <tuple>
#include "tupleutility.hh"

namespace Dune
{
namespace Fem
{
  /**
   *  \brief This class makes use of the "Covariant Return Types" feature in C++ for a std::tuple<>.
   *
   *  The following code will not compile, unless DerivedObjectType is really derived from BaseObjectType.
   *  \code{.cpp}
   *    typedef WhatEver BaseObjectType;
   *    typedef WhatElse DerivedObjectType;
   *
   *    struct Base
   *    {
   *      virtual BaseObjectType* getObject()
   *      {
   *        return new BaseObjectType();
   *      }
   *    };
   *
   *    struct Derived
   *      : public Base
   *    {
   *      virtual DerivedObjectType* getObject()
   *      {
   *        return new DerivedObjectType();
   *      }
   *    };
   *  \endcode
   *
   *  With this background, we want to allow to "enlarge" existing tuple types in derived classes
   *  through a small wrapper.
   *
   *  Access to real tuple type should be done by using the dereferencing operator.
   *
   *  Here is an easy example for a pointer of tuples:
   *
   *  \code{.cpp}
   *    typedef std::tuple< int* > BaseObjectType;
   *    typedef std::tupel< int*, double* > DerivedObjectType;
   *
   *    struct Base
   *    {
   *      int intData_;
   *
   *      Base()
   *       : intData_( 42 )
   *      {}
   *
   *      virtual BaseObjectType* getObject()
   *      {
   *        return new BaseObjectType( std::make_tuple( intData_ ) );
   *      }
   *    };
   *
   *    struct Derived
   *      : public Base
   *    {
   *      double doubleData_;
   *
   *      Derived()
   *        : Base(),
   *          doubleData_( 3.14 )
   *      {}
   *
   *      virtual DerivedObjectType* getObject()
   *      {
   *        return new DerivedObjectType( *getObject(), std::make_tuple( doubleData ) );
   *      }
   *    };
   *  \endcode
   *
   *  \note The first parameter is always the Type of the base class. Do not mix it up!
   */
  template< class ... TupleArgs >
  class CovariantTuple;

  /**
   * \brief Specialization for the base classes using CovariantTuple.
   *
   * \tparam Tuple A std::tuple which should be wrapped.
   */
  template< class Tuple >
  class CovariantTuple< Tuple >
    : public Tuple
  {
    static_assert( is_tuple< Tuple >::value, "CovariantTuple is expecting a tuple" );

  public:
    //! type of the wrapped tuple
    typedef Tuple type;

    CovariantTuple( type t )
      : Tuple( t ),
        t_( this )
    {}

    /**
     * \brief Get access to the real tuple
     */
    type& operator* ()
    {
      return *t_;
    }

    /**
     * \brief Get access to the real tuple
     */
    type* operator-> ()
    {
      return &t_;
    }

  protected:
    type* t_;
  };

  /**
   * \brief Specialization for the derived classes which wants to enlarge an existing method of the base class.
   *
   * \tparam TupleBase type of the tuple returned by a base class
   * \tparam Tuples variadic number of tuples which should be added to the tuple of the base class to
   *         extend a virtual method in a derived class.
   */
  template< class TupleBase, class... Tuples >
  class CovariantTuple< TupleBase, Tuples... >
    : public CovariantTuple< TupleBase >
  {
    typedef CovariantTuple< TupleBase > BaseType;

  public:
    //! type of the concatenated, wrapped tuple
    typedef tuple_concat_t< TupleBase, Tuples... > type;

    /**
     *  \brief Constructor
     *
     *  \param t Tuple for the base class
     *  \param tuples arbritary number of other tuples which should be concatenated.
     */
    CovariantTuple( TupleBase t, Tuples... tuples )
      : BaseType( t ),
        t_( tuple_cat( *BaseType::t_, tuples... ) )
    {}

    /**
     * \brief Get access to the real, concatenated tuple
     */
    type& operator* ()
    {
      return t_;
    }

    /**
     * \brief Get access to the real, concatenated tuple
     */
    type* operator-> ()
    {
      return &t_;
    }

  private:
    type t_;

  };

}
}
#endif
