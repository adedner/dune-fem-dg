#ifndef DUNE_COVARIANT_TUPLE_HH
#define DUNE_COVARIANT_TUPLE_HH

#include <tuple>
#include "tupleutility.hh"

 template< class ... TupleArgs >
  class CovariantTuple;


  template< class Tuple >
  class CovariantTuple< Tuple >
    : public Tuple
  {
    static_assert( is_tuple< Tuple >::value, "CovariantTuple is expecting a tuple" );

  public:
    typedef Tuple type;

    CovariantTuple( type t )
      : Tuple( t ),
        t_( this )
    {}

    type& operator* ()
    {
      return *t_;
    }

    type* operator-> ()
    {
      return &t_;
    }

  protected:
    type* t_;
  };


  template< class TupleBase, class... Tuples >
  class CovariantTuple< TupleBase, Tuples... >
    : public CovariantTuple< TupleBase >
  {
    typedef CovariantTuple< TupleBase > BaseType;

  public:
    typedef typename tuple_concat< TupleBase, Tuples... >::type type;

    CovariantTuple( TupleBase t, Tuples... tuples )
      : BaseType( t ),
        t_( tuple_cat( *BaseType::t_, tuples... ) )
    {}

    type& operator* ()
    {
      return t_;
    }

    type* operator-> ()
    {
      return &t_;
    }

  private:
    type t_;

  };


#endif
