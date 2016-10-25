#ifndef DUNE_TUPLE_CONCAT_HH
#define DUNE_TUPLE_CONCAT_HH

#include <tuple>
#include <type_traits>
#include <utility>

namespace Dune
{
namespace Fem
{

  /**
   * \brief concats tuples and flattens the result
   *
   * \tparam Args variadic list of tuples
   *
   * Usage
   * \code{.cpp}
   *   typedef std::tuple< double > FirstTuple;
   *   typedef std::tuple< int >    SecondTuple;
   *   typedef std::tuple< bool >   ThirdTuple;
   *
   *   typedef typename tuple_concat< FirstTuple, SecondTuple, ThirdTuple >::type ConcatTuple;
   *   //type of ConcatTuple is std::tuple< double, int, bool >
   * \endcode
   */
  template < class... Args >
  class tuple_concat
  {
    template < class... Args2 >
    static auto apply(Args2&& ... args2) -> decltype(std::tuple_cat(args2...))
    {
      return std::tuple_cat(args2...);
    }
    public:
    //! type of a concatenated, flattened tuple
    typedef decltype(apply(std::declval<Args>()...)) type;
  };


  /**
   * \brief Extracts the first element out of a variadic list of Arguments
   *
   * \tparam Args a variadic lists of elements.
   */
  template < class... Args >
  struct tuple_head
  {
    static_assert( sizeof ... ( Args ) > 0,
                   "tuple_head<> needs at least one template argument" );
    typedef typename std::tuple_element< 0, std::tuple< Args... > >::type type;
  };


  /**
   * \brief typedef check if a type is a std::tuple<>
   *
   * usage
   * \code{.cpp}
   *   typedef
   *   static const int =
   * \endcode
   *
   * \tparam T type to be checked.
   */
  template< class T >
  struct is_tuple
  {
    static const int value = false;
  };

  template< class... TElem >
  struct is_tuple< std::tuple< TElem... > >
  {
    static const int value = true;
  };



  template< class FullTupleImp, class IndexSequenceImp >
  struct tuple_reducer;

  template< class FullTupleImp, std::size_t... i >
  struct tuple_reducer< FullTupleImp, std::index_sequence< i... > >
  {
    typedef std::tuple< typename std::tuple_element< i, FullTupleImp >::type... > type;

    static type apply( const FullTupleImp& tuple )
    {
      return std::make_tuple( std::get< i >( tuple )... );
    }

  };

  template< class TupleImp >
  struct drop_tuple;

  template< class Tuple >
  struct drop_tuple< std::tuple< Tuple > >
  {
    static_assert( std::is_same<std::tuple< Tuple >,std::tuple< Tuple > >::value ,
                   "We do not unwrap a simple tuple to avoid unexpected behaviour!" );
  };

  template< class... TupleArg >
  struct drop_tuple< std::tuple< std::tuple< TupleArg... > > >
  {
    typedef std::tuple< TupleArg... > type;
  };
}
}

#endif
