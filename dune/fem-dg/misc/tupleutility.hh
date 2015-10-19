#ifndef DUNE_TUPLE_CONCAT_HH
#define DUNE_TUPLE_CONCAT_HH

#include <tuple>
#include <type_traits>

template < class... Args >
class tuple_concat
{
  template < class... Args2 >
  static auto apply(Args2&& ... args2) -> decltype(std::tuple_cat(args2...))
  {
    return std::tuple_cat(args2...);
  }
  public:
  typedef decltype(apply(std::declval<Args>()...)) type;
};

template < class... Args >
struct tuple_head
{
  static_assert( sizeof ... ( Args ) > 0,
                 "tuple_head<> needs at least one template argument" );
  typedef typename std::tuple_element< 0, std::tuple< Args... > >::type type;
};


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

template< class Tuple, class T >
struct has_type
{ static_assert( is_tuple<Tuple>::value, "First template parameter is not of type std::tuple<...>" ); };

template< class T >
struct has_type< std::tuple<>, T > : std::false_type {};

template< typename U, typename... Ts, class T>
struct has_type< std::tuple<U, Ts...>, T > : has_type<T, std::tuple< Ts... > > {};

template< class... Ts, class T>
struct has_type< std::tuple<T, Ts...>, T > : std::true_type {};


template< class T >
struct wrap_tuple
{
  typedef typename std::conditional< is_tuple< T >::value,
                                     T,
                                     std::tuple<T>
                                      >::type                 type;
};

//
//template< class T1, class T2 >
//class tuple_check
//{
//  static_assert( false, "no tuple inserted" );
//};
//
//template< class T1, class T2 >
//class tuple_check
//{
//  static_assert( is_tuple< T1 >::value, "First Parameter is no std::tuple" );
//  static_assert( is_tuple< T2 >::value, "Second Parameter is no std::tuple" );
//
//
//  public:
//  static_assert( std::tuple_size< T1 >::value == std::tuple_size< T2 >::value,
//                 "size of tuple does not fit" );
//
//  static_assert( ForLoop< T1, 0, std::tuple_size< T1 >::value >
//                 "ddd" );
//
//};

#endif
