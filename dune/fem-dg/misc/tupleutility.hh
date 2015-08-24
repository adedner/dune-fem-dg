#ifndef DUNE_TUPLE_CONCAT_HH
#define DUNE_TUPLE_CONCAT_HH

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



#endif
