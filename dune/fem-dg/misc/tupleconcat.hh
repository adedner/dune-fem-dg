#ifndef DUNE_TUPLE_CONCAT_HH
#define DUNE_TUPLE_CONCAT_HH

#include <tuple>

template <class... Args>
auto func(Args&& ... args) -> decltype(std::tuple_cat(args...))
{ return std::tuple_cat(args...); }

template <class... Args>
struct tuple_concat {
    typedef decltype(func(std::declval<Args>()...)) type;
};


#endif
