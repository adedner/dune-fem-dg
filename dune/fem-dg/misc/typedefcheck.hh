#ifndef DUNE_TYPEDEF_CHECK_HH
#define DUNE_TYPEDEF_CHECK_HH


#define CHECK_TYPEDEF_EXISTS( tname )                                                                \
    template< class T, bool exists >                                                                \
    struct tname##_Helper                                                                             \
    { typedef void type; };                                                                         \
    template< class T  >                                                                            \
    struct tname##_Helper< T, true >                                                                  \
    { typedef typename T::tname type; };                                                           \
    template< class T >                                                                             \
    struct tname##s {                                                                                \
      template <typename TT>                                                                        \
      static auto apply(TT const&) -> decltype( typename TT::tname(), std::true_type()) {          \
        return std::true_type();                                                                    \
      }                                                                                             \
      static std::false_type apply(...) { return std::false_type(); }                               \
      typedef typename tname##_Helper< T, decltype( apply( std::declval<T>() ) )::value >::type type; \
    };                                                                                              \

#endif
