#ifndef DUNE_TYPEDEF_CHECK_HH
#define DUNE_TYPEDEF_CHECK_HH


#define CHECK_TYPEDEF_EXISTS( tname )                                                                   \
    template< class T, bool exists >                                                                    \
    struct tname##_Helper                                                                               \
    { typedef void type; };                                                                             \
    template< class T  >                                                                                \
    struct tname##_Helper< T, true >                                                                    \
    { typedef typename T::tname type; };                                                                \
    template< class T >                                                                                 \
    struct tname##s {                                                                                   \
      template <typename TT>                                                                            \
      static auto apply(TT const&) -> decltype( std::declval<typename TT::tname>(), std::true_type()) { \
        return std::true_type();                                                                        \
      }                                                                                                 \
      static std::false_type apply(...) { return std::false_type(); }                                   \
      typedef typename tname##_Helper< T, decltype( apply( std::declval<T>() ) )::value >::type type;   \
    };                                                                                                  \




//// SFINAE test
//template<class T>
//class has_type_type
//{
//    typedef char one;
//    typedef long two;
//
//    template <class C> static one test( decltype(T::type*) ) ;
//    template <typename C> static two test(...);
//
//public:
//    enum { value = sizeof(test<T>(0)) == sizeof(char) };
//};
//
//
//
//template< typename T, typename Check = void > struct HasXYZ
//{ static const bool value = false; };
//
//template< typename T > struct HasXYZ< T, typename T::XYZ >
//{ static const bool value = true; };
//
//
//template <typename T, typename Check = void>
//struct A
//{ static const bool value = false; };
//
//template <typename T>
//struct A<T,void>
//{ static const bool value = true; };


#endif
