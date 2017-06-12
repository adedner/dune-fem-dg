#ifndef DUNE_FUNCTION_CHECK_HH
#define DUNE_FUNCTION_CHECK_HH

#define CHECK_FUNCTION_EXISTS( funcName )                                 \
  template <typename T>                                                   \
  class has_member_function_##funcName                                             \
  {                                                                       \
      typedef char one;                                                   \
      typedef long two;                                                   \
      template <typename C> static one test( decltype(&C::funcName) );    \
      template <typename C> static two test(...);                         \
  public:                                                                 \
      enum { value = sizeof(test<T>(0)) == sizeof(char) };                \
  };

#endif
