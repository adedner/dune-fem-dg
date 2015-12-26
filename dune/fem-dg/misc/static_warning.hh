#ifndef DUNE_FEM_STATIC_WARNING_HH
#define DUNE_FEM_STATIC_WARNING_HH

//solution taken and modified from stackoverflow.com by "Managu"
//uses deprecation warning attribute

#if defined(__GNUC__)
#define DEPRECATE(foo, msg) foo __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
#define DEPRECATE(foo, msg) __declspec(deprecated(msg)) foo
#endif

#if defined(__GNUC__) || defined(_MSC_VER)

#define PP_CAT(x,y) PP_CAT1(x,y)
#define PP_CAT1(x,y) x##y

namespace detailed
{
    struct true_type {};
    struct false_type {};
    template <int test> struct converter : public true_type {};
    template <> struct converter<0> : public false_type {};
}

#define STATIC_WARNING_IMPL(cond, msg, id) \
struct PP_CAT(static_warning,id) { \
  DEPRECATE(void _(::detailed::false_type const& ),msg) {}; \
  void _(::detailed::true_type const& ) {}; \
  PP_CAT(static_warning,id)() {_(::detailed::converter<(cond)>());} \
}

#define STATIC_WARNING(cond, msg, id) \
    STATIC_WARNING_IMPL(cond, msg,id) PP_CAT(_localvar_,id)

#else
#define STATIC_WARNING(cond, msg, id)
#endif

/**
 * \def static_warning(cond,msg)
 * Same as static_assert(), but generates a warning instead of an error.
 *
 * \note This macro uses the deprecated attribute to generate warnings.
 * Thus the -Wno-deprecated-declaration compiler option has to be turned off.
 *
 * Example usage:
 * \code{.cpp}
 * #line 1
 * static_warning(1==2, "Failed with 1 and 2");
 * static_warning(1<2, "Succeeded with 1 and 2");
 *
 * struct Foo
 * {
 *   static_warning(2==3, "2 and 3: oops");
 *   static_warning(2<3, "2 and 3 worked");
 * };
 *
 * void func()
 * {
 *   static_warning(3==4, "Not so good on 3 and 4");
 *   static_warning(3<4, "3 and 4, check");
 * }
 *
 * template <typename T> struct wrap
 * {
 *   typedef T type;
 *   static_warning(4==5, "Bad with 4 and 5");
 *   static_warning(4<5, "Good on 4 and 5");
 * };
 *
 * template struct wrap<int>;
 * \endcode
 *
 * \note Only supported for gcc...
 */
#define static_warning(cond, msg) \
    STATIC_WARNING(cond, msg, __COUNTER__)


#endif
