#ifndef DUNE_FEM_STATIC_WARNING_HH
#define DUNE_FEM_STATIC_WARNING_HH

//solution taken from stackoverflow.com by "Managu"
//uses deprecation warning attribute

#if defined(__GNUC__)
#define DEPRECATE(foo, msg) foo __attribute__((deprecated(msg)))
#warning "gnu"
#elif defined(_MSC_VER)
#define DEPRECATE(foo, msg) __declspec(deprecated(msg)) foo
#warning "msc"
#else
#error This compiler is not supported
#endif

#define PP_CAT(x,y) PP_CAT1(x,y)
#define PP_CAT1(x,y) x##y

namespace detailed
{
    struct true_type {};
    struct false_type {};
    template <int test> struct converter : public true_type {};
    template <> struct converter<0> : public false_type {};
}

#define STATIC_WARNING(cond, msg) \
struct PP_CAT(static_warning,__LINE__) { \
  DEPRECATE(void _(::detailed::false_type const& ),msg) {}; \
  void _(::detailed::true_type const& ) {}; \
  PP_CAT(static_warning,__LINE__)() {_(::detailed::converter<(cond)>());} \
}

// Note: using STATIC_WARNING_TEMPLATE changes the meaning of a program in a small way.
// It introduces a member/variable declaration.  This means at least one byte of space
// in each structure/class instantiation.  STATIC_WARNING should be preferred in any
// non-template situation.
//  'token' must be a program-wide unique identifier.
#define STATIC_WARNING_TEMPLATE(token, cond, msg) \
    STATIC_WARNING(cond, msg) PP_CAT(PP_CAT(_localvar_, token),__LINE__)

//Example usage:
//
//#line 1
//STATIC_WARNING(1==2, "Failed with 1 and 2");
//STATIC_WARNING(1<2, "Succeeded with 1 and 2");
//
//struct Foo
//{
//  STATIC_WARNING(2==3, "2 and 3: oops");
//  STATIC_WARNING(2<3, "2 and 3 worked");
//};
//
//void func()
//{
//  STATIC_WARNING(3==4, "Not so good on 3 and 4");
//  STATIC_WARNING(3<4, "3 and 4, check");
//}
//
//template <typename T> struct wrap
//{
//  typedef T type;
//  STATIC_WARNING(4==5, "Bad with 4 and 5");
//  STATIC_WARNING(4<5, "Good on 4 and 5");
//  STATIC_WARNING_TEMPLATE(WRAP_WARNING1, 4==5, "A template warning");
//};


#endif
