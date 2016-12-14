#ifndef DUNE_FEM_INTEGRAL_CONSTANT_HH
#define DUNE_FEM_INTEGRAL_CONSTANT_HH

#include <type_traits>
#include <utility>
#include <string>

namespace Dune
{
namespace Fem
{
  // Global Type Definitions
  // -----------------------

  /**
   * \brief shortcut for std::integral_constant<> used by container classes.
   *
   * \ingroup Container
   */
  template< unsigned long int i>
  using _index = std::integral_constant< unsigned long int, i >;


  template< bool b >
  using _bool = std::integral_constant< bool, b >;

  //! shortcut for std::integral_constant<bool,false>
  typedef _bool<false> __f;

  //! shortcut for std::integral_constant<bool,true>
  typedef _bool<true> __t;

  /**
   * \brief shortcut for std::integer_sequence<> used by container classes.
   *
   * \ingroup Container
   */
  template< unsigned long int... i>
  using _indices = std::integer_sequence< unsigned long int, i... >;

  //! shortcut for std::integral_constant<>
  typedef _index< 0> __0;
  //! shortcut for std::integral_constant<>
  typedef _index< 1> __1;
  //! shortcut for std::integral_constant<>
  typedef _index< 2> __2;
  //! shortcut for std::integral_constant<>
  typedef _index< 3> __3;
  //! shortcut for std::integral_constant<>
  typedef _index< 4> __4;
  //! shortcut for std::integral_constant<>
  typedef _index< 5> __5;
  //! shortcut for std::integral_constant<>
  typedef _index< 6> __6;
  //! shortcut for std::integral_constant<>
  typedef _index< 7> __7;
  //! shortcut for std::integral_constant<>
  typedef _index< 8> __8;
  //! shortcut for std::integral_constant<>
  typedef _index< 9> __9;
  //! shortcut for std::integral_constant<>
  typedef _index<10> __10;
  //! shortcut for std::integral_constant<>
  typedef _index<11> __11;
  //! shortcut for std::integral_constant<>
  typedef _index<12> __12;
  //! shortcut for std::integral_constant<>
  typedef _index<13> __13;
  //! shortcut for std::integral_constant<>
  typedef _index<14> __14;
  //! shortcut for std::integral_constant<>
  typedef _index<15> __15;
  //! shortcut for std::integral_constant<>
  typedef _index<16> __16;
  //! shortcut for std::integral_constant<>
  typedef _index<17> __17;
  //! shortcut for std::integral_constant<>
  typedef _index<18> __18;
  //! shortcut for std::integral_constant<>
  typedef _index<19> __19;
  //! shortcut for std::integral_constant<>
  typedef _index<20> __20;
  //! shortcut for std::integral_constant<>
  typedef _index<21> __21;
  //! shortcut for std::integral_constant<>
  typedef _index<22> __22;
  //! shortcut for std::integral_constant<>
  typedef _index<23> __23;
  //! shortcut for std::integral_constant<>
  typedef _index<24> __24;
  //! shortcut for std::integral_constant<>
  typedef _index<25> __25;
  //! shortcut for std::integral_constant<>
  typedef _index<26> __26;
  //! shortcut for std::integral_constant<>
  typedef _index<27> __27;
  //! shortcut for std::integral_constant<>
  typedef _index<28> __28;
  //! shortcut for std::integral_constant<>
  typedef _index<29> __29;
  //! shortcut for std::integral_constant<>
  typedef _index<30> __30;
  //! shortcut for std::integral_constant<>
  typedef _index<31> __31;
  //! shortcut for std::integral_constant<>
  typedef _index<32> __32;
  //! shortcut for std::integral_constant<>
  typedef _index<33> __33;
  //! shortcut for std::integral_constant<>
  typedef _index<34> __34;
  //! shortcut for std::integral_constant<>
  typedef _index<35> __35;
  //! shortcut for std::integral_constant<>
  typedef _index<36> __36;
  //! shortcut for std::integral_constant<>
  typedef _index<37> __37;
  //! shortcut for std::integral_constant<>
  typedef _index<38> __38;
  //! shortcut for std::integral_constant<>
  typedef _index<39> __39;
  //! shortcut for std::integral_constant<>
  typedef _index<40> __40;
  //! shortcut for std::integral_constant<>
  typedef _index<41> __41;
  //! shortcut for std::integral_constant<>
  typedef _index<42> __42;
  //! shortcut for std::integral_constant<>
  typedef _index<43> __43;
  //! shortcut for std::integral_constant<>
  typedef _index<44> __44;
  //! shortcut for std::integral_constant<>
  typedef _index<45> __45;
  //! shortcut for std::integral_constant<>
  typedef _index<46> __46;
  //! shortcut for std::integral_constant<>
  typedef _index<47> __47;
  //! shortcut for std::integral_constant<>
  typedef _index<48> __48;
  //! shortcut for std::integral_constant<>
  typedef _index<49> __49;


  __0  _0;
  __1  _1;
  __2  _2;
  __3  _3;
  __4  _4;
  __5  _5;
  __6  _6;
  __7  _7;
  __8  _8;
  __9  _9;
  __10 _10;
  __11 _11;
  __12 _12;
  __13 _13;
  __14 _14;
  __15 _15;
  __16 _16;
  __17 _17;
  __18 _18;
  __19 _19;
  __20 _20;
  __21 _21;
  __22 _22;
  __23 _23;
  __24 _24;
  __25 _25;
  __26 _26;
  __27 _27;
  __28 _28;
  __29 _29;
  __30 _30;
  __31 _31;
  __32 _32;
  __33 _33;
  __34 _34;
  __35 _35;
  __36 _36;
  __37 _37;
  __38 _38;
  __39 _39;
  __40 _40;
  __41 _41;
  __42 _42;
  __43 _43;
  __44 _44;
  __45 _45;
  __46 _46;
  __47 _47;
  __48 _48;
  __49 _49;

  //some easy print routines
  template< unsigned long int i >
  static std::string print( _index<i> index )
  {
    return std::to_string( i );
  }

  template< unsigned long int i, unsigned long int... is >
  static std::string print( _index<i> index, _index<is>... indices )
  {
    return std::to_string( i ) + " | " + print( indices... );
  }

  template< unsigned long int... is >
  static std::string print( std::tuple< _index<is>... > )
  {
    return print( _index<is>()... );
  }

  template< unsigned long int... is, unsigned long int... js >
  static std::string print( std::tuple< _index<is>... >, std::tuple< _index<js>... > )
  {
    return print( _index<is>()..., _index<js>()... );
  }
}
}

#endif // #ifndef DUNE_FEM_DG_SIMULATOR_HH
