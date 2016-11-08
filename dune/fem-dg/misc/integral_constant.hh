#ifndef DUNE_FEM_INTEGRAL_CONSTANT_HH
#define DUNE_FEM_INTEGRAL_CONSTANT_HH

namespace Dune
{
namespace Fem
{
  // Global Type Definitions
  // -----------------------
  template< unsigned long int i>
  using _index = std::integral_constant< unsigned long int, i >;

  template< unsigned long int... i>
  using _indices = std::integer_sequence< unsigned long int, i... >;

  typedef _index< 0> __0;
  typedef _index< 1> __1;
  typedef _index< 2> __2;
  typedef _index< 3> __3;
  typedef _index< 4> __4;
  typedef _index< 5> __5;
  typedef _index< 6> __6;
  typedef _index< 7> __7;
  typedef _index< 8> __8;
  typedef _index< 9> __9;
  typedef _index<10> __10;
  typedef _index<11> __11;
  typedef _index<12> __12;
  typedef _index<13> __13;
  typedef _index<14> __14;
  typedef _index<15> __15;
  typedef _index<16> __16;
  typedef _index<17> __17;
  typedef _index<18> __18;
  typedef _index<19> __19;


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

}
}

#endif // #ifndef DUNE_FEM_DG_SIMULATOR_HH
