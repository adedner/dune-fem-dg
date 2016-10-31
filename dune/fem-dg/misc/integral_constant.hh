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

  _index< 0> _0;
  _index< 1> _1;
  _index< 2> _2;
  _index< 3> _3;
  _index< 4> _4;
  _index< 5> _5;
  _index< 6> _6;
  _index< 7> _7;
  _index< 8> _8;
  _index< 9> _9;
  _index<10> _10;
  _index<11> _11;
  _index<12> _12;
  _index<13> _13;
  _index<14> _14;
  _index<15> _15;
  _index<16> _16;
  _index<17> _17;
  _index<18> _18;
  _index<19> _19;




}
}

#endif // #ifndef DUNE_FEM_DG_SIMULATOR_HH
