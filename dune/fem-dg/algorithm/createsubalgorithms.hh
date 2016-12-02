#ifndef DUNE_FEMDG_COUPLEDALGORITHMS_HH
#define DUNE_FEMDG_COUPLEDALGORITHMS_HH

#include <tuple>
#include <type_traits>
#include <memory>
#include <utility>
#include <iostream>
#include <dune/fem-dg/misc/integral_constant.hh>

namespace Dune
{
namespace Fem
{

  /**
   *  \brief Creates a tuple of uncoupled sub algorithms
   *
   *  \tparam SubAlgorithmImp
   */
  template< class... SubAlgorithmsImp >
  class CreateSubAlgorithms
  {
  public:
    typedef std::tuple< std::shared_ptr< SubAlgorithmsImp > ... >    SubAlgorithmTupleType;
  private:
    static_assert( std::tuple_size<SubAlgorithmTupleType>::value > 0, "We need at least one Sub-Algorithm" );

    template< int i >
    using Element = typename std::tuple_element< i, SubAlgorithmTupleType >::type::element_type;

    template< unsigned long int ...i, class GlobalContainerImp  >
    static SubAlgorithmTupleType apply ( std::index_sequence< i... >, std::shared_ptr< GlobalContainerImp > cont )
    {
      return std::make_tuple( std::make_shared<Element<i> >( cont->sub( _index<i>() ) )... );
    }
  public:

    template< class GlobalContainerImp >
    static SubAlgorithmTupleType apply ( std::shared_ptr< GlobalContainerImp > cont  )
    {
      return apply( typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >(), cont );
    }
  };



} // namespace Fem
} // namespace Dune

#endif
