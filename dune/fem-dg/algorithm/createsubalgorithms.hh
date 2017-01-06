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
    typedef typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value > sequence;

    static_assert( std::tuple_size<SubAlgorithmTupleType>::value > 0, "We need at least one Sub-Algorithm" );

    template< int i >
    using SubAlgorithm = typename std::tuple_element_t< i, SubAlgorithmTupleType >::element_type;

    template< unsigned long int ...i, class GlobalContainerImp  >
    static decltype(auto) apply ( std::index_sequence< i... >, const std::shared_ptr< GlobalContainerImp >& cont )
    {
      return std::make_tuple( std::make_shared<SubAlgorithm<i> >( cont->sub( _index<i>() ), cont->extra( _index<i>() ) )... );
    }


    template< class GridImp >
    static decltype(auto) gridItem ( const GridImp& grid )
    {
      return const_cast<typename std::remove_const<GridImp>::type& >( grid );
    }

    template< unsigned long int ...i, class GlobalContainerImp  >
    static decltype(auto) grids ( std::index_sequence< i... >, const std::shared_ptr< GlobalContainerImp >& cont )
    {
      return std::tie( gridItem( (*(cont->sub(_index<i>() )))(_0)->solution()->gridPart().grid() )... );
    }
  public:

    template< class GlobalContainerImp >
    static decltype(auto) apply ( std::shared_ptr< GlobalContainerImp > cont )
    {
      return apply( sequence(), cont );
    }

    template< class GlobalContainerImp >
    static decltype(auto) grids ( std::shared_ptr< GlobalContainerImp > cont )
    {
      return grids( sequence(), cont );
    }
  };



} // namespace Fem
} // namespace Dune

#endif
