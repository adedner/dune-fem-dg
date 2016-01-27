#ifndef DUNE_FEMDG_COUPLEDALGORITHMS_HH
#define DUNE_FEMDG_COUPLEDALGORITHMS_HH

namespace Dune
{
namespace Fem
{


  /**
   *  \brief Creates a tuple of uncoupled sub algorithms
   *
   *  \note This is the default implementation which should be used
   *  for uncoupled algorithms.
   *
   *  \tparam SubAlgorithmTupleImp
   *  \tparam GridImp
   */
  template< class SubAlgorithmTupleImp, class GridImp >
  class UncoupledSubAlgorithms
  {
    typedef SubAlgorithmTupleImp  SubAlgorithmTupleType;
    typedef GridImp               GridType;

    template< int i >
    using Element = typename std::remove_pointer< typename std::tuple_element< i, SubAlgorithmTupleType >::type >::type;


    template< std::size_t i >
    static Element<i>* createSubAlgorithm( GridType& grid )
    {
      static typename Element<i>::ContainerType container( grid );
      return new Element<i>( grid, container );
    }

    template< std::size_t ...i >
    static SubAlgorithmTupleType apply ( Std::index_sequence< i... >, GridType &grid )
    {
      return std::make_tuple( createSubAlgorithm<i>( grid )... );
    }

  public:
    /**
     * \brief Creates a tuple of sub algorithms
     *
     * \param grid the grid which is needed to construct the algorithms
     * \return The constructed tuple of algorithms
     */
    static SubAlgorithmTupleType apply ( GridType &grid )
    {
      return apply( typename Std::make_index_sequence_impl< std::tuple_size< SubAlgorithmTupleType >::value >::type(), grid );
    }
  };


} // namespace Fem
} // namespace Dune

#endif
