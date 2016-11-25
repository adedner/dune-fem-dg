#ifndef DUNE_FEMDG_CONTAINER_HH
#define DUNE_FEMDG_CONTAINER_HH

#include <memory>


namespace Dune
{
namespace Fem
{

  template< class... >
  struct SubOrderSelect;

  template< class OrderHead, class... Args >
  struct SubOrderSelect< std::tuple<OrderHead>,std::tuple<Args...> >
  {
    typedef std::tuple< typename std::tuple_element<OrderHead::value,std::tuple<Args...> >::type > type;
  };

  template< class OrderHead, class... OrderArgs, class... Args >
  struct SubOrderSelect< std::tuple<OrderHead, OrderArgs...>,std::tuple<Args...> >
  {
    typedef typename tuple_concat< std::tuple< typename std::tuple_element<OrderHead::value,std::tuple<Args...> >::type >,
                                   typename SubOrderSelect< std::tuple<OrderArgs...>, std::tuple<Args...> >::type >::type
                                                                                                                     type;
  };



  /**
   * \brief Defines a global container
   *
   * \ingroup Container
   *
   * Before reading this summary, you should have read \ref Container.
   *
   * ![A global container](container7.png)
   */
  template< class Item2TupleImp, class Item1TupleImp, class SubOrderRowImp, class SubOrderColImp, class... DiscreteFunctions >
  class GlobalContainer
  {
    typedef TwoArgContainer< ArgContainerArgWrapper< Item2TupleImp >::template _t2Inv, ArgContainerArgWrapper< Item1TupleImp >::template _t1,
                             std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > > ContainerType;

    static_assert( is_tuple< SubOrderRowImp >::value, "SubOrderRowImp has to be a std::tuple<std::tuple<>...>" );
  public:

    /**
     * \brief constructor. Delegates everything
     */
    template< class GridImp, class ... Args>
    GlobalContainer( const std::shared_ptr< GridImp >& gridptr, Args&&... args )
    : cont_( gridptr, std::forward<Args>(args)... )
    {}

    /**
     * \brief returns the i's container
     */
    template< unsigned long int i >
    decltype(auto) sub( _index<i> index )
    {
      static_assert( std::tuple_size< SubOrderRowImp >::value > i,
                     "SubOrderRowImp does not contain the necessary sub structure information.\
                      SubOrderRowImp has to be a std::tuple containing i std::tuple's!" );

      static const SubOrderRowImp order;

      std::cout << "###CREATE: sub container " << print( index )
                << " from global container containing elements " << print( std::get<i>(order) )  << std::endl;
      return cont_( std::get<i>(order), std::get<i>(order) );
    }

    const std::string name() const
    {
      return std::string("global");
    }
  private:
    ContainerType cont_;
  };




  /**
   * \brief Defines a global container.
   *
   * \ingroup Container
   */
  template< class Item2TupleImp, class Item1TupleImp, class SubOrderRowImp, class SubOrderColImp, class... DiscreteFunctions >
  class CombinedGlobalContainer;

  /**
   * \brief Defines a global container.
   *
   * \ingroup Container
   *
   * Before reading this summary, you should have read \ref Container.
   *
   * Global containers has got a method `sub<i>()` which returns the
   * container for the `i`'s sub algorithm.
   *
   * In comparison to the class GlobalContainer this class is the more general one.
   *
   * The template parameters have to have the following structure to work properly:
   * - The first argument has to be a `std::tuple<>` of the same template structure
   *   used for TwoArgContainers.
   * - The first argument has to be a `std::tuple<>` of the same template structure
   *   used for OneArgContainers.
   * - The third argument has to be a `std::tuple<>` and describing the argument
   *   structure of the row arguments.
   * - The forth argument has to be a `std::tuple<>` and describing the argument
   *   structure of the column arguments.
   * - The fifth and further template arguments has to be the arguments for the templates.
   *
   * \warning Each of the first four template classes should contain the same number of elements.
   *
   * The following picture gives an overview, how the class is created and how
   * the template parameters has to be choosen.
   *
   * ![A combined global container](container6.png)
   */
  template< class... Item2TupleImp, class... Item1TupleImp, class... SubOrderRowImp, class... SubOrderColImp, class... DiscreteFunctions >
  class CombinedGlobalContainer< std::tuple< Item2TupleImp... >,
                                 std::tuple< Item1TupleImp... >,
                                 std::tuple< SubOrderRowImp... >,
                                 std::tuple< SubOrderColImp... >,
                                 DiscreteFunctions...>
  {
    //std::tuple< std::tuple< std::tuple< Head2...> >, Item2TupleImp... >

    static const int size1 = std::tuple_size< std::tuple< Item1TupleImp... > >::value;
    static const int size2 = std::tuple_size< std::tuple< Item2TupleImp... > >::value;
    static_assert( size1 == size2, "Item2TupleImp and Item1TupleImp has to contain the same numbers of elements." );

    template< unsigned long int i >
    using ContainerItem = TwoArgContainer< ArgContainerArgWrapper< typename std::tuple_element<i, std::tuple<Item2TupleImp...> >::type >::template _t2Inv,
                                           ArgContainerArgWrapper< typename std::tuple_element<i, std::tuple<Item1TupleImp...> >::type >::template _t1,
                                           typename SubOrderSelect< typename std::tuple_element<i,std::tuple<SubOrderRowImp...> >::type, std::tuple< DiscreteFunctions...> >::type,
                                           typename SubOrderSelect< typename std::tuple_element<i,std::tuple<SubOrderColImp...> >::type, std::tuple< DiscreteFunctions...> >::type >;

    template< unsigned long int i >
    using SharedContainerItem = std::shared_ptr< ContainerItem<i> >;

    typedef typename tuple_copy_t< size1, SharedContainerItem >::type ContainerType;

    static_assert( is_tuple< std::tuple<SubOrderRowImp...> >::value, "SubOrderRowImp has to be a std::tuple<std::tuple<>...>" );
  protected:
    static std::make_integer_sequence< unsigned long int, size1 > sequence;

    template< unsigned long int i, class... Args >
    static decltype(auto) createItem( Args&&... args )
    {
      return std::make_shared<ContainerItem<i> >( args... );
    }

    template< unsigned long int ...i, class SameObject>
    static decltype(auto) createContainer( _indices<i...>, std::shared_ptr< SameObject > obj, const std::string name )
    {
      return std::make_tuple( createItem<i>( obj, name )... );
    }
    template< unsigned long int ...i >
    static decltype(auto) createContainer( _indices<i...>, const std::string name )
    {
      return std::make_tuple( createItem<i>( name )... );
    }

  public:

    /**
     * \brief constructor. Delegates everything
     */
    template< class GridImp, class ... Args>
    explicit CombinedGlobalContainer( const std::shared_ptr< GridImp >& gridptr, Args&&... args )
    : cont_( createContainer( sequence, gridptr, "" /*std::forward<Args>(args)...*/ ) )
    {}

    /**
     * \brief returns the i's container
     */
    template< unsigned long int i >
    decltype(auto) sub( _index<i> index )
    {
      static_assert( std::tuple_size< std::tuple<SubOrderRowImp...> >::value > i,
                     "SubOrderRowImp does not contain the necessary sub structure information.\
                      SubOrderRowImp has to be a std::tuple containing i std::tuple's!" );

      static const std::tuple<SubOrderRowImp...> rowOrder;
      static const std::tuple<SubOrderColImp...> colOrder;

      const auto& cont = std::get<i>( cont_ );

      std::cout << "###CREATE: sub container " << print( index )
                << " from combined global container containing elements " << print( std::get<i>(rowOrder) )<< " x " << print( std::get<i>(colOrder) )  << std::endl;
      return (*cont)( std::get<i>(rowOrder), std::get<i>(colOrder) );
    }

    const std::string name() const
    {
      return std::string("global");
    }
  private:
    ContainerType cont_;
  };

}
}
#endif // FEMHOWTO_STEPPER_HH
