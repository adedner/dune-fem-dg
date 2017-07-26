#ifndef DUNE_FEMDG_CONTAINER_HH
#define DUNE_FEMDG_CONTAINER_HH

#include <memory>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/algorithm/sub/containers.hh>

namespace Dune
{
namespace Fem
{

  template< class... >
  struct SubOrderSelect;

  template< class OrderHead, class... Args >
  struct SubOrderSelect< std::tuple<OrderHead>,std::tuple<Args...> >
  {
    typedef std::tuple< std::tuple_element_t<OrderHead::value,std::tuple<Args...> > > type;
  };

  template< class OrderHead, class... OrderArgs, class... Args >
  struct SubOrderSelect< std::tuple<OrderHead, OrderArgs...>,std::tuple<Args...> >
  {
    typedef tuple_concat_t< std::tuple< std::tuple_element_t<OrderHead::value,std::tuple<Args...> > >,
                            typename SubOrderSelect< std::tuple<OrderArgs...>, std::tuple<Args...> >::type > type;
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
  template< class Item2TupleImp, class Item1TupleImp, class SubOrderRowImp, class SubOrderColImp, class ExtraArg, class... DiscreteFunctions >
  class GlobalContainer
  {
    typedef TwoArgContainer< ArgContainerArgWrapper< Item2TupleImp >::template _t2Inv,
                             ArgContainerArgWrapper< Item1TupleImp >::template _t1,
                             ArgContainerArgWrapper< Item1TupleImp >::template _t1,
                             std::tuple< DiscreteFunctions... >,
                             std::tuple< DiscreteFunctions... > >                      ContainerType;

    static_assert( is_tuple< SubOrderRowImp >::value, "SubOrderRowImp has to be a std::tuple<std::tuple<>...>" );
  public:

    /**
     * \brief constructor. Delegates everything
     */
    template< class... SameObjects >
    explicit GlobalContainer( const std::string name, const std::shared_ptr< SameObjects >& ... obj )
    : cont_( std::get<0>(std::tie( obj... ) ), name )
    {}

    /**
     * \brief returns the i's container
     */
    template< unsigned long int i >
    decltype(auto) sub( _index<i> index ) const
    {
      static_assert( std::tuple_size< SubOrderRowImp >::value > i,
                     "SubOrderRowImp does not contain the necessary sub structure information.\
                      SubOrderRowImp has to be a std::tuple containing i std::tuple's!" );

      static const SubOrderRowImp rowOrder;
      static const SubOrderColImp colOrder;

      std::cout << "###CREATE: sub container " << print( index )
                << " from global container containing elements " << print( std::get<i>(rowOrder) )  << std::endl;
      return cont_( std::get<i>(rowOrder), std::get<i>(colOrder) );
    }

    template< unsigned long int i >
    decltype(auto) extra( _index<i> index ) const
    {
      return ExtraArg::init( *this );
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
  template< class Item2TupleImp, class Item1TupleImp, class SubOrderRowImp, class SubOrderColImp, class ExtraArg, class... DiscreteFunctions >
  class CombinedGlobalContainer
  {
    static_assert( static_fail< Item2TupleImp >::value, "Please check your template arguments" );
  };

  /**
   * \brief Defines a global container. Partial specialization
   *
   * \ingroup Container
   *
   * \copydoc CombinedGlobalContainer
   */
  template< class... Item2TupleImp, class... Item1TupleImp, class... SubOrderRowImp, class... SubOrderColImp, class... ExtraArgs, class... DiscreteFunctions >
  class CombinedGlobalContainer< std::tuple< Item2TupleImp... >,
                                 std::tuple< Item1TupleImp... >,
                                 std::tuple< SubOrderRowImp... >,
                                 std::tuple< SubOrderColImp... >,
                                 std::tuple< ExtraArgs... >,
                                 DiscreteFunctions...>
  {
    static const int size1 = std::tuple_size< std::tuple< Item1TupleImp... > >::value;
    static const int size2 = std::tuple_size< std::tuple< Item2TupleImp... > >::value;
    static_assert( size1 == size2, "Item2TupleImp and Item1TupleImp has to contain the same numbers of elements." );

    template< unsigned long int i >
    using ContainerItem = TwoArgContainer< ArgContainerArgWrapper< std::tuple_element_t<i, std::tuple<Item2TupleImp...> > >::template _t2Inv,
                                           ArgContainerArgWrapper< std::tuple_element_t<i, std::tuple<Item1TupleImp...> > >::template _t1,
                                           ArgContainerArgWrapper< std::tuple_element_t<i, std::tuple<Item1TupleImp...> > >::template _t1,
                                           typename SubOrderSelect< std::tuple_element_t<i,std::tuple<SubOrderRowImp...> >, std::tuple< DiscreteFunctions...> >::type,
                                           typename SubOrderSelect< std::tuple_element_t<i,std::tuple<SubOrderColImp...> >, std::tuple< DiscreteFunctions...> >::type >;

    template< unsigned long int i >
    using SharedContainerItem = std::shared_ptr< ContainerItem<i> >;

    typedef typename tuple_copy_t< size1, SharedContainerItem >::type ContainerType;

    static_assert( is_tuple< std::tuple<SubOrderRowImp...> >::value, "SubOrderRowImp has to be a std::tuple<std::tuple<>...>" );
  protected:
    static_assert( size1 >= 0, "invalid integer_sequence: Throw this additional assertion here because \
                                gcc-6 won't stop compiling for a very long time... :(" );

    static std::make_integer_sequence< unsigned long int, size1 > sequence;

    template< unsigned long int i, class SameObject >
    static decltype(auto) createItem( const std::string name, const std::shared_ptr< SameObject>& obj )
    {
      return std::make_shared<ContainerItem<i> >( obj, name );
    }

    template< unsigned long int ...i, class... SameObject>
    static decltype(auto) createContainer( _indices<i...>, const std::string name, std::tuple< std::shared_ptr< SameObject >... > obj)
    {
      return std::make_tuple( createItem<i>( name, std::get<i>(obj) )... );
    }
  public:

    /**
     * \brief constructor. Delegates everything
     */
    template< class... SameObjects >
    explicit CombinedGlobalContainer( const std::string name, const std::shared_ptr< SameObjects >& ... obj )
    : cont_( createContainer( sequence, name, std::make_tuple( obj... ) ) )
    {}

    /**
     * \brief returns the i's container
     */
    template< unsigned long int i >
    decltype(auto) sub( _index<i> index ) const
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

    template< unsigned long int i >
    decltype(auto) extra( _index<i> index ) const
    {
      typedef std::tuple_element_t<i,std::tuple< ExtraArgs... > > ExtraArgType;
      return ExtraArgType::init( *this );
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
   *
   */
  template< class Container, class SubOrderRowImp, class SubOrderColImp, class ExtraArgs, class... DiscreteFunctions >
  class NewCombinedGlobalContainer
  {
    static_assert( static_fail< Container >::value, "Please check your template arguments" );
  };

  /**
   * \brief Defines a global container. Partial specialization
   *
   * \ingroup Container
   *
   * \copydoc NewCombinedGlobalContainer
   */
  template< class... Containers, class... SubOrderRowImp, class... SubOrderColImp, class... ExtraArgs, class... DiscreteFunctions >
  class NewCombinedGlobalContainer< std::tuple< Containers... >,
                                    std::tuple< SubOrderRowImp... >,
                                    std::tuple< SubOrderColImp... >,
                                    std::tuple< ExtraArgs... >,
                                    DiscreteFunctions...>
  {
    typedef std::tuple< typename Containers::RowOneArgType... > RowItem1TupleImp;
    typedef std::tuple< typename Containers::ColOneArgType... > ColItem1TupleImp;
    typedef std::tuple< typename Containers::TwoArgType... >    Item2TupleImp;


    static const int size = std::tuple_size< RowItem1TupleImp >::value;

    template< unsigned long int i >
    using Item2TupleItem = std::tuple_element_t<i,Item2TupleImp>;

    template< unsigned long int i >
    using RowItem1TupleItem = std::tuple_element_t<i,RowItem1TupleImp>;

    template< unsigned long int i >
    using ColItem1TupleItem = std::tuple_element_t<i,ColItem1TupleImp>;


    template< unsigned long int i >
    using ContainerItem = TwoArgContainer< Item2TupleItem<i>::template _t2,
                                           RowItem1TupleItem<i>::template _t1,
                                           ColItem1TupleItem<i>::template _t1,
                                           typename SubOrderSelect< std::tuple_element_t<i,std::tuple<SubOrderRowImp...> >, std::tuple< DiscreteFunctions...> >::type,
                                           typename SubOrderSelect< std::tuple_element_t<i,std::tuple<SubOrderColImp...> >, std::tuple< DiscreteFunctions...> >::type >;

    template< unsigned long int i >
    using SharedContainerItem = std::shared_ptr< ContainerItem<i> >;

    typedef typename tuple_copy_t< size, SharedContainerItem >::type ContainerType;

    static_assert( is_tuple< std::tuple<SubOrderRowImp...> >::value, "SubOrderRowImp has to be a std::tuple<std::tuple<>...>" );
  protected:
    static_assert( size >= 0, "invalid integer_sequence: Throw this additional assertion here because \
                               gcc-6 won't stop compiling for a very long time... :(" );
    static std::make_integer_sequence< unsigned long int, size > sequence;

    template< unsigned long int i, class SameObject >
    static decltype(auto) createItem( const std::string name, const std::shared_ptr< SameObject>& obj )
    {
      return std::make_shared<ContainerItem<i> >( obj, name );
    }

    template< unsigned long int ...i, class... SameObject>
    static decltype(auto) createContainer( _indices<i...>, const std::string name, std::tuple< std::shared_ptr< SameObject >... > obj)
    {
      return std::make_tuple( createItem<i>( name, std::get<i>(obj) )... );
    }
  public:

    /**
     * \brief constructor. Delegates everything
     */
    template< class... SameObjects >
    explicit NewCombinedGlobalContainer( const std::string name, const std::shared_ptr< SameObjects >& ... obj )
    : cont_( createContainer( sequence, name, std::make_tuple( obj... ) ) )
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

    template< unsigned long int i >
    decltype(auto) extra( _index<i> index )
    {
      typedef std::tuple_element_t<i,std::tuple< ExtraArgs... > > ExtraArgType;
      return ExtraArgType::init( *this );
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
