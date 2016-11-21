#ifndef DUNE_FEMDG_CONTAINER_HH
#define DUNE_FEMDG_CONTAINER_HH

#include <memory>
#include <dune/common/shared_ptr.hh>


namespace Dune
{
namespace Fem
{


  template< class DiscreteFunctionImp >
  struct ContainerItem
  {

    typedef DiscreteFunctionImp                                      DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType         GridPartType;
    typedef typename GridPartType::GridType                          GridType;

    //ContainerItem( const std::string& name )
    //: name_( name ),
    //  gridPart_( nullptr ),
    //  space_( nullptr ),
    //  df_( nullptr )
    //{}

    ContainerItem( const std::string& name, const std::shared_ptr<GridType >& grid )
    : name_( name ),
      grid_( grid ),
      gridPart_( std::make_shared< GridPartType >( *grid_ ) ),
      space_( std::make_shared< DiscreteFunctionSpaceType >( *gridPart_ ) ),
      df_( std::make_shared< DiscreteFunctionType >( name, *space_ ) )
    {}

    //be carefull
    ContainerItem( const std::string& name, const std::shared_ptr<GridPartType>& gridPart )
    : name_( name ),
      grid_( nullptr ),
      gridPart_( gridPart ),
      space_( std::make_shared< DiscreteFunctionSpaceType >( *gridPart_ ) ),
      df_( std::make_shared< DiscreteFunctionType >( name, *space_ ) )
    {}

    //be carefull
    ContainerItem( const std::string& name, const std::shared_ptr<DiscreteFunctionSpaceType>& space )
    : name_( name ),
      grid_( nullptr ),
      gridPart_( nullptr ),
      space_( space ),
      df_( std::make_shared< DiscreteFunctionType >( name, space ) )
    {}

    ContainerItem( const std::string& name, const ContainerItem& cont )
    : name_( name ),
      grid_( cont.grid_ ),
      gridPart_( cont.gridPart_ ),
      space_( cont.space_ ),
      df_( std::make_shared< DiscreteFunctionType >( name, *space_ ) )
    {
    }

    std::shared_ptr< DiscreteFunctionType > shared() const
    {
      return df_;
    }

  protected:
    const std::string                            name_;
    std::shared_ptr< GridType >                  grid_;
    std::shared_ptr< GridPartType >              gridPart_;
    std::shared_ptr< DiscreteFunctionSpaceType > space_;
    std::shared_ptr< DiscreteFunctionType >      df_;
  };

  //forward declaration
  template< template<class,int...> class, class Arg>
  struct OneArgContainer
  {
    static_assert( static_fail< Arg >::value, "second argument has to be a tuple" );
  };

  //forward declaration
  template< template<class,class,int...> class, template<class,int...> class, class Arg, class >
  class TwoArgContainer
  {
    static_assert( static_fail< Arg >::value, "third/fourth argument has to be a tuple" );
  };


  //only for tuples
  template< template< class, int... > class OneArgImp, class... Args >
  struct OneArgContainer< OneArgImp, std::tuple< Args... > >
  {
    typedef std::tuple< Args... >                ArgTupleType;
    typedef typename VectorPacker< OneArgImp, ArgTupleType >::shared_type Item1TupleType;

  public:
    template< unsigned long int i >
    using Item1 = typename std::tuple_element< i, Item1TupleType>::type::element_type;

    static const int numArgs = std::tuple_size< ArgTupleType >::value;

  protected:
    static const int size = std::tuple_size< Item1TupleType >::value;
    static std::make_integer_sequence< unsigned long int, size > sequence;

    ////// Creation
    template< unsigned long int i, class SameObject >
    static decltype(auto) createItem1( std::shared_ptr< SameObject > obj, const std::string name )
    {
      std::cout << "###CREATE: item1 ('" << name << "'): " << print( _index<i>() ) << std::endl;
      return std::make_shared<Item1<i> >( obj, name );
    }
    template< unsigned long int i >
    static decltype(auto) createItem1( const std::string name )
    {
      std::cout << "###CREATE: item1 ('" << name << "'): " << print( _index<i>() ) << std::endl;
      return std::make_shared< Item1<i> >( name );
    }
    template< unsigned long int ...i, class SameObject>
    static decltype(auto) createContainer( _indices<i...>, std::shared_ptr< SameObject > obj, const std::string name )
    {
      return std::make_tuple( createItem1<i>( obj, name )... );
    }
    template< unsigned long int ...i >
    static decltype(auto) createContainer( _indices<i...>, const std::string name )
    {
      return std::make_tuple( createItem1<i>( name )... );
    }

    ////// Copy
    template< class Item1Tuple, unsigned long int ...i >
    static decltype(auto) copyContainer( const Item1Tuple& item1, std::tuple< _index<i>...  > )
    {
      std::cout << "###COPY: item1 ('-'): " << print( _index<i>()... ) << std::endl;
      return std::make_tuple( std::get<i>( item1 )... );
    }

    //Be your own (template) friend
    template< template<class,int...> class, class >
    friend class OneArgContainer;

    template< template<class,class,int...> class, template<class,int...> class, class, class >
    friend class TwoArgContainer;

  public:
    // copy, for internal use only
    OneArgContainer( const Item1TupleType& item )
    : item1_( item )
    {}

    // owning container
    template< class SameObject >
    OneArgContainer( std::shared_ptr< SameObject > obj, const std::string name = "" )
    : item1_( createContainer( sequence, obj, name ) )
    {}

    // non owning container, for coupling
    OneArgContainer( const std::string name = "" )
    : item1_( createContainer( sequence, name ) )
    {}

    // item access
    template< unsigned long int i >
    decltype(auto) operator() ( _index<i> index )
    {
      const auto& res = std::get<i>( item1_ );
      //std::cout << "###ACCESS: item1 ('" << res->solution()->name() << "')" << print( index ) << std::endl;
      return std::get<i>( item1_ );
    }

    // sub container
    template< unsigned long int... i >
    decltype(auto) operator() ( std::tuple< _index<i>... > index )
    {
      typedef std::tuple< typename std::tuple_element<i,ArgTupleType>::type...> NewArgTupleType;
      typedef OneArgContainer< OneArgImp, NewArgTupleType > SubContainerType;

      return std::make_shared< SubContainerType >( copyContainer(item1_, index ) );
    }
  protected:
    Item1TupleType item1_;
  };



  template< template<class,class,int...> class TwoArgImp, template<class,int...> class OneArgImp, class... RowArgs, class... ColArgs >
  class TwoArgContainer< TwoArgImp, OneArgImp, std::tuple< RowArgs... >, std::tuple< ColArgs... > >
    : public OneArgContainer< OneArgImp, std::tuple< RowArgs... > >
  {
    typedef OneArgContainer< OneArgImp, std::tuple< RowArgs... > > BaseType;
    typedef OneArgContainer< OneArgImp, std::tuple< ColArgs... > > FakeColBaseType;

    typedef std::tuple< RowArgs... >        RowArgTupleType;
    typedef std::tuple< ColArgs... >        ColArgTupleType;

    typedef typename MatrixPacker< TwoArgImp, RowArgTupleType, ColArgTupleType >::shared_type  Item2TupleType;
    typedef typename VectorPacker< OneArgImp, RowArgTupleType >::shared_type                   Item1TupleType;

    typedef Item1TupleType                                                                     RowItem1TupleType;
    typedef typename VectorPacker< OneArgImp, ColArgTupleType >::shared_type                   ColItem1TupleType;
  public:

    template< unsigned long int i, unsigned long int j >
    using Item2 = typename std::tuple_element< j, typename std::tuple_element< i, Item2TupleType>::type >::type::element_type;

    using BaseType::operator();

    static const int numRowArgs = std::tuple_size< RowArgTupleType >::value;
    static const int numColArgs = std::tuple_size< ColArgTupleType >::value;

  protected:
    using BaseType::item1_;

    static const int size = BaseType::size;
    using BaseType::sequence;

    ///// Creation
    template< unsigned long int i, unsigned long int j >
    decltype(auto) createItem2( const std::string name )
    {
      std::cout << "###CREATE: item2 ('" << name << "'): " << print( _index<i>(), _index<j>() ) << std::endl;
      return std::make_shared<Item2<i,j> >( BaseType::operator()( _index<i>() ),
                                            BaseType::operator()( _index<j>() ),
                                            name );
    }
    template< unsigned long int i, unsigned long int ...j >
    decltype(auto) createContainerRow( _indices<j...>, const std::string name )
    {
      return std::make_tuple( createItem2<i,j>( name )... );
    }
    template< unsigned long int ...i, unsigned long int ...j >
    decltype(auto) createContainer( _indices<i...> row, _indices<j...> col, const std::string name )
    {
      return std::make_tuple( createContainerRow<i>( col, name )... );
    }




   // template< int i, unsigned long int... i1, unsigned long int... i2 >
   // constexpr unsigned long int entry1(std::tuple< std::pair< _index<i1>, _index<i2> >... >)
   // {
   //   typedef std::tuple< std::pair< _index<i1>, _index<i2> >... > TuplePair;
   //   return typename std::tuple_element<i, typename TuplePair::first_type >::type::value;
   // }

   // template< int i, unsigned long int... i1, unsigned long int... i2 >
   // constexpr unsigned long int entry2(std::tuple< std::pair< _index<i1>, _index<i2> >... >)
   // {
   //   typedef std::tuple< std::pair< _index<i1>, _index<i2> >... > TuplePair;
   //   return typename std::tuple_element<i, typename TuplePair::second_type >::type::value;
   // }




   // template< unsigned long int i, unsigned long int j, unsigned long int... j1, unsigned long int... j2 >
   // decltype(auto) copyContainerRowSubSub( std::tuple< std::pair< _index<j1>, _index<j2> >... > tupPair)
   // {
   //   typedef std::tuple< std::pair< _index<j1>, _index<j2> >... > TuplePair;

   //   static unsigned long int jj = std::get<j>( tupPair );

   //   return std::make_tuple( std::get<j>( std::get<i>( item2_ ) )... );
   // }

   // template< unsigned long int i, unsigned long int... j, unsigned long int... j1, unsigned long int... j2 >
   // decltype(auto) copyContainerRowSub( _indices<j...>)
   // {
   //   return std::make_tuple( std::get<j>( std::get<i>( item2_ ) )... );
   // }

   // template< unsigned long int ii, unsigned long int ii2, unsigned long int... j1, unsigned long int... j2 >
   // decltype(auto) copyContainerRow( std::tuple< std::pair< _index<j1>, _index<j2> >... > jPair )
   // {
   //   typedef std::tuple< std::pair< _index<j1>, _index<j2> >... > TupleImp;
   // static const int tupIndSize = std::tuple_size< TupleImp >::value;
   //   static std::make_integer_sequence< unsigned long int, tupIndSize > seq;
   //   return copyContainerRowSub<ii>( seq );
   // }

   // //current work
   // template< unsigned long int i, unsigned long int j, unsigned long int... i1, unsigned long int... i2, unsigned long int... j1, unsigned long int... j2 >
   // decltype(auto) copyContainerRow( std::tuple< std::pair< _index<i1>, _index<i2> >... > iTup,
   //                                  std::tuple< std::pair< _index<j1>, _index<j2> >... > jTup )
   // {
   //   typedef std::tuple< RowArgs... >        RowArgTupleImp;
   //   typedef std::tuple< ColArgs... >        ColArgTupleImp;

   //   typedef typename MatrixPacker< TwoArgImp, RowArgTupleImp, ColArgTupleImp >::type  RowItem2TupleImp;
   //   typedef typename MatrixPacker< TwoArgImp, RowArgTupleImp, ColArgTupleImp >::type  ColItem2TupleImp;

   //   template< unsigned long int i, unsigned long int j >
   //   using RowItem2 = typename std::tuple_element< j, typename std::tuple_element< i, RowItem2TupleImp>::type >::type:

   //   template< unsigned long int i, unsigned long int j >
   //   using ColItem2 = typename std::tuple_element< j, typename std::tuple_element< i, ColItem2TupleImp>::type >::type:

   //   typedef std::tuple< std::pair< _index<i1>, _index<i2> >... > ITuplePair;
   //   typedef std::tuple< std::pair< _index<j1>, _index<j2> >... > JTuplePair;
   //   //number of subcontainer
   //   static unsigned long int iSub = typename std::tuple_element<i, typename ITuplePair::first_type >::type::value;
   //   static unsigned long int jSub = typename std::tuple_element<j, typename JTuplePair::first_type >::type::value;
   //   //
   //   static unsigned long int ii = typename std::tuple_element<i, typename ITuplePair::second_type >::type::value;
   //   static unsigned long int jj = typename std::tuple_element<j, typename JTuplePair::second_type >::type::value;

   //   //Item2????
   //   return std::make_shared<Item2<iSub,jSub> >( std::get<iSub>( item2 )( _index<ii>() ),
   //                                               std::get<jSub>( item2 )( _index<jj>() ),
   //                                               "dummy" );
   // }



   // template< unsigned long int... i1, unsigned long int... i2, unsigned long int... j1, unsigned long int... j2 >
   // decltype(auto) copyContainer(std::tuple< std::pair< _index<i1>, _index<i2> >... > iTup,
   //                              std::tuple< std::pair< _index<j1>, _index<j2> >... > jTup)
   // {
   //   typedef std::tuple< std::pair< _index<i1>, _index<i2> >... > IPairType;
   //   typedef std::tuple< std::pair< _index<j1>, _index<j2> >... > JPairType;
   //   //number sub container
   //   static unsigned long int ii = entry1<i>( iTup );

   //   static const int iPairSize = std::tuple_size< IPairType >::value;
   //   static std::make_integer_sequence< unsigned long int, iPairSize > seq;

   //   //number of sub element in subcontainer
   //   static unsigned long int ii2 = entry2<i>( iTup );
   //   return std::make_tuple( copyContainerRow<ii,ii2>( col )... );
   // }


    ///// Copy
    template< unsigned long int i, unsigned long int ...j >
    decltype(auto) copyContainerRow( std::tuple< _index<j>... > )
    {
      std::cout << "###COPY: item2 ('-'): " << print( std::make_tuple( _index<j>()... ) ) << std::endl;
      return std::make_tuple( std::get<j>( std::get<i>( item2_ ) )... );
    }
    template< unsigned long int ...i, unsigned long int ...j >
    decltype(auto) copyContainer( std::tuple< _index<i>... > row, std::tuple< _index<j>... > col )
    {
      return std::make_tuple( copyContainerRow<i>( col )... );
    }

    //Be your own (template) friend
    template< template<class,class,int...> class, template<class,int...> class, class, class >
    friend class TwoArgContainer;

    template< template<class,int...> class, class >
    friend class OneArgContainer;

    //todo: protected would be much better...
    //protected:
  public:
    // copy, for internal use only
    TwoArgContainer( const RowItem1TupleType& rowItem1,
                     const ColItem1TupleType& colItem1,
                     const Item2TupleType& item2 )
    : BaseType( rowItem1 ),
      item2_( item2 ),
      colItem1_( colItem1 )
    {}
  public:
    template< class SameObject >
    TwoArgContainer( std::shared_ptr< SameObject > obj, const std::string name = "" )
    : BaseType( obj, name ),
      item2_( createContainer( sequence, sequence, name + "matrix" ) ),
      colItem1_( FakeColBaseType::createContainer( sequence, obj, name ) )
    {}

    TwoArgContainer( const std::string name = "" )
    : BaseType( name ),
      item2_( createContainer( sequence, sequence, name + "matrix" ) ),
      colItem1_( FakeColBaseType::createContainer( sequence, name ) )
    {}


    // item acess
    template< unsigned long int i, unsigned long int j >
    decltype(auto) operator() ( _index<i> row, _index<j> col  )
    {
      const auto& res = std::get<j>( std::get<i>( item2_ ) );
      //std::cout << "###ACCESS: item2 ('-')" << print( row, col ) << std::endl;
      return res;
    }

    // sub Container
    template< unsigned long int... i, unsigned long int... j >
    decltype(auto) operator() ( std::tuple< _index<i>... > row, std::tuple< _index<j>... > col )
    {
      typedef std::tuple< typename std::tuple_element<i,RowArgTupleType>::type...> NewRowArgTupleType;
      typedef std::tuple< typename std::tuple_element<j,ColArgTupleType>::type...> NewColArgTupleType;

      typedef OneArgContainer< OneArgImp, NewRowArgTupleType > RowOneArgContainerType;
      typedef OneArgContainer< OneArgImp, NewColArgTupleType > ColOneArgContainerType;

      typedef TwoArgContainer< TwoArgImp, OneArgImp, NewRowArgTupleType, NewColArgTupleType > SubContainerType;

      std::cout << "###CREATE: local sub container " << print( row, col ) << std::endl;
      return std::make_shared< SubContainerType >( RowOneArgContainerType::copyContainer( item1_, row ),
                                                   ColOneArgContainerType::copyContainer( colItem1_, col ),
                                                   copyContainer( row, col ) );
    }
  protected:
    Item2TupleType          item2_;
    ColItem1TupleType       colItem1_;
  };


  //empty item1 oder item2 container
  template< class... MatrixContainerImp >
  struct EmptyContainerItem
  {
    template< class ... Args>
    EmptyContainerItem( Args&&... args )
    {}
  };

  namespace details
  {
    template< unsigned long int, class... >
    struct int_const_tuple;

    template< unsigned long int id, class Arg >
    struct int_const_tuple< id, Arg >
    {
      typedef std::tuple< _index<id> > type;
    };

    template< unsigned long int id, class Arg, class... Args >
    struct int_const_tuple< id, Arg, Args... >
    {
      typedef typename tuple_concat< std::tuple< _index<id> >, typename int_const_tuple<id+1,Args...>::type >::type type;
    };
  }

  //generate tuple containing a consecutive number on integral_constant<>
  template< class Arg, class... Args >
  struct index_tuple
  {
    typedef typename details::int_const_tuple< 0, Arg, Args... >::type type;
  };


  template< class... TupleMatrices >
  struct tuple_matrix_combiner;


  template< class TupleMatrix >
  struct tuple_matrix_combiner< TupleMatrix >
  {
    typedef TupleMatrix type;
  };

  template< class TupleMatrix1, class TupleMatrix2 >
  struct tuple_matrix_combiner< TupleMatrix1, TupleMatrix2 >
  {
    typedef typename tuple_matrix_combine< TupleMatrix1, TupleMatrix2, _t<EmptyContainerItem> >::type type;
  };

  template< class TupleMatrix1, class TupleMatrix2, class TupleMatrix3, class... TupleMatrixArgs >
  struct tuple_matrix_combiner< TupleMatrix1, TupleMatrix2, TupleMatrix3, TupleMatrixArgs... >
  {
    typedef typename tuple_matrix_combiner< typename tuple_matrix_combine< TupleMatrix1, TupleMatrix2, _t<EmptyContainerItem> >::type, TupleMatrix3, TupleMatrixArgs... >::type type;
  };




  //TODO this is not a sub-container, it is a global container...
  //Put this class to another header...
  template< class Item2TupleImp, class Item1TupleImp, class SubOrderRowImp, class SubOrderColImp, class... DiscreteFunctions >
  class GlobalContainer
  {
    typedef TwoArgContainer< _template< Item2TupleImp >::template _t2Inv, _template< Item1TupleImp >::template _t1,
                             std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > > ContainerType;

    static_assert( is_tuple< SubOrderRowImp >::value, "SubOrderRowImp has to be a std::tuple<std::tuple<>...>" );
  public:

    // constructor: do not touch, delegate everything
    template< class GridImp, class ... Args>
    GlobalContainer( const std::shared_ptr< GridImp >& gridptr, Args&&... args )
    : cont_( gridptr, std::forward<Args>(args)... )
    {}

    // sub container
    template< unsigned long int i >
    decltype(auto) sub( _index<i> index )
    {
      static_assert( std::tuple_size< SubOrderRowImp >::value > i,
                     "SubOrderRowImp does not contain the necessary sub structure information.\
                      SubOrderRowImp has to be a std::tuple containing i std::tuple's!" );
      //default implementation:
      //static const typename index_tuple< DiscreteFunctions... >::type order;

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


  //Forward declaration
  template< class Item2TupleImp, class Item1TupleImp, class SubOrderRowImp, class SubOrderColImp, class... DiscreteFunctions >
  class CombinedGlobalContainer;

  template< class... Item2TupleImp, class... Item1TupleImp, class... SubOrderRowImp, class... SubOrderColImp, class... DiscreteFunctions >
  class CombinedGlobalContainer< std::tuple< Item2TupleImp... >,
                                 std::tuple< Item1TupleImp... >,
                                 std::tuple< SubOrderRowImp... >,
                                 std::tuple< SubOrderColImp... >,
                                 DiscreteFunctions...>
  {
    static const int size1 = std::tuple_size< std::tuple< Item1TupleImp... > >::value;
    static const int size2 = std::tuple_size< std::tuple< Item2TupleImp... > >::value;
    static_assert( size1 == size2, "Item2TupleImp and Item1TupleImp has to contain the same numbers of elements." );

    template< unsigned long int i >
    using ContainerItem = TwoArgContainer< _template< typename std::tuple_element<i, std::tuple<Item2TupleImp...> >::type >::template _t2Inv,
                                           _template< typename std::tuple_element<i, std::tuple<Item1TupleImp...> >::type >::template _t1,
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

    // constructor: do not touch, delegate everything
    template< class GridImp, class ... Args>
    explicit CombinedGlobalContainer( const std::shared_ptr< GridImp >& gridptr, Args&&... args )
    : cont_( createContainer( sequence, gridptr, "" /*std::forward<Args>(args)...*/ ) )
    {}

    // sub container
    template< unsigned long int i >
    decltype(auto) sub( _index<i> index )
    {
      static_assert( std::tuple_size< std::tuple<SubOrderRowImp...> >::value > i,
                     "SubOrderRowImp does not contain the necessary sub structure information.\
                      SubOrderRowImp has to be a std::tuple containing i std::tuple's!" );
      //default implementation:
      //static const typename index_tuple< DiscreteFunctions... >::type order;

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
