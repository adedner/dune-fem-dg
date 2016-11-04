#ifndef DUNE_FEMDG_CONTAINER_HH
#define DUNE_FEMDG_CONTAINER_HH

#include <memory>



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

    ContainerItem( const std::string& name )
    : name_( name ),
      gridPart_( nullptr ),
      space_( nullptr ),
      df_( nullptr )
    {}

    ContainerItem( const std::string& name, GridType& grid)
    : name_( name ),
      gridPart_( std::make_unique< GridPartType >( grid ) ),
      space_( std::make_unique< DiscreteFunctionSpaceType >( *gridPart_ ) ),
      df_( std::make_shared< DiscreteFunctionType >( name, *space_ ) )
    {}

    ContainerItem( const std::string& name, const GridPartType& gridPart )
    : name_( name ),
      gridPart_( nullptr ),
      space_( std::make_unique< DiscreteFunctionSpaceType >( *gridPart_ ) ),
      df_( std::make_shared< DiscreteFunctionType >( name, *space_ ) )
    {}

    ContainerItem( const std::string& name, const DiscreteFunctionSpaceType& space )
    : name_( name ),
      gridPart_( nullptr ),
      space_( nullptr ),
      df_( std::make_shared< DiscreteFunctionType >( name, space ) )
    {}

    ContainerItem( const std::string& name, const ContainerItem& cont )
    : name_( name ),
      gridPart_( nullptr ),
      space_( nullptr ),
      df_( std::make_shared< DiscreteFunctionType >( name, cont.space() ) )
    {}

    const GridType& grid() const
    {
      assert( df_ );
      return df_->space().gridPart().grid();
    }
    GridType& grid()
    {
      assert( df_ );
      return df_->space().gridPart().grid();
    }

    const GridPartType& gridPart() const
    {
      assert( df_ );
      return df_->space().gridPart();
    }
    GridPartType& gridPart()
    {
      assert( df_ );
      return df_->space().gridPart();
    }

    const DiscreteFunctionSpaceType& space() const
    {
      assert( df_ );
      return df_->space();
    }
    DiscreteFunctionSpaceType& space()
    {
      assert( df_ );
      return const_cast< DiscreteFunctionSpaceType& >( df_->space() );
    }

    const DiscreteFunctionType& operator()()
    {
      assert( df_ );
      return *df_;
    }

    std::shared_ptr< DiscreteFunctionType > shared()
    {
      return df_;
    }


    template< class WrongArg >
    ContainerItem& operator=( WrongArg& arg )
    {
      static_assert( std::is_same< WrongArg, void >::value, "types does not match" );
      return *this;
    }

    ContainerItem& operator=( const GridType& grid)
    {
      gridPart_ = std::make_unique< GridPartType >( grid );
      space_ = std::make_unique< DiscreteFunctionSpaceType >( *gridPart_ );
      df_ = std::make_shared< DiscreteFunctionType >( name_, *space_ );
      return *this;
    }

    ContainerItem& operator=( GridType& grid)
    {
      gridPart_ = std::make_unique< GridPartType >( grid );
      space_ = std::make_unique< DiscreteFunctionSpaceType >( *gridPart_ );
      df_ = std::make_shared< DiscreteFunctionType >( name_, *space_ );
      return *this;
    }

    ContainerItem& operator=( const GridPartType& gridPart )
    {
      space_ = std::make_unique< DiscreteFunctionSpaceType >( gridPart );
      df_ = std::make_shared< DiscreteFunctionType >( name_, *space_ );
      return *this;
    }
    ContainerItem& operator=( GridPartType& gridPart )
    {
      space_ = std::make_unique< DiscreteFunctionSpaceType >( gridPart );
      df_ = std::make_shared< DiscreteFunctionType >( name_, *space_ );
      return *this;
    }

    ContainerItem& operator=( const DiscreteFunctionSpaceType& space )
    {
      df_ = std::make_shared< DiscreteFunctionType >( name_, space );
      return *this;
    }
    ContainerItem& operator=( DiscreteFunctionSpaceType& space )
    {
      df_ = std::make_shared< DiscreteFunctionType >( name_, space );
      return *this;
    }

    ContainerItem& operator=( DiscreteFunctionType& df )
    {
      df_ = std::make_shared< DiscreteFunctionType >( name_, df.space() );
      return *this;
    }
    ContainerItem& operator=( const DiscreteFunctionType& df )
    {
      df_ = std::make_shared< DiscreteFunctionType >( name_, df.space() );
      return *this;
    }

    ContainerItem& operator=( ContainerItem& cont )
    {
      df_ = cont.shared();
      return *this;
    }
    ContainerItem& operator=( const ContainerItem& cont )
    {
      df_ = cont.shared();
      return *this;
    }

  private:
    const std::string                            name_;
    std::unique_ptr< GridPartType >              gridPart_;
    std::unique_ptr< DiscreteFunctionSpaceType > space_;
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

  protected:
    static const int size = std::tuple_size< Item1TupleType >::value;
    static std::make_integer_sequence< unsigned long int, size > sequence;

    ////// Creation
    template< unsigned long int i, class SameObject >
    static std::shared_ptr< Item1<i> > createItem1( SameObject& obj, const std::string name )
    {
      return std::make_shared<Item1<i> >( obj, name );
    }
    template< unsigned long int i >
    static std::shared_ptr< Item1<i> > createItem1( const std::string name )
    {
      return std::make_shared< Item1<i> >( name );
    }
    template< unsigned long int ...i, class SameObject>
    static Item1TupleType createContainer( _indices<i...>, SameObject& obj, const std::string name )
    {
      return std::make_tuple( createItem1<i>( obj, name )... );
    }
    template< unsigned long int ...i >
    static Item1TupleType createContainer( _indices<i...>, const std::string name )
    {
      return std::make_tuple( createItem1<i>( name )... );
    }

    ////// Copy
    template< class Item1Tuple, unsigned long int ...i >
    static Item1TupleType copyContainer( const Item1Tuple& item1, std::tuple< _index<i>...  > )
    {
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
    OneArgContainer( SameObject& obj, const std::string name = "" )
    : item1_( createContainer( sequence, obj, name ) )
    {}

    // non owning container, for coupling
    OneArgContainer( const std::string name = "" )
    : item1_( createContainer( sequence, name ) )
    {}

    // item access
    template< unsigned long int i >
    std::shared_ptr< Item1<i> > operator() ( _index<i> index )
    {
      return std::get<i>( item1_ );
    }

    // sub container
    template< unsigned long int... i >
    std::shared_ptr< OneArgContainer< OneArgImp, std::tuple< typename std::tuple_element<i,ArgTupleType>::type...> > >
    operator() ( std::tuple< _index<i>... > index )
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

  protected:
    using BaseType::item1_;

    static const int size = BaseType::size;
    using BaseType::sequence;

    ///// Creation
    template< unsigned long int i, unsigned long int j >
    std::shared_ptr< Item2<i,j> > createItem2( const std::string name )
    {
      return std::make_shared<Item2<i,j> >( BaseType::operator()( _index<i>() ),
                                            BaseType::operator()( _index<j>() ),
                                            name );
    }
    template< unsigned long int i, unsigned long int ...j >
    std::tuple< std::shared_ptr< Item2<i,j> > ... >
    createContainerRow( _indices<j...>, const std::string name )
    {
      return std::make_tuple( createItem2<i,j>( name )... );
    }
    template< unsigned long int ...i, unsigned long int ...j >
    auto createContainer( _indices<i...> row, _indices<j...> col, const std::string name )
      -> decltype( std::make_tuple( createContainerRow<i>( col, name )... ) )

    {
      return std::make_tuple( createContainerRow<i>( col, name )... );
    }


    ///// Copy
    template< unsigned long int i, unsigned long int ...j >
    std::tuple< std::shared_ptr< Item2<i,j> > ... >
    copyContainerRow( std::tuple< _index<j>... > )
    {
      return std::make_tuple( std::get<j>( std::get<i>( item2_ ) )... );
    }
    template< unsigned long int ...i, unsigned long int ...j >
    auto copyContainer( std::tuple< _index<i>... > row, std::tuple< _index<j>... > col )
      -> decltype( std::make_tuple( copyContainerRow<i>( col )... ) )
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
    TwoArgContainer( SameObject& obj, const std::string name = "" )
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
    std::shared_ptr< Item2< i, j > > operator() ( _index<i> row, _index<j> col  )
    {
      return std::get<j>( std::get<i>( item2_ ) );
    }

    // sub Container
    template< unsigned long int... i, unsigned long int... j >
    std::shared_ptr< TwoArgContainer< TwoArgImp, OneArgImp, std::tuple< typename std::tuple_element<i,RowArgTupleType>::type...>,
                                                            std::tuple< typename std::tuple_element<j,ColArgTupleType>::type...> > >
    operator() ( std::tuple< _index<i>... > row, std::tuple< _index<j>... > col )
    {
      typedef std::tuple< typename std::tuple_element<i,RowArgTupleType>::type...> NewRowArgTupleType;
      typedef std::tuple< typename std::tuple_element<j,ColArgTupleType>::type...> NewColArgTupleType;

      typedef OneArgContainer< OneArgImp, NewRowArgTupleType > RowOneArgContainerType;
      typedef OneArgContainer< OneArgImp, NewColArgTupleType > ColOneArgContainerType;

      typedef TwoArgContainer< TwoArgImp, OneArgImp, NewRowArgTupleType, NewColArgTupleType > SubContainerType;

      return std::make_shared< SubContainerType >( RowOneArgContainerType::copyContainer( item1_, row ),
                                                   ColOneArgContainerType::copyContainer( colItem1_, col ),
                                                   copyContainer( row, col ) );
    }
  protected:
    Item2TupleType          item2_;
    ColItem1TupleType       colItem1_;
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


  //TODO this is not a sub-container, it is a global container...
  //Put this class to another header...
  template< class Item2TupleImp, class Item1TupleImp, class SubOrderImp, class... DiscreteFunctions >
  class GlobalContainer
    : public TwoArgContainer< _template< Item2TupleImp >::template _t2, _template< Item1TupleImp >::template _t1,
                              std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > >
  {
    typedef TwoArgContainer< _template< Item2TupleImp >::template _t2, _template< Item1TupleImp >::template _t1,
                             std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > > BaseType;

    static_assert( is_tuple< SubOrderImp >::value, "SubOrderImp has to be a std::tuple<std::tuple<>...>" );
  public:
    using BaseType::operator();

    // constructor: do not touch, delegate everything
    template< class ... Args>
    GlobalContainer( Args&&... args )
    : BaseType( args... )
    {}

    // sub container
    template< unsigned long int i >
    decltype(auto) sub( _index<i> index )
    {
      static_assert( std::tuple_size< SubOrderImp >::value > i,
                     "SubOrderImp does not contain the necessary sub structure information." );
      //default implementation:
      //static const typename index_tuple< DiscreteFunctions... >::type order;
      static const SubOrderImp order;
      return BaseType::operator()( std::get<i>(order), std::get<i>(order) );
    }
  };

}
}
#endif // FEMHOWTO_STEPPER_HH
