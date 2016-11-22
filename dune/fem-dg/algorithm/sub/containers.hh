#ifndef DUNE_FEMDG_SUB_CONTAINER_HH
#define DUNE_FEMDG_SUB_CONTAINER_HH

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


  /**
   * \brief A placeholder for empty item1 or item2 container.
   */
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



}
}
#endif // FEMHOWTO_STEPPER_HH
