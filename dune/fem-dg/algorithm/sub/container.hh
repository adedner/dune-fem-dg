#ifndef DUNE_FEMDG_CONTAINER_HH
#define DUNE_FEMDG_CONTAINER_HH

#include <memory>



namespace Dune
{
namespace Fem
{

  //helper class to generic pack a matrix and vector
  namespace
  {

    template< template<class,int...> class, class, int >
    struct VectorPack;

    //Normal
    template<template<class> class OneArgImp, class... Rows, int row >
    struct VectorPack< OneArgImp, std::tuple< Rows... >, row >
    {
      typedef std::tuple< std::shared_ptr< OneArgImp< Rows > >...  > shared_type;
      typedef std::tuple<                  OneArgImp< Rows >  ...  >        type;
    };

    //Generic, i.e. each entry can be indivually choosen via <int> argument
    template<template<class,int> class OneArgImp, class RowHead, int row >
    struct VectorPack< OneArgImp, std::tuple< RowHead >, row >
    {
      typedef std::tuple< std::shared_ptr< typename OneArgImp< RowHead, row >::type > > shared_type;
      typedef std::tuple<                  typename OneArgImp< RowHead, row >::type   >        type;
    };

    template<template<class,int> class OneArgImp, class RowHead, class... Rows, int row >
    struct VectorPack< OneArgImp, std::tuple< RowHead, Rows... >, row >
    {
    private:
      typedef std::tuple< std::shared_ptr< typename OneArgImp< RowHead, row >::type > >                                SharedTupleHead;
      typedef std::tuple<                  typename OneArgImp< RowHead, row >::type   >                                      TupleHead;
      typedef typename drop_tuple< typename VectorPack< OneArgImp, std::tuple< Rows >..., row+1 >::shared_type >::type SharedTupleEnd;
      typedef typename drop_tuple< typename VectorPack< OneArgImp, std::tuple< Rows >..., row+1 >::type        >::type       TupleEnd;
    public:
      typedef std::tuple< SharedTupleHead, SharedTupleEnd > shared_type;
      typedef std::tuple< TupleHead,       TupleEnd       >        type;

    };

    template< template<class,int...> class OneArgImp, class Rows >
    struct VectorPacker
    {
      typedef typename VectorPack< OneArgImp, Rows, 0 >::type               type;
      typedef typename VectorPack< OneArgImp, Rows, 0 >::shared_type shared_type;
    };




    template< template<class,class,int...> class PairImp, class Rows, class Cols, int row, int col >
    struct MatrixPack;

    //Normal
    template<template<class,class> class PairImp, class RowHead, class... Cols, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< Cols... >, row, col >
    {
      //typedef std::tuple< std::tuple<                  typename PairImp::template type< RowHead, Cols > ... > >        type;
      typedef std::tuple< std::tuple<                  typename PairImp< RowHead, Cols >::type  ... > >        type;
      typedef std::tuple< std::tuple< std::shared_ptr< typename PairImp< RowHead, Cols >::type >... > > shared_type;
    };

    template<template<class,class> class PairImp, class RowHead, class RowHead2, class... Rows, class... Cols, int row, int col>
    struct MatrixPack< PairImp, std::tuple< RowHead, RowHead2, Rows... >, std::tuple< Cols... >, row, col >
    {
    private:
      typedef std::tuple< RowHead2, Rows... > RowTuple;
      typedef std::tuple< Cols... > ColTuple;

      typedef std::tuple<                  typename PairImp< RowHead, Cols >::type  ... >       TupleHead;
      typedef std::tuple< std::shared_ptr< typename PairImp< RowHead, Cols >::type >... > SharedTupleHead;
      typedef typename drop_tuple< typename MatrixPack< PairImp, RowTuple, ColTuple, row, col >::type        >::type       TupleEnd;
      typedef typename drop_tuple< typename MatrixPack< PairImp, RowTuple, ColTuple, row, col >::shared_type >::type SharedTupleEnd;
    public:
      typedef std::tuple< TupleHead,       TupleEnd       > type;
      typedef std::tuple< SharedTupleHead, SharedTupleEnd > shared_type;
    };

    //Generic, i.e. each entry can be indivually choosen via <int,int> argument
    template<template<class,class,int,int> class PairImp, class RowHead, class ColHead, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< ColHead >, row, col >
    {
      typedef std::tuple< std::tuple<                  typename PairImp< RowHead, ColHead, row, col >::type >   >  type;
      typedef std::tuple< std::tuple< std::shared_ptr< typename PairImp< RowHead, ColHead, row, col >::type > > >  shared_type;
    };

    template<template<class,class,int,int> class PairImp, class RowHead, class ColHead, class... Cols, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< ColHead, Cols... >, row, col >
    {
    private:
      typedef                  typename PairImp< RowHead, ColHead, row, col >::type   TupleHead;
      typedef std::shared_ptr< typename PairImp< RowHead, ColHead, row, col >::type > SharedTupleHead;
      typedef typename drop_tuple< typename MatrixPack< PairImp, std::tuple< RowHead>, std::tuple< Cols...>, row, col+1 >::type        >::type       TupleEnd;
      typedef typename drop_tuple< typename MatrixPack< PairImp, std::tuple< RowHead>, std::tuple< Cols...>, row, col+1 >::shared_type >::type SharedTupleEnd;
    public:
      typedef std::tuple< std::tuple< TupleHead,       TupleEnd       > >         type;
      typedef std::tuple< std::tuple< SharedTupleHead, SharedTupleEnd > >  shared_type;
    };

    template<template<class,class,int,int> class PairImp, class RowHead, class RowHead2, class... Rows, class... Cols, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead, RowHead2, Rows... >, std::tuple< Cols... >, row, col >
    {
    private:
      typedef std::tuple< RowHead2, Rows... > RowTuple;
      typedef std::tuple< Cols... > ColTuple;

      typedef std::tuple<                  typename PairImp< RowHead, Cols, row, col >::type  ... >                        TupleHead;
      typedef std::tuple< std::shared_ptr< typename PairImp< RowHead, Cols, row, col >::type >... >                  SharedTupleHead;
      typedef typename drop_tuple< typename MatrixPack< PairImp, RowTuple, ColTuple, row+1, 0 >::type        >::type       TupleEnd;
      typedef typename drop_tuple< typename MatrixPack< PairImp, RowTuple, ColTuple, row+1, 0 >::shared_type >::type SharedTupleEnd;
    public:
      typedef std::tuple< std::tuple< TupleHead,       TupleEnd       > >         type;
      typedef std::tuple< std::tuple< SharedTupleHead, SharedTupleEnd > >  shared_type;
    };

    template< template<class,class,int...> class TwoArgImp, class Rows, class Cols >
    struct MatrixPacker
    {
      typedef typename MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::type               type;
      typedef typename MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::shared_type shared_type;
    };

  }

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
  template< template<class,int...> class, class >
  struct OneArgContainer;

  //forward declaration
  template< template<class,class,int...> class, template<class,int...> class, class, class >
  class TwoArgContainer;


  //only for tuples
  template< template< class > class OneArgImp, class... Args >
  struct OneArgContainer< OneArgImp, std::tuple< Args... > >
  {
    typedef std::tuple< Args... >                ArgTupleType;

    typedef typename VectorPacker< OneArgImp, ArgTupleType >::shared_type Item1TupleType;

  public:

    template< unsigned long int i >
    using DiscreteFunction = typename std::tuple_element< i, ArgTupleType>::type;

    template< unsigned long int i >
    using Item1 = typename std::tuple_element< i, Item1TupleType>::type::element_type;

  protected:
    template< unsigned long int... i >
    using SubContainer = OneArgContainer< OneArgImp, Item1<i>... >;

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
    static Item1TupleType createContainer( std::integer_sequence< unsigned long int, i... >, SameObject& obj, const std::string name )
    {
      return std::make_tuple( createItem1<i>( obj, name )... );
    }
    template< unsigned long int ...i >
    static Item1TupleType createContainer( std::integer_sequence< unsigned long int, i... >, const std::string name )
    {
      return std::make_tuple( createItem1<i>( name )... );
    }

    ////// Copy
    template< class Item1Tuple, unsigned long int ...i >
    static Item1TupleType copyContainer( const Item1Tuple& item1, std::tuple< std::integral_constant< unsigned long int, i... > > )//std::integer_sequence< unsigned long int, i... > )
    {
      return std::make_tuple( std::get<i>( item1 )... );
    }

    //Be your own (template) friend
    template< template<class,int...> class, class >
    friend class OneArgContainer;

    template< template<class,class,int...> class, template<class,int...> class, class, class >
    friend class TwoArgContainer;

    // copy, for internal use only
    OneArgContainer( const Item1TupleType& item )
    : item1_( item )
    {}
  public:

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
    std::shared_ptr< Item1<i> > operator() ( std::integral_constant<unsigned long int, i> index )
    {
      return std::get<i>( item1_ );
    }

    // sub container
    template< unsigned long int... i >
    std::shared_ptr< SubContainer< i...> >
    operator() ( std::tuple< std::integral_constant<unsigned long int, i>... > index )
    {
      return std::make_shared< SubContainer< i...> >( copyContainer(item1_, index ) );
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

    using BaseType::item1_;

  protected:
    static const int size = BaseType::size;
    using BaseType::sequence;

    ///// Creation
    template< unsigned long int i, unsigned long int j >
    std::shared_ptr< Item2<i,j> > createItem2( const std::string name )
    {
      return std::make_shared<Item2<i,j> >( BaseType::operator()( std::integral_constant< unsigned long int, i>() ),
                                            BaseType::operator()( std::integral_constant< unsigned long int, j>() ),
                                            name );
    }
    template< unsigned long int i, unsigned long int ...j >
    std::tuple< std::shared_ptr< Item2<i,j> > ... >
    createContainerRow( std::integer_sequence< unsigned long int, j... >,
                        const std::string name )
    {
      return std::make_tuple( createItem2<i,j>( name )... );
    }
    template< unsigned long int ...i, unsigned long int ...j >
    auto createContainer( std::integer_sequence< unsigned long int, i... > row,
                          std::integer_sequence< unsigned long int, j... > col,
                          const std::string name )
      -> decltype( std::make_tuple( createContainerRow<i>( col, name )... ) )

    {
      return std::make_tuple( createContainerRow<i>( col, name )... );
    }


    ///// Copy
    template< unsigned long int i, unsigned long int ...j >
    std::tuple< std::shared_ptr< Item2<i,j> > ... >
    copyContainerRow( std::tuple< std::integral_constant< unsigned long int, j... > > )
    {
      return std::make_tuple( std::get<j>( std::get<i>( item2_ ) )... );
    }
    template< unsigned long int ...i, unsigned long int ...j >
    auto copyContainer( std::tuple< std::integral_constant< unsigned long int, i... > > row,
                        std::tuple< std::integral_constant< unsigned long int, j... > > col )
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
    std::shared_ptr< Item2< i, j > > operator() ( std::integral_constant<unsigned long int, i> row,
                                                  std::integral_constant<unsigned long int, j> col  )
    {
      return std::get<j>( std::get<i>( item2_ ) );
    }

    // sub Container
    template< unsigned long int... i, unsigned long int... j >
    std::shared_ptr< TwoArgContainer< TwoArgImp, OneArgImp, std::tuple< typename std::tuple_element<i,RowArgTupleType>::type...>,
                                                     std::tuple< typename std::tuple_element<j,ColArgTupleType>::type...> > >
    operator() ( std::tuple< std::integral_constant<unsigned long int, i>... > row,
                 std::tuple< std::integral_constant<unsigned long int, j>... > col )
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

  template< template<class> class Obj, template<class,class> class InnerObj  >
  struct TemplateContainer
  {
    template< class R, class C >
    struct Object
    {
      typedef Obj< InnerObj< R,C > > type;
    };
  };

  template< template<class> class Obj, template<int,int> class InnerObj >
  struct GeneralTemplateContainer
  {
    template< class R, class C, int r, int c >
    struct Object
    {
      typedef Obj< typename InnerObj< r, c >::template type<R,C> > type;
    };
  };

}
}
#endif // FEMHOWTO_STEPPER_HH
