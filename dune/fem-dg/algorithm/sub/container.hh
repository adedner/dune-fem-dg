#ifndef DUNE_FEMDG_CONTAINER_HH
#define DUNE_FEMDG_CONTAINER_HH

#include <memory>



namespace Dune
{
namespace Fem
{

  //helper class to generic pack a matrix
  namespace
  {
    template< template<class,class> class PairImp, class Rows, class Cols >
    struct MatrixPack;

    template<template<class,class> class PairImp, class RowHead, class... Cols >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< Cols... > >
    {
      typedef std::tuple< std::tuple< std::shared_ptr< PairImp< RowHead, Cols > > ... > > type;
    };

    template<template<class,class> class PairImp, class RowHead, class RowHead2, class... Rows, class... Cols >
    struct MatrixPack< PairImp, std::tuple< RowHead, RowHead2, Rows... >, std::tuple< Cols... > >
    {
    private:
      typedef std::tuple< RowHead2, Rows... > RowTuple;
      typedef std::tuple< Cols... > ColTuple;

      typedef std::tuple< std::shared_ptr< PairImp< RowHead, Cols > >... > TupleHead;
      typedef typename drop_tuple< typename MatrixPack< PairImp, RowTuple, ColTuple >::type >::type TupleEnd;
    public:
      typedef std::tuple< TupleHead, TupleEnd > type;
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

  template< class >
  struct OneArgContainer;


  //only for tuples
  template< class... Args >
  struct OneArgContainer< std::tuple< Args... > >
  {
    typedef std::tuple< std::shared_ptr< Args >... >               Item1TupleType;
    typedef std::tuple< typename Args::Object ... >                ObjectTupleType;

  public:

    template< unsigned long int i >
    using DiscreteFunction = typename std::tuple_element< i, ObjectTupleType>::type;

    template< unsigned long int i >
    using Item1 = typename std::tuple_element< i, Item1TupleType>::type::element_type;

  protected:
    template< unsigned long int... i >
    using SubContainer = OneArgContainer< Item1<i>... >;

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
    template< unsigned long int ...i >
    Item1TupleType copyContainer( std::integer_sequence< unsigned long int, i... > )
    {
      return std::make_tuple( std::get<i>( item1_ )... );
    }

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
      return std::make_shared< SubContainer< i...> >( copyContainer( index ) );
    }
  protected:
    Item1TupleType item1_;
  };





}
}
#endif // FEMHOWTO_STEPPER_HH
