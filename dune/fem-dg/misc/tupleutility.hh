#ifndef DUNE_TUPLE_CONCAT_HH
#define DUNE_TUPLE_CONCAT_HH

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/fem/operator/linear/spoperator.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief concats tuples and flattens the result
   *
   * \tparam Args variadic list of tuples
   *
   * Usage
   * \code{.cpp}
   *   typedef std::tuple< double > FirstTuple;
   *   typedef std::tuple< int >    SecondTuple;
   *   typedef std::tuple< bool >   ThirdTuple;
   *
   *   typedef typename tuple_concat< FirstTuple, SecondTuple, ThirdTuple >::type ConcatTuple;
   *   //type of ConcatTuple is std::tuple< double, int, bool >
   * \endcode
   */
  template < class... Args >
  class tuple_concat
  {
    template < class... Args2 >
    static auto apply(Args2&& ... args2) -> decltype(std::tuple_cat(args2...))
    {
      return std::tuple_cat(args2...);
    }
    public:
    //! type of a concatenated, flattened tuple
    typedef decltype(apply(std::declval<Args>()...)) type;
  };


  /**
   * \brief Extracts the first element out of a variadic list of Arguments
   *
   * \tparam Args a variadic lists of elements.
   */
  template < class... Args >
  struct tuple_head
  {
    static_assert( sizeof ... ( Args ) > 0,
                   "tuple_head<> needs at least one template argument" );
    typedef typename std::tuple_element< 0, std::tuple< Args... > >::type type;
  };


  /**
   * \brief typedef check if a type is a std::tuple<>
   *
   * usage
   * \code{.cpp}
   *   typedef
   *   static const int =
   * \endcode
   *
   * \tparam T type to be checked.
   */
  template< class T >
  struct is_tuple
  {
    static const int value = false;
  };

  template< class... TElem >
  struct is_tuple< std::tuple< TElem... > >
  {
    static const int value = true;
  };



  template< class FullTupleImp, class IndexSequenceImp >
  struct tuple_reducer;

  template< class FullTupleImp, std::size_t... i >
  struct tuple_reducer< FullTupleImp, std::index_sequence< i... > >
  {
    typedef std::tuple< typename std::tuple_element< i, FullTupleImp >::type... > type;

    static type apply( const FullTupleImp& tuple )
    {
      return std::make_tuple( std::get< i >( tuple )... );
    }

  };

  //helper struct to write a static_assert() which will always fail
  //note:
  //static_assert( false, "error" );
  //will always be false, especially in partial specializations of a class/struct
  //use
  //static_assert( static_fail< YourClassTemplate >, "error" );
  //instead.
  template<class... T>
  struct static_fail : std::false_type{};

  template<template<class...> class... T>
  struct static_fail_t : std::false_type{};



  template< class TupleImp >
  struct drop_tuple;

  template< class Tuple >
  struct drop_tuple< std::tuple< Tuple > >
  {
    static_assert( static_fail<Tuple>::value ,
                   "We do not unwrap a simple tuple to avoid unexpected behaviour!" );
  };

  template< class... TupleArg >
  struct drop_tuple< std::tuple< std::tuple< TupleArg... > > >
  {
    typedef std::tuple< TupleArg... > type;
  };


  ////hot C++17 stuff, replaces all non-type template overloads
  //template< auto >                  constexpr int template_level(){ return -1; }

  template<bool>                    constexpr int template_level(){ return -1; }
  template<char>                    constexpr int template_level(){ return -1; }
  template<signed char>             constexpr int template_level(){ return -1; }
  template<unsigned char>           constexpr int template_level(){ return -1; }
  template<short int>               constexpr int template_level(){ return -1; }
  template<int>                     constexpr int template_level(){ return -1; }
  template<long int>                constexpr int template_level(){ return -1; }
  template<long long int>           constexpr int template_level(){ return -1; }
  template<unsigned short int>      constexpr int template_level(){ return -1; }
  template<unsigned int>            constexpr int template_level(){ return -1; }
  template<unsigned long int>       constexpr int template_level(){ return -1; }
  template<unsigned long long int>  constexpr int template_level(){ return -1; }

  //simple template
  template<class>                   constexpr int template_level(){ return 0; }
  //template template class
  template<template<class...>class> constexpr int template_level(){ return 1; }
  //template template template class
  template<template<template<class...>class...>class> constexpr int template_level(){ return 2; }


  // is_same for templates...
  template< template<class... > class, template<class...> class >
  struct is_same_template : std::false_type{};

  template< template<class... > class A >
  struct is_same_template<A,A> : std::true_type{};



  //helper class to generic pack a matrix and vector
  namespace details
  {

    //forward declaration
    template< template<class,int...> class, class, int >
    struct VectorPack;
    //forward declaration
    template< template<class,class,int...> class PairImp, class Rows, class Cols, int row, int col >
    struct MatrixPack;


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
      typedef std::tuple< std::shared_ptr< typename OneArgImp< RowHead, row >::type > >   SharedTupleHead;
      typedef std::tuple<                  typename OneArgImp< RowHead, row >::type   >         TupleHead;
      typedef typename VectorPack< OneArgImp, std::tuple< Rows >..., row+1 >::shared_type  SharedTupleEnd;
      typedef typename VectorPack< OneArgImp, std::tuple< Rows >..., row+1 >::type               TupleEnd;
    public:
      typedef typename tuple_concat< SharedTupleHead, SharedTupleEnd >::type shared_type;
      typedef typename tuple_concat< TupleHead,       TupleEnd       >::type        type;

    };


    ////Normal
    //generate one row
    template<template<class,class> class PairImp, class RowHead, class... Cols, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< Cols... >, row, col >
    {
      //typedef std::tuple< std::tuple<                  typename PairImp::template type< RowHead, Cols > ... > >        type;
      typedef std::tuple< std::tuple<                  typename PairImp< RowHead, Cols >::type  ... > >        type;
      typedef std::tuple< std::tuple< std::shared_ptr< typename PairImp< RowHead, Cols >::type >... > > shared_type;
    };

    //generate whole matrix
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

    ////Generic, i.e. each entry can be indivually choosen via <int,int> argument
    //generate last col element of a row
    template<template<class,class,int,int> class PairImp, class RowHead, class ColHead, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< ColHead >, row, col >
    {
      typedef std::tuple< std::tuple<                  typename PairImp< RowHead, ColHead, row, col >::type >   >  type;
      typedef std::tuple< std::tuple< std::shared_ptr< typename PairImp< RowHead, ColHead, row, col >::type > > >  shared_type;
    };

    //generate one row
    template<template<class,class,int,int> class PairImp, class RowHead, class ColHead, class... Cols, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< ColHead, Cols... >, row, col >
    {
    private:
      typedef std::tuple<                  typename PairImp< RowHead, ColHead, row, col >::type   > TupleHead;
      typedef std::tuple< std::shared_ptr< typename PairImp< RowHead, ColHead, row, col >::type > > SharedTupleHead;
      typedef typename drop_tuple< typename MatrixPack< PairImp, std::tuple< RowHead>, std::tuple< Cols...>, row, col+1 >::type        >::type       TupleEnd;
      typedef typename drop_tuple< typename MatrixPack< PairImp, std::tuple< RowHead>, std::tuple< Cols...>, row, col+1 >::shared_type >::type SharedTupleEnd;
    public:
      typedef std::tuple< typename tuple_concat< TupleHead,       TupleEnd       >::type >        type;
      typedef std::tuple< typename tuple_concat< SharedTupleHead, SharedTupleEnd >::type > shared_type;
    };

    //generate whole matrix
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
      typedef std::tuple< TupleHead,       TupleEnd       >         type;
      typedef std::tuple< SharedTupleHead, SharedTupleEnd >  shared_type;
    };
  }


  template< template<class,int...> class OneArgImp, class Rows >
  struct VectorPacker
  {
    typedef typename details::VectorPack< OneArgImp, Rows, 0 >::type               type;
    typedef typename details::VectorPack< OneArgImp, Rows, 0 >::shared_type shared_type;

    static const int size = std::tuple_size< Rows >::value;
  };

  template< template<class,class,int...> class TwoArgImp, class Rows, class Cols >
  struct MatrixPacker
  {
    typedef typename details::MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::type               type;
    typedef typename details::MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::shared_type shared_type;

    static const int rows = std::tuple_size< Rows >::value;
    static const int col = std::tuple_size< Rows >::value;
  };


  //forward declaration
  template< class TupleMatrix >
  struct tuple_matrix;

  template< class... >
  struct tuple_min_rows;

  template< class Row >
  struct tuple_min_rows< Row >
  {
    static const unsigned long int value = std::tuple_size<Row>::value;
  };

  template< class Row1, class Row2, class... Args >
  struct tuple_min_rows< Row1, Row2, Args... >
  {
    static const unsigned long int value = std::min( std::tuple_size<Row1>::value, tuple_min_rows< Row2, Args... >::value );
  };

  template< class... >
  struct tuple_max_rows;

  template< class Row >
  struct tuple_max_rows< Row >
  {
    static const unsigned long int value = std::tuple_size<Row>::value;
  };

  template< class Row1, class Row2, class... Args >
  struct tuple_max_rows< Row1, Row2, Args... >
  {
    static const unsigned long int value = std::max( std::tuple_size<Row1>::value, tuple_max_rows< Row2, Args... >::value );
  };

  template< class... Args >
  struct tuple_matrix< std::tuple<Args...> >
  {

    static const int rows = std::tuple_size< std::tuple<Args...> >::value;
    static const int cols = tuple_max_rows< Args... >::value;

    //check whether is rectangular
    static const bool isMatrix = (tuple_max_rows< Args... >::value == tuple_min_rows< Args... >::value);

  };


  //forward declaration
  template< int size, class TupleElement >
  struct tuple_copy
  {
    typedef typename tuple_concat< std::tuple< TupleElement >, typename tuple_copy<size-1,TupleElement>::type >::type type;
  };

  template< class TupleElement >
  struct tuple_copy<1,TupleElement>
  {
    typedef std::tuple<TupleElement> type;
  };

  template< class TupleElement >
  struct tuple_copy<0,TupleElement>
  {
    typedef std::tuple<> type;
  };

  //forward declaration
  template< int size, template<unsigned long int> class TupleElement >
  struct tuple_copy_t
  {
    typedef typename tuple_concat< typename tuple_copy_t<size-1,TupleElement>::type, std::tuple< TupleElement<size-1> > >::type type;
  };

  template< template<unsigned long int> class TupleElement >
  struct tuple_copy_t<1,TupleElement>
  {
    typedef std::tuple< TupleElement<0> > type;
  };

  template< template<unsigned long int> class TupleElement >
  struct tuple_copy_t<0,TupleElement>
  {
    typedef std::tuple<> type;
  };

  //forward declaration
  template< class TupleMatrix1, class TupleMatrix2, class Default >
  struct tuple_matrix_combine;

  template< class TupleRow1, int cols2, class Default/*= _t< EmptyContainerItem >*/ >
  struct tuple_matrix_combine_row1
  {
    typedef typename tuple_concat< TupleRow1,
                                   typename tuple_copy<cols2,Default>::type>::type type;
  };

  template< class TupleRow2, int cols1, class Default/*= _t< EmptyContainerItem >*/ >
  struct tuple_matrix_combine_row2
  {
    typedef typename tuple_concat< typename tuple_copy<cols1,Default>::type,
                                   TupleRow2 >::type type;
  };

  template< class... TupleRowArgs1, class... TupleRowArgs2, class Default >
  struct tuple_matrix_combine<std::tuple<TupleRowArgs1...>, std::tuple<TupleRowArgs2... >, Default >
  {
    static const int col1 = tuple_matrix<std::tuple<TupleRowArgs1...> >::cols;
    static const int col2 = tuple_matrix<std::tuple<TupleRowArgs2...> >::cols;

    typedef std::tuple< typename tuple_matrix_combine_row1< TupleRowArgs1, col2, Default >::type... > Tuple1;
    typedef std::tuple< typename tuple_matrix_combine_row2< TupleRowArgs2, col1, Default >::type... > Tuple2;

    typedef typename tuple_concat< Tuple1, Tuple2 >::type type;


    //static const int rows1 = tuple_matrix< TupleMatrix1 >::rows;
    //static const int cols1 = tuple_matrix< TupleMatrix1 >::cols;
    //static const bool isMatrix1 = tuple_matrix< TupleMatrix1 >::isMatrix;

    //static const int rows2 = tuple_matrix< TupleMatrix2 >::rows;
    //static const int cols2 = tuple_matrix< TupleMatrix2 >::cols;
    //static const bool isMatrix2 = tuple_matrix< TupleMatrix2 >::isMatrix;

    //static_assert( isMatrix1 && isMatrix2, "Sorry, only rectangular matrices, i.e. real matrices, allowed!" );

  };


  namespace details
  {
    //simple struct to produce nested TupleRows
    template< template<class...> class... >
    struct Nested;

    template< template<class...> class ArgHead >
    struct Nested< ArgHead >
    {
      template<class... I>
      struct _t{ typedef ArgHead<I...> type; };
    };

    template< template<class> class ArgHead, template<class...> class ArgHead2, template<class...> class... Args >
    struct Nested< ArgHead, ArgHead2, Args... >
    {
      template<class... I>
      struct _t{ typedef ArgHead<typename Nested<ArgHead2,Args...>::template _t<I...>::type > type; };
    };
  }


  // create a template tuple via
  // std::tuple< _t< YourTemplates... > >
  template< template< class... >class ... Args >
  using _t = details::Nested< Args... >;

  //expects a tuple of tuple
  template< class TupleImp >
  struct _template
  {
    template<class R,class C,int r,int c>
    struct _t2{ typedef typename std::tuple_element<c,typename std::tuple_element<r,TupleImp>::type>::type::template _t<R,C>::type type; };

    template<class R,class C,int r,int c>
    struct _t2Inv{ typedef typename std::tuple_element<c,typename std::tuple_element<r,TupleImp>::type>::type::template _t<C,R>::type type; };

    template<class R,int r>
    struct _t1{ /*static_assert( r<std::tuple_size<TupleImp>::value, "selected tuple element does not exist." );*/typedef typename std::tuple_element<r,TupleImp>::type::template _t<R>::type type; };
  };

  //forward declaration, expects a template which is valid for all
  template< template<class...> class... >
  struct _template2;

  template< template<class...> class ArgHead >
  struct _template2< ArgHead >
  {
    template<class R,class C,int r,int c>
    struct _t2{ typedef ArgHead< R,C > type; };

    template<class R,class C,int r,int c>
    struct _t2Inv{ typedef ArgHead< C,R > type; };

    template<class R,int r>
    struct _t1{ typedef ArgHead<R> type; };
  };

  template< template<class...> class ArgHead, template<class...> class ArgHead2, template<class...> class... Args >
  struct _template2< ArgHead, ArgHead2, Args... >
  {
    template<class R,class C,int r,int c>
    struct _t2{ typedef ArgHead< typename _template2<ArgHead2,Args... >::template _t2< R,C,r,c >::type > type; };

    template<class R,int r>
    struct _t1{ typedef ArgHead< typename _template2<ArgHead2,Args...>::template _t1<R,r>::type > type; };
  };

}
}

#endif
