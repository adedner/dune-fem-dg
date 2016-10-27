#ifndef DUNE_TUPLE_CONCAT_HH
#define DUNE_TUPLE_CONCAT_HH

#include <tuple>
#include <type_traits>
#include <utility>

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

  template< class TupleImp >
  struct drop_tuple;

  template< class Tuple >
  struct drop_tuple< std::tuple< Tuple > >
  {
    static_assert( std::is_same<std::tuple< Tuple >,std::tuple< Tuple > >::value ,
                   "We do not unwrap a simple tuple to avoid unexpected behaviour!" );
  };

  template< class... TupleArg >
  struct drop_tuple< std::tuple< std::tuple< TupleArg... > > >
  {
    typedef std::tuple< TupleArg... > type;
  };


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

    template< template<class,class,int...> class TwoArgImp, class Rows, class Cols >
    struct MatrixPacker
    {
      typedef typename MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::type               type;
      typedef typename MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::shared_type shared_type;
    };

  }


  //helper class to store a variadic couple template parameters
  namespace{

    //forward declarations
    template<template<class...>class...>
    struct TemplateTupleRow;

    //forward declarations
    template<class...>
    struct TemplateTuple;

    //Row, last element
    template<template<class...>class ArgHead >
    struct TemplateTupleRow<ArgHead>
    {
      template<class... I>
      struct _t
      {
        typedef std::tuple<ArgHead<I...> > Tuple;

        template< int i >
        using type = typename std::tuple_element<i,Tuple>::type;
      };
    };

    //Row, general
    template< template<class...>class ArgHead,
              template<class...>class ArgHead2,
              template<class...>class... ArgEnd >
    struct TemplateTupleRow< ArgHead, ArgHead2, ArgEnd... >
    {
      template<class... I>
      struct _t
      {
        typedef typename tuple_concat< std::tuple< ArgHead<I...> >,
                                       typename TemplateTupleRow< ArgHead2,ArgEnd... >::template _t<I...>::Tuple
                                      >::type                                                              Tuple;

        template< int i >
        using type = typename std::tuple_element<i,Tuple>::type;
      };
    };

    //[Row x] Cols, last element
    template< template<class,class>class... ArgRowHead >
    struct TemplateTuple< TemplateTupleRow< ArgRowHead...> >
    {
      template<class Row,class Col>
      struct _t
      {
        typedef std::tuple< typename TemplateTupleRow< ArgRowHead...>::template _t< Row, Col>::Tuple > Tuple;

        template< int r, int c >
        using type = typename std::tuple_element<r,typename std::tuple_element<c,Tuple>::type>::type;
      };
    };

    //[Row x] Cols, general
    template< template<class,class>class... ArgRowHead,
              class TemplateArgRowHead,
              class... TemplateArgRowEnd >
    struct TemplateTuple< TemplateTupleRow< ArgRowHead...>, TemplateArgRowHead, TemplateArgRowEnd... >
    {
      template<class Row,class Col>
      struct _t
      {
        typedef std::tuple< typename TemplateTupleRow< ArgRowHead...>::template _t< Row, Col>::Tuple,
                            typename drop_tuple< typename TemplateTuple< TemplateArgRowHead, TemplateArgRowEnd... >::template _t<Row,Col>::Tuple >::type >
                                                                                                                Tuple;

        template< int r, int c >
        using type = typename std::tuple_element<r,typename std::tuple_element<c,Tuple>::type>::type;
      };
    };

  }

  //No (to much) advanced variadic template template magic here
  //because it would not allow an easy usable template structure...

  //forward declaration
  template< template<class...> class... >
  struct template_unique;

  template< template<class...> class... >
  struct template_general;


  ////
  //matrix and vector (1/2 arg type)
  template< template<class...> class InnerObj >
  struct template_unique< InnerObj >
  {
    template<class... I>
    struct _t
    {
      typedef InnerObj<I...>  type;
    };
  };

  template< template<class> class OuterObj,
            template<class...> class InnerObj >
  struct template_unique< OuterObj, InnerObj >
  {
    template<class... I>
    struct _t
    {
      typedef OuterObj< InnerObj<I...> > type;
    };
  };

  template< template<class> class OutestObj,
            template<class> class OuterObj,
            template<class...> class InnerObj >
  struct template_unique< OutestObj, OuterObj, InnerObj >
  {
    template<class... I>
    struct _t
    {
      typedef OutestObj< OuterObj< InnerObj<I...> > > type;
    };
  };


  ////
  //matrix (2 arg type)
  template< template<class,class> class InnerObj >
  struct template_general< InnerObj >

  {
    template<class R,class C,int r,int c>
    struct _t
    {
      typedef typename InnerObj<R,C>::template type<r,c> type;
    };
  };
  template< template<class> class Obj,
            template<class,class> class InnerObj >
  struct template_general< Obj, InnerObj >
  {
    template<class R,class C,int r,int c>
    struct _t
    {
      typedef Obj< typename InnerObj<R,C>::template type<r,c> > type;
    };
  };
  template< template<class> class OuterObj,
            template<class> class Obj,
            template<class,class> class InnerObj >
  struct template_general< OuterObj, Obj, InnerObj >
  {
    template<class R,class C,int r,int c>
    struct _t
    {
      typedef OuterObj< Obj< typename InnerObj<R,C>::template type<r,c> > > type;
    };
  };

  //vector (1 arg type)
  template< template<class> class InnerObj >
  struct template_general< InnerObj >

  {
    template<class R,int r>
    struct _t
    {
      typedef typename InnerObj<R>::template type<r> type;
    };
  };
  template< template<class> class Obj,
            template<class> class InnerObj >
  struct template_general< Obj, InnerObj >
  {
    template<class R,int r>
    struct _t
    {
      typedef Obj< typename InnerObj<R>::template type<r> > type;
    };
  };
  template< template<class> class OuterObj,
            template<class> class Obj,
            template<class> class InnerObj >
  struct template_general< OuterObj, Obj, InnerObj >
  {
    template<class R,int r>
    struct _t
    {
      typedef Obj< typename InnerObj<R>::template type<r> > type;
    };
  };





}
}

#endif
