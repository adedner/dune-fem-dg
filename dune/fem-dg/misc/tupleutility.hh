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

  /**
   * \brief helper struct to write a static_assert() which will always fail
   *
   * Note that, e.g.
   *
   * \code
   * static_assert( false, "error" );
   * \ȩndcode
   *
   * will always be false, especially in partial specializations of a class/struct.
   * Thus, a small workaround is needed.
   *
   * Use
   *
   * \code
   * static_assert( static_fail< YourClassTemplate >, "error" );
   * \endcode
   *
   * instead.
   */
  template<class... T>
  struct static_fail : std::false_type{};

  template<template<class...> class... T>
  struct static_fail_t : std::false_type{};



  /**
   * \brief helper class to drop a double nested std::tuple.
   */
  template< class Tuple >
  struct drop_tuple
  {
    static_assert( static_fail<Tuple>::value,
                   "Argument has to be a std::tuple< std::tuple<...> >!" );
  };

  /**
   * \brief partial specialization.
   */
  template< class Tuple >
  struct drop_tuple< std::tuple< Tuple > >
  {
    static_assert( static_fail<Tuple>::value ,
                   "We do not unwrap a simple (non nested) std::tuple to avoid unexpected behaviour!" );
  };

  /**
   * \brief partial specialization of helper class to drop a double nested std::tuple.
   *
   * This means, that the following statement is true:
   * \code
   * static const bool res = std::is_same< std::tuple< std::tuple< int, double > >, std::tuple<int, double > >::value;
   * \endcode
   */
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


  /**
   * \brief std::is_same only allows to compare types of classes. This is a template template version.
   */
  template< template<class... > class, template<class...> class >
  struct is_same_template : std::false_type{};

  /**
   * \brief partial specialization, same type;
   */
  template< template<class... > class A >
  struct is_same_template<A,A> : std::true_type{};


  namespace details
  {

    /**
     * \brief helper class used by VectorPacker class.
     */
    template< template<class,int...> class, class, int >
    struct VectorPack;

    /**
     * \brief helper class used by MatrixPacker class.
     */
    template< template<class,class,int...> class PairImp, class Rows, class Cols, int row, int col >
    struct MatrixPack;


    /**
     * \brief partial specialization for unique* template parameters.
     */
    template<template<class> class OneArgImp, class... Rows, int row >
    struct VectorPack< OneArgImp, std::tuple< Rows... >, row >
    {
      typedef std::tuple< std::shared_ptr< OneArgImp< Rows > >...  > shared_type;
      typedef std::tuple<                  OneArgImp< Rows >  ...  >        type;
    };

    /**
     * \brief partial specialization for *non unique* template parameters, last element.
     */
    template<template<class,int> class OneArgImp, class RowHead, int row >
    struct VectorPack< OneArgImp, std::tuple< RowHead >, row >
    {
      typedef std::tuple< std::shared_ptr< typename OneArgImp< RowHead, row >::type > > shared_type;
      typedef std::tuple<                  typename OneArgImp< RowHead, row >::type   >        type;
    };

    /**
     * \brief partial specialization for *non unique* template parameters.
     */
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


    /**
     * \brief partial specialization for unique* template parameters, last row.
     */
    template<template<class,class> class PairImp, class RowHead, class... Cols, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< Cols... >, row, col >
    {
      typedef std::tuple< std::tuple<                  typename PairImp< RowHead, Cols >::type  ... > >        type;
      typedef std::tuple< std::tuple< std::shared_ptr< typename PairImp< RowHead, Cols >::type >... > > shared_type;
    };

    /**
     * \brief partial specialization for unique* template parameters.
     */
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

    /**
     * \brief partial specialization for *non unique* template parameters, last element in last row.
     */
    template<template<class,class,int,int> class PairImp, class RowHead, class ColHead, int row, int col >
    struct MatrixPack< PairImp, std::tuple< RowHead >, std::tuple< ColHead >, row, col >
    {
      typedef std::tuple< std::tuple<                  typename PairImp< RowHead, ColHead, row, col >::type >   >  type;
      typedef std::tuple< std::tuple< std::shared_ptr< typename PairImp< RowHead, ColHead, row, col >::type > > >  shared_type;
    };

    /**
     * \brief partial specialization for *non unique* template parameters, last row.
     */
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

    /**
     * \brief partial specialization for *non unique* template parameters.
     */
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


  /*
   * \brief Helper class to create a std::tuple<> of classes which are described by a templated class.
   *
   * The templated class has to contain a typedef called 'type' and contain appropiate
   * template arguments. Despite the fact, that the first template argument accepted by `VectorPacker`'
   * allows one type argument and arbritary* `int` arguments, only two versions are valid
   * (one zero and one `int` arguments). Others implementation will (hopefully) fail.
   *
   * We will give a detailed explanation of both versions:
   *
   *
   * Creating a tuple of unique (wrapped) elements
   * ---------------------------------------------
   *
   * The simplest version is, neglecting all `int` arguments.
   * Then, we can just define a simple class
   * \code
   * template< class Arg >
   * struct OneUniqueArg
   * {
   *   //The following codeline is the simplest typedef, but not
   *   //the way we want to use this class ('cause it is boring...)
   *   //typedef int type;
   *
   *   //more fancy stuff: Yes, we want to wrap an abritary class
   *   //(std::vector, here) around the template argument.
   *   typedef std::vector< Arg > type;
   * };
   * \endcode
   *
   * Furthermore, let us define a tuple of some arbitrary classes.
   *
   * \code
   *   typedef std::tuple< int, char, std::string > ArgTupleType;
   * \endcode
   *
   * Using the typedefs
   *
   * \code
   *   typedef typename VectorPacker< OneUniqueArg, ArgTupleType >::type         TupleType;
   *   typedef typename VectorPacker< OneUniqueArg, ArgTupleType >::shared_type  SharedTupleType;
   * \endcode
   *
   * the following code should be acceptable
   *
   * \code
   * typedef std::tuple< std::vector< int >,
   *                     std::vector< char >,
   *                     std::vector< std::string > >                       TupleType2;
   * typedef std::tuple< std::shared_ptr< std::vector< int > >,
   *                     std::shared_ptr< std::vector< char > >,
   *                     std::shared_ptr< std::vector< std::string > > >    SharedTupleType2;
   * static_assert( std::is_same< TupleType, TupleType2 >::value, "Upps!" );
   * static_assert( std::is_same< SharedTupleType, SharedTupleType2 >::value, "Upps!" );
   * \endcode
   *
   * The `shared_type` is mainly a convenience typedef.
   *
   * Maybe, this already looks a little bit fancy, but we want to go
   * a step deeper.
   *
   *
   * Creating a tuple of individual (wrapped) elements
   * -------------------------------------------------
   *
   * Now, we want to exchange the `std::vector<>` of the last section by a
   * `std::list<>`, but: Only* for the second* tuple element.
   *
   * This means, we want to define
   * \code
   * std::tuple< std::vector< int >,
   *             std::list< char >,
   *             std::vector< std::string > >;
   * \endcode
   *
   * This is the point, where the `int` argument of the template class comes into play:
   *
   * \code
   * template< class Arg, int >
   * struct OneNonUniqueArg;
   * \endcode
   *
   * Via partial specialization of classes, it is now possible to define different
   * cases, i.e.
   *
   * \code
   * // first element: std::vector
   * template< class Arg >
   * struct OneNonUniqueArg< Arg, 0 >
   * { typedef std::vector< Arg > type; };
   *
   * // second element: std::list
   * template< class Arg >
   * struct OneNonUniqueArg< Arg, 1 >
   * { typedef std::list< Arg > type; };
   *
   * // third element: std::vector
   * template< class Arg >
   * struct OneNonUniqueArg< Arg, 2 >
   * { typedef std::vector< Arg > type;};
   * \endcode
   *
   * Now, using the typedefs
   *
   * \code
   *   typedef typename VectorPacker< OneNonUniqueArg, ArgTupleType >::type         TupleType3;
   *   typedef typename VectorPacker< OneNonUniqueArg, ArgTupleType >::shared_type  SharedTupleType3;
   * \endcode
   *
   * the following code should be acceptable
   *
   * \code
   * typedef std::tuple< std::vector< int >,
   *                     std::list< char >,
   *                     std::vector< std::string > >                       TupleType4;
   * typedef std::tuple< std::shared_ptr< std::vector< int > >,
   *                     std::shared_ptr< std::list< char > >,
   *                     std::shared_ptr< std::vector< std::string > > >    SharedTupleType4;
   * static_assert( std::is_same< TupleType3, TupleType4 >::value, "Upps again!" );
   * static_assert( std::is_same< SharedTupleType3, SharedTupleType4 >::value, "Upps again!" );
   * \endcode
   *
   * \warning We need this helper class for generic definition of complex type structures.
   * Using this class directly in the way described above does not really make sense,
   * because (as you might have noticed) it may be much to complex.
   * To escape this template hell, you should also have a look at
   * the template alias _t (which is a shortcut for the class details::Nested) and
   * the _template and _template2 wrapper classes.
   *
   *
   * \tparam OneArgImp template class describing
   * \tparam Rows std::tuple of arguments
   *
   */
  template< template<class,int...> class OneArgImp, class Rows >
  struct VectorPacker
  {
    static_assert( is_tuple<Rows>::value, "The template argument 'Rows' has to be a std::tuple<>!" );
    typedef typename details::VectorPack< OneArgImp, Rows, 0 >::type               type;
    typedef typename details::VectorPack< OneArgImp, Rows, 0 >::shared_type shared_type;

    static const int size = std::tuple_size< Rows >::value;
  };

  /*
   * \brief Helper class to create a std::tuple< std::tuple<>... > of classes which are
   * described by a templated class.
   *
   * In general, this class is a generalization of the class VectorPacker with the only
   * difference that we have got two arguments and the resulting type is a std::tuple<> of
   * some std::tuple<>.
   *
   * The templated class `TwoArgImp` has to contain a typedef called 'type' and contain appropiate
   * template arguments. Despite the fact, that the first template argument accepted by `MatrixPacker`'
   * allows two type arguments and arbritary* `int` arguments, only two versions are valid
   * (one zero and two `int` arguments). Others implementation will (hopefully) fail.
   *
   * We will give a detailed explanation of both versions.
   *
   * \note Since this is very similar to the VectorPacker class, we will drop the std::shared_ptr<>
   * case, here.
   *
   *
   * Creating a tuple of unique (wrapped) elements
   * ---------------------------------------------
   *
   * The simplest version is, neglecting all `int` arguments.
   * Then, we can just define a simple class
   * \code
   * template< class Row, class Col >
   * struct TwoUniqueArg
   * {
   *   //more fancy stuff: Yes, we want to wrap an abritary class
   *   //(std::pair, here) around the template argument.
   *   typedef std::pair< Row, Col > type;
   * };
   * \endcode
   *
   * Furthermore, let us define a tuple of some arbitrary classes.
   *
   * \code
   *   typedef std::tuple< int, char, std::string > ArgTupleType;
   * \endcode
   *
   * Using the typedefs
   *
   * \code
   *   typedef typename MatrixPacker< TwoUniqueArg, ArgTupleType, ArgTupleType >::type         TupleType;
   * \endcode
   *
   * the following code should be acceptable
   *
   * \code
   * using std::string;
   * typedef std::tuple< std::tuple< std::pair< int,    int >, std::pair< int,    char >, std::pair< int,    string > >,
   *                     std::tuple< std::pair< char,   int >, std::pair< char,   char >, std::pair< char,   string > >,
   *                     std::tuple< std::pair< string, int >, std::pair< string, char >, std::pair< string, string > > >
   *                                                                              TupleType2;
   * static_assert( std::is_same< TupleType, TupleType2 >::value, "Upps! again" );
   * \endcode
   *
   * Now, we want to go a step deeper.
   *
   *
   * Creating a tuple of individual (wrapped) elements
   * -------------------------------------------------
   *
   * Now, we want to exchange the `std::pair<>` of the last section by a
   * `std::map<>`, but: Only* for the tuple element in the *second row* and *second column*.
   *
   * This means, we want to define
   * \code
   * using std::string;
   * std::tuple< std::tuple< std::pair< int,    int >, std::pair< int,    char >, std::pair< int,    string > >,
   *             std::tuple< std::pair< char,   int >, std::map < char,   char >, std::pair< char,   string > >,
   *             std::tuple< std::pair< string, int >, std::pair< string, char >, std::pair< string, string > > >
   * \endcode
   *
   * This is the point, where the `int` arguments of the template class comes into play:
   *
   * \code
   * template< class Row, class Col, int, int >
   * struct TwoNonUniqueArg;
   * \endcode
   *
   * Via partial specialization of classes, it is now possible to define different
   * cases, i.e.
   *
   * \code
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 0, 0 >
   * { typedef std::pair< Row, Col > type; };
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 0, 1 >
   * { typedef std::pair< Row, Col > type; };
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 0, 2 >
   * { typedef std::pair< Row, Col > type;};
   *
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 1, 0 >
   * { typedef std::pair< Row, Col > type; };
   * // -----> here comes the std::map! <------
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 1, 1 >
   * { typedef std::map < Row, Col > type; };
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 1, 2 >
   * { typedef std::pair< Row, Col > type;};
   *
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 2, 0 >
   * { typedef std::pair< Row, Col > type; };
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 2, 1 >
   * { typedef std::pair< Row, Col > type; };
   * template< class Row, class Col >
   * struct TwoNonUniqueArg< Row, Col, 2, 2 >
   * { typedef std::pair< Row, Col > type;};
   * \endcode
   *
   * Now, using the typedefs
   *
   * \code
   *   typedef typename MatrixPacker< TwoNonUniqueArg, ArgTupleType >::type         TupleType3;
   * \endcode
   *
   * the following code should be acceptable
   *
   * \code
   * typedef std::tuple< std::tuple< std::pair< int,    int >, std::pair< int,    char >, std::pair< int,    string > >,
   *                     std::tuple< std::pair< char,   int >, std::map < char,   char >, std::pair< char,   string > >,
   *                     std::tuple< std::pair< string, int >, std::pair< string, char >, std::pair< string, string > > >
   *                                                                                TupleType4;
   * static_assert( std::is_same< TupleType3, TupleType4 >::value, "Upps again!" );
   * \endcode
   *
   * \warning We need this helper class for generic definition of complex type structures.
   * Using this class directly in the way described above does not really make sense,
   * because (as you might have noticed) it may be much to complex.
   * To escape this template hell, you should also have a look at
   * the template alias _t (which is a shortcut for the class details::Nested) and
   * the _template and _template2 wrapper classes.
   *
   *
   * \tparam TwoArgImp template class describing
   * \tparam Rows std::tuple of row arguments
   * \tparam Cols std::tuple of col arguments
   *
   */
  template< template<class,class,int...> class TwoArgImp, class Rows, class Cols >
  struct MatrixPacker
  {
    static_assert( is_tuple<Rows>::value, "The template argument 'Rows' has to be a std::tuple<>!" );
    static_assert( is_tuple<Cols>::value, "The template argument 'Cols' has to be a std::tuple<>!" );
    typedef typename details::MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::type               type;
    typedef typename details::MatrixPack< TwoArgImp, Rows, Cols, 0, 0 >::shared_type shared_type;

    static const int rows = std::tuple_size< Rows >::value;
    static const int col = std::tuple_size< Rows >::value;
  };


  namespace details
  {
    /**
     * \brief helper class to
     *
     * See template alias _t for detailed explanation;
     */
    template< template<class...> class... >
    struct Nested;

    /**
     * \brief partical specialization for one template argument.
     */
    template< template<class...> class ArgHead >
    struct Nested< ArgHead >
    {
      template<class... I>
      struct _t{ typedef ArgHead<I...> type; };
    };

    /**
     * \brief partical specialization for more than one template argument.
     */
    template< template<class> class ArgHead, template<class...> class ArgHead2, template<class...> class... Args >
    struct Nested< ArgHead, ArgHead2, Args... >
    {
      template<class... I>
      struct _t{ typedef ArgHead<typename Nested<ArgHead2,Args...>::template _t<I...>::type > type; };
    };
  }


  /**
   * \brief This alias template is an important helper to define a std::tuple
   * (or any other classes expecting classes and not templates as template argument) of templates.
   *
   * Usually, it is not possible to store e.g. a std::tuple of template classes, i.e.
   *
   * \code
   * typedef std::tuple< std::vector<int>, std::list<char> > Type;
   * \endcode
   *
   * is possible, but
   * \code
   * typedef std::tuple< std::vector, std::list > Type;
   * \endcode
   *
   * will fail.
   *
   * To escape this problem, it is possible to wrap this template:
   *
   * \code
   * typedef std::tuple< _t< std::vector >, _t< std::list > > Type;
   * \endcode
   *
   * To access the original template, use (e.g. for the first element of the tuple)
   * \code
   * typedef typename std::tuple_element<0, Type >::type FirstElementType;
   * template< class... Args >
   * using MyTemplate = FirstElementType::template _t< Args... >;
   * \endcode
   *
   * Nesting
   * -------
   *
   * Another feature is the usage of nested template: In order to use this feature simply add
   * more than one template argument to the type alias _t.
   *
   * Writing
   *
   * \code
   * typedef _t< std::shared_ptr, std::vector > NestedType;
   * \endcode
   *
   * and
   *
   * \code
   * template< class... Args >
   * using MyTemplate = NestedType::template _t< Args... >;
   * \endcode
   *
   * the following command line will compile
   *
   * \code
   * static_assert( std::is_same< MyTemplate<int>, std::shared::ptr<std::vector<int> > >::value, "Upps" );
   * \endcode
   *
   * \warning Note, that all outer templates have to have one template argument,
   * i.e. the nesting has to be well defined afterwards.
   *
   */
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


}
}

#endif
