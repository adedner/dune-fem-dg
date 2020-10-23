#ifndef DUNE_FEMDG_SUB_CONTAINER_HH
#define DUNE_FEMDG_SUB_CONTAINER_HH

#include <memory>
#include <dune/fem-dg/misc/tupleutility.hh>
#include <dune/fem-dg/misc/uniquefunctionname.hh>

namespace Dune
{
namespace Fem
{

  /** \addtogroup Container
   * @{
   *
   * This short essay explains, why a general container class is needed.
   *
   * Why Algorithms need Containers
   * ------------------------------
   *
   * What is a container class?
   *
   * Each algorithm (see \ref Algorithms) has got some data.
   * This data can be collected in a class, we will call from now on
   * a container class.
   * Container classes for sub algorithms are simply called 'container';
   * container classes for algorithms are called 'global container'. The latter
   * one can be seen as simple tuple of containers.
   *
   * A major problem is that data from a container is a highly (sub) algorithm dependent
   * data structure because only the (sub) algorithm itself knows which data is needed.
   *
   * As long as the data has not to be shared between different (sub) algorithms, everything
   * works fine and a different container types between (sub) algorithm are not relevant.
   *
   * In order to share data in a generic way, we need the knowledge of the data (and its structure)
   * of the (sub) algorithm on a global level. Therefore, it is beneficial to create this data
   * on a global level and pass it to the constructor. Non-shared data can (of course)
   * be constructed inside the algorithm.
   *
   * A draft of such an algorithm could be:
   *
   * \code
   * struct MySubAlgorithm
   * {
   *   typedef MySubAlgorithmContainerType ContainerType;
   *
   *   MySubAlgorithm( const ContainerType& container )
   *    : container_( container )
   *   {}
   *
   * private:
   *   const ContainerType2& container_;
   * };
   * \endcode
   *
   * This is cool, but not generic, since we always have to construct a container
   * on a global level of type `MySubAlgorithm::ContainerType`.
   *
   * In order to circumvent this problem, we can use templated constructors.
   *
   * So, an improved draft of an algorithm could be:
   *
   * \code
   * struct MySubAlgorithm2
   * {
   *   typedef MySubAlgorithmContainerType ContainerType;
   *
   *   template< class ContainerImp >
   *   MySubAlgorithm2( const ContainerImp& container )
   *    : container_( container )
   *   {}
   *
   * private:
   *   const ContainerType& container_;
   * };
   * \endcode
   *
   *
   * But there are still two issues:
   *
   * The first issue is that `ContainerImp` has to be convertable to
   * `ContainerType`.
   *
   * And the second issue is: Are we really interested in the container class itself
   * or are we interested in the data contained therein? In particular, we are only interested
   * in the data that is stored in the container class.
   *
   * The final draft is:
   *
   * \code
   * struct MySubAlgorithm3
   * {
   *   typedef MySubAlgorithmContainerType       ContainerType;
   *   typedef typename ContainerType::DataType1 DataType1;
   *   typedef typename ContainerType::DataType2 DataType2;
   *   typedef typename ContainerType::DataType3 DataType3;
   *
   *
   *   template< class ContainerImp >
   *   MySubAlgorithm3( const ContainerImp& container )
   *   : data1_( container.data1() ),
   *     data2_( container.data2() ),
   *     data3_( container.data3() )
   *   {}
   *
   * private:
   *   std::shared_ptr< DataType1 > data1_;
   *   std::shared_ptr< DataType2 > data2_;
   *   std::shared_ptr< DataType3 > data3_;
   *   //Note: do not store container directly!
   * };
   * \endcode
   *
   * Here, we have used `std::shared_ptr`.
   * Of course, the methods `data1()`, `data2()` and `data3()` have
   * to return a `std::shared_ptr` of the corresponding data type.
   *
   * In a nutshell:
   * * Containers contain data the algorithm needs.
   * * We do not want to store the container itself inside the algorithm.
   * * Constructors of algorithms has got a templated container argument.
   * * The user has to ensure that the container contains the correct data.
   *
   * The last item yields to the question on how to ensure receive the correct data.
   *
   * Our solution here is to provide a general container base class
   * (see OneArgContainer and TwoArgContainer).
   *
   * In general, we cannot totally prevent the user to shoot oneself in his foot
   * and receive nasty compiler errors.
   *
   * But, we can make life a little bit more simple and stick to the following design goals:
   * - provide as much flexibility as necessary,
   * - keep to simple rules.
   *
   * Furthermore, it would be helpful to
   * - reduce the length of types (i.e. avoid redundant information:
   *   use templates instead of full classes),
   * - throw helpful static assertions whenever possible.
   *
   * The price we have to pay is 'some' template code under the hood.
   *
   * Structure of Containers
   * ----------------------
   *
   * Dealing with the solution of partial differential equations
   * and systems of partial differential equations, each system can be described
   * as a general operator \f$ \mathcal{L} \f$ which maps a tuple of data
   * to another tuple of data.
   *
   * Without loss of generality let us suppose that \f$ \mathcal{L} \f$ is a discrete
   * space operator and the 'tuple of data' are discrete functions.
   *
   * This mapping can be written as a system, which is visualized in the following
   * picture.
   *
   * ![A general structure of systems](container1.png)
   *
   * Let us first collect some observations
   * - There are blue squares containing one argument, and some squares contain two arguments.
   * - Without its arguments each blue square corresponds to a template in C++: This means, that
   *   we have got the same unique template for many squares.
   * - It makes sense to store all elements containing one argument in one container
   *   (`OneArgContainer`) and all elements containing two arguments in another container
   *   (`TwoArgContainer`).
   * - Different elements may have the same arguments: We have to wrap them in another
   *   class (red classes)
   * - The arguments \f$ V_0,\ldots,V_n\f$ and \f$W_0,\ldots,W_m\f$ are repeated a dozen of
   *   times: We only want them to store once!
   *
   * The final goal is to collect all square in the general `OneArgContainer` and `TwoArgContainer`,
   * respectively.
   *
   * This is already an interesting design, but we want to go a little
   * bit further and allow each single element (each blue square) to
   * be individually described, i.e. of different type.
   *
   * ![A more general structure of systems](container2.png)
   *
   * Got it? We want define all templates in a individual, short way.
   *
   * The class `TwoArgContainer` needs four templates:
   * - one for all templates with two arguments,
   * - one for all templates with one arguments (`TwoArgContainer` inherits from `OneArgContainer`),
   * - one containing a std::tuple of all row arguments,
   * - one containing a std::tuple of all column arguments.
   *
   * Similar, the `OneArgContainer` class needs only two templates:
   * - one for all templates with one arguments,
   * - one containing a std::tuple of all arguments.
   *
   * We will explain later on, how to really create a struct which can store
   * several templates.
   *
   * Next, we will mention all basic methods each container should have.
   *
   *
   * Basic methods of Containers
   * ---------------------------
   *
   * The `OneArgContainer` and the `TwoArgContainer` class both overload the bracket operator
   * with two different versions:
   * - element access via a `std::integral_constant<unsigned long int,i>` argument,
   * - creation of a sub container containing only specified rows and columns of the original
   *   container. This is specified via a `std::tuple` of `std::integral_constant<unsigned long int,i>`
   *
   * Of course, `OneArgContainer` needs one argument and `TwoArgContainer` needs two arguments
   * (one for row(s) and one for column(s)).
   *
   *
   * \note For `std::integral_constant<unsigned long int,i>` abbreviations are provided:
   * Let `i` be an integer number (smaller than 50), then type is called `__i` and the
   * corresponding variable `_i`.
   *
   * The functionality of the bracket operator is shown in the following example.
   *
   * ![How to create sub containers and to access elements](container3.png)
   *
   * Once again, note that `TwoArgContainer` inherits from `OneArgContainer`.
   *
   * Creation of Container
   * ---------------------
   *
   * There are different ways to create the first template argument of the `OneArgContainer`:
   * - For unique templates:
   *   + you can directly use a template,
   *   + use `ArgContainerArgWrapperUnique`: plug the template therein.
   * - For non unique templates:
   *   + use `ArgContainerArgWrapper`: define a `std::tuple< _t<...>... >`
   *   and plug this struct therein.
   *
   * The struct `_t` is a type alias and allows to store a template inside a `std::tuple`.
   *
   * The following figure will explain the creation of the container more detailed:
   *
   * ![How to create the first template argument of the OneArgContainer](container5.png)
   *
   * For the class `TwoArgContainer` we use the similar structure with the difference
   * that for non unique templates a `std::tuple<std::tuple< _t<...>... >...>` is needed.
   *
   * ![How to create the first template argument of the TwoArgContainer](container4.png)
   *
   * That's it!
   *
   * Global Container and Combined Global Containers
   * -----------------------------------------------
   *
   * The classes `GlobalContainer` and `CombinedGlobalContainer` are helper classes
   * to easily create a global container combining all (sub) containers.
   *
   * These classes only contains a `sub<i>()` method which extracts the `i`s
   * container.
   *
   * Here are two pictures visualizing both classes
   *
   * ![CombinedContainer](container7.png)
   *
   * ![GlobalCombinedContainer](container6.png)
   *
   *
  ** @}*/



  /**
   * \brief Wrapper class to create a discrete function depending on grid, grid part or space.
   *
   * \note This is a convenience class. Do we really need it?
   *
   * \ingroup Container
   */
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
  template< template<class,class,int...> class, template<class,int...> class, template<class,int...> class, class Arg, class >
  class TwoArgContainer
  {
    static_assert( static_fail< Arg >::value, "fourth/fifth argument has to be a tuple" );
  };

  /**
   * \brief Base class for all container classes taking one argument.
   *
   * \ingroup Container
   *
   * Before start reading, you should have an overview: \ref Container.
   *
   * ![How to create sub containers and to access elements](container3_1.png)
   *
   * There are different ways to create the first template argument.
   *
   * ![How to create the first template argument of the OneArgContainer](container5.png)
   *
   *
   */
  template< template< class, int... > class OneArgImp, class... Args >
  struct OneArgContainer< OneArgImp, std::tuple< Args... > >
  {
    typedef std::tuple< Args... >                ArgTupleType;
    typedef typename VectorPacker< OneArgImp, ArgTupleType >::shared_type Item1TupleType;

  public:

    typedef OneArgContainerStore< OneArgImp > OneArgType;

    template< unsigned long int i >
    using Item1 = typename std::tuple_element_t< i, Item1TupleType>::element_type;

    static const int numArgs = std::tuple_size< ArgTupleType >::value;

  protected:
    static const int size = std::tuple_size< Item1TupleType >::value;
    static_assert( size >= 0, "invalid integer_sequence: Throw this additional assertion here because \
                               gcc-6 won't stop compiling for a very long time... :(" );

    typedef std::make_integer_sequence< unsigned long int, size > SequenceType;

    ////// Creation
    template< unsigned long int i, class SameObject >
    static decltype(auto) createItem1( std::shared_ptr< SameObject > obj, const std::string name )
    {
      std::cout << "###CREATE: item1 ('" << name << "'): " << print( _index<i>() ) << std::endl;
      return std::make_shared<Item1<i> >( obj, name );
    }
    template< unsigned long int ...i, class SameObject>
    static decltype(auto) createContainer( _indices<i...>, std::shared_ptr< SameObject > obj, const std::string name )
    {
      return std::make_tuple( createItem1<i>( obj, name )... );
    }

    ////// Copy
    template< class Item1Tuple, unsigned long int ...i >
    static decltype(auto) copyContainer( const Item1Tuple& item1, std::tuple< _index<i>...  > )
    {
      std::cout << "###COPY: item1 ('-'): " << print( _index<i>()... ) << std::endl;
      return std::make_tuple( std::get<i>( item1 )... );
    }

    //Be your own (template) friend
    template< template<class,class,int...> class, template<class,int...> class, template<class,int...> class, class, class >
    friend class TwoArgContainer;

    template< template<class,int...> class, class >
    friend struct OneArgContainer;
  public:
    /**
     * \brief constructor, for internal use only.
     */
    OneArgContainer( const Item1TupleType& item )
    : item1_( item )
    {}

    /**
     * \brief constructor
     *
     * \param obj A common object to create all container elements.
     * \param name just a simple name for this container.
     */
    template< class SameObject >
    OneArgContainer( std::shared_ptr< SameObject > obj, const std::string name = "" )
    : item1_( createContainer( SequenceType(), obj, name ) )
    {}

    /**
     * \brief Allows access to the i's element.
     *
     * \param index the i's element we want to access
     */
    template< unsigned long int i >
    decltype(auto) operator() ( _index<i> index ) const
    {
      //const auto& res = std::get<i>( item1_ );
      //std::cout << "###ACCESS: item1 ('" << res->solution()->name() << "')" << print( index ) << std::endl;
      return std::get<i>( item1_ );
    }

    /**
     * \brief This more advanced method allows to create sub containers out of an already existing container.
     *
     * \param index std::tuple of elements which should be contained in the new sub container.
     *
     * \warning Just for the case, that the sub container only contains one index:
     * Do not forget to wrap a std::tuple<> around the std::integral_constant<>.
     * Otherwise you may accidently call the wrong method...
     */
    template< unsigned long int... i >
    decltype(auto) operator() ( std::tuple< _index<i>... > index ) const
    {
      typedef std::tuple< std::tuple_element_t<i,ArgTupleType>...> NewArgTupleType;

      typedef std::tuple< _index<i>... > Maps;

      typedef MappedOneArgContainer< OneArgImp, Maps >            MappedOneArgImp;

      typedef OneArgContainer< MappedOneArgImp::template _t1, NewArgTupleType > SubContainerType;

      return std::make_shared< SubContainerType >( copyContainer(item1_, index ) );
    }
  protected:
    Item1TupleType item1_;
  };



  /**
   * \brief Base class for all container classes taking two argument.
   *
   * \ingroup Container
   *
   * Before start reading, you should have an overview: \ref Container.
   *
   * ![How to create sub containers and to access elements](container3_2.png)
   *
   * There are different ways to create the first template argument.
   * The second argument is similarly constructed, see `OneArgContainer`.
   *
   * ![How to create the first template argument of the TwoArgContainer](container4.png)
   *
   *
   */
  template< template<class,class,int...> class TwoArgImp,
            template<class,int...> class RowOneArgImp,
            template<class,int...> class ColOneArgImp,
            class... RowArgs,
            class... ColArgs >
  class TwoArgContainer< TwoArgImp, RowOneArgImp, ColOneArgImp, std::tuple< RowArgs... >, std::tuple< ColArgs... > >
    : public OneArgContainer< ColOneArgImp, std::tuple< RowArgs... > >
  {
    typedef OneArgContainer< RowOneArgImp, std::tuple< RowArgs... > > BaseType;
    typedef OneArgContainer< ColOneArgImp, std::tuple< ColArgs... > > FakeColBaseType;

    typedef std::tuple< RowArgs... >        RowArgTupleType;
    typedef std::tuple< ColArgs... >        ColArgTupleType;

    typedef typename MatrixPacker< TwoArgImp, RowArgTupleType, ColArgTupleType >::shared_type  Item2TupleType;
    typedef typename VectorPacker< RowOneArgImp, RowArgTupleType >::shared_type                Item1TupleType;

    typedef Item1TupleType                                                                     RowItem1TupleType;
    typedef typename VectorPacker< ColOneArgImp, ColArgTupleType >::shared_type                ColItem1TupleType;
  public:

    typedef OneArgContainerStore< RowOneArgImp > RowOneArgType;
    typedef OneArgContainerStore< ColOneArgImp > ColOneArgType;
    typedef TwoArgContainerStore< TwoArgImp >    TwoArgType;

    template< unsigned long int i, unsigned long int j >
    using Item2 = typename std::tuple_element_t< j, std::tuple_element_t< i, Item2TupleType> >::element_type;

    using BaseType::operator();

    static const int numRowArgs = std::tuple_size< RowArgTupleType >::value;
    static const int numColArgs = std::tuple_size< ColArgTupleType >::value;

  protected:
    using BaseType::item1_;
    typedef typename BaseType::SequenceType SequenceType;

    static const int size = BaseType::size;

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
    decltype(auto) copyContainerRow( std::tuple< _index<j>... > ) const
    {
      std::cout << "###COPY: item2 ('-'): " << print( std::make_tuple( _index<j>()... ) ) << std::endl;
      return std::make_tuple( std::get<j>( std::get<i>( item2_ ) )... );
    }
    template< unsigned long int ...i, unsigned long int ...j >
    decltype(auto) copyContainer( std::tuple< _index<i>... > row, std::tuple< _index<j>... > col ) const
    {
      return std::make_tuple( copyContainerRow<i>( col )... );
    }

    //Be your own (template) friend
    template< template<class,class,int...> class, template<class,int...> class, template<class,int...> class, class, class >
    friend class TwoArgContainer;

    template< template<class,int...> class, class >
    friend struct OneArgContainer;

  public:
    /**
     * \brief constructor, for internal use only.
     */
    TwoArgContainer( const RowItem1TupleType& rowItem1,
                     const ColItem1TupleType& colItem1,
                     const Item2TupleType& item2 )
    : BaseType( rowItem1 ),
      item2_( item2 ),
      colItem1_( colItem1 )
    {}

    /**
     * \brief constructor.
     *
     * \param obj A common object to create all container elements.
     * \param name just a simple name for this container.
     */
    template< class SameObject >
    TwoArgContainer( std::shared_ptr< SameObject > obj, const std::string name = "" )
    : BaseType( obj, name ),
      item2_( createContainer( SequenceType(), SequenceType(), name + "matrix" ) ),
      colItem1_( FakeColBaseType::createContainer( SequenceType(), obj, name ) )
    {}


    /**
     * \brief Allows access to the element in the i's row and j's column.
     *
     * \param row the i's row we want to access
     * \param col the j's col we want to access
     */
    template< unsigned long int i, unsigned long int j >
    decltype(auto) operator() ( _index<i> row, _index<j> col ) const
    {
      const auto& res = std::get<j>( std::get<i>( item2_ ) );
      //std::cout << "###ACCESS: item2 ('-')" << print( row, col ) << std::endl;
      return res;
    }

    /**
     * \brief This more advanced method allows to create sub containers out of an already existing container.
     *
     * \param row std::tuple of rows which should be contained in the new sub container.
     * \param col std::tuple of columns which should be contained in the new sub container.
     *
     * \warning Just for the case, that the sub container only contains one row or one line:
     * Do not forget to wrap a std::tuple<> around the std::integral_constant<>.
     * Otherwise you may accidently call the wrong method...
     */
    template< unsigned long int... i, unsigned long int... j >
    decltype(auto) operator() ( std::tuple< _index<i>... > row, std::tuple< _index<j>... > col ) const
    {
      typedef std::tuple< std::tuple_element_t<i,RowArgTupleType>...>   NewRowArgTupleType;
      typedef std::tuple< std::tuple_element_t<j,ColArgTupleType>...>   NewColArgTupleType;

      typedef std::tuple< _index<i>... >                                RowMaps;
      typedef std::tuple< _index<j>... >                                ColMaps;

      typedef MappedOneArgContainer< RowOneArgImp, RowMaps >            MappedRowOneArgImp;
      typedef MappedOneArgContainer< ColOneArgImp, ColMaps >            MappedColOneArgImp;

      typedef MappedTwoArgContainer< TwoArgImp, RowMaps, ColMaps >      MappedTwoArgImp;

      typedef TwoArgContainer< MappedTwoArgImp::template _t2,
                               MappedRowOneArgImp::template _t1,
                               MappedColOneArgImp::template _t1,
                               NewRowArgTupleType, NewColArgTupleType>  SubContainerType;

      std::cout << "###CREATE: local sub container " << print( row, col ) << std::endl;
      return std::make_shared< SubContainerType >( BaseType::copyContainer( item1_, row ),
                                                   FakeColBaseType::copyContainer( colItem1_, col ),
                                                   copyContainer( row, col ) );
    }
  protected:
    Item2TupleType          item2_;
    ColItem1TupleType       colItem1_;
  };





  template< class Arg >
  struct ExtraArgContainer
  {
    static_assert( static_fail<Arg>::value, "wrong template args." );
  };

  /**
   * \brief Base class for all container classes taking one argument.
   *
   * \ingroup Container
   *
   * Before start reading, you should have an overview: \ref Container.
   */
  template< class... Args >
  struct ExtraArgContainer< std::tuple< Args... > >
  {
    typedef std::tuple< Args... >                    ArgTupleType;
  public:
    template< unsigned long int i >
    using Item = std::tuple_element_t<i,ArgTupleType>;

    static const int size = std::tuple_size< ArgTupleType >::value;

  protected:
    ////// Copy
    template< class ItemTuple, unsigned long int ...i >
    static decltype(auto) copyContainer( const ItemTuple& item, std::tuple< _index<i>...  > )
    {
      return std::make_tuple( std::get<i>( item )... );
    }

    //template< class... >
    //friend class ExtraArgContainer;
  public:
    /**
     * \brief constructor.
     */
    ExtraArgContainer( const ArgTupleType& item )
    : item_( item )
    {}

    /**
     * \brief Allows access to the i's element.
     *
     * \param index the i's element we want to access
     */
    template< unsigned long int i >
    decltype(auto) operator() ( _index<i> index ) const
    {
      static_assert( i<size, "index element does not exist!" );
      return std::get<i>( item_ );
    }

    /**
     * \brief This more advanced method allows to create sub containers out of an already existing container.
     *
     * \param index std::tuple of elements which should be contained in the new sub container.
     *
     * \warning Just for the case, that the sub container only contains one index:
     * Do not forget to wrap a std::tuple<> around the std::integral_constant<>.
     * Otherwise you may accidently call the wrong method...
     */
    template< unsigned long int... i >
    decltype(auto) operator() ( std::tuple< _index<i>... > index ) const
    {
      typedef ExtraArgContainer< std::tuple< Item<i>... > > SubExtraContainerType;

      return std::make_shared< SubExtraContainerType >( copyContainer(item_, index ) );
    }
  protected:
    ArgTupleType item_;
  };





  /**
   * \brief Base class for all container classes taking one argument.
   *
   * \ingroup Container
   *
   * Before start reading, you should have an overview: \ref Container.
   */
  template< class... Args >
  struct ExtraArg
  {
    typedef std::tuple< Args... >                    ArgTupleType;
  public:
    template< unsigned long int i >
    using Access = std::tuple_element_t< i, ArgTupleType>;

  protected:
  public:
    static const int size = std::tuple_size< ArgTupleType >::value;
    static_assert( size >= 0, "invalid integer_sequence: Throw this additional assertion here because \
                               gcc-6 won't stop compiling for a very long time... :(" );
    typedef std::make_integer_sequence< unsigned long int, size > SequenceType;

    ////// Creation
    template< unsigned long int i, class ContainerImp >
    static decltype(auto) createItem( const ContainerImp& cont )
    {
      //typedef Access<i> AccessType;
      //return value is/should be shared_ptr...
      return Access<i>::apply( cont );
    }
    template< unsigned long int ...i, class ContainerImp >
    static decltype(auto) createContainer( _indices<i...>, const ContainerImp& cont )
    {
      return std::make_tuple( createItem<i>( cont )... );
    }

    template< class... >
    friend struct ExtraArg;
  public:
    //no default constructor
    ExtraArg() = delete;

    template< class ContainerImp >
    static decltype(auto) init( const ContainerImp& cont )
    {
      typedef decltype(std::declval<ExtraArg<Args...> >().createContainer( SequenceType(), cont ) ) ReturnType;

      return std::make_shared<ExtraArgContainer<ReturnType> >( createContainer( SequenceType(), cont ) );
    }
  };






  /**
   * \brief A placeholder for empty item1 or item2 container.
   *
   * \ingroup Container
   */
  template< class... >
  struct EmptyContainerItem
  {
    template< class ... Args>
    EmptyContainerItem( Args&&... args )
    {}
  };

  /**
   * \brief Container only storing the solution.
   *
   * \ingroup Container
   */
  template <class DiscreteFunctionImp >
  struct SolutionContainerItem
  {
    using CItem = ContainerItem< DiscreteFunctionImp >;
  public:
    using DiscreteFunction = DiscreteFunctionImp;

    // owning container
    template< class SameObject >
    SolutionContainerItem( const std::shared_ptr<SameObject>& obj, const std::string name = "" )
    : stringId_( FunctionIDGenerator::instance().nextId() ),
      solution_(      std::make_shared< CItem >( name + "u" + stringId_, obj ) )
    {}

    //solution
    std::shared_ptr< DiscreteFunction > solution() const
    {
      return solution_->shared();
    }
  private:
    const std::string        stringId_;
    std::shared_ptr< CItem > solution_;
  };




  //Simple Select classes
  //----------------------

  template< class Id >
  struct SolutionSelect
  {
    template< class ContainerImp >
    static decltype(auto) apply( const std::shared_ptr<ContainerImp>& cont )
    {
      return (*cont)( Id() )->solution();
    }
  };

  template< class Id >
  struct ExactSolutionSelect
  {
    template< class ContainerImp >
    static decltype(auto) apply( const std::shared_ptr<ContainerImp>& cont )
    {
      return (*cont)( Id() )->exactSolution();
    }
  };

  template< class RowId, class ColId >
  struct MatrixSelect
  {
    template< class ContainerImp >
    static decltype(auto) apply( const std::shared_ptr<ContainerImp>& cont )
    {
      return (*cont)( RowId(), ColId() )->matrix();
    }
  };


  //For IODataTuple
  template< class Id >
  struct IOSolutionSelect
  {
    template< class ContainerImp >
    static decltype(auto) apply( const std::shared_ptr<ContainerImp>& cont )
    {
      return (*cont)( Id() )->solution().get();
    }
  };
  template< class Id >
  struct IOExactSolutionSelect
  {
    template< class ContainerImp >
    static decltype(auto) apply( const std::shared_ptr<ContainerImp>& cont )
    {
      return (*cont)( Id() )->exactSolution().get();
    }
  };

  //Container Access
  //----------------------

  template< template<class> class Select, class... >
  struct OneArgAccess
  {
    static_assert( static_fail_t< Select >::value, "failed" );
  };

  //local version
  template< template<class> class Select, unsigned long int local >
  struct OneArgAccess< Select, _index<local> >
  {
    template< class ContainerImp >
    static decltype(auto) apply( const ContainerImp& cont )
    {
      return Select< _index<local> >::apply( cont );
    }
  };

  //global version
  template< template<class> class Select, unsigned long int global, unsigned long int local >
  struct OneArgAccess< Select, _index<global>, _index<local> >
  {
    template< class ContainerImp >
    static decltype(auto) apply( const ContainerImp& cont )
    {
      return Select< _index<local> >::apply( cont.sub( _index<global>() ) );
    }
  };

  template< template<class> class Select, class... Args >
  using _e = OneArgAccess<Select,Args...>;

  template< class... Args >
  using _ee = OneArgAccess<SolutionSelect,Args...>;

  template< template<class,class> class Select, class... >
  struct TwoArgAccess;

  //local version
  template< template<class,class> class Select, unsigned long int row, unsigned long int col >
  struct TwoArgAccess< Select, _index<row>, _index<col> >
  {
    template< class ContainerImp >
    static decltype(auto) apply( const ContainerImp& cont )
    {
      return Select< _index<row>, _index<col> >::apply( cont );
    }
  };

  //global version
  template< template<class,class> class Select, unsigned long int global, unsigned long int row, unsigned long int col >
  struct TwoArgAccess< Select, _index<global>, _index<row>, _index<col> >
  {
    template< class ContainerImp >
    static decltype(auto) apply( const ContainerImp& cont )
    {
      return Select< _index<row>, _index<col> >::apply( cont.sub( _index<global>() ) );
    }
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
