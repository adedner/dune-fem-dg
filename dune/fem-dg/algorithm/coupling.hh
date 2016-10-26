#ifndef DUNE_FEMDG_COUPLEDALGORITHMS_HH
#define DUNE_FEMDG_COUPLEDALGORITHMS_HH

#include <tuple>
#include <type_traits>
#include <memory>
#include <utility>

namespace Dune
{
namespace Fem
{

  /**
   *  \brief Creates a tuple of uncoupled sub algorithms
   *
   *  \note This is the default implementation which should be used
   *  for uncoupled algorithms.
   *
   *  \tparam GridImp
   *  \tparam SubAlgorithmImp
   */
  template< class GridImp, class... SubAlgorithmsImp >
  class UncoupledSubAlgorithms
  {
    typedef GridImp               GridType;

    template <class Object>
    class AlgDeleter
    {
    public:
      typedef Object ObjectType;

      //! \brief empty constructor.
      AlgDeleter() : obj_( nullptr )
      {}

      //! \brief Constructor taking the object that needs to deleted.
      explicit AlgDeleter( Object* obj )
        : obj_( obj )
      {}

      //! \brief Delete an object and the additional one.
      template<class Ptr>
      void operator()(Ptr* ptr)
      {
          delete ptr;
          delete obj_;
      }
    protected:
      Object* obj_;
    };
  public:
    typedef std::tuple< std::shared_ptr< SubAlgorithmsImp > ... >    SubAlgorithmTupleType;

    template< int i >
    using ElementPtr = typename std::tuple_element< i, SubAlgorithmTupleType >::type;

    template< int i >
    using Element = typename ElementPtr<i>::element_type;

    template< int i >
    using Deleter = AlgDeleter< typename Element<i>::ContainerType >;

    template< int i >
    using Container = typename Deleter<i>::ObjectType;

  private:

    template< std::size_t i >
    static ElementPtr<i> createSubAlgorithm( GridType& grid )
    {
      Container<i>* container = new Container<i>( grid );

      //std::tuple<int,double,std::string> t{1,-0.5,"tuple"};
      //Dune::Fem::for_each(container,
      //                    [](const auto& entry, auto I)
      //                    {
      //                      std::cout << "new container: name: " << entry->template solution<I>()->name() << ", size: " << entry->template solution<I>()->size() << std::endl;
      //                    });
      //std::cout << "new container: name: " << sol.name() << ", size: " << sol.size() << std::endl;

      //ElementPtr<i> ptr( new Element<i>( grid, *container ), Deleter<i>( container ) );

      //grids = ElementPtr<i>::initializeGrid();
      //gridParts = ElementPtr<i>::initializeGridParts();
      //spaces = ElementPtr<i>::initializeSpaces();

      std::shared_ptr< Container<i> > container = std::make_shared< Container<i> >( grid );
      ElementPtr<i> ptr = std::make_shared<Element<i> >( grid, container );


      auto sol = (*container)(_0)->solution();
      std::cout << "new container: name: " << sol->name() << ", size: " << sol->size() << std::endl;
      std::cout << "new container: name: " << ptr->solution().name() << ", size: " << ptr->solution().size() << std::endl;
      return ptr;
    }

    template< std::size_t ...i >
    static SubAlgorithmTupleType apply ( std::index_sequence< i... >, GridType &grid )
    {
      return std::make_tuple( createSubAlgorithm<i>( grid )... );
    }

  public:
    /**
     * \brief Creates a tuple of sub algorithms
     *
     * \param grid the grid which is needed to construct the algorithms
     * \return The constructed tuple of algorithms
     */
    static SubAlgorithmTupleType apply ( GridType &grid )
    {
      return apply( typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >(), grid );
    }

    template< class GlobalContainerImp >
    static SubAlgorithmTupleType apply ( GridType &grid, std::shared_ptr< GlobalContainerImp > cont  )
    {
      return apply( typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >(), grid, cont );
    }
  };


//  template< class GridImp, class GlobalContainerImp, class... SubAlgorithmsImp >
//  class ThetaSchemeCoupling
//  {
//    typedef GridImp               GridType;
//
//    template <class Object>
//    class AlgDeleter
//    {
//    public:
//      typedef Object ObjectType;
//
//      //! \brief empty constructor.
//      AlgDeleter() : obj_( nullptr )
//      {}
//
//      //! \brief Constructor taking the object that needs to deleted.
//      explicit AlgDeleter( Object* obj )
//        : obj_( obj )
//      {}
//
//      //! \brief Delete an object and the additional one.
//      template<class Ptr>
//      void operator()(Ptr* ptr)
//      {
//          delete ptr;
//          delete obj_;
//      }
//    protected:
//      Object* obj_;
//    };
//  public:
//    typedef GlobalContainerImp    GlobalContainerType;
//
//    typedef std::tuple< std::shared_ptr< SubAlgorithmsImp > ... >    SubAlgorithmTupleType;
//
//    typedef std::tuple< std::shared_ptr< typename SubAlgorithmsImp::ContainerType > ... >    SubContainerTupleType;
//
//    template< int i >
//    using ElementPtr = typename std::tuple_element< i, SubAlgorithmTupleType >::type;
//
//    template< int i >
//    using Element = typename ElementPtr<i>::element_type;
//
//    template< int i >
//    using Deleter = AlgDeleter< typename Element<i>::ContainerType >;
//
//    template< int i >
//    using Container = typename Deleter<i>::ObjectType;
//
//  private:
//
//    template< std::size_t i >
//    static ElementPtr<i> createSubAlgorithm( GridType& grid )
//    {
//      Container<i>* container = new Container<i>( grid );
//
//      //std::shared_ptr< Container<i> > container = std::make_shared< Container<i> >( grid );
//
//#if 0
//      template< class ... SubAlgorithms >
//      class GlobalContainer
//      {
//        typedef std::tuple< SubAlgorithms::ContainerType ... > ContainerTuple;
//
//
//        template< int i >
//        using Container = typename std::tuple_element<0, ContainerTuple >::type;
//
//        typedef Container< 0 >::DiscreteFunctionType::DiscreteFunctionSpaceType::GridPartType::GridType GridType;
//
//      public:
//        GlobalContainer( GridType& grid )
//        : containers_( std::make_tuple( std::make_unique< SubAlgorithms >( grid, "" ) ) )
//        {}
//
//        template< int i >
//        struct Sub
//        {};
//
//        template<>
//        struct Sub<0>
//        {
//          std::shared_ptr< ContainerType<0> > container() const
//          {
//            return std::make_shared< ContainerType<0> > ( ...);
//          }
//        };
//
//        template< int i >
//        const Cot
//
//
//      private:
//        ContainerTuple& containers_;
//      };
//
//
//
//      container->template solution<0>() = globalContainer->template solution<12>();
//      container->template solution<1>() = globalContainer->template solution<3>();
//#endif
//      auto sol = container->template solution<0>();
//
//      std::cout << "new container: name: " << sol->name() << ", size: " << sol->size() << std::endl;
//      ElementPtr<i> ptr( new Element<i>( grid, *container ), Deleter<i>( container ) );
//      return ptr;
//    }
//
//    template< std::size_t i >
//    static std::shared_ptr< Container<i> > createContainers( GridType& grid )
//    {
//      return std::make_shared< Container<i> >( grid );
//    }
//
//    static SubContainerTupleType createContainerTuple( GridType& grid )
//    {
//      auto c0 = createContainers<0>( grid );
//      auto c1 = createContainers<1>( grid, c0.template sub<0>() );
//      auto c2 = createContainers<2>( grid, c0 );
//      auto tuple = std::make_tuple( c0, c1, c2 );
//
//      return std::make_tuple( createContainers<0>( grid ), createContainers<1>( grid ), createContainers<2>( grid ) );
//    }
//
//    //static std::tuple< std::shared_ptr< > > sharedObjects( GridType& grid )
//    //{
//    //  return std::make_tuple( std::make_shared<GridPartType<0> >( grid ) );
//    //}
//
//    template< std::size_t ...i >
//    static SubAlgorithmTupleType apply ( std::index_sequence< i... >, GridType &grid )
//    {
//      return std::make_tuple( createSubAlgorithm<i>( grid )... );
//    }
//
//  public:
//    /**
//     * \brief Creates a tuple of sub algorithms
//     *
//     * \param grid the grid which is needed to construct the algorithms
//     * \return The constructed tuple of algorithms
//     */
//    static SubAlgorithmTupleType apply ( GridType &grid )
//    {
//      return apply( typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >(), grid );
//    }
//  };




} // namespace Fem
} // namespace Dune

#endif
