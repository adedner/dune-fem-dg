#ifndef DUNE_FEMDG_COUPLEDALGORITHMS_HH
#define DUNE_FEMDG_COUPLEDALGORITHMS_HH

#include <tuple>
#include <type_traits>
#include <memory>
#include <utility>
#include <iostream>
#include <dune/fem-dg/misc/integral_constant.hh>

namespace Dune
{
namespace Fem
{

  /**
   *  \brief Creates a tuple of uncoupled sub algorithms
   *
   *  \tparam SubAlgorithmImp
   */
  template< class... SubAlgorithmsImp >
  class CreateSubAlgorithms
  {
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

    static_assert( std::tuple_size<SubAlgorithmTupleType>::value > 0, "We need at least one Sub-Algorithm" );

    typedef typename std::tuple_element<0,SubAlgorithmTupleType>::type::element_type::GridType GridType;

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
      std::shared_ptr< Container<i> > container = std::make_shared< Container<i> >( grid );
      ElementPtr<i> ptr = std::make_shared<Element<i> >( grid, container );

      auto sol = (*container)(_0)->solution();
      std::cout << "new container: name: " << sol->name() << ", size: " << sol->size() << std::endl;
      std::cout << "new container: name: " << ptr->solution().name() << ", size: " << ptr->solution().size() << std::endl;
      return ptr;
    }


    template< unsigned long int i, class GlobalContainerImp  >
    static decltype(auto) createSubAlgorithm( std::shared_ptr< GlobalContainerImp > cont )
    {
      static _index<i> idx;
      return std::make_shared<Element<i> >( cont->sub( idx ) );
    }


    template< std::size_t ...i >
    static SubAlgorithmTupleType apply ( std::index_sequence< i... >, GridType &grid )
    {
      return std::make_tuple( createSubAlgorithm<i>( grid )... );
    }

    template< unsigned long int ...i, class GlobalContainerImp  >
    static SubAlgorithmTupleType apply ( std::index_sequence< i... >, std::shared_ptr< GlobalContainerImp > cont )
    {
      return std::make_tuple( createSubAlgorithm<i>( cont )... );
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
    static SubAlgorithmTupleType apply ( std::shared_ptr< GlobalContainerImp > cont  )
    {
      return apply( typename std::make_index_sequence< std::tuple_size< SubAlgorithmTupleType >::value >(), cont );
    }
  };



} // namespace Fem
} // namespace Dune

#endif
