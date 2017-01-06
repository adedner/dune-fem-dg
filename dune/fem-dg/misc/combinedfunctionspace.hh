#ifndef DUNE_FEMDG_COMBINEDFUNCTIONSPACE_HH
#define DUNE_FEMDG_COMBINEDFUNCTIONSPACE_HH

#include <dune/fem/space/common/functionspace.hh>
#include <dune/common/fmatrix.hh>
#include <dune/fem/common/utility.hh>

namespace Dune
{
namespace Fem
{

  // Forward declaration of
  // base class for vector valued function spaces.
  template< class... FunctionSpaceImp >
  class CombinedFunctionSpace;

  //! \brief Traits class for vector function spaces
  template< class FunctionSpaceHead, class... FunctionSpaceTail >
  struct CombinedVectorSpaceTraits
  {
    /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainFieldType */
    typedef typename FunctionSpaceHead::DomainField   DomainFieldType;
    /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeFieldType */
    typedef typename FunctionSpaceHead::RangeField    RangeFieldType;


    //! size of domian space
    static const int dimDomain = FunctionSpaceHead::dimDomain;

    //! size of range space
    static constexpr const int dimRange = Std::sum( static_cast< int >( FunctionSpaceHead::dimRange ),
                                                    static_cast< int >( FunctionSpaceTail::dimRange ) ...  );

    /** \copydoc Dune::Fem::FunctionSpaceInterface::DomainType */
    typedef FieldVector< DomainFieldType, dimDomain >   DomainType;

    /** \copydoc Dune::Fem::FunctionSpaceInterface::RangeType */
    typedef FieldVector< RangeFieldType, dimRange>      RangeType;

    /** \brief linear mapping type */
    typedef FieldMatrix< RangeFieldType, dimRange, dimDomain > LinearMappingType;

    /** \brief scalar function space type */
    typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1>
      ScalarFunctionSpaceType;

  };

  /** @ingroup CombinedFunctionSpace
      \brief A vector valued function space.

      FunctionSpace defines what the types of the domain vector
      space and the range vector space for a function are.
  */
  template< class... FunctionSpaceImp >
  class CombinedFunctionSpace
  : public FunctionSpaceInterface< CombinedVectorSpaceTraits< FunctionSpaceImp... > >
  {
    typedef CombinedFunctionSpace< FunctionSpaceImp... > ThisType;
    typedef std::tuple< FunctionSpaceImp... >            FunctionSpaceTupleType;

  public:
    template< int i >
    struct SubFunctionSpace
    {
      typedef typename std::tuple_element_t< i, FunctionSpaceTupleType >::DomainFieldType         DomainFieldType;
      typedef typename std::tuple_element_t< i, FunctionSpaceTupleType >::RangeFieldType          RangeFieldType;
      static const int dimDomain = std::tuple_element_t< i, FunctionSpaceTupleType >::dimDomain;
      static const int dimRange = std::tuple_element_t< i, FunctionSpaceTupleType >::dimRange;
      typedef typename std::tuple_element_t< i, FunctionSpaceTupleType >::DomainType              DomainType;
      typedef typename std::tuple_element_t< i, FunctionSpaceTupleType >::RangeType               RangeType;
      typedef typename std::tuple_element_t< i, FunctionSpaceTupleType >::JacobianRangeType       JacobianRangeType;
      typedef typename std::tuple_element_t< i, FunctionSpaceTupleType >::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef typename std::tuple_element_t< i, FunctionSpaceTupleType >::HessianRangeType        HessianRangeType;
    };
    typedef ThisType FunctionSpaceType;
  };
} // namespace Fem

} // namespace Dune

#endif
