#ifndef ASSEMBLER_TRAITS_HH
#define ASSEMBLER_TRAITS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune
{
namespace Fem
{
  // traits for the operator passes
  template< class ModelImp,
            class MatrixContainerImp,
            class AdvectionFluxImp,
            class DiffusionFluxImp,
            class DomainDiscreteFunctionImp,
            class RangeDiscreteFunctionImp = DomainDiscreteFunctionImp,
            class ExtraParameterTupleImp = std::tuple<> >
  struct DefaultAssemblerTraits
  {
    typedef typename DomainDiscreteFunctionImp::GridPartType             GridPartType;
    typedef typename GridPartType::GridType                              GridType;

    typedef ModelImp                                                     ModelType;
    typedef AdvectionFluxImp                                             AdvectionFluxType;
    typedef DiffusionFluxImp                                             DiffusionFluxType;

    static const int polynomialOrder = (DomainDiscreteFunctionImp::DiscreteFunctionSpaceType::polynomialOrder==-1) ?
                                        0 : DomainDiscreteFunctionImp::DiscreteFunctionSpaceType::polynomialOrder;

    typedef DomainDiscreteFunctionImp                                    DomainDiscreteFunctionType;
    typedef RangeDiscreteFunctionImp                                     RangeDiscreteFunctionType;
    typedef MatrixContainerImp                                           MatrixContainerType;

    typedef Fem::CachingQuadrature< GridPartType, 0 >                    VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 >                    FaceQuadratureType;

    typedef ExtraParameterTupleImp                                       ExtraParameterTupleType;
  };

  //TODO more general structure for use with multitypeblockmatrix
  template< class DiscreteFunctionTupleImp, template<class,class > class MatrixContainerImp >
  struct AssemblerTraitsList
  {
    static const int size = std::tuple_size< DiscreteFunctionTupleImp >::value;

    template< std::size_t i, std::size_t j >
    using DomainDiscreteFunction = std::tuple_element_t< j, DiscreteFunctionTupleImp >;

    template< std::size_t i, std::size_t j >
    using RangeDiscreteFunction = std::tuple_element_t< i, DiscreteFunctionTupleImp >;

    template< std::size_t i, std::size_t j >
    using Matrix = MatrixContainerImp< typename DomainDiscreteFunction<i,j>::DiscreteFunctionSpaceType,
                                       typename RangeDiscreteFunction<i,j>::DiscreteFunctionSpaceType >;

    template< class DomainDF, class RangeDF >
    using MatrixContainer = MatrixContainerImp< typename DomainDF::DiscreteFunctionSpaceType, typename RangeDF::DiscreteFunctionSpaceType >;

    ////todo improve
    //template< std::size_t i, std::size_t j >
    //using MatrixContainerType = SparseRowLinearOperator< DomainDiscreteFunctionType<i>, RangeDiscreteFunctionType<j> >

    ////todo improve
    //template< std::size_t i, std::size_t j >
    //using OperatorType = std::tuple_element_t< i, OperatorTupleImp >;

#if HAVE_ISTL
    ////todo allow multitype block matrix from istl!
    //
    //typedef MultiTypeBlockVector<Matrix<FieldMatrix<double,3,3> >, Matrix<FieldMatrix<double,3,1> > > RowType0;
    //typedef MultiTypeBlockVector<Matrix<FieldMatrix<double,1,3> >, Matrix<FieldMatrix<double,1,1> > > RowType1;

    //std::integral_constant<int, 0> _0;
    //std::integral_constant<int, 1> _1;

    //MultiTypeBlockMatrix<RowType0,RowType1> multiMatrix;
#endif
  };

}
}

#endif
