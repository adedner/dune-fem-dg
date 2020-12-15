#ifndef DUNE_FEM_DG_ESTIMATORBASE_HH
#define DUNE_FEM_DG_ESTIMATORBASE_HH

//- Dune-fem includes
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/function/localfunction/const.hh>

#include <dune/fem-dg/operator/adaptation/adaptation.hh>

namespace Dune
{
namespace Fem
{
  // EstimatorBase
  // -------------

  /** \class EstimatorBase
   *  \brief estimates and marks grid elements for refinement/coarsening
   *
   *  Base class for estimating and marking grid elements for refinement and
   *  coarsening. The EstimatorBase marks the grid elements by default neither
   *  for refinement nor for coarsening.
   *
   *  \tparam DiscreteFunction Discrete function type
   */
  template< class DiscreteFunction >
  class EstimatorBase
  {
    typedef EstimatorBase< DiscreteFunction >                         ThisType;

  public:
    typedef DiscreteFunction DiscreteFunctionType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
    typedef ConstLocalFunction< DiscreteFunctionType >                  ConstLocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: DomainFieldType     DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType      RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainType          DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeType           RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType   JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType        GridPartType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType        IteratorType;

    typedef typename GridPartType :: GridType                         GridType;
    typedef typename GridPartType :: IndexSetType                     IndexSetType;
    typedef typename GridPartType :: IntersectionIteratorType         IntersectionIteratorType;

    typedef typename GridPartType :: IntersectionType                 IntersectionType;

    typedef typename GridPartType :: template Codim< 0 > :: EntityType        ElementType;
    typedef typename GridType :: template Codim< 0 > :: Entity                GridElementType;
    typedef typename ElementType::Geometry                             GeometryType;
    static const int dimension = GridType :: dimension;

    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >           ElementQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 1 >           FaceQuadratureType;
    typedef Dune::FieldMatrix<double,dimension,dimension>             JacobianInverseType;
    typedef std :: vector< double >                                   IndicatorType;

  protected:
    const DiscreteFunctionSpaceType& dfSpace_;
    GridPartType& gridPart_;
    const IndexSetType& indexSet_;
    GridType &grid_;
    IndicatorType indicator_;

  public:
    //! Constructor
    explicit EstimatorBase ( const DiscreteFunctionSpaceType& space )
    : dfSpace_( space ),
      gridPart_( dfSpace_.gridPart() ),
      indexSet_( gridPart_.indexSet() ),
      grid_( gridPart_.grid() ),
      indicator_() // create empty vector, resize before use
    {
    }

    void clear( IndicatorType& indicator )
    {
      typedef typename IndicatorType :: iterator IteratorType;
      const IteratorType end = indicator.end();
      for( IteratorType it = indicator.begin(); it != end; ++it )
        *it = 0.0;
    }

    virtual void clear()
    {
      clear( indicator_ );
    }


    virtual void estimateLocal( const ElementType& it )
    {
      DUNE_THROW( Dune::NotImplemented, "Use your own local estimator !" );
    }

    virtual void estimate()
    {
      clear( indicator_ );
      for( const auto& en : elements( dfSpace_.gridPart() ) )
        estimateLocal( en );
    }


    /** \brief estimates and marks grid cells for refining/coarsening
     */
    virtual void estimateAndMark( const DiscreteFunctionType& )
    {}


    virtual void markLocal( const ElementType& entity )
    {}


    //! mark elements
    void mark()
    {
      // loop over all elements
      for( const auto& en : elements( dfSpace_.gridPart() ) )
        markLocal( en );
    }

  };

}
}

#endif // #ifndef ESTIMATOR_HH
