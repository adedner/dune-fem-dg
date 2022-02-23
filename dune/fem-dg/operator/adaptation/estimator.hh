#ifndef METSTROEM_ESTIMATOR_HH
#define METSTROEM_ESTIMATOR_HH

//- Dune-fem includes
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>

namespace Dune
{
namespace Fem
{
  // Estimator
  // ---------

  /** \class Estimator
   *  \brief A class for estimating and marking grid elements for refinement and coarsening
   *
   *  This is a indicator based on a 'gradient' of the numerical solution.
   *  It evaluates maximum of the numerical solution in a grid element, and compares
   *  this value with the corresponding neighbor's value.
   *  The estimate is based on two indicators, e.g. one indicator could be
   *  the density and the second vertical velocity.
   *
   *  \tparam DiscreteFunction Discrete function type
   */
  template< class DiscreteFunction, class Problem >
  class Estimator :
    public EstimatorBase< DiscreteFunction >,
    public ComputeMinMaxVolume
  {
    typedef Estimator< DiscreteFunction, Problem >              ThisType;
    typedef EstimatorBase< DiscreteFunction >                   BaseType;
  public:
    typedef Problem                                             ProblemType;
    typedef DiscreteFunction                                    DiscreteFunctionType;

    typedef typename BaseType :: DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;
    typedef typename BaseType :: ConstLocalFunctionType         ConstLocalFunctionType;
    typedef typename BaseType :: DomainFieldType                DomainFieldType;
    typedef typename BaseType :: RangeFieldType                 RangeFieldType;
    typedef typename BaseType :: DomainType                     DomainType;
    typedef typename BaseType :: RangeType                      RangeType;
    typedef typename BaseType :: JacobianRangeType              JacobianRangeType;
    typedef typename BaseType :: GridPartType                   GridPartType;
    typedef typename BaseType :: IteratorType                   IteratorType;

    typedef typename BaseType :: GridType                       GridType;
    typedef typename BaseType :: IndexSetType                   IndexSetType;
    typedef typename BaseType :: IntersectionIteratorType       IntersectionIteratorType;
    typedef typename BaseType :: IntersectionType               IntersectionType;
    typedef typename BaseType :: ElementType                    ElementType;
    typedef typename BaseType :: GridElementType                GridElementType;
    typedef typename BaseType :: GeometryType                   GeometryType;

    static const int dimension = GridType :: dimension;

    typedef typename BaseType :: ElementQuadratureType          ElementQuadratureType;
    typedef typename BaseType :: FaceQuadratureType             FaceQuadratureType;
    typedef typename BaseType :: JacobianInverseType            JacobianInverseType;
    typedef typename BaseType :: IndicatorType                  IndicatorType;

    typedef Dune::ReferenceElements
      < DomainFieldType, dimension >                            ReferenceElementContainerType;

    typedef ComputeMinMaxVolume  ComputeMinMaxVolumeType;
  public:
    using BaseType :: mark;
    using BaseType :: clear;
    using BaseType :: estimateLocal ;
    using BaseType :: estimate ;

    using ComputeMinMaxVolumeType :: coarsestVolume ;
    using ComputeMinMaxVolumeType :: finestVolume ;
    using ComputeMinMaxVolumeType :: computeGlobalMinMax ;

  protected:
    using BaseType :: dfSpace_;
    using BaseType :: gridPart_;
    using BaseType :: indexSet_;
    using BaseType :: grid_;
    using BaseType :: indicator_;
    IndicatorType*    indicator2Ptr_;
    const double refineTolerance_;
    const double coarseTolerance_;
    double ind1MaxDiff_;
    double ind2MaxDiff_;
    const int neighborRefLevel_;
    size_t numberOfElements_ ;

  protected:
    const ProblemType& problem_;

  protected:
    //! \brief calculates the coordinates of the barycenter for given grid entity
    const DomainType& localBarycenterPoint( const ElementType& entity ) const
    {
      const auto referenceElement
                = ReferenceElementContainerType::general( entity.type() );
      return referenceElement.position( 0, 0 );
    }

    /** \brief caculates a first quantity (e.g. density, pot. temperature) whose
     *    gradient will be tracked for the grid adaptation
     *
     *  \param[in] entity Grid entity
     *  \param[in] x Point in local coordinates w.r.t. \a entity
     *  \param[in] u Values of the numerical solution in \a x
     */
    virtual double indicator1( const ElementType& entity,
                               const DomainType& x, const RangeType& u ) const
    {
      const DomainType& xgl = entity.geometry().global( x );
      return problem_.indicator1( xgl, u );
    }

    /** \brief caculates a second quantity (e.g. density, pot. temperature) whose
     *    gradient will be tracked for the grid adaptation
     *
     *  \param[in] entity Grid entity
     *  \param[in] x Point in local coordinates w.r.t. \a entity
     *  \param[in] u Values of the numerical solution in \a x
     */
    virtual double indicator2( const ElementType& entity,
                               const DomainType& x, const RangeType& u ) const
    {
      const DomainType& xgl = entity.geometry().global( x );
      return problem_.indicator2( xgl, u );
    }

    void markNeighborsForRefinement( const ElementType& entity, const int level ) const
    {
      if (level <= 0)
        return;

      // also mark all neighbors of the actual entity for refinement
      for (const auto& intersection : intersections(gridPart_, entity) )
      {
        if( intersection.neighbor() )
        {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
          const GridElementType &neighbor = intersection.outside();
#else
          typename GridElementType::EntityPointer outside = intersection.outside();
          const GridElementType& neighbor = Dune :: Fem :: gridEntity( *outside );
#endif
          // only do the following when the neighbor is not a ghost entity
          if( neighbor.partitionType() != Dune::GhostEntity )
          {
            if ( (neighbor.geometry().volume() > finestVolume()) || (! neighbor.isRegular()) )
            {
              // mark for refinement
              grid_.mark( 1, neighbor );
            }

            // mark further neighbors
            markNeighborsForRefinement( neighbor, level-1 );
          }
        }
      }
    }

  public:
    //! \brief Constructor
    explicit Estimator ( const DiscreteFunctionSpaceType &space,
                         const Problem& problem,
                         const AdaptationParameters& param = AdaptationParameters() )
    : BaseType( space ),
      ComputeMinMaxVolumeType( space.gridPart(),
                               param.coarsestLevel( Dune::DGFGridInfo<GridType>::refineStepsForHalf() ),
                               param.finestLevel( Dune::DGFGridInfo<GridType>::refineStepsForHalf() ) ),
      indicator2Ptr_( 0 ),
      refineTolerance_( param.refinementTolerance() ),
      coarseTolerance_( param.coarsenTolerance() ),
      neighborRefLevel_( param.neighborRefLevel() ),
      numberOfElements_( 0 ),
      problem_( problem )
    {
      if( problem.twoIndicators() )
      {
        indicator2Ptr_ = new IndicatorType( indexSet_.size( 0 ) );
      }
    }

    //! clear the stored indicators
    void clear()
    {
      BaseType::clear();
      if( indicator2Ptr_ )
        clear( *indicator2Ptr_ );

      numberOfElements_ = 0;
    }

    //! return number of leaf element for this process (global comm still needed)
    size_t numberOfElements () const
    {
      return numberOfElements_;
    }

    /** \brief calculates the maximum of the differences between the values in the center of
     *  the current grid entity and its neighbors
     *
     *  \param[in] entity Grid entity
     *  \param[out] ind1Min Minimal difference of the first indicator quantity (e.g. density)
     *    values in the center of \a entity and its neighbor
     *  \param[out] ind1Max Maximal difference of the first indicator quantity (e.g. density)
     *    values in the center of \a entity and its neighbor
     *  \param[out] ind2Min Minimal difference of the second indicator quantity (e.g. density)
     *    values in the center of \a entity and its neighbor
     *  \param[out] ind2Min Maximal difference of the second indicator quantity (e.g. density)
     *    values in the center of \a entity and its neighbor
     *
     *  \note \a indicator1_ and \a indicator2_ are assigned its correspond values
     *    for \a entity and its neighbor
     */
    void estimateLocal( const ConstLocalFunctionType& lf,
                        ConstLocalFunctionType& lfnb,
                        const ElementType& entity, double& ind1Min, double& ind1Max,
                        double& ind2Min, double& ind2Max )
    {
      const int enIdx = indexSet_.index( entity );

      RangeType val( 0. );
      RangeType valnb( 0. );

      // get local function on the element
      const int quadOrder = ( lf.order()==0 ? 1 : lf.order() );
      ElementQuadratureType quad( entity, quadOrder );

      // get max and min of the indicator quantity
      double ind1LocMax = -1E100;
      double ind1LocMin =  1E100;
      double ind2LocMax = -1E100;
      double ind2LocMin =  1E100;
      for( const auto qp : quad )
      {
        DomainType xEn = qp.position();
        lf.evaluate( xEn, val );
        const double ind1 = indicator1( entity, xEn, val );
        ind1LocMax = std::max( ind1LocMax, ind1 );
        ind1LocMin = std::min( ind1LocMin, ind1 );

        double ind2 = 0.;
        if( indicator2Ptr_ )
        {
          ind2 = indicator2( entity, xEn, val );
          ind2LocMax = std::max( ind2LocMax, ind2 );
          ind2LocMin = std::min( ind2LocMin, ind2 );
        }
      }
      ind1Max = std::max( ind1Max, ind1LocMax );
      ind1Min = std::min( ind1Min, ind1LocMin );
      if( indicator2Ptr_ )
      {
        ind2Max = std::max( ind2Max, ind2LocMax );
        ind2Min = std::min( ind2Min, ind2LocMin );
      }

      // iterate over neighbors
      for (const auto& intersection : intersections(gridPart_, entity) )
      {
        if( intersection.neighbor() )
        {
          // access neighbor
          const ElementType neighbor = intersection.outside();
          const int nbIdx = indexSet_.index( neighbor );

          // handle face from one side only
          if ( neighbor.partitionType() == Dune::GhostEntity ||
                entity.level() > neighbor.level() ||
              (entity.level() == neighbor.level() && enIdx < nbIdx) )
          {
            // get local function on the neighbor element
            auto guard = bindGuard( lfnb, neighbor );
            ElementQuadratureType quadNeigh( entity, quadOrder );

            // get max and min of the indicator quantity in the neighbor
            double ind1nbLocMax = -1E100;
            double ind1nbLocMin =  1E100;
            double ind2nbLocMax = -1E100;
            double ind2nbLocMin =  1E100;

            for( const auto qp : quadNeigh )
            {
              DomainType xNe = qp.position();
              lfnb.evaluate( xNe, valnb );
              const double ind1nb = indicator1( neighbor, xNe, valnb );
              ind1nbLocMax = std::max( ind1nbLocMax, ind1nb );
              ind1nbLocMin = std::min( ind1nbLocMin, ind1nb );

              if( indicator2Ptr_ )
              {
                const double ind2nb = indicator2( neighbor, xNe, valnb );
                ind2nbLocMax = std::max( ind2nbLocMax, ind2nb );
                ind2nbLocMin = std::min( ind2nbLocMin, ind2nb );
              }
            }

            // calcualte local difference and assign indicators
            double ind1LocalDiff = std::abs( ind1LocMax - ind1nbLocMin );
            ind1LocalDiff = std::max( ind1LocalDiff, std::abs( ind1LocMin - ind1nbLocMax ) );
            indicator_[enIdx] = std::max( indicator_[enIdx], ind1LocalDiff );
            indicator_[nbIdx] = std::max( indicator_[nbIdx], ind1LocalDiff );
            if( indicator2Ptr_ )
            {
              IndicatorType& indicator2_ = *indicator2Ptr_;
              double ind2LocalDiff = std::abs( ind2LocMax - ind2nbLocMin );
              ind2LocalDiff = std::max( ind2LocalDiff, std::abs( ind2LocMin - ind2nbLocMax ) );
              indicator2_[enIdx] = std::max( indicator2_[enIdx], ind2LocalDiff );
              indicator2_[nbIdx] = std::max( indicator2_[nbIdx], ind2LocalDiff );
            }
          }
        } // neighbor
      } // iteratore through neighbors

    }


    //! \brief calculates indicators
    void estimate( const DiscreteFunctionType& uh )
    {
      indicator_.resize( indexSet_.size( 0 ) );
      clear();

      double indMax[ 2 ] = { std::numeric_limits< double > :: min (), std::numeric_limits< double > :: min () };
      double indMin[ 2 ] = { std::numeric_limits< double > :: max (), std::numeric_limits< double > :: max () };
      ind2MaxDiff_   = 0;

      ConstLocalFunctionType lf( uh );
      ConstLocalFunctionType lfnb( uh );

      numberOfElements_ = 0 ;
      for( const auto& en : elements( dfSpace_.gridPart() ) )
      {
        auto guard = bindGuard( lf, en );

        // do local estimation
        estimateLocal( lf, lfnb, en, indMin[ 0 ], indMax[ 0 ], indMin[ 1 ], indMax[ 1 ] );
        // count number of elements
        ++ numberOfElements_ ;
      }

      // global max and min
      computeGlobalMinMax( gridPart_.grid().comm(), 2, indMax, indMin );

      // return global max differences
      ind1MaxDiff_ = indMax[ 0 ] - indMin [ 0 ];
      if( indicator2Ptr_ )
        ind2MaxDiff_ = indMax[ 1 ] - indMin[ 1 ];
    }


    /** \brief marks an grid entity for refinement/coarsening
     */
    virtual void markLocal( const ElementType& entity )
    {
      // do not mark ghost elements
      if( entity.partitionType() == Dune::GhostEntity ) return ;

      // get local error indicator
      const int entityId = indexSet_.index(entity);

      double localIndicator1 = indicator_[ entityId ];
      const double locRefTol1 = refineTolerance_ * ind1MaxDiff_;
      const double locCoarTol1 = coarseTolerance_ * ind1MaxDiff_;

      typedef typename ElementType :: Geometry Geometry ;
      const Geometry& geometry = entity.geometry();
      // check if element is allowed to be refined by the problem settings
      const DomainType& xEn = geometry.center();
      const double volume   = geometry.volume();

      const bool problemAllows = problem_.allowsRefinement( xEn );
      bool toBeRefined = (localIndicator1 > locRefTol1) && problemAllows;
      bool toBeCoarsend = (localIndicator1 < locCoarTol1);

      if( indicator2Ptr_ )
      {
        const IndicatorType& indicator2_ = *indicator2Ptr_;
        double localIndicator2 = indicator2_[ entityId ];
        const double locRefTol2 =  refineTolerance_ * ind2MaxDiff_;
        const double locCoarTol2 =  coarseTolerance_ * ind2MaxDiff_;
        toBeRefined = (toBeRefined || ((localIndicator2 > locRefTol2) && problemAllows));
        toBeCoarsend = (toBeCoarsend && (localIndicator2 < locCoarTol2));
      }

      // get grid's entity object form marking
      const GridElementType& element = Dune :: Fem :: gridEntity( entity );

      // allow refinement if volume of element is still bigger then smallest volume allowed
      if ( toBeRefined && (volume > finestVolume()) )
      {
        // mark for refinement
        grid_.mark( 1, element );

        // also mark distant neighbors of the actual entity for refinement
        // markNeighborsForRefinement( entity, neighborRefLevel_ );
      }
      else if ( toBeCoarsend && (volume < coarsestVolume()) )
      {
        // mark for coarsening
        grid_.mark( -1, element );
      }
      else
      {
        // don't do anything
        grid_.mark( 0, element );
      }
    }


    //! \brief estimate and mark grid entities for refinement/coarsening
    void estimateAndMark( const DiscreteFunctionType& uh )
    {
      estimate( uh );
      mark();
    }

  };

}
}
#endif // #ifndef ESTIMATOR_HH
