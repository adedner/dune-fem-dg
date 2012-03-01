#ifndef METSTROEM_ESTIMATOR_HH
#define METSTROEM_ESTIMATOR_HH

//- Dune-fem includes 
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>


// Estimator
// ---------

/** \class Estimator
 *  \brief A class for estimating and marking grid elements for refinement and coarsening
 *
 *  This is a indicator based on a 'gradient' of the numerical solution.
 *  It evaluates numerical solution in an element's barycenter, and compares
 *  this value with the corresponding neighbor's value.
 *  The estimate is based on two indicators, e.g. one indicator could be 
 *  the density and the second vertical velocity.
 *  
 *  \note One may expect that this indicator will not give good results for higher order 
 *    methods as its results are solely based on a single value within one grid entity.
 *
 *  \tparam DiscreteFunction Discrete function type
 */
template< class DiscreteFunction, class Problem >
class Estimator : public EstimatorBase< DiscreteFunction >
{
  typedef Estimator< DiscreteFunction, Problem >              ThisType;
  typedef EstimatorBase< DiscreteFunction >                   BaseType;
public:
  typedef Problem                                             ProblemType;
  typedef DiscreteFunction                                    DiscreteFunctionType;

  typedef typename BaseType :: DiscreteFunctionSpaceType      DiscreteFunctionSpaceType;
  typedef typename BaseType :: LocalFunctionType              LocalFunctionType;

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
  typedef typename BaseType :: ElementPointerType             ElementPointerType;
  typedef typename BaseType :: GridElementType                GridElementType;
  typedef typename BaseType :: GeometryType                   GeometryType;

  static const int dimension = GridType :: dimension;

  typedef typename BaseType :: ElementQuadratureType          ElementQuadratureType;
  typedef typename BaseType :: FaceQuadratureType             FaceQuadratureType;
  typedef typename BaseType :: JacobianInverseType            JacobianInverseType;
  typedef typename BaseType :: IndicatorType                  IndicatorType;
  typedef Dune::GenericReferenceElement
    < DomainFieldType, dimension >                            ReferenceElementType;
  typedef Dune::GenericReferenceElements
    < DomainFieldType, dimension >                            ReferenceElementContainerType;

public:
  using BaseType :: mark;
  using BaseType :: clear;

private:
  using BaseType :: uh_;
  using BaseType :: dfSpace_;
  using BaseType :: gridPart_;
  using BaseType :: indexSet_;
  using BaseType :: grid_;
  using BaseType :: indicator_;
  IndicatorType*    indicator2Ptr_;
  const double refineTolerance_;
  const double coarseTolerance_;
  const int finestLevel_;
  const int coarsestLevel_;
  double ind1MaxDiff_;
  double ind2MaxDiff_;
  const int neighborRefLevel_;

protected:
  const ProblemType& problem_;

protected:
  //! \brief calculates the coordinates of the barycenter for given grid entity
  const DomainType& localBarycenterPoint( const ElementType& entity ) const
  {
    const ReferenceElementType& referenceElement
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
    const IntersectionIteratorType nbend = gridPart_.iend( entity ); 
    for (IntersectionIteratorType nb = gridPart_.ibegin( entity ); 
         nb != nbend; ++nb)
    {
      const IntersectionType& intersection = *nb ;
      if( intersection.neighbor() )
      {
        ElementPointerType outside = intersection.outside();
        const GridElementType& neighbor = Dune :: Fem :: gridEntity( *outside ); 
        // only do the following when the neighbor is not a ghost entity 
        if( neighbor.partitionType() != GhostEntity ) 
        { 
          if ( (neighbor.level() < finestLevel_) || (! neighbor.isRegular()) )
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
  explicit Estimator ( const DiscreteFunctionType &uh, 
                       const Problem& problem,
                       const Dune::AdaptationParameters& param = Dune::AdaptationParameters() )
  : BaseType( uh ),
    indicator2Ptr_( 0 ),
    refineTolerance_( param.refinementTolerance() ),
    coarseTolerance_( param.coarsenTolerance() ),
    finestLevel_( param.finestLevel( Dune::DGFGridInfo<GridType>::refineStepsForHalf() ) ),
    coarsestLevel_( param.coarsestLevel( Dune::DGFGridInfo<GridType>::refineStepsForHalf() ) ),
    neighborRefLevel_( param.neighborRefLevel() ),
    problem_( problem )
  {
    if( problem.twoIndicators() )
    {
      indicator2Ptr_ = new IndicatorType( indexSet_.size( 0 ) );
    }
  }


  void clear()
  {
    BaseType::clear();
    if( indicator2Ptr_ )
      clear( *indicator2Ptr_ );
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
  void estimateLocal( const ElementType& entity, double& ind1Min, double& ind1Max,
                      double& ind2Min, double& ind2Max )
  { 
    const int enIdx = indexSet_.index( entity ); 
    
    RangeType val( 0. );
    RangeType valnb( 0. );

    // get local function on the element
    LocalFunctionType lf = uh_.localFunction( entity );
    const int quadOrder = ( lf.order()==0 ? 1 : lf.order() );
    ElementQuadratureType quad( entity, quadOrder );
    const int numQuad = quad.nop();

    for( int qp=0; qp<numQuad; ++qp )
    {
      DomainType xEn = quad.point( qp );
      lf.evaluate( xEn, val );
      const double ind1 = indicator1( entity, xEn, val );
      ind1Max = std::max( ind1Max, ind1 );
      ind1Min = std::min( ind1Min, ind1 );

      double ind2 = 0.;
      if( indicator2Ptr_ )
      {
        ind2 = indicator2( entity, xEn, val );
        ind2Max = std::max( ind2Max, ind2 );
        ind2Min = std::min( ind2Min, ind2 );
      }

      // iterate over neighbors
      const IntersectionIteratorType nbend = gridPart_.iend( entity ); 
      for (IntersectionIteratorType nb = gridPart_.ibegin( entity ); 
           nb != nbend; ++nb)
      {
        const IntersectionType& intersection = *nb ;
        if( intersection.neighbor() )
        {
          // access neighbor
          ElementPointerType outside = intersection.outside();
          const ElementType& neighbor = *outside; 
          const int nbIdx = indexSet_.index( neighbor );

          // handle face from one side only
          if ( neighbor.partitionType() == GhostEntity || 
                entity.level() > neighbor.level() ||
              (entity.level() == neighbor.level() && enIdx < nbIdx) )
          {
            // get local function on the neighbor element
            LocalFunctionType lfnb = uh_.localFunction( neighbor );
            ElementQuadratureType quadNeigh( entity, quadOrder );
            const int numQuadNe = quad.nop();

            // run over neighbor quadrature points
            for( int qpNe=0; qpNe<numQuadNe; ++qpNe )
            {
              DomainType xNe = quadNeigh.point( qpNe );
              lfnb.evaluate( xNe, valnb );
              const double ind1nb = indicator1( neighbor, xNe, valnb );
              double ind1LocalDiff = std::abs( ind1 - ind1nb );
              indicator_[enIdx] = std::max( indicator_[enIdx], ind1LocalDiff );
              indicator_[nbIdx] = std::max( indicator_[nbIdx], ind1LocalDiff );

              if( indicator2Ptr_ )
              {
                IndicatorType& indicator2_ = *indicator2Ptr_;
                const double ind2nb = indicator2( neighbor, xNe, valnb );
                double ind2LocalDiff = std::abs( ind2 - ind2nb );
                indicator2_[enIdx] = std::max( indicator2_[enIdx], ind2LocalDiff );
                indicator2_[nbIdx] = std::max( indicator2_[nbIdx], ind2LocalDiff );
              }
            }
          }
        }
      }
    } 
  }


  //! \brief calculates indicators
  void estimate()
  {
    indicator_.resize( indexSet_.size( 0 )); 
    clear();

    double ind1Max = -1E100;
    double ind1Min =  1E100;
    double ind2Max = -1E100;
    double ind2Min =  1E100;
    ind2MaxDiff_ = 0.;

    const IteratorType end = dfSpace_.end();
    for( IteratorType it = dfSpace_.begin(); it != end; ++it )
      estimateLocal( *it, ind1Min, ind1Max, ind2Min, ind2Max );

    // global max 
    {
      double commBuff[ 2 ] = { ind1Max, ind2Max };

      // get global max of indMax 
      gridPart_.grid().comm().max( &commBuff[ 0 ], 2 );

      ind1Max = commBuff[ 0 ];
      ind2Max = commBuff[ 1 ];
    }

    // global min 
    {
      double commBuff[ 2 ] = { ind1Min, ind2Min };

      // get global min of indMin 
      gridPart_.grid().comm().min( &commBuff[ 0 ], 2 );

      ind1Min = commBuff[ 0 ];
      ind2Min = commBuff[ 1 ];
    }

    // return global max differences
    ind1MaxDiff_ = ind1Max - ind1Min;
    if( indicator2Ptr_ )
      ind2MaxDiff_ = ind2Max - ind2Min;
  }


  /** \brief marks an grid entity for refinement/coarsening
   */
  virtual void markLocal( const ElementType& entity )
  {
    // get local error indicator 
    const int entityId = indexSet_.index(entity);

    double localIndicator1 = indicator_[ entityId ]; 
    const double locRefTol1 = refineTolerance_ * ind1MaxDiff_;
    const double locCoarTol1 = coarseTolerance_ * ind1MaxDiff_;

    // check if element is allowed to be refined by the problem settings
    const DomainType& xEn = entity.geometry().center();
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

    const GridElementType& element = Dune :: Fem :: gridEntity( entity );

    // do not mark ghost elements 
    if( element.partitionType() == GhostEntity ) return ;

    const int elemLevel = element.level();

    if ( toBeRefined && (elemLevel < finestLevel_) )
    {
      // mark for refinement 
      grid_.mark( 1, element );

      // also mark distant neighbors of the actual entity for refinement
      // markNeighborsForRefinement( entity, neighborRefLevel_ );
    }
    else if ( toBeCoarsend && (elemLevel > coarsestLevel_) )
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
  void estimateAndMark()
  {
    estimate();
    mark();
  }
  
};


#endif // #ifndef ESTIMATOR_HH 
