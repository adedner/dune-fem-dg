#ifndef DUNE_FEM_LIMITERDISCRETEMODEL_HH
#define DUNE_FEM_LIMITERDISCRETEMODEL_HH

#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>

namespace Dune {

namespace Fem {

  ///////////////////////////////////////////////////////////////////
  //
  //  --StandardLimiterDiscreteModel
  //
  ///////////////////////////////////////////////////////////////////

  // **********************************************
  // **********************************************
  // **********************************************
  template <class GlobalPassTraitsImp, class Model, int passId = -1 >
  class StandardLimiterDiscreteModel;

  template <class GlobalPassTraitsImp, class Model, int passId >
  struct LimiterTraits
    : public LimiterDefaultTraits<GlobalPassTraitsImp,Model,passId>
  {
    typedef StandardLimiterDiscreteModel<GlobalPassTraitsImp,Model,passId> DGDiscreteModelType;
  };

  template <class GlobalPassTraitsImp, class Model, int passId >
  class StandardLimiterDiscreteModel :
    public LimiterDefaultDiscreteModel<GlobalPassTraitsImp, Model , passId >
  {
    typedef LimiterDefaultDiscreteModel<GlobalPassTraitsImp, Model, passId > BaseType;

    // These type definitions allow a convenient access to arguments of pass.
    integral_constant< int, passId > uVar;
  public:
    typedef LimiterTraits<GlobalPassTraitsImp,Model, passId > Traits;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits:: DomainType DomainType;
    typedef typename Traits:: LocalDomainType LocalDomainType;
    typedef typename Traits:: DomainFieldType DomainFieldType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename GridPartType::IntersectionIteratorType :: Intersection IntersectionType;
    typedef typename GridPartType::template Codim<0>::EntityType        EntityType;
    typedef typename GridPartType::template Codim<0>::EntityPointerType EntityPointerType;
    typedef typename GridType::template Codim<0>::Entity                GridEntityType;

    enum { dimGrid = GridType :: dimension };

    // type of surface domain type
    typedef FieldVector< DomainFieldType, dimGrid - 1 > FaceLocalDomainType;

    typedef typename Traits::DestinationType DestinationType;

    // Indicator for Limiter
    typedef typename Traits::IndicatorType   IndicatorType ;

    enum { dimRange = Traits :: dimRange };
    enum { evaluateJacobian = false };

  protected:
    using BaseType :: model_;
    using BaseType :: inside ;
    using BaseType :: outside ;

    IndicatorType* indicator_;

    double refTol_;
    double crsTol_;
    int finLevel_;
    int crsLevel_;
    const bool shockIndicatorAdaptivty_;

  public:
    //! constructor
    StandardLimiterDiscreteModel(const Model& mod,
                                 const int polOrd,
                                 const AdaptationParameters& param = AdaptationParameters() )
      : BaseType(mod),
        indicator_( 0 ),
        refTol_( -1 ),
        crsTol_( -1 ),
        finLevel_( 0 ),
        crsLevel_( 0 ),
        shockIndicatorAdaptivty_( param.shockIndicator() )
    {
      if( shockIndicatorAdaptivty_ )
      {
        refTol_   = param.refinementTolerance();
        crsTol_   = param.coarsenTolerance();
        finLevel_ = param.finestLevel( DGFGridInfo<GridType>::refineStepsForHalf() );
        crsLevel_ = param.coarsestLevel( DGFGridInfo<GridType>::refineStepsForHalf() );
      }
    }

    void setIndicator(IndicatorType* indicator)
    {
      indicator_ = indicator;
    }

    void setEntity(const EntityType& en)
    {
      BaseType :: setEntity ( en );
    }

    void indicatorMax() const
    {
    }

    enum { Coarsen = -1 , None = 0, Refine = 1 };

    //! mark element for refinement or coarsening
    void adaptation(GridType& grid,
                    const EntityType& entity,
                    const RangeType& shockIndicator,
                    const RangeType& adaptIndicator) const
    {
      const double val = adaptIndicator[0];

      // set indicator
      if( indicator_ )
      {
        typedef typename IndicatorType :: LocalFunctionType  LocalFunctionType;
        LocalFunctionType lf = indicator_->localFunction( entity );
        assert( lf.numDofs() == 3 );
        lf[0] = shockIndicator[0];
        lf[1] = adaptIndicator[0];
        lf[2] = val; // store real adapt indicator
      }

      // only mark for adaptation if this type of adaptation is enabled
      if( ! shockIndicatorAdaptivty_ ) return ;

      // get real grid entity from possibly wrapped entity
      const GridEntityType& gridEntity = Fem :: gridEntity( entity );

      // get refinement marker
      int refinement = grid.getMark( gridEntity );

      // if element already is marked for refinement then do nothing
      if( refinement >= Refine ) return ;

      // use component 1 (max of adapt indicator)
      {
        const int level = gridEntity.level();
        if( (( val > refTol_ ) || (val < 0)) && level < finLevel_)
        {
          refinement = Refine;
        }
        else if ( (val >= 0) && (val < crsTol_)  && level > crsLevel_ )
        {
          if( refinement == None )
          {
            refinement = Coarsen;
          }
        }
      }

      // set new refinement marker
      grid.mark( refinement, gridEntity );
    }

    template <class FaceQuadratureImp,
              class ArgumentTuple,
              class JacobianTuple>
    double numericalFlux(const IntersectionType& it,
                         const double time,
                         const FaceQuadratureImp& innerQuad,
                         const FaceQuadratureImp& outerQuad,
                         const int quadPoint,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         const JacobianTuple& jacLeft,
                         const JacobianTuple& jacRight,
                         RangeType& shockIndicator,
                         RangeType& adaptIndicator,
                         JacobianRangeType&,
                         JacobianRangeType& ) const
    {
      const FaceLocalDomainType& x = innerQuad.localPoint( quadPoint );

      if (! physical(inside(), innerQuad.point( quadPoint ), uLeft[ uVar ] ) ||
          ! physical(outside(), outerQuad.point( quadPoint ), uRight[ uVar ] ) )
      {
        adaptIndicator = shockIndicator = 1e10;
        return -1.;
      }
      else
      {
        // evaluate adaptation indicator
        model_.adaptationIndicator(it, time, x, uLeft[ uVar ], uRight[ uVar ], adaptIndicator );
        model_.jump( it, time, x, uLeft[ uVar ], uRight[ uVar ], shockIndicator );
        return 1.;
      }
    }

    //! returns difference between internal value and boundary
    //! value
    template <class FaceQuadratureImp,
              class ArgumentTuple,
              class JacobianTuple>
    double boundaryFlux(const IntersectionType& it,
                        const double time,
                        const FaceQuadratureImp& innerQuad,
                        const int quadPoint,
                        const ArgumentTuple& uLeft,
                        const JacobianTuple& jacLeft,
                        RangeType& adaptIndicator,
                        JacobianRangeType& gDiffLeft ) const
    {
      const FaceLocalDomainType& x = innerQuad.localPoint( quadPoint );

      RangeType uRight;

      // evaluate boundary value
      model_.boundaryValue(it, time, x, uLeft[ uVar ], uRight );

      if (! physical(inside(), innerQuad.point( quadPoint ), uLeft[ uVar ] ) ||
          ! physical(inside(), innerQuad.point( quadPoint ), uRight ) )
      {
        adaptIndicator = 1e10;
        return -1.;
      }
      else
      {
        model_.jump( it, time, x, uLeft[ uVar ], uRight, adaptIndicator);
        return 1.;
      }
    }

    //! returns difference between internal value and boundary
    //! value
    inline void boundaryValue(const IntersectionType& it,
                              const double time,
                              const FaceLocalDomainType& x,
                              const RangeType& uLeft,
                              RangeType& uRight) const
    {
      model_.boundaryValue(it, time, x, uLeft, uRight);
    }

    /** \brief returns true if model provides physical check */
    bool hasPhysical() const { return model_.hasPhysical(); }

    /** \brief check physical values */
    template <class ArgumentTuple>
    bool checkPhysical( const EntityType& entity,
                        const LocalDomainType& xLocal,
                        const ArgumentTuple& u ) const
    {
      return physical( entity, xLocal, u[ uVar ] );
    }

    /** \brief check physical values */
    bool physical( const EntityType& entity,
                   const LocalDomainType& xLocal,
                   const RangeType& u ) const
    {
      return model_.physical( entity, xLocal, u );
    }

    /** \brief adjust average values, e.g. transform to primitive or something similar */
    void adjustAverageValue( const EntityType& entity,
                             const LocalDomainType& xLocal,
                             RangeType& u ) const
    {
      model_.adjustAverageValue( entity, xLocal, u );
    }
  };

} // end namespace Fem
} // end namespace Dune
#endif
