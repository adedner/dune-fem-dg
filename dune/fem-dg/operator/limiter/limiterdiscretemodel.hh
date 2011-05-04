#ifndef DUNE_FEM_LIMITERDISCRETEMODEL_HH
#define DUNE_FEM_LIMITERDISCRETEMODEL_HH

#include <dune/fem/space/dgspace.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/operator/limiter/limitpass.hh>

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
  public:
    typedef LimiterTraits<GlobalPassTraitsImp,Model, passId > Traits;

    typedef typename Traits:: DomainType DomainType;
    typedef FieldVector<typename 
          DomainType :: field_type, Traits::dimDomain-1> FaceDomainType;
    
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType :: Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;

    typedef typename Traits::DestinationType DestinationType;

    enum { dimRange = Traits :: dimRange };

  protected:
    using BaseType :: model_;

    double refTol_;
    double crsTol_;
    int finLevel_;
    int crsLevel_;

  public:
    //! constructor 
    StandardLimiterDiscreteModel(const Model& mod, 
                         const int polOrd ) 
      : BaseType(mod), 
        refTol_(1), crsTol_(0.1), finLevel_(0), crsLevel_(0)
    {
      refTol_   = Parameter :: getValue("RefinementTolerance", refTol_ );
      crsTol_   = Parameter :: getValue("CoarseningTolerance", crsTol_ );
      finLevel_ = Parameter :: getValue("FinestLevel", finLevel_ );
      crsLevel_ = Parameter :: getValue("CoarsestLevel", crsLevel_ );
    }

    template < class IndicatorType >  
    void setIndicator(IndicatorType* ind) 
    {
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
                    const EntityType& en,
                    const RangeType& shockIndicator, 
                    const RangeType& adaptIndicator) const 
    {
      const double val = adaptIndicator[0];

      // get refinement marker 
      int refinement = grid.getMark( en );
      
      // if element already is marked for refinement then do nothing 
      if( refinement >= Refine ) return ;
      
      // use component 1 (max of adapt indicator)
      {
        const int level = en.level();
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
      grid.mark( refinement, en );
    }

    template <class ArgumentTuple>
    double numericalFlux(const IntersectionType& it,
                         const double time, 
                         const FaceDomainType& x,
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         RangeType& shockIndicator,
                         RangeType& adaptIndicator) const
    {
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);

      if (!physical(argULeft) || !physical(argURight)) 
      {
        adaptIndicator = shockIndicator = 1e10;
        return -1.;
      } 
      else 
      {
        // evaluate adaptation indicator 
        model_.adaptationIndicator(it, x, argULeft, argURight, adaptIndicator );
        model_.jump(argULeft, argURight, shockIndicator );
        return 1.;
      }
    }

    //! returns difference between internal value and boundary 
    //! value 
    template <class ArgumentTuple>
    double boundaryFlux(const IntersectionType& it,
                        const double time, 
                        const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& adaptIndicator) const
    {
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      RangeType uRight; 

      // evaluate boundary value 
      model_.boundaryValue(it, time,x,argULeft,uRight);
      
      if (!physical(argULeft) || !physical(uRight) ) 
      {
        adaptIndicator = 1e10;
        return -1.;
      }
      else 
      {
        const UType& argURight = uRight;
        model_.jump(argULeft, argURight, adaptIndicator);
        return 1.;
      }
    }
    
    //! returns difference between internal value and boundary 
    //! value 
    inline void boundaryValue(const IntersectionType& it,
                              const double time, 
                              const FaceDomainType& x,
                              const RangeType& uLeft,
                              RangeType& uRight) const
    {
      model_.boundaryValue(it, time, x, uLeft, uRight);
    }

    /** \brief returns true if model provides physical check */
    bool hasPhysical() const { return model_.hasPhysical(); } 

    /** \brief check physical values (called from StandardLimiterDiscreteModelCaller) */
    template <class ArgumentTuple>
    bool checkPhysical(const ArgumentTuple& u) const
    {
      // take first component 
      return physical( Element<0> :: get( u ) );
    }

    /** \brief check physical values */
    template <class ArgumentTuple>
    bool checkPhysical( const EntityType& en,
                   const DomainType& x,
                   const ArgumentTuple& u ) const
    { 
      return physical( Element<0> :: get( u ) );
    }

    /** \brief check physical values */
    bool physical( const EntityType& en,
                   const DomainType& x,
                   const RangeType& u ) const
    { 
      return physical( u );
    }

    /** \brief check physical values */
    bool physical(const RangeType& u) const
    {
      return model_.physical( u );
    }
  };

} // end namespace Fem 
} // end namespace Dune 
#endif
