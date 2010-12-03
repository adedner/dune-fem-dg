#ifndef DUNE_FEM_LIMITER_HH
#define DUNE_FEM_LIMITER_HH 

#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/fvspace.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/limiter/limitpass.hh>

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
        model_.jump(argULeft, uRight, adaptIndicator);
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
    bool physical(const RangeType& u) const
    {
      return model_.physical( u );
    }
  };

  template <class Model, class DiscreteFunctionType> 
  class LimitedReconstruction
  {
    // MethodOrderTraits
    template <class GridPartImp, class ftype, int dimRange,int polOrd>
    class PassTraits 
    {
    public:
      typedef ftype        FieldType;
      typedef GridPartImp  GridPartType ;
      typedef typename GridPartType :: GridType GridType;
      typedef typename GridType :: ctype ctype;
      enum { dimDomain = GridType :: dimensionworld };

      typedef CachingQuadrature<GridPartType,0> VolumeQuadratureType;
      typedef CachingQuadrature<GridPartType,1> FaceQuadratureType;

      typedef FunctionSpace< ctype, FieldType, dimDomain, dimRange> FunctionSpaceType; 
      typedef DiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, 
                 polOrd> DiscreteFunctionSpaceType;

      typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
      typedef DestinationType DiscreteFunctionType;

      // Indicator for Limiter 
      typedef FunctionSpace< ctype, FieldType, dimDomain, 3> FVFunctionSpaceType;
      typedef FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, SimpleStorage> IndicatorSpaceType;
      typedef AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;
    };

  public:  
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    typedef typename GridPartType :: GridType   GridType;

    typedef typename GridType::template Codim<0>::Entity EntityType;

    // make space exchangeable 
    typedef DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 1 >
      ReconstructionSpaceType;

    // type of storage for reconstruction 
    typedef AdaptiveDiscreteFunction< ReconstructionSpaceType >  ReconstructionType ; 

    //! define local function, this could also be a temporary local function 
    typedef typename ReconstructionType :: LocalFunctionType  LocalFunctionType;

    // pass ids 
    enum PassIdType{ u = 0 , limitPass = 1 };
    enum { dimRange = DiscreteFunctionSpaceType :: dimRange };

    typedef PassTraits<Model, dimRange, 1> PassTraitsType;
    typedef StandardLimiterDiscreteModel<PassTraitsType, Model, u > DiscreteModel1Type;

  public:
    LimitedReconstruction( const DiscreteFunctionSpaceType& space )
      : space_( space.gridPart() ) 
      , reconstruction_( "LimitedReconstruction" , space_ );
    {} 

    //! calculate internal reconstruction 
    void update( const DiscreteFunctionType& arg ) 
    {
      // apply limit pass in any case
      limitPass_.enableFirstCall();
      // calculate reconstruction 
      limitPass_( arg, reconstruction_ );
    }
    
    //! return local reconstruction 
    LocalFunctionType localFunction( const EntityType& entity ) 
    {
      return reconstruction_.localFunction( entity );
    }

    //! return local reconstruction 
    const LocalFunctionType localFunction( const EntityType& entity ) const 
    {
      return reconstruction_.localFunction( entity );
    }

  private:
    LimitedReconstruction( const LimitedReconstruction& );

  protected:
    ReconstructionSpaceType space_;
    ReconstructionType reconstruction_;
  };

} // end namespace Fem 

} // end namespace Dune 

#endif // DUNE_FEM_LIMITER_HH
