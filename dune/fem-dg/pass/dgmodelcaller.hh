#ifndef DUNE_FEM_DG_DISCRETEMODELCALLER_HH
#define DUNE_FEM_DG_DISCRETEMODELCALLER_HH

#include <utility>
#include <vector>
#include <memory>

#include <dune/common/fvector.hh> 

#include <dune/fem/version.hh>

#include <dune/fem/pass/callerutility.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>

#include <dune/fem/pass/modelcallerdefault.hh> 

namespace Dune {

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class CDGDiscreteModelCaller 
    : public Fem::DiscreteModelCallerDefault< DiscreteModelImp, ArgumentImp, SelectorImp >
  {
    typedef Fem::DiscreteModelCallerDefault< DiscreteModelImp, ArgumentImp, SelectorImp > BaseType;
  public:
    typedef DiscreteModelImp DiscreteModelType;  
    typedef typename BaseType :: EntityType              EntityType;
    typedef typename BaseType :: JacobianRangeType       JacobianRangeType;
    typedef typename BaseType :: RangeType               RangeType;
    typedef typename BaseType :: JacobianRangeTupleType  JacobianRangeTupleType;
    typedef typename BaseType :: RangeTupleType          RangeTupleType;
    typedef typename BaseType :: LocalFunctionTupleType  LocalFunctionTupleType;
    typedef typename BaseType :: VolumeQuadratureType    VolumeQuadratureType;
    typedef typename BaseType :: FaceQuadratureType      FaceQuadratureType;
    typedef typename BaseType :: Intersection            Intersection;
    typedef typename BaseType :: MassFactorType          MassFactorType;
    typedef typename BaseType :: RangeCreator            RangeCreator;
    typedef typename BaseType :: JacobianCreator         JacobianCreator;

    typedef std::vector< RangeTupleType > RangeTupleVectorType ;
    typedef std::vector< JacobianRangeTupleType > JacobianRangeTupleVectorType ;

    enum { evaluateJacobian = DiscreteModelType :: evaluateJacobian };

    using BaseType::evaluateQuadrature;
    using BaseType::evaluateQuad;

    CDGDiscreteModelCaller(DiscreteModelType& problem) :
      BaseType( ),
      problem_( problem )
#ifndef NDEBUG
      , quadInnerId_( 0 )
      , quadOuterId_( 0 )
      , quadId_( 0 )
#endif
    {
      // we need at least size 1 
      vecValuesEn_.push_back( RangeTupleType( RangeCreator::apply() ) );
      vecValuesNb_.push_back( RangeTupleType( RangeCreator::apply() ) );

      valuesVec_.push_back( RangeTupleType( RangeCreator::apply() ) );
      jacobiansVec_.push_back( JacobianRangeTupleType( JacobianCreator::apply() )  );
      
      vecJacobiansEn_.push_back( JacobianRangeTupleType( JacobianCreator::apply() )  );
      vecJacobiansNb_.push_back( JacobianRangeTupleType( JacobianCreator::apply() )  );
    }

    void setEntity ( const EntityType &entity )
    {
      BaseType :: setEntity( entity );
      problem_.setEntity( entity );
    }

    template <class QuadratureType>
    void setEntity ( const EntityType &entity,
                     const QuadratureType& quad)
    {
      setEntity( entity );

      // evaluate all local functions for whole quadrature 
      resizeAndEvaluate( quad, data_->localFunctionsSelf(), valuesVec_,
          jacobiansVec_ , problem_.hasSource () || evaluateJacobian ); 
      // only when we have a source or diffusion term we need the jacobians 

#ifndef NDEBUG 
      quadId_ = quad.id();
#endif
    }

    void setNeighbor( const EntityType &neighbor )
    {
      BaseType :: setNeighbor( neighbor );
      problem_.setNeighbor( neighbor );
    }

    const JacobianRangeTupleType& 
    jacobianValue(const JacobianRangeTupleVectorType& jacobiansVec,
                  const int quadPoint) const
    {
      assert( ( evaluateJacobian ) ? (int) jacobiansVec.size() > quadPoint : true ); 
      return ( evaluateJacobian ) ? jacobiansVec[ quadPoint ] : jacobiansVec[ 0 ];
    }

    void analyticalFlux(const EntityType& en, 
                        const VolumeQuadratureType& quad, 
                        const int quadPoint,
                        JacobianRangeType& res) 
    {
      // make sure we get the right quadrature 
      assert( quadId_ == quad.id() );
      assert( (int) valuesVec_.size() > quadPoint );

      // for CDG: forward calculated jacobians to analyticalFlux in user's discrete model
      problem_.analyticalFlux(en, time_, quad.point(quadPoint), valuesVec_[ quadPoint ],
                              jacobianValue(jacobiansVec_, quadPoint) , res);
    }

    double analyticalFluxAndSource(const EntityType& en, 
                                   const VolumeQuadratureType& quad, 
                                   const int quadPoint,
                                   JacobianRangeType& fluxRes, 
                                   RangeType& sourceRes) 
    {
      // make sure we git the right quadrature 
      assert( quadId_ == quad.id() );
      assert( (int) valuesVec_.size() > quadPoint );
      
      // for CDG: forward calculated jacobians to analyticalFlux in user's discrete model
      problem_.analyticalFlux(en, time_, quad.point(quadPoint), valuesVec_[ quadPoint ],
                              jacobianValue(jacobiansVec_, quadPoint), fluxRes);

      // return time step restriction for source term (zero if none)
      return 
        problem_.source(en, time_, quad.point(quadPoint), valuesVec_[ quadPoint ],
                        jacobianValue(jacobiansVec_, quadPoint), sourceRes);
    }

    // Ensure: entities set correctly before call
    template <class QuadratureImp> 
    void initializeIntersection( const EntityType &neighbor,
                                 const Intersection& intersection,
                                 const QuadratureImp& quadInner,
                                 const QuadratureImp& quadOuter)
    {
      assert( intersection.neighbor () );

      setNeighbor( neighbor );

      // only when we have diffusion term jacobians are needed 
      resizeAndEvaluate( quadInner, data_->localFunctionsSelf(),  vecValuesEn_, 
                         vecJacobiansEn_,  evaluateJacobian );
      resizeAndEvaluate( quadOuter, data_->localFunctionsNeigh(), vecValuesNb_,
                         vecJacobiansNb_ , evaluateJacobian );

#ifndef NDEBUG 
      quadInnerId_ = quadInner.id();
      quadOuterId_ = quadOuter.id();
#endif

      // for CDG Pass 
      problem_.initializeIntersection( intersection, 
                                       time_,
                                       quadInner, 
                                       quadOuter, 
                                       vecValuesEn_,
                                       vecValuesNb_ );
    }

    // Ensure: entities set correctly before call
    template <class QuadratureImp> 
    void initializeBoundary(const Intersection& intersection,
                            const QuadratureImp& quadInner)
    {
      assert( intersection.boundary () );

      // only when we have diffusion term jacobians are needed 
      resizeAndEvaluate( quadInner, data_->localFunctionsSelf(), vecValuesEn_, 
                         vecJacobiansEn_, evaluateJacobian );
#ifndef NDEBUG 
      quadInnerId_ = quadInner.id();
#endif
      // for CDG Pass 
      problem_.initializeBoundary( intersection, 
                                   time_,
                                   quadInner, 
                                   vecValuesEn_);
    }

    // Ensure: entities set correctly before call
    template <class QuadratureType>
    double numericalFlux(const Intersection& intersection,
                         const QuadratureType& quadInner, 
                         const QuadratureType& quadOuter, 
                         const int quadPoint,
                         RangeType& resEn, 
                         RangeType& resNeigh,
                         JacobianRangeType& resEnGrad, 
                         JacobianRangeType& resNeighGrad) 
    {
      assert( vecValuesEn_.size() >= quadInner.nop() ); 
      assert( quadInnerId_ == quadInner.id() );
      assert( vecValuesNb_.size() >= quadInner.nop() ); 
      assert( quadOuterId_ == quadOuter.id() );

      return problem_.numericalFlux(intersection, time_, 
                                    quadInner,
                                    quadOuter,
                                    quadPoint,
                                    vecValuesEn_[ quadPoint ], 
                                    vecValuesNb_[ quadPoint ],
                                    jacobianValue(vecJacobiansEn_, quadPoint),
                                    jacobianValue(vecJacobiansNb_, quadPoint),
                                    resEn, 
                                    resNeigh,
                                    resEnGrad, 
                                    resNeighGrad);
    }

    double boundaryFlux(const Intersection& intersection,
                        const FaceQuadratureType& quad,
                        const int quadPoint,
                        RangeType& boundaryFlux,
                        JacobianRangeType& boundaryGradFlux ) 
    {
      assert( vecValuesEn_.size() >= quad.nop() ); 
      assert( quadInnerId_ == quad.id() );

      return problem_.boundaryFlux(intersection, time_, 
                                   quad, quadPoint,
                                   vecValuesEn_[ quadPoint ], 
                                   jacobianValue(vecJacobiansEn_, quadPoint),
                                   boundaryFlux, boundaryGradFlux );
    }

    //! return true when a mass matrix has to be build  
    bool hasMass () const  { return problem_.hasMass(); }

    //! evaluate mass matrix factor 
    void mass(const EntityType& en,
              const VolumeQuadratureType& quad,
              const int quadPoint,
              MassFactorType& m) const
    {
      assert( quadId_ == quad.id() );
      assert( (int) valuesVec_.size() > quadPoint );
      
      // call problem implementation 
      problem_.mass(en, time_, quad.point(quadPoint),
                    valuesVec_[ quadPoint ],
                    m);
    }
    
  protected:  
    template <class QuadratureImp> 
    void resizeAndEvaluate(const QuadratureImp& faceQuad,
                           LocalFunctionTupleType& lfs,
                           RangeTupleVectorType& vecValues,
                           JacobianRangeTupleVectorType& vecJacobians,
                           const bool evalJacobians ) 
    {
#if DUNE_VERSION_NEWER_REV(DUNE_FEM,1,1,0)
      evaluateQuadrature( faceQuad, lfs, vecValues );
      if ( evalJacobians )
      {
        evaluateQuadrature( faceQuad, lfs, vecJacobians );
      }
#else 
      const size_t quadNop = faceQuad.nop();
      assert( vecValues.size() == vecJacobians.size() );
      if( vecValues.size() < quadNop ) 
      {
        while( vecValues.size() < quadNop ) 
        {
          vecValues.push_back( RangeTupleType( RangeCreator::apply() ) );
          vecJacobians.push_back( JacobianRangeTupleType( JacobianCreator::apply() ) );
        }
      }

      for(size_t quadPoint = 0; quadPoint < quadNop; ++quadPoint )
      {
        evaluateQuad(faceQuad, quadPoint, lfs, vecValues[ quadPoint ]);
      }

      if ( evalJacobians ) 
      {
        for(size_t quadPoint = 0; quadPoint < quadNop; ++quadPoint )
        {
          evaluateJacobianQuad( faceQuad, quadPoint, lfs, vecJacobians[ quadPoint ] );
        }
      }
#endif
    }

  private:
    CDGDiscreteModelCaller(const CDGDiscreteModelCaller&);
    CDGDiscreteModelCaller& operator=(const CDGDiscreteModelCaller&);

  protected:
    DiscreteModelType& problem_;

  protected:
    using BaseType::data_;
    using BaseType::valuesEn_;
    using BaseType::valuesNeigh_;
    using BaseType::jacobians_;
    using BaseType::time_;

    RangeTupleVectorType vecValuesEn_;
    RangeTupleVectorType vecValuesNb_;

    RangeTupleVectorType valuesVec_;
    JacobianRangeTupleVectorType jacobiansVec_;

    // for jacobians on entity use variable from BaseType
    JacobianRangeTupleVectorType vecJacobiansEn_;
    // need also variable for jacobians on neighbor 
    JacobianRangeTupleVectorType vecJacobiansNb_;

#ifndef NDEBUG
    size_t quadInnerId_;
    size_t quadOuterId_;
    size_t quadId_;
#endif

  };

}

#endif
