#ifndef DUNE_FEM_DG_PRIMALDISCRETEMODELS_HH
#define DUNE_FEM_DG_PRIMALDISCRETEMODELS_HH

// Dune includes
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/operator/dg/discretemodelcommon.hh>
#include <dune/fem-dg/operator/fluxes/dgprimalfluxes.hh>

namespace Dune {

  // AdvectionDiffusionDGPrimalModel
  //--------------------------------

  //! \brief forward declaration
  template <class Model, class NumFlux, 
            DGDiffusionFluxIdentifier diffFluxId,
            int polOrd, int passUId,
            bool returnAdvectionPart, bool returnDiffsuionPart >
  class AdvectionDiffusionDGPrimalModel;


  // AdvectionDiffusionDGPrimalTraits
  //---------------------------------

  template <class Model, class NumFlux,
            DGDiffusionFluxIdentifier diffFluxId,
            int polOrd, int passUId,
            bool returnAdvectionPart, bool returnDiffusionPart >
  struct AdvectionDiffusionDGPrimalTraits 
   : public AdvectionTraits< Model, NumFlux, polOrd, passUId, -1, returnAdvectionPart> 
  {
    typedef AdvectionDiffusionDGPrimalModel
      < Model, NumFlux, diffFluxId, polOrd, passUId, 
        returnAdvectionPart, returnDiffusionPart >                   DGDiscreteModelType;
  };


  // AdvectionDiffusionDGPrimalModel
  //--------------------------------

  template< class Model, 
            class NumFlux, 
            DGDiffusionFluxIdentifier diffFluxId,
            int polOrd, 
            int passUId,
            bool returnAdvectionPart, bool returnDiffusionPart >
  class AdvectionDiffusionDGPrimalModel :
    public AdvectionModel< Model, NumFlux, polOrd, passUId, -1, returnAdvectionPart > 
  {
  public:
    typedef AdvectionDiffusionDGPrimalTraits
      < Model, NumFlux, diffFluxId, polOrd, passUId, returnAdvectionPart, returnDiffusionPart >
                                                           Traits; /*@LST0E@*/

    typedef AdvectionModel< Model, NumFlux, polOrd, 
                            passUId, -1, returnAdvectionPart>    BaseType;

    using BaseType :: uVar; 
    using BaseType :: inside;
    using BaseType :: outside;
    using BaseType :: model_;
    using BaseType :: uBnd_;

  public:
    enum { dimDomain = Traits :: dimDomain };
    enum { dimRange  = Traits :: dimRange };

    enum { advection = returnAdvectionPart }; // true if advection is enabled 
    enum { diffusion = returnDiffusionPart }; // this should be disabled for LDG 

    typedef FieldVector< double, dimDomain >               DomainType;
    typedef FieldVector< double, dimDomain-1 >             FaceDomainType;

#if defined TESTOPERATOR
#warning NO MASSOPERATOR APPLIED
    static const bool ApplyInverseMassOperator = false;
#else
    static const bool ApplyInverseMassOperator = true;
#endif

    typedef typename Traits :: GridPartType                            GridPartType;
    typedef typename Traits :: GridType                                GridType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename GridType :: template Codim< 0 > :: Entity         EntityType;
    typedef typename GridType :: template Codim< 0 > :: EntityPointer  EntityPointerType;
    typedef typename Traits :: RangeFieldType                          RangeFieldType;
    typedef typename Traits :: DomainFieldType                         DomainFieldType;
    typedef typename Traits :: RangeType                               RangeType;
    typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

    typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // type of diffusion flux implementation 
    typedef DGPrimalDiffusionFlux< DiscreteFunctionSpaceType, Model, diffFluxId > DiffusionFluxType;

    enum { evaluateJacobian = DiffusionFluxType :: evaluateJacobian  }; // we need to evaluate jacobians here
  public:
    /**
     * @brief constructor
     */
    AdvectionDiffusionDGPrimalModel(const Model& mod,
                                    const NumFlux& numf,
                                    const DiffusionFluxType& diffflux )
      : BaseType( mod, numf ),
        diffFlux_( diffflux )
    {
    }

    //! copy constructor (for thread parallel progs mainly)
    AdvectionDiffusionDGPrimalModel(const AdvectionDiffusionDGPrimalModel& other ) 
      : BaseType( other ),
        diffFlux_( other.diffFlux_ )
    {
    }

    inline bool hasSource() const 
    {                 
      return ( model_.hasNonStiffSource() || model_.hasStiffSource() );
    } /*@\label{dm:hasSource}@*/

    inline bool hasFlux() const { 
      return (advection || diffusion) && model_.hasFlux(); }     /*@LST0E@*/
                                                  /*@LST0S@*/
    /**
     * @brief analytical flux function$
     */
    template <class ArgumentTuple, class JacobianTuple >
    double source( const EntityType& en,
                   const double time, 
                   const DomainType& x,
                   const ArgumentTuple& u, 
                   const JacobianTuple& jac, 
                   RangeType& s ) const
    {
      s = 0;

      double dtEst = std::numeric_limits< double > :: max();

      if (diffusion && model_.hasStiffSource() )
      {
        const double dtStiff = model_.stiffSource( en, time, x, u[uVar], jac[uVar], s );
        dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
        maxDiffTimeStep_ = std::max( dtStiff, maxDiffTimeStep_ );
      }


      if (advection && model_.hasNonStiffSource() )
      {
        RangeType sNonStiff (0);
        const double dtNon = model_.nonStiffSource( en, time, x, u[uVar], jac[uVar], sNonStiff );

        // add to source 
        s += sNonStiff;

        dtEst = ( dtNon > 0 ) ? std::min( dtEst, dtNon ) : dtEst;
        maxAdvTimeStep_  = std::max( dtNon, maxAdvTimeStep_ );
      }

      // return the fastest wave from source terms
      return dtEst ;
    }

    void switchUpwind() const 
    { 
      // reset max time steps 
      BaseType :: switchUpwind() ;
      // switch upwind direction if necessary 
      diffFlux_.switchUpwind();
    }

    template <class QuadratureImp, class ArgumentTupleVector > 
    void initializeIntersection(const Intersection& it,
                                const double time,
                                const QuadratureImp& quadInner, 
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec) 
    {
      if( diffusion ) 
      {
        typedef DGFluxTupleToVectorConverter< ArgumentTupleVector, passUId, 0> ConverterType;
        // call diffusion flux 
        diffFlux_.initializeIntersection( it, inside(), outside(),
                                          time, quadInner, quadOuter,
                                          ConverterType( uLeftVec ), 
                                          ConverterType( uRightVec ) );
      }
    }

    template <class QuadratureImp, class ArgumentTupleVector > 
    void initializeBoundary(const Intersection& it,
                            const double time,
                            const QuadratureImp& quadInner, 
                            const ArgumentTupleVector& uLeftVec)
    {
      if( diffusion ) 
      {
        typedef DGFluxTupleToVectorConverter< ArgumentTupleVector, passUId, 0> ConverterType;
        if( diffFlux_.hasLifting() ) 
        {
          const bool hasBoundaryValue = model_.hasBoundaryValue(it, time, quadInner.localPoint( 0 ));

          // turns uLeftVec into a fake 
          // vector of RangeTypes by storing the passUId 
          ConverterType uLeft ( uLeftVec ) ;

          const size_t quadNop = quadInner.nop();
          if( uBndVec_.size() < quadNop ) uBndVec_.resize( quadNop );

          for(size_t qp = 0; qp < quadNop; ++qp)
          {
            assert( hasBoundaryValue == 
                model_.hasBoundaryValue(it, time, quadInner.localPoint( qp )) );

            if( hasBoundaryValue ) 
              model_.boundaryValue(it, time, quadInner.localPoint( qp ), 
                                   uLeft[ qp ], uBndVec_[ qp ] );
            else
              // do something bad to uBndVec as it shouldn't be used
              uBndVec_[ qp ] = std::numeric_limits< double >::quiet_NaN();
          }

          // call diffusion flux 
          diffFlux_.initializeBoundary( it, inside(), 
                                        time, quadInner, 
                                        uLeft, uBndVec_ );
        }
      }
    } 

  public:
    /**
     * @brief flux function on interfaces between cells
     *
     * @param it intersection
     * @param time current time given by TimeProvider
     * @param x coordinate of required evaluation local to \c it
     * @param uLeft DOF evaluation on this side of \c it
     * @param uRight DOF evaluation on the other side of \c it
     * @param gLeft result for this side of \c it
     * @param gRight result for the other side of \c it
     * @return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class QuadratureImp,
              class ArgumentTuple, 
              class JacobianTuple >          /*@LST0S@*/
    double numericalFlux(const Intersection& it,
                         const double time,
                         const QuadratureImp& faceQuadInner,
                         const QuadratureImp& faceQuadOuter,
                         const int quadPoint, 
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         const JacobianTuple& jacLeft,
                         const JacobianTuple& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      ///////////////////////////////////////////////////////////////
      // for SIPG the total numerical flux for multiplication with phi
      // is given with
      //    gLeft = numflux(f(u))*n - {G(u)grad(u)}*n + C11/h mu [u]*n
      // and for STSIPG
      //    gLeft = numflux(f(u))*n - {G(u)grad(u)}*n + C11/h {G(u)}[u]*n
      // The total numerical flux for multiplication with grad(phi) is
      //    gDiffLeft = -0.5*G(u-)[u]
      // where h = min(|e+|,|e-|) / |it|.
      //
      // INFO: SIPG gives a suboptimal EOC for k=0,2,
      //       STSIPG gives an optimal results for EOC
      //       see R.Hartman,P.Houston JCP2008
      ///////////////////////////////////////////////////////////////
      
      /*****************************
       * Advection                 *
       ****************************/
      double wave = BaseType :: 
        numericalFlux( it, time, faceQuadInner, faceQuadOuter,
                       quadPoint, uLeft, uRight, jacLeft, jacRight, 
                       gLeft, gRight, gDiffLeft, gDiffRight );

      double diffTimeStep = 0.0;
      if( diffusion ) 
      {
        RangeType dLeft, dRight; 

        diffTimeStep = 
          diffFlux_.numericalFlux(it, *this, 
                                  time, faceQuadInner, faceQuadOuter, quadPoint,
                                  uLeft[ uVar ], uRight[ uVar ], 
                                  jacLeft[ uVar ], jacRight[ uVar ], 
                                  dLeft, dRight,
                                  gDiffLeft, gDiffRight);

        gLeft  += dLeft;
        gRight += dRight;
      }
      else
      {
        gDiffLeft = 0;
        gDiffRight = 0;
      }

      maxAdvTimeStep_  = std::max( wave, maxAdvTimeStep_ );
      maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );

      return std::max( wave, diffTimeStep );
    }


    /**
     * @brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class QuadratureImp, 
              class ArgumentTuple, class JacobianTuple>
    double boundaryFlux(const Intersection& it,
                        const double time, 
                        const QuadratureImp& faceQuadInner,
                        const int quadPoint,
                        const ArgumentTuple& uLeft,
                        const JacobianTuple& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
    {
      /****************************/
      /* Advection                *
       ****************************/
      double wave = BaseType :: 
        boundaryFlux( it, time, faceQuadInner, quadPoint,
                      uLeft, jacLeft, gLeft, gDiffLeft );
                                  
      /****************************/
      /* Diffusion                 *
       ****************************/
      double diffTimeStep = 0.0;

      const bool hasBoundaryValue = 
        model_.hasBoundaryValue( it, time, faceQuadInner.localPoint(0) );

      if( diffusion && hasBoundaryValue ) 
      {
        // diffusion boundary flux for Dirichlet boundaries 
        RangeType dLeft;
        diffTimeStep = 
          diffFlux_.boundaryFlux(it, 
                                 *this, 
                                 time, faceQuadInner, quadPoint,
                                 uLeft[ uVar ], uBnd_, // is set during call of  BaseType::boundaryFlux
                                 jacLeft[ uVar ], 
                                 dLeft,
                                 gDiffLeft);
        gLeft += dLeft;
      }
      else if ( diffusion ) 
      {
        RangeType diffBndFlux;
        model_.diffusionBoundaryFlux( it, time, faceQuadInner.localPoint(quadPoint),
                                      uLeft[uVar], jacLeft[uVar], diffBndFlux );
        gLeft += diffBndFlux;
      }
      else 
        gDiffLeft = 0;

      maxAdvTimeStep_  = std::max( wave, maxAdvTimeStep_ );
      maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );

      wave = std::max( wave, diffTimeStep );
      return wave;
    }
                                                  /*@LST0S@*/
    /**
     * @brief analytical flux function$
     */
    template <class ArgumentTuple, class JacobianTuple >
    void analyticalFlux( const EntityType& en,
                         const double time, 
                         const DomainType& x,
                         const ArgumentTuple& u, 
                         const JacobianTuple& jac, 
                         JacobianRangeType& f ) const
    {
      /*****************************
       * Advection                 *
       ****************************/
      BaseType :: analyticalFlux( en, time, x, u, jac, f );

      /*****************************
       * Diffusion                 *
       ****************************/
      if( diffusion && hasFlux() ) 
      {
        JacobianRangeType diffmatrix;
        model_.diffusion(en, time, x, u[ uVar ], jac[ uVar ], diffmatrix);
        // primal case 
        f -= diffmatrix;
      }
    }

    //! return reference to diffusion flux 
    DiffusionFluxType& diffusionFlux() const { return diffFlux_; }

  protected:
    // store flux here for thread parallel programs 
    // since the diffusion flux contains a lot of 
    // temporary class variables  
    mutable DiffusionFluxType diffFlux_;
    // storage for boundary values 
    std::vector< RangeType > uBndVec_;

    // defined in AdvectionModel 
    using BaseType :: maxAdvTimeStep_ ;
    using BaseType :: maxDiffTimeStep_ ;
  };                                              /*@LST0E@*/

} // namespace Dune

#endif
