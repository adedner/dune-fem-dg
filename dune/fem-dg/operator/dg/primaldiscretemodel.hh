#ifndef DUNE_FEM_DG_PRIMALDISCRETEMODELS_HH
#define DUNE_FEM_DG_PRIMALDISCRETEMODELS_HH

// Dune-Fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/operator/dg/discretemodelcommon.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/fluxes.hh>

#include <dune/fem-dg/operator/fluxes/diffusion/fluxes.hh>

namespace Dune
{
namespace Fem
{

  // AdvectionDiffusionDGPrimalModelBase
  //------------------------------------

  /**
   *  \brief discrete model for advection diffusion operator
   *
   *  \ingroup DiscreteModels
   *  \ingroup PassBased
   */
  template< class TraitsImp >
  class AdvectionDiffusionDGPrimalModelBase :
    public TraitsImp::AdvectionModelType
  {
  public:
    typedef TraitsImp                               Traits;
    typedef typename Traits::AdvectionModelType     BaseType;

    using BaseType::uVar;
    using BaseType::inside;
    using BaseType::outside;
    using BaseType::model_;
    using BaseType::uBnd_;

  public:
    static const bool advection = Traits::advection; // true if advection is enabled
    static const bool diffusion = Traits::diffusion; // this should be disabled for LDG

    enum { passUId = Traits::passUId };

    typedef typename BaseType::DomainType         DomainType;

    typedef typename BaseType::ModelType          ModelType ;
    typedef typename BaseType::AdvectionFluxType  AdvectionFluxType ;

#if defined TESTOPERATOR
#warning NO MASSOPERATOR APPLIED
    enum { ApplyInverseMassOperator = false };
#else
    enum { ApplyInverseMassOperator = true };
#endif

    typedef typename BaseType::GridPartType                            GridPartType;
    typedef typename BaseType::GridType                                GridType;
    typedef typename BaseType::IntersectionIteratorType                IntersectionIteratorType;
    typedef typename BaseType::IntersectionType                        IntersectionType;
    typedef typename BaseType::EntityType                              EntityType;
    typedef typename BaseType::RangeFieldType                          RangeFieldType;
    typedef typename BaseType::DomainFieldType                         DomainFieldType;
    typedef typename BaseType::RangeType                               RangeType;
    typedef typename BaseType::JacobianRangeType                       JacobianRangeType;

    typedef typename BaseType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // type of diffusion flux implementation
    typedef typename Traits::DiffusionFluxType                           DiffusionFluxType;

    enum { evaluateJacobian = DiffusionFluxType::evaluateJacobian }; // we need to evaluate jacobians here
  public:
    /**
     * \brief constructor
     */
    AdvectionDiffusionDGPrimalModelBase(const ModelType& mod,
                                        const AdvectionFluxType& numf,
                                        const DiffusionFluxType& diffflux )
      : BaseType( mod, numf ),
        diffFlux_( diffflux )
    {
    }

    //! copy constructor (for thread parallel progs mainly)
    AdvectionDiffusionDGPrimalModelBase(const AdvectionDiffusionDGPrimalModelBase& other )
      : BaseType( other ),
        diffFlux_( other.diffFlux_ )
    {
    }

    //! return true if source term is enabled (this might be overloaded in a derived class)
    inline bool hasSource() const
    {
      if( advection == diffusion )
        return model_.hasStiffSource() || model_.hasNonStiffSource() ;
      else if ( advection && ! diffusion )
        return model_.hasNonStiffSource() ;
      else
        return model_.hasStiffSource() ;
    }

    inline bool hasFlux() const
    {
      return (advection || diffusion) && model_.hasFlux();
    }

    /**
     * \brief analytical flux function
     */
    template <class LocalEvaluation>
    double source( const LocalEvaluation& local,
                   RangeType& s ) const
    {
      s = 0;

      double dtEst = std::numeric_limits< double >::max();

      if ( model_.hasStiffSource() )
      {
        const double dtStiff = model_.stiffSource( local, local.values()[uVar], local.jacobians()[uVar], s );
        dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
      }

      if ( model_.hasNonStiffSource() )
      {
        RangeType sNonStiff (0);
        const double dtNon = model_.nonStiffSource( local, local.values()[uVar], local.jacobians()[uVar], sNonStiff );

        // add to source
        s += sNonStiff;

        dtEst = ( dtNon > 0 ) ? std::min( dtEst, dtNon ) : dtEst;
      }

      // return the fastest wave from source terms
      return dtEst;
    }

    void switchUpwind() const
    {
      // reset max time steps
      BaseType::switchUpwind() ;
      // switch upwind direction if necessary
      diffFlux_.switchUpwind();
    }

    template <class LocalEvaluationVec >
    void initializeIntersection( const LocalEvaluationVec& left,
                                 const LocalEvaluationVec& right )
    {
      if( diffusion )
      {
        // call diffusion flux
        diffFlux_.initializeIntersection( left[uVar], right[uVar], left[uVar].values(), right[uVar].values() );
      }
    }

    template <class LocalEvaluationVec >
    void initializeBoundary(const LocalEvaluationVec& local )
    {
      if( diffusion )
      {
        if( diffFlux_.hasLifting() )
        {
          // use unsigned int, otherwise operator [] on LocalEvaluation has
          // different meaning
          unsigned int qp = 0;
          const bool hasBoundaryValue = model_.hasBoundaryValue( local[ qp ] );
          const bool hasRobinBoundaryValue = model_.hasRobinBoundaryValue( local[ qp ] );

          const unsigned int quadNop = local.quadrature().nop();
          if( uBndVec_.size() < quadNop )
            uBndVec_.resize( quadNop );

          for(; qp < quadNop; ++qp)
          {
            assert( hasBoundaryValue == model_.hasBoundaryValue( local[qp] ) );
            assert( hasRobinBoundaryValue == model_.hasRobinBoundaryValue( local[qp] ) );

            if( hasBoundaryValue || hasRobinBoundaryValue )
              model_.boundaryValue(local[qp], local[qp].values()[uVar], uBndVec_[qp] );
            else
              // do something bad to uBndVec as it shouldn't be used
              uBndVec_[qp] = std::numeric_limits< double >::quiet_NaN();
          }

          // call diffusion flux
          diffFlux_.initializeBoundary( local[uVar], local[uVar].values(), uBndVec_ );
        }
      }
    }

  public:
    /**
     * \brief flux function on interfaces between cells
     *
     * \param left local evaluation
     * \param right local evaluation
     * \param gLeft result for this side of \c it
     * \param gRight result for the other side of \c it
     * \return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      /*****************************
       * Advection                 *
       ****************************/
      double advectionWaveSpeed = BaseType ::
        numericalFlux( left, right,
                       gLeft, gRight, gDiffLeft, gDiffRight );

      double diffusionWaveSpeed = 0.0;
      if( diffusion )
      {
        RangeType dLeft, dRight;

        diffusionWaveSpeed =
          diffFlux_.numericalFlux(left, right,
                                  left.values()[uVar], right.values()[uVar],
                                  left.jacobians()[uVar], right.jacobians()[uVar],
                                  dLeft, dRight,
                                  gDiffLeft, gDiffRight);

        gLeft  += dLeft;
        gRight += dRight;
      }
      else
      {
        gDiffLeft  = 0;
        gDiffRight = 0;
      }

      return advectionWaveSpeed + diffusionWaveSpeed;
    }


    /**
     * \brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      /****************************/
      /* Advection                *
       ****************************/
      double advectionWaveSpeed = BaseType::boundaryFlux( left, gLeft, gDiffLeft );

      /****************************/
      /* Diffusion                 *
       ****************************/
      double diffusionWaveSpeed = 0.0;

      const bool hasBoundaryValue = model_.hasBoundaryValue( left );
      const bool hasRobinBoundaryValue = model_.hasRobinBoundaryValue( left );


      if( diffusion )
      {
        if( hasBoundaryValue || hasRobinBoundaryValue)
        {
          // diffusion boundary flux for Dirichlet boundaries
          RangeType dLeft ( 0 );
          diffusionWaveSpeed =
            diffFlux_.boundaryFlux(left[uVar],
                                   left[uVar].values(), uBnd_, // is set during call of  BaseType::boundaryFlux
                                   left[uVar].jacobians(),
                                   dLeft,
                                   gDiffLeft);
          gLeft += dLeft;
        }
        if( !hasBoundaryValue )
        {
          RangeType diffBndFlux ( 0 );
          model_.diffusionBoundaryFlux( left,
                                        left.values()[uVar], left.jacobians()[uVar], diffBndFlux );
          gLeft += diffBndFlux;
        }

      }
      else
        gDiffLeft = 0;

      return advectionWaveSpeed + diffusionWaveSpeed;
    }

    /**
     * \brief analytical flux function
     */
    template <class LocalEvaluation>
    void analyticalFlux( const LocalEvaluation& local,
                         JacobianRangeType& f ) const
    {
      /*****************************
       * Advection                 *
       ****************************/
      BaseType::analyticalFlux( local, f );

      /*****************************
       * Diffusion                 *
       ****************************/
      if( diffusion && hasFlux() )
      {
        JacobianRangeType diffmatrix;
        model_.diffusion( local,
                          local.values()[ uVar ], local.jacobians()[ uVar ], diffmatrix);
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
  };



  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  // AdvectionDiffusionDGPrimalModel
  //--------------------------------

  //! \brief forward declaration
  template <class OpTraits,
            bool enableAdvection, bool enableDiffusion,
            int pUId, int... passIds >
  class AdvectionDiffusionDGPrimalModel;


  // AdvectionDiffusionDGPrimalTraits
  //---------------------------------

  template <class OpTraits,
            bool enableAdvection, bool enableDiffusion,
            int pUId, int... passIds >
  struct AdvectionDiffusionDGPrimalTraits
   : public AdvectionTraits< OpTraits, enableAdvection, pUId >
  {
    typedef typename OpTraits::DiffusionFluxType               DiffusionFluxType;

    enum { advection = enableAdvection };
    enum { diffusion = enableDiffusion };

    enum { passUId = pUId };

    // type of base class
    typedef AdvectionModel< OpTraits, enableAdvection, passUId >     AdvectionModelType ;

    // type of my discrete model
    typedef AdvectionDiffusionDGPrimalModel< OpTraits, enableAdvection, enableDiffusion, passUId > DGDiscreteModelType;
  };

  // AdvectionDiffusionDGPrimalModel
  //--------------------------------
  /**
   *  \brief discrete model for advection diffusion operator
   *
   *  \ingroup DiscreteModels
   *  \ingroup PassBased
   */
  template< class OpTraits,
            bool enableAdvection, bool enableDiffusion,
            int passUId, int... passIds >
  class AdvectionDiffusionDGPrimalModel :
    public AdvectionDiffusionDGPrimalModelBase<
              AdvectionDiffusionDGPrimalTraits < OpTraits, enableAdvection, enableDiffusion, passUId > >
  {
  public:
    typedef AdvectionDiffusionDGPrimalTraits
      < OpTraits, enableAdvection, enableDiffusion, passUId >   Traits;

    typedef AdvectionDiffusionDGPrimalModelBase< Traits > BaseType ;
    typedef typename BaseType::DiffusionFluxType  DiffusionFluxType ;
    typedef typename BaseType::ModelType          ModelType ;
    typedef typename BaseType::AdvectionFluxType  AdvectionFluxType ;

    // constructor: do not touch/delegate everything
    template< class ... Args>
    AdvectionDiffusionDGPrimalModel( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}
  };


  ///////////////////////////////////////////////////////////////////////
  //
  // AdvectionDiffusionDGPrimalTraits
  //
  //////////////////////////////////////////////////////////////////////

  // AdvectionDiffusionDGPrimalModel
  //--------------------------------

  //! \brief forward declaration
  template <class OpTraits,
            bool enableAdvection, bool enableDiffusion,
            int passUId, int... passIds >
  class AdaptiveAdvectionDiffusionDGPrimalModel;


  // AdvectionDiffusionDGPrimalTraits
  //---------------------------------
  /**
   *  \brief discrete model for advection diffusion operator
   *
   *  \ingroup DiscreteModels
   *  \ingroup PassBased
   */
  template <class OpTraits,
            bool enableAdvection, bool enableDiffusion,
            int passUId, int... passIds >
  struct AdaptiveAdvectionDiffusionDGPrimalTraits
   : public AdvectionDiffusionDGPrimalTraits
        < OpTraits, enableAdvection, enableDiffusion, passUId, passIds... >
  {
    // type of base class
    typedef AdaptiveAdvectionModel< OpTraits, enableAdvection, passUId >     AdvectionModelType ;

    // type of my discrete model
    typedef AdaptiveAdvectionDiffusionDGPrimalModel
      < OpTraits, enableAdvection, enableDiffusion, passUId >                 DGDiscreteModelType;
  };

  // AdaptiveAdvectionDiffusionDGPrimalModel
  //----------------------------------------
  /**
   *  \brief discrete model for advection diffusion operator
   *
   *  \ingroup DiscreteModels
   *  \ingroup PassBased
   */
  template< class OpTraits,
            bool enableAdvection, bool enableDiffusion,
            int passUId, int... passIds >
  class AdaptiveAdvectionDiffusionDGPrimalModel :
    public AdvectionDiffusionDGPrimalModelBase<
       AdaptiveAdvectionDiffusionDGPrimalTraits< OpTraits, enableAdvection, enableDiffusion, passUId > >
  {
  public:
    typedef AdaptiveAdvectionDiffusionDGPrimalTraits
      < OpTraits, enableAdvection, enableDiffusion, passUId >   Traits;

    typedef AdvectionDiffusionDGPrimalModelBase< Traits >  BaseType;

    using BaseType::uVar;
    using BaseType::inside;
    using BaseType::outside;
    using BaseType::model_;
    using BaseType::uBnd_;

  public:
    static const bool advection = enableAdvection; // true if advection is enabled
    static const bool diffusion = enableDiffusion; // this should be disabled for LDG

    typedef typename BaseType::IntersectionType    IntersectionType;
    typedef typename BaseType::EntityType          EntityType;
    typedef typename BaseType::RangeType           RangeType;
    typedef typename BaseType::JacobianRangeType   JacobianRangeType;
    typedef typename BaseType::DomainType          DomainType;

    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // type of diffusion flux implementation
    typedef typename BaseType::AdvectionFluxType  AdvectionFluxType;
    typedef typename BaseType::DiffusionFluxType  DiffusionFluxType;
    typedef typename BaseType::ModelType          ModelType;

    enum { evaluateJacobian = DiffusionFluxType::evaluateJacobian  }; // we need to evaluate jacobians here
  public:
    // constructor: do not touch/delegate everything
    template< class ... Args>
    AdaptiveAdvectionDiffusionDGPrimalModel( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

  public:
    /**
     * \brief flux function on interfaces between cells
     *
     * \param[in]  left local evaluation context of inside cell
     * \param[in]  right local evaluation context of outside cell
     * \param[out] uLeft DOF evaluation on this side of \c it
     * \param[out] uRight DOF evaluation on the other side of \c it
     * \param[out] gLeft result for this side of \c it
     * \param[out] gRight result for the other side of \c it
     * \return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      return 0.0;
#if 0
      // do not take into account part of the domain which is not
      // of computation significance (i.e. a damping layer)
      if( ! model_.allowsRefinement( left ) )
        return 0.;

      // call numerical flux of base type
      const double ldt = BaseType::numericalFlux( it, time, faceQuadInner, faceQuadOuter, quadPoint,
                                              uLeft, uRight, jacLeft, jacRight,
                                              gLeft, gRight, gDiffLeft, gDiffRight );

      if( diffusion && adaptation_ )
      {
        JacobianRangeType diffmatrix;
        RangeType error;
        DomainType normal = it.integrationOuterNormal( faceQuadInner.localPoint( quadPoint ) );

        // G(u-)grad(u-) for multiplication with phi
        // call on inside
        model_.diffusion( inside(), time, faceQuadInner.point( quadPoint ),
                          uLeft[uVar], jacLeft[uVar], diffmatrix );

        diffmatrix.mv( normal, error );

        // G(u+)grad(u+) for multiplication with phi
        // call on outside
        model_.diffusion( outside(), time, faceQuadOuter.point( quadPoint ),
                          uRight[uVar], jacRight[uVar], diffmatrix );

        diffmatrix.mmv( normal, error );

        error *= 0.5;

        // calculate grid width
        double weight = weight_ *
          (0.5 * ( this->enVolume() + this->nbVolume() ) / normal.two_norm() );

        // add error to indicator
        enIndicator_.add( error, weight );
        nbIndicator_.addChecked( error, weight );
      }

      return ldt ;
#endif
    }


    /**
     * \brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      const double ldt = BaseType::boundaryFlux( left, gLeft, gDiffLeft );

      if( diffusion && adaptation_ )
      {
        // implement error indicator here
      }

      return ldt ;
    }

  protected:
    // defined in AdvectionModel
    using BaseType::adaptation_ ;
    using BaseType::enIndicator_;
    using BaseType::nbIndicator_;
    using BaseType::weight_ ;
  };

} // end of namespace
} // end of namespace

#endif
