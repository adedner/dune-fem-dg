#ifndef DUNE_FEM_DG_LDGAVERAGEFLUX_HH
#define DUNE_FEM_DG_LDGAVERAGEFLUX_HH

// Dune-Fem includes
#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief A diffusion flux for the LDG scheme.
   *
   * \ingroup DiffusionFluxes
   */
  template <class DiscreteFunctionSpaceImp,
            class ModelImp,
            class FluxParameterImp >
  class LDGAverageDiffusionFlux :
    public LDGDiffusionFluxBase< DiscreteFunctionSpaceImp, ModelImp, FluxParameterImp >
  {
    typedef LDGDiffusionFluxBase< DiscreteFunctionSpaceImp, ModelImp, FluxParameterImp > BaseType;
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    enum { dimRange  = DiscreteFunctionSpaceType :: dimRange };

    typedef typename DiscreteFunctionSpaceType :: DomainType           DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType       RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeType            RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType    JacobianRangeType;

    typedef FieldVector< DomainFieldType, dimDomain-1 > FaceDomainType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType         GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename GridPartType :: GridType                          GridType;
    typedef typename DiscreteFunctionSpaceType :: EntityType           EntityType;
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType :: polynomialOrder };

    typedef typename BaseType :: ParameterType  ParameterType;

    // type of gradient space
    typedef typename DiscreteFunctionSpaceType ::
        template ToNewDimRange< dimGradRange > :: Type   DiscreteGradientSpaceType;

    typedef typename DiscreteGradientSpaceType :: RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType :: JacobianRangeType GradientJacobianType;

    // jacobians of the functions do not have to be evaluated for this flux
    enum { evaluateJacobian = false };

  protected:
    using BaseType::determineDirection;
    using BaseType::model_;
    using BaseType::cflDiffinv_;
    using BaseType::numericalFlux ;
    using BaseType::parameter ;


  public:
    /**
     * \brief constructor
     */
    LDGAverageDiffusionFlux(GridPartType& gridPart,
                            const ModelImp& mod,
                            const ParameterType& param ) :
      BaseType( gridPart, mod, param ),
      penalty_( parameter().penalty() ),
      // Set CFL number for penalty term (compare diffusion in first pass)
      penaltyTerm_( std::abs(  penalty_ ) > 0 )
    {
      if( Fem::Parameter::verbose () )
      {
        std::cout << "LDGAverageDiffusionFlux: penalty = " << penalty_ << std::endl;
      }
    }

    LDGAverageDiffusionFlux(const LDGAverageDiffusionFlux& other)
      : BaseType( other ),
        penalty_( other.penalty_ ),
        penaltyTerm_( other.penaltyTerm_ )
    {}

    //! returns true if lifting has to be calculated
    const bool hasLifting () const { return false; }

  protected:
    enum { realLDG = true };
    double theta( const Intersection& intersection ) const
    {
      if( realLDG )
      {
        // LDG thete is 1 or 0
        if( determineDirection( intersection ) )
          return 1.0;
        else
          return 0.0;
      }
      else
        // Average fluxes
        return 0.5;
    }

  public:
    /**
     * \brief flux function on interfaces between cells
     *
     * \param left local evaluation
     * \param right local evaluation
     * \param uLeft DOF evaluation on this side of \c intersection
     * \param uRight DOF evaluation on the other side of \c intersection
     * \param gLeft result for this side of \c intersection
     * \param gRight result for the other side of \c intersection
     * \return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class LocalEvaluation>
    double gradientNumericalFlux(const LocalEvaluation& left,
                                 const LocalEvaluation& right,
                                 const RangeType& uLeft,
                                 const RangeType& uRight,
                                 GradientRangeType& gLeft,
                                 GradientRangeType& gRight,
                                 GradientJacobianType& gDiffLeft,
                                 GradientJacobianType& gDiffRight) const
    {
      const DomainType normal = left.intersection().integrationOuterNormal( left.localPosition() );

      // get factor for each side
      const double thetaLeft  = theta( left.intersection() );
      const double thetaRight = 1.0 - thetaLeft;

      GradientJacobianType diffmatrix;

      double diffTimeStep = 0.0;

      if( thetaLeft > 0 )
      {
        //TODO use diffusionTimeStep
        diffTimeStep =
          /* central differences (might be suboptimal) */
          model_.diffusion(left, uLeft, diffmatrix );

        diffmatrix.mv(normal, gLeft );
        gLeft *= thetaLeft ;
      }
      else
        gLeft = 0;

      if( thetaRight > 0 )
      {
        const double diffStepRight =
          model_.diffusion( right, uRight, diffmatrix );

        diffmatrix.mv(normal, gRight);

        // add to flux
        gLeft.axpy( thetaRight, gRight );

        diffTimeStep = std::max( diffTimeStep, diffStepRight );
      }

      // copy flux
      gRight = gLeft;

#ifndef NDEBUG
      gDiffLeft = 0;
      gDiffRight = 0;
#endif

      // upper bound for the next time step length
      return diffTimeStep * cflDiffinv_;
    }


    template <class LocalEvaluation>
    double gradientBoundaryFlux(const LocalEvaluation& left,
                                const RangeType& uLeft,
                                const RangeType& uBnd,
                                GradientRangeType& gLeft,
                                GradientJacobianType& gDiffLeft) const
    {
      const DomainType normal = left.intersection().integrationOuterNormal( left.localPosition() );

      // get factor for each side
      const double thetaLeft  = theta( left.intersection() );
      const double thetaRight = 1.0 - thetaLeft;

      GradientJacobianType diffmatrix;

      // calculate uVal
      RangeType uVal( 0 );

      if( thetaLeft > 0 )
        uVal.axpy( thetaLeft , uLeft );
      if( thetaRight > 0 )
        uVal.axpy( thetaRight, uBnd );

      //TODO use diffusionTimeStep
      const double diffTimeStep =
            model_.diffusion(left, uVal, diffmatrix );

      diffmatrix.mv(normal, gLeft);

#ifndef NDEBUG
      gDiffLeft = 0;
#endif

      return diffTimeStep * cflDiffinv_;
    }


    /**
     * \brief flux function on interfaces between cells
     *
     * \param left local evaluation
     * \param right local evaluation
     * \param uLeft DOF evaluation on this side of \c intersection
     * \param uRight DOF evaluation on the other side of \c intersection
     * \param gLeft result for this side of \c intersection
     * \param gRight result for the other side of \c intersection
     * \return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const GradientRangeType& sigmaLeft,
                         const GradientRangeType& sigmaRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft, // not used here (only for primal passes)
                         JacobianRangeType& gDiffRight )
    {
      const DomainType normal = left.intersection().integrationOuterNormal( left.localPosition() );

     /**********************************
      * Diffusion sigma Flux (Pass 2)  *
      **********************************/
      JacobianRangeType diffmatrix;

      //std::cout << "uleft = " << uLeft << "  sigLeft " << sigmaLeft << std::endl;

      /* Central differences */
      const double diffTimeLeft =
        model_.diffusion( left, uLeft, sigmaLeft, diffmatrix);

      RangeType diffflux;
      diffmatrix.mv(normal, diffflux);

      const double diffTimeRight =
        model_.diffusion( right, uRight, sigmaRight, diffmatrix);
      diffmatrix.umv(normal, diffflux);
      diffflux *= 0.5;

      double diffTimeStep = std::max( diffTimeLeft, diffTimeRight );

      // add penalty term ( enVolume() is available since we derive from
      //    DiscreteModelDefaultWithInsideOutside)
      const double factor = penalty_ * diffTimeStep ;

      RangeType jump( uLeft );
      jump -= uRight;
      diffflux.axpy(factor, jump);

      gLeft  = diffflux;
      gRight = diffflux;

#ifndef NDEBUG
      gDiffLeft = 0;
      gDiffRight = 0;
#endif

      // timestep restict to diffusion timestep
      // WARNING: reconsider this
      diffTimeStep *= cflDiffinv_;
      return diffTimeStep;
    }


    /**
     * \brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        const RangeType& uLeft,
                        const RangeType& uRight,
                        const GradientRangeType& sigmaLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft )
    {
      // get local point
      const DomainType normal = left.intersection().integrationOuterNormal( left.localPosition() );

      /****************************/
      /* Diffusion (Pass 2)       */
      /****************************/
      JacobianRangeType diffmatrix;
      //TODO use diffusionTimeStep
      double diffTimeStep =
        model_.diffusion(left, uLeft, sigmaLeft, diffmatrix);

      diffmatrix.mv(normal, gLeft);

      // add penalty term
      const double factor = penalty_ * diffTimeStep;

      RangeType jump( uLeft );
      jump -= uRight;
      gLeft.axpy(factor,jump);

#ifndef NDEBUG
      gDiffLeft = 0;
#endif

      diffTimeStep *= cflDiffinv_;
      return diffTimeStep;
    }
  protected:
    const double penalty_;
    const bool penaltyTerm_;

  }; // end LDGAverageDiffusionFlux

} // end namespace
} // end namespace
#endif
