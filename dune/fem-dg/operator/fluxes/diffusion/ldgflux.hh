#ifndef DUNE_FEM_DG_LDGAVERAGEFLUX_HH
#define DUNE_FEM_DG_LDGAVERAGEFLUX_HH

// local includes
#include "fluxbase.hh"

namespace Dune
{
namespace Fem
{

  /**
   * \brief Implementation of Bassi-Rebay 1 and Local DG fluxes.
   *
   * \ingroup DiffusionFluxes
   */
  template <class DiscreteFunctionSpaceImp,
            class ModelImp,
            class FluxParameterImp>
  class LDGDiffusionFluxImpl :
    public LDGDiffusionFluxBase< DiscreteFunctionSpaceImp, ModelImp, FluxParameterImp >
  {
    typedef LDGDiffusionFluxBase< DiscreteFunctionSpaceImp, ModelImp, FluxParameterImp > BaseType;
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    static const int dimDomain = DiscreteFunctionSpaceType :: dimDomain;
    static const int dimRange  = DiscreteFunctionSpaceType :: dimRange;

    typedef typename DiscreteFunctionSpaceType :: DomainType           DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType       RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeType            RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType    JacobianRangeType;

    typedef FieldVector< DomainFieldType, dimDomain-1 > FaceDomainType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType         GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename GridPartType :: IntersectionType                  Intersection;
    typedef typename GridPartType :: GridType                          GridType;
    typedef typename DiscreteFunctionSpaceType :: EntityType           EntityType;
    static const int dimGradRange = dimDomain * dimRange;
    static const int polOrd = DiscreteFunctionSpaceType :: polynomialOrder;

    typedef typename BaseType :: ParameterType  ParameterType;

    // type of gradient space
    typedef typename DiscreteFunctionSpaceType ::
        template ToNewDimRange< dimGradRange > :: Type   DiscreteGradientSpaceType;

    typedef typename DiscreteGradientSpaceType :: RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType :: JacobianRangeType GradientJacobianType;

    // jacobians of the functions do not have to be evaluated for this flux
    static const bool evaluateJacobian = false;

  protected:
    typedef typename BaseType::IdEnum         EnumType;

    using BaseType::determineDirection;
    using BaseType::model_;
    using BaseType::cflDiffinv_;
    using BaseType::numericalFlux ;
    using BaseType::parameter ;
    using BaseType::nonconformingFactor_;
    using BaseType::dimensionFactor_;

  public:
    /**
     * \brief constructor
     */
    LDGDiffusionFluxImpl(GridPartType& gridPart,
                         const ModelImp& mod,
                         const ParameterType& param,
                         const EnumType staticMethod ) :
      BaseType( gridPart, mod, param ),
      method_( staticMethod == EnumType::local ? param.getMethod() : staticMethod ),
      penalty_( parameter().penalty() ),
      // Set CFL number for penalty term (compare diffusion in first pass)
      penaltyTerm_( std::abs(  penalty_ ) > 0 )
    {
      if( Fem::Parameter::verbose () )
      {
        std::cout << "LDGDiffusionFluxImpl: penalty = " << penalty_ << std::endl;
      }
    }

    LDGDiffusionFluxImpl(const LDGDiffusionFluxImpl& other)
      : BaseType( other ),
        method_( other.method_ ),
        penalty_( other.penalty_ ),
        penaltyTerm_( other.penaltyTerm_ )
    {}

    //! returns true if lifting has to be calculated
    const bool hasLifting () const { return false; }

  protected:
    static const bool realLDG = true;
    double theta( const DomainType& normal ) const
    {
      if( method_ == EnumType::ldg )
      {
        // LDG thete is 1 or 0
        if( determineDirection( normal ) )
          return 1.0;
        else
          return 0.0;
      }
      else
      {
        assert( method_ == EnumType::br1 );
        // Average fluxes (Bassi-Rebay 1)
        return 0.5;
      }
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
      const double thetaLeft  = theta( normal );
      const double thetaRight = 1.0 - thetaLeft;

      GradientJacobianType diffmatrix;

      double diffTimeStep = 0.0;

      // select left or right or average
      if( thetaLeft > 0 )
      {
        //TODO use diffusionTimeStep
        diffTimeStep = 0;

        /* central differences (might be suboptimal) */
        model_.jacobian(left, uLeft, diffmatrix );

        diffmatrix.mv(normal, gLeft );
        gLeft *= thetaLeft ;
      }
      else
        gLeft = 0;

      if( thetaRight > 0 )
      {
        //const double diffStepRight = 0;

        // right jacobian
        model_.jacobian( right, uRight, diffmatrix );

        diffmatrix.mv(normal, gRight);

        // add to flux
        gLeft.axpy( thetaRight, gRight );

        // diffTimeStep = std::max( diffTimeStep, diffStepRight );
      }

      // copy flux
      gRight = gLeft;

#ifndef NDEBUG
      gDiffLeft = 0;
      gDiffRight = 0;
#endif

      // upper bound for the next time step length
      return diffTimeStep; // * cflDiffinv_;
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
      const double thetaLeft  = theta( normal );
      const double thetaRight = 1.0 - thetaLeft;

      GradientJacobianType diffmatrix;

      // calculate uVal
      RangeType uVal( 0 );

      // select u from uLeft and uRight
      if( thetaLeft > 0 )
        uVal.axpy( thetaLeft , uLeft );
      if( thetaRight > 0 )
        uVal.axpy( thetaRight, uBnd );

      // compute jacobian of u
      model_.jacobian(left, uVal, diffmatrix );

      diffmatrix.mv(normal, gLeft);

#ifndef NDEBUG
      gDiffLeft = 0;
#endif

      //TODO use diffusionTimeStep
      const double diffTimeStep = 0; // * cflDiffinv_;
      return diffTimeStep;
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
                         const JacobianRangeType& sigmaLeft,
                         const JacobianRangeType& sigmaRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft, // not used here (only for primal passes)
                         JacobianRangeType& gDiffRight )
    {
      const DomainType normal = left.intersection().integrationOuterNormal( left.localPosition() );
      const double faceLengthSqr = normal.two_norm2();

      /**********************************
       * Diffusion sigma Flux (Pass 2)  *
       **********************************/
      JacobianRangeType diffmatrix;

      // diffusion for uLeft
      model_.diffusion( left, uLeft, sigmaLeft, diffmatrix);

      RangeType diffflux;
      diffmatrix.mv(normal, diffflux);

      // diffusion for uRight
      model_.diffusion( right, uRight, sigmaRight, diffmatrix);
      // add to diffflux
      diffmatrix.umv(normal, diffflux);
      diffflux *= 0.5;

      const double faceVolumeEstimate = dimensionFactor_ *
        (left.intersection().conforming() ? faceLengthSqr
          : (nonconformingFactor_ * faceLengthSqr));

      const double diffTimeLeft =
        model_.diffusionTimeStep( left, faceVolumeEstimate, uLeft );

      const double diffTimeRight =
        model_.diffusionTimeStep( right, faceVolumeEstimate, uRight );

      double diffTimeStep = std::max( diffTimeLeft, diffTimeRight );

      // add penalty factor
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
                        const JacobianRangeType& sigmaLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft )
    {
      // get local point
      const DomainType normal = left.intersection().integrationOuterNormal( left.localPosition() );
      const double faceLengthSqr = normal.two_norm2();

      /****************************/
      /* Diffusion (Pass 2)       */
      /****************************/
      JacobianRangeType diffmatrix;

      model_.diffusion(left, uLeft, sigmaLeft, diffmatrix);

      diffmatrix.mv(normal, gLeft);

      const double faceVolumeEstimate = dimensionFactor_ * faceLengthSqr;

      double diffTimeStep =
        model_.diffusionTimeStep( left, faceVolumeEstimate, uLeft );

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
    const EnumType method_;
    const double penalty_;
    const bool penaltyTerm_;

  }; // end LDGDiffusionFluxImpl

} // end namespace
} // end namespace
#endif
