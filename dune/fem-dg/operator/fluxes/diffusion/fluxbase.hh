#ifndef DUNE_FEM_DG_DIFFUSIONFLUXES_HH
#define DUNE_FEM_DG_DIFFUSIONFLUXES_HH

#include <type_traits>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>

#include "parameters.hh"

namespace Dune
{
namespace Fem
{



  /**
   * \brief Base class for all diffusion fluxes.
   *
   * \ingroup DiffusionFluxes
   */
  template <class DiscreteFunctionSpaceImp,
            class Model,
            class FluxParameterImp >
  class DGDiffusionFluxBase
  {
  public:
    static const bool evaluateJacobian = false;
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    typedef FluxParameterImp                                         ParameterType;
    typedef typename ParameterType::IdEnum                           IdEnum;

    static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;
    static const int dimRange  = DiscreteFunctionSpaceType::dimRange;

    typedef typename DiscreteFunctionSpaceType::DomainType           DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType       RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeType            RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType    JacobianRangeType;

    typedef FieldVector< DomainFieldType, dimDomain-1 >              FaceDomainType;

    typedef typename DiscreteFunctionSpaceType::GridPartType         GridPartType;
    typedef typename GridPartType::IntersectionIteratorType          IntersectionIterator;
    typedef typename GridPartType::IntersectionType                  IntersectionType;
    typedef typename GridPartType::GridType                          GridType;
    typedef typename DiscreteFunctionSpaceType::EntityType           EntityType;
    static const int dimGradRange = dimDomain * dimRange;
    static const int polOrd = DiscreteFunctionSpaceType::polynomialOrder;

    // type of gradient space
    typedef typename DiscreteFunctionSpaceType::
        template ToNewDimRange< dimGradRange >::Type                 DiscreteGradientSpaceType;

    typedef typename DiscreteGradientSpaceType::RangeType            GradientType;

    typedef Fem::TemporaryLocalFunction< DiscreteGradientSpaceType > LiftingFunctionType;

    typedef Fem::CachingQuadrature< GridPartType, 0>                 VolumeQuadratureType ;

    DomainType upwindDefault() const
    {
      DomainType upwind ( M_PI );
      // set upwind to some strange numbers
      if( dimDomain > 1 ) upwind[1] = M_LN2 ;
      if( dimDomain > 2 ) upwind[2] = M_E ;
      return upwind ;
    }
  public:
    /**
     * \brief constructor
     */
    DGDiffusionFluxBase( const Model& mod,
                         const bool initUpwind,
                         const ParameterType& parameters )
    : model_(mod),
      param_( parameters ),
      upwind_( upwindDefault() ),
      // Set CFL number for penalty term (compare diffusion in first pass)
      cflDiffinv_(8. * (polOrd+1.) ),
      dimensionFactor_( 2.0 * ( dimDomain ) ),
      nonconformingFactor_( 2.0 * ( dimDomain - 1 ) )
    {
      if( initUpwind )
      {
        parameter().upwind( upwind_ );
        if( Fem::Parameter :: verbose() )
          std::cout << "Using upwind = " << upwind_ << std::endl;
      }
    }

    //! copy constructor
    DGDiffusionFluxBase( const DGDiffusionFluxBase& other )
    : model_( other.model_ ),
      param_( other.parameter() ),
      upwind_( other.upwind_ ),
      cflDiffinv_( other.cflDiffinv_ ),
      dimensionFactor_( other.dimensionFactor_ ),
      nonconformingFactor_( other.nonconformingFactor_ )
    {
    }

    //! returns true if lifting has to be calculated
    //const bool hasLifting () const { return false; }

  protected:
    bool determineDirection( const bool areaSwitch, const double enVolume, const double nbVolume,
                             const IntersectionType& intersection ) const
    {
      if (areaSwitch && std::abs( enVolume - nbVolume ) > 1e-8)
      {
        // we need en=K^- and nb=K^+ such that c = (|K^-|/|K^+|) <= 1
        // outside is K^- if outside volume is smaller
        return ( enVolume > nbVolume );
      }
      else
      {
        Fem::Quadrature< DomainFieldType, dimDomain - 1 > quad( intersection.type() , 0 );
        const DomainType normal = intersection.outerNormal( quad.point( 0 ) );
        return determineDirection( normal );
      }
    }

    bool determineDirection( const DomainType& normal ) const
    {
       return ( normal * upwind_ ) < 0 ;
    }

  public:
    void switchUpwind()
    {
      upwind_ *= -1.0;
    }

    bool hasLifting() const { return false; }

    void diffusionFluxName( std::ostream& out ) const
    {
      out << ParameterType::methodNames( param_.getMethod() );
    }

    void diffusionFluxPenalty( std::ostream& out ) const
    {
      out << param_.penalty();
    }

    void diffusionFluxLiftFactor( std::ostream& out ) const
    {
    }


    template <class LocalEvaluation, class ArgumentTupleVector >
    void initializeIntersection(const LocalEvaluation& left,
                                const LocalEvaluation& right,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec)
    {
    }

    template <class LocalEvaluation, class ArgumentTupleVector>
    void initializeBoundary(const LocalEvaluation& local,
                            const ArgumentTupleVector& uLeftVec,
                            const std::vector< RangeType >& uRight)
    {
    }

    //template <class QuadratureImp, class ArgumentTupleVector >
    //void initializeIntersection(const IntersectionType& intersection,
    //                            const EntityType& inside,
    //                            const EntityType& outside,
    //                            const double time,
    //                            const QuadratureImp& quadInner,
    //                            const QuadratureImp& quadOuter,
    //                            const ArgumentTupleVector& uLeftVec,
    //                            const ArgumentTupleVector& uRightVec)
    //{
    //}

    //template <class QuadratureImp, class ArgumentTupleVector>
    //void initializeBoundary(const IntersectionType& intersection,
    //                        const EntityType& entity,
    //                        const double time,
    //                        const QuadratureImp& quadInner,
    //                        const ArgumentTupleVector& uLeftVec,
    //                        const std::vector< RangeType >& uRight)
    //{
    //}

  public:
    /**
     * \brief flux function on interfaces between cells for advection and diffusion
     *
     * \param[in] right local evaluation
     * \param[in] left local evaluation
     * \param[in] uLeft DOF evaluation on this side of \c it
     * \param[in] uRight DOF evaluation on the other side of \c it
     * \param[out] gLeft num. flux projected on normal on this side
     * \           of \c it for multiplication with \f$ \phi \f$
     * \param[out] gRight advection flux projected on normal for the other side
     * \           of \c it for multiplication with \f$ \phi \f$
     * \param[out] gDiffLeft num. flux projected on normal on this side
     * \           of \c it for multiplication with \f$ \nabla\phi \f$
     * \param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * \return wave speed estimate (multiplied with the integration element of the intersection),
     *              to estimate the time step |T|/wave.
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const JacobianRangeType& jacLeft,
                         const JacobianRangeType& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight )
    {
      assert( false );
      abort();
      return 0.0;
    }


    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        const RangeType& uLeft,
                        const RangeType& uRight,
                        const JacobianRangeType& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft )
    {
      assert( false );
      abort();
      return 0.0;
    }

    const Model &model () const { return model_; }
    const ParameterType& parameter() const { return param_; }

  protected:
    const Model   &model_;
    const ParameterType& param_;
    DomainType upwind_;
    const double cflDiffinv_;
    const double dimensionFactor_;
    const double nonconformingFactor_;
  }; // end DGPrimalDiffusionFlux


  /**
   * \brief A diffusion flux for the LDG scheme.
   *
   * \ingroup DiffusionFluxes
   */
  template <class DiscreteFunctionSpaceImp,
            class ModelImp,
            class FluxParameterImp >
  class LDGDiffusionFluxBase :
    public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, ModelImp, FluxParameterImp >
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, ModelImp, FluxParameterImp > BaseType;
  public:
    typedef DiscreteFunctionSpaceImp                DiscreteFunctionSpaceType;

    static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;
    static const int dimRange  = DiscreteFunctionSpaceType::dimRange;

    typedef typename BaseType::DomainType               DomainType;
    typedef typename BaseType::RangeFieldType           RangeFieldType;
    typedef typename BaseType::DomainFieldType          DomainFieldType;
    typedef typename BaseType::RangeType                RangeType;
    typedef typename BaseType::JacobianRangeType        JacobianRangeType;

    typedef typename BaseType::FaceDomainType           FaceDomainType;

    typedef typename BaseType::GridPartType             GridPartType;
    typedef typename BaseType::IntersectionIterator     IntersectionIterator;
    typedef typename BaseType::IntersectionType         IntersectionType;
    typedef typename BaseType::GridType                 GridType;
    typedef typename BaseType::EntityType               EntityType;
    static const int dimGradRange = dimDomain * dimRange;
    static const int polOrd = DiscreteFunctionSpaceType::polynomialOrder;

    typedef typename BaseType::ParameterType             ParameterType;

    // type of gradient space
    typedef typename DiscreteFunctionSpaceType::
      template ToNewDimRange< dimGradRange >::Type       DiscreteGradientSpaceType;

    typedef typename DiscreteGradientSpaceType::RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType::JacobianRangeType GradientJacobianType;

    // jacobians of the functions do not have to be evaluated for this flux
    static const bool evaluateJacobian = false;

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
     LDGDiffusionFluxBase(GridPartType& gridPart,
                          const ModelImp& mod,
                          const ParameterType& param ) :
      BaseType( mod, true, param )
    {}

    //! copy constructor
    LDGDiffusionFluxBase( const LDGDiffusionFluxBase& other )
    : BaseType( other )
    {}

    //! returns true if lifting has to be calculated
    const bool hasLifting () const { return false; }

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
      assert( false );
      abort();
      return 0.0;
    }


    template <class LocalEvaluation>
    double gradientBoundaryFlux(const LocalEvaluation& left,
                                const RangeType& uLeft,
                                const RangeType& uBnd,
                                GradientRangeType& gLeft,
                                GradientJacobianType& gDiffLeft) const
    {
      assert( false );
      abort();
      return 0.0;
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
      assert( false );
      abort();
      return 0.0;
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
      assert( false );
      abort();
      return 0.0;
    }
  }; // end LDGAverageDiffusionFlux


} // end namespace
} // end namespace
#endif
