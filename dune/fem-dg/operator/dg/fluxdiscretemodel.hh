#ifndef DUNE_FEM_DG_FLUXDISCRETEMODEL_HH
#define DUNE_FEM_DG_FLUXDISCRETEMODEL_HH

// Dune-Fem includes
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/operator/dg/discretemodelcommon.hh>
#include <dune/fem-dg/operator/fluxes/ldgflux.hh>

//*************************************************************
namespace Dune
{
namespace Fem
{


  // GradientModel
  //--------------

  template <class Traits, int passUId>
  class GradientModel;


  // GradientTraits
  //---------------

  template <class Traits, int passUId >
  struct GradientTraits : public Traits
  {
    typedef GradientModel< Traits, passUId >         DGDiscreteModelType;
  };


  // GradientModel
  //--------------

  template < class OpTraits, int passUId>
  class GradientModel :
    public Fem::DGDiscreteModelDefaultWithInsideOutside
      < GradientTraits< OpTraits, passUId>, passUId >
  {
    typedef Fem::DGDiscreteModelDefaultWithInsideOutside
              < GradientTraits< OpTraits, passUId >,  passUId > BaseType;

    using BaseType :: inside;
    using BaseType :: outside;

    // This type definition allows a convenient access to arguments of passes.
    integral_constant< int, passUId > uVar;

  public:
    typedef GradientTraits< OpTraits, passUId >        Traits;
    typedef typename Traits :: ModelType                             ModelType;
    typedef typename Traits :: FluxType                              NumFluxType;

    typedef typename Traits :: DomainType                            DomainType;
    typedef typename Traits :: FaceDomainType                        FaceDomainType;
    typedef typename Traits :: RangeType                             RangeType;
    typedef typename Traits :: GridType                              GridType;
    typedef typename Traits :: JacobianRangeType                     JacobianRangeType;
    typedef typename Traits :: GridPartType
              :: IntersectionIteratorType                            IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection            Intersection;
    typedef typename BaseType :: EntityType                          EntityType;

    enum { evaluateJacobian = NumFluxType :: evaluateJacobian };

    // necessary for TESTOPERATOR
    // not sure how it works for dual operators!
    enum { ApplyInverseMassOperator = true };

  public:
    /**
     * \brief constructor
     */
    GradientModel(const ModelType& mod,
                  const NumFluxType& numf) :
      model_( mod ),
      gradientFlux_( numf ),
      cflDiffinv_( 2.0 * ( Traits::polynomialOrder + 1) )
    {
      #if defined TESTOPERATOR
        std::cerr <<"didn't test how to use TESTOPERATOR with dual formulation";
        abort();
      #endif
    }

    void setTime ( double time ) { const_cast< ModelType& >( model_ ).setTime( time ); }

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }

    template< class LocalEvaluation >
    inline double source( const LocalEvaluation&,
                          RangeType& ) const
    {
      return 0.;
    }

    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeIntersection(const Intersection& it,
                                const double time,
                                const QuadratureImp& quadInner,
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec)
    {
    }

    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeBoundary(const Intersection& it,
                            const double time,
                            const QuadratureImp& quadInner,
                            const ArgumentTupleVector& uLeftVec)
    {
    }

    //! dummy method
    void switchUpwind() const {}

    /**
     * \brief flux function on interfaces between cells for advection and diffusion
     *
     * @param[in] it intersection
     * @param[in] time current time given by TimeProvider
     * @param[in] x coordinate of required evaluation local to \c it
     * @param[in] uLeft DOF evaluation on this side of \c it
     * @param[in] uRight DOF evaluation on the other side of \c it
     * @param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * @param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * @note For dual operators we have \c gDiffLeft = 0 and \c gDiffRight = 0.
     *
     * \return wave speed estimate (multiplied with the integration element of the intersection),
     *              to estimate the time step |T|/wave.
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      return gradientFlux_.gradientNumericalFlux(
                left.intersection(), left.entity(), right.entity(), left.time(),
                left.quadrature(), right.quadrature(), left.index(),
                left.values()[ uVar ], right.values()[ uVar ],
                gLeft, gRight, gDiffLeft, gDiffRight);
    }

    /**
     * \brief method required by LocalDGPass
     */
    template <class LocalEvaluation>
    void analyticalFlux(const LocalEvaluation& local,
                        JacobianRangeType& f)
    {
      model_.jacobian( local.entity(), local.time(), local.point(), local.values()[ uVar ], f);
    }

    /**
     * \brief same as numericalFlux() but for the boundary
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      return boundaryFluxImpl( left.intersection(), left.time(),
                               left.quadrature(), left.index(),
                               left.values()[ uVar ], gLeft, gDiffLeft );
    }

  protected:
    template <class QuadratureImp,
              class UType>
    double boundaryFluxImpl(const Intersection& it,
                        double time,
                        const QuadratureImp& faceQuadInner,
                        const int quadPoint,
                        const UType& uLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );

      UType uRight;

      if( model_.hasBoundaryValue(it,time,x) )
      {
        model_.boundaryValue(it, time, x, uLeft, uRight);
      }
      else
      {
        uRight = uLeft;
      }

      return gradientFlux_.gradientBoundaryFlux(it, inside(),
                             time, faceQuadInner, quadPoint,
                             uLeft,
                             uRight,
                             gLeft,
                             gDiffLeft );
    }

  private:
    const ModelType&   model_;
    const NumFluxType& gradientFlux_;
    const double cflDiffinv_;
  };


  // AdvectionDiffusionLDGModel
  //---------------------------

  template< class Traits,
            int passUId, int passGradId,
            bool advectionPartExists, bool diffusionPartExists >
  class AdvectionDiffusionLDGModel;


  // AdvectionDiffusionLDGTraits
  //----------------------------

  template <class Traits,
            int passUId, int passGradId,
            bool advectionPartExists, bool diffusionPartExists >
  struct AdvectionDiffusionLDGTraits
  : public AdvectionTraits
          < Traits, passUId, passGradId, advectionPartExists>

  {
    typedef AdvectionDiffusionLDGModel
      < Traits, passUId, passGradId,
        advectionPartExists, diffusionPartExists >                   DGDiscreteModelType;
  };


  // AdvectionDiffusionLDGModel
  //---------------------------

  template< class OpTraits,
            int passUId, int passGradId,
            bool advectionPartExists, bool diffusionPartExists >
  class AdvectionDiffusionLDGModel :
    public AdvectionModel< OpTraits, passUId, passGradId, advectionPartExists >
  {
  public:
    typedef AdvectionDiffusionLDGTraits
      < OpTraits, passUId, passGradId,
          advectionPartExists, diffusionPartExists >  Traits;

    typedef AdvectionModel< OpTraits, passUId, passGradId, advectionPartExists>    BaseType;

    using BaseType :: inside ;
    using BaseType :: outside ;
    using BaseType :: model_;
    using BaseType :: uBnd_;

    // These type definitions allow a convenient access to arguments of pass.
    using BaseType :: uVar;

    integral_constant< int, passGradId> sigmaVar;

  public:
    enum { dimDomain = Traits :: dimDomain };
    enum { dimRange  = Traits :: dimRange };

    enum { advection = advectionPartExists };
    enum { diffusion = diffusionPartExists };

    typedef typename Traits :: DomainType         DomainType;
    typedef typename Traits :: FaceDomainType     FaceDomainType;

#if defined TESTOPERATOR
    enum { ApplyInverseMassOperator = false };
#else
    enum { ApplyInverseMassOperator = true };
#endif

    typedef typename Traits :: GridPartType                            GridPartType;
    typedef typename Traits :: GridType                                GridType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename BaseType :: EntityType                            EntityType;
    typedef typename Traits :: RangeFieldType                          RangeFieldType;
    typedef typename Traits :: DomainFieldType                         DomainFieldType;
    typedef typename Traits :: RangeType                               RangeType;
    typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

    typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename Traits :: ModelType ModelType;
    typedef typename Traits :: FluxType  AdvectionFluxType;
    typedef LDGDiffusionFlux< DiscreteFunctionSpaceType, ModelType > DiffusionFluxType;
    enum { evaluateJacobian = false };
  public:
    /**
     * \brief constructor
     */
    AdvectionDiffusionLDGModel(const ModelType& mod,
                               const AdvectionFluxType& numf,
                               DiffusionFluxType& diffflux)
      : BaseType( mod, numf ),
        diffFlux_( diffflux ),
        penalty_( 1.0 ),
        cflDiffinv_( 8.0 * ( Traits::polynomialOrder + 1) )
    {}

    bool hasSource() const
    {
      return ( model_.hasNonStiffSource() || model_.hasStiffSource() );
    }

    bool hasFlux() const { return advection || diffusion; };

    /**
     * \brief analytical flux function$
     */
    template <class LocalEvaluation>
    double source( const LocalEvaluation& local,
                   RangeType& s ) const
    {
      s = 0;

      double dtEst = std::numeric_limits< double > :: max();

      typedef typename DiffusionFluxType :: GradientRangeType GradientRangeType;
      Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > uJac( local.values()[ sigmaVar ] );

      if (diffusion)
      {
        const double dtStiff =
          model_.stiffSource( local.entity(), local.time(), local.point(), local.values()[uVar], uJac, s );
        dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
      }

      if (advection)
      {
        RangeType sNonStiff(0);
        const double dtNon =
          model_.nonStiffSource( local.entity(), local.time(), local.point(), local.values()[uVar], uJac, sNonStiff );

        s += sNonStiff;

        dtEst = ( dtNon > 0 ) ? std::min( dtEst, dtNon ) : dtEst;
      }

      // return the fastest wave from source terms
      return dtEst;
    }

    void switchUpwind() const
    {
      BaseType :: switchUpwind();
      diffFlux_.switchUpwind();
    }

  public:
    /**
     * \brief flux function on interfaces between cells for advection and diffusion
     *
     * \param[in]  left local evaluation context of inside cell
     * \param[in]  right local evaluation context of outside cell
     * \param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * \param[out] gRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \phi \f$
     * \param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * \param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * \note For dual operators we have \c gDiffLeft = 0 and \c gDiffRight = 0.
     *
     * \return wave speed estimate (multiplied with the integration element of the intersection),
     *              to estimate the time step |T|/wave.
     */
    template< class LocalEvaluation >
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight )
    {
      // advection

      const double wave = BaseType ::
        numericalFlux( left, right,
                       gLeft, gRight, gDiffLeft, gDiffRight );

      // diffusion

      double diffTimeStep = 0.0;
      if( diffusion )
      {
        RangeType dLeft, dRight;
        diffTimeStep =
          diffFlux_.numericalFlux(left.intersection(), *this,
                                  left.time(), left.quadrature(), right.quadrature(), left.index(),
                                  left.values()[ uVar ], right.values()[ uVar ],
                                  left.values() [ sigmaVar ], right.values()[ sigmaVar ],
                                  dLeft, dRight,
                                  gDiffLeft, gDiffRight);

        gLeft  += dLeft;
        gRight += dRight;
      }

      gDiffLeft  = 0;
      gDiffRight = 0;

      return std::max( wave, diffTimeStep );
    }


    /**
     * \brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      // advection

      const double wave = BaseType ::
        boundaryFlux( left, gLeft, gDiffLeft );

      // diffusion

      double diffTimeStep = 0.0;

      bool hasBoundaryValue =
        model_.hasBoundaryValue( left.intersection(), left.time(), left.localPoint() );

      if( diffusion && hasBoundaryValue )
      {
        // diffusion boundary flux for Dirichlet boundaries
        RangeType dLeft ( 0 );
        diffTimeStep = diffFlux_.boundaryFlux( left.intersection(),
                               *this,
                               left.time(), left.quadrature(), left.index(),
                               left.values()[ uVar ], uBnd_, // is set during call of  BaseType::boundaryFlux
                               left.values()[ sigmaVar ],
                               dLeft,
                               gDiffLeft);
        gLeft += dLeft;
      }
      else if ( diffusion )
      {
        RangeType diffBndFlux ( 0 );

        typedef typename DiffusionFluxType :: GradientRangeType GradientRangeType;
        Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > uJac( left.values()[ sigmaVar ] );

        model_.diffusionBoundaryFlux( left.intersection(), left.time(), left.localPoint(),
                                      left.values()[uVar], uJac, diffBndFlux );
        gLeft += diffBndFlux;
      }
      else
        gDiffLeft = 0;

      return std::max( wave, diffTimeStep );
    }

    /**
     * \brief analytical flux function$
     */
    template <class LocalEvaluation>
    void analyticalFlux( const LocalEvaluation& local,
                         JacobianRangeType& f ) const
    {
      // advection

      BaseType :: analyticalFlux( local, f );

      // diffusion

      if( diffusion )
      {
        JacobianRangeType diffmatrix;

        typedef typename DiffusionFluxType :: GradientRangeType GradientRangeType;
        Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > uJac( local.values()[ sigmaVar ] );

        model_.diffusion( local.entity(), local.time(), local.point(), local.values()[ uVar ], uJac, diffmatrix);
        // ldg case
        f += diffmatrix;
      }
    }
  protected:
    DiffusionFluxType& diffFlux_;
    const double penalty_;
    const double cflDiffinv_;
  };

}
}
#endif
