#ifndef DUNE_FEM_DG_FLUXDISCRETEMODEL_HH
#define DUNE_FEM_DG_FLUXDISCRETEMODEL_HH

#include <type_traits>

// Dune-Fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/operator/dg/discretemodelcommon.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/ldgflux.hh>

//*************************************************************
namespace Dune
{
namespace Fem
{


  // GradientModel
  //--------------

  template <class Traits, int passUId, int... passIds >
  class GradientModel;


  // GradientTraits
  //---------------

  template <class Traits, int... passIds >
  struct GradientTraits : public Traits
  {
    typedef GradientModel< Traits, passIds... >     DGDiscreteModelType;
  };


  // GradientModel
  //--------------

  template < class OpTraits, int passUId, int... passIds >
  class GradientModel :
    public Fem::DGDiscreteModelDefaultWithInsideOutside
      < GradientTraits< OpTraits, passUId, passIds...>, passUId, passIds... >
  {
    typedef Fem::DGDiscreteModelDefaultWithInsideOutside
              < GradientTraits< OpTraits, passUId, passIds... >, passUId, passIds... > BaseType;

    using BaseType::inside;
    using BaseType::outside;

    // This type definition allows a convenient access to arguments of passes.
    std::integral_constant< int, passUId > uVar;

  public:
    typedef GradientTraits< OpTraits, passUId, passIds... >        Traits;
    typedef typename Traits::ModelType                             ModelType;
    typedef typename Traits::DiffusionFluxType                     DiffusionFluxType;

    typedef typename Traits::DiscreteFunctionSpaceType::FunctionSpaceType::DomainType        DomainType;
    typedef typename Traits::DiscreteFunctionSpaceType::FunctionSpaceType::RangeType         RangeType;
    typedef typename Traits::DiscreteFunctionSpaceType::FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridType                              GridType;
    typedef typename Traits::GridPartType                          GridPartType;
    typedef typename GridPartType::IntersectionIteratorType        IntersectionIterator;
    typedef typename GridPartType::IntersectionType                Intersection;
    typedef typename BaseType::EntityType                          EntityType;

    enum { evaluateJacobian = DiffusionFluxType::evaluateJacobian };

    // necessary for TESTOPERATOR
    // not sure how it works for dual operators!
    enum { ApplyInverseMassOperator = true };

  public:
    /**
     * \brief constructor
     */
    GradientModel(const ModelType& mod,
                  const DiffusionFluxType& diffFlux) :
      model_( mod ),
      gradientFlux_( diffFlux ),
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
     * \param[in] left local evaluation
     * \param[in] right local evaluation
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
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      return gradientFlux_.gradientNumericalFlux( left, right,
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
      model_.jacobian( local, local.values()[ uVar ], f);
    }

    /**
     * \brief same as numericalFlux() but for the boundary
     */
    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const
    {
      return boundaryFluxImpl( left, left.values()[ uVar ], gLeft, gDiffLeft );
    }

  protected:
    template <class LocalEvaluation, class UType >
    double boundaryFluxImpl( const LocalEvaluation& left,
                             const UType& uLeft,
                             RangeType& gLeft,
                             JacobianRangeType& gDiffLeft ) const
    {
      const DomainType normal = left.intersection().integrationOuterNormal( left.localPosition() );

      UType uRight;

      if( model_.hasBoundaryValue( left ) )
        model_.boundaryValue( left, uLeft, uRight);
      else
        uRight = uLeft;

      return gradientFlux_.gradientBoundaryFlux( left, uLeft, uRight, gLeft, gDiffLeft );
    }

  private:
    const ModelType&   model_;
    const DiffusionFluxType& gradientFlux_;
    const double cflDiffinv_;
  };


  // AdvectionDiffusionLDGModel
  //---------------------------

  template< class Traits,
            bool enableAdvection, bool enableDiffusion,
            int passUId, int passGradId, int... passIds >
  class AdvectionDiffusionLDGModel;


  // AdvectionDiffusionLDGTraits
  //----------------------------

  template <class Traits,
            bool enableAdvection, bool enableDiffusion,
            int... passIds >
  struct AdvectionDiffusionLDGTraits
  : public AdvectionTraits< Traits, enableAdvection, passIds... >
  {
    typedef AdvectionDiffusionLDGModel< Traits, enableAdvection, enableDiffusion, passIds... >
                                                  DGDiscreteModelType;
  };


  // AdvectionDiffusionLDGModel
  //---------------------------

  template< class OpTraits,
            bool enableAdvection, bool enableDiffusion,
            int passUId, int passGradId, int... passIds >
  class AdvectionDiffusionLDGModel :
    public AdvectionModel< OpTraits, enableAdvection, passUId, passGradId, passIds... >
  {
  public:
    typedef AdvectionDiffusionLDGTraits< OpTraits, enableAdvection, enableDiffusion,
                                         passUId, passGradId, passIds... > Traits;

    typedef AdvectionModel< OpTraits, enableAdvection, passUId, passGradId, passIds... >    BaseType;

    using BaseType::inside ;
    using BaseType::outside ;
    using BaseType::model_;
    using BaseType::uBnd_;

    // These type definitions allow a convenient access to arguments of pass.
    using BaseType::uVar;

    std::integral_constant< int, passGradId> sigmaVar;

  public:

    typedef typename Traits::ModelType                       ModelType;

    enum { dimRange  = ModelType::dimRange };
    enum { dimDomain = ModelType::Traits::dimDomain };

    enum { advection = enableAdvection };
    enum { diffusion = enableDiffusion };

    typedef typename BaseType::DomainType                    DomainType;
    typedef typename BaseType::RangeFieldType                RangeFieldType;
    typedef typename BaseType::DomainFieldType               DomainFieldType;
    typedef typename BaseType::RangeType                     RangeType;
    typedef typename BaseType::JacobianRangeType             JacobianRangeType;

#if defined TESTOPERATOR
    enum { ApplyInverseMassOperator = false };
#else
    enum { ApplyInverseMassOperator = true };
#endif

    typedef typename Traits::GridPartType                     GridPartType;
    typedef typename Traits::GridType                         GridType;
    typedef typename GridPartType::IntersectionIteratorType   IntersectionIterator;
    typedef typename GridPartType::IntersectionType           Intersection;
    typedef typename BaseType::EntityType                     EntityType;
    typedef typename Traits::DiscreteFunctionSpaceType        DiscreteFunctionSpaceType;

    typedef typename Traits::AdvectionFluxType                AdvectionFluxType;
    typedef typename Traits::DiffusionFluxType                DiffusionFluxType;
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
          model_.stiffSource( local, local.values()[uVar], uJac, s );
        dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
      }

      if (advection)
      {
        RangeType sNonStiff(0);
        const double dtNon =
          model_.nonStiffSource( local, local.values()[uVar], uJac, sNonStiff );

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
        typedef typename DiffusionFluxType::GradientRangeType GradientRangeType;
        Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > jacLeft ( left.values()[ sigmaVar ] );
        Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > jacRight( right.values()[ sigmaVar ] );

        diffTimeStep =
          diffFlux_.numericalFlux(left, right,
                                  left.values()[ uVar ], right.values()[ uVar ],
                                  jacLeft, jacRight,
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

      bool hasBoundaryValue = model_.hasBoundaryValue( left );
      bool hasRobinBoundaryValue = model_.hasRobinBoundaryValue( left );

      if( diffusion )
      {
        if( hasBoundaryValue || hasRobinBoundaryValue )
        {
          // diffusion boundary flux for Dirichlet boundaries
          RangeType dLeft ( 0 );
          typedef typename DiffusionFluxType::GradientRangeType GradientRangeType;
          Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > uJac( left.values()[ sigmaVar ] );

          diffTimeStep = diffFlux_.boundaryFlux( left,
                                                 left.values()[ uVar ],
                                                 uBnd_, // is set during call of  BaseType::boundaryFlux
                                                 uJac,
                                                 dLeft,
                                                 gDiffLeft);
          gLeft += dLeft;
        }
        if( !hasBoundaryValue )
        {
          RangeType diffBndFlux ( 0 );

          typedef typename DiffusionFluxType::GradientRangeType GradientRangeType;
          Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > uJac( left.values()[ sigmaVar ] );

          model_.diffusionBoundaryFlux( left, left.values()[uVar], uJac, diffBndFlux );
          gLeft += diffBndFlux;
        }
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

        typedef typename DiffusionFluxType::GradientRangeType GradientRangeType;
        Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > uJac( local.values()[ sigmaVar ] );

        model_.diffusion( local, local.values()[ uVar ], uJac, diffmatrix);
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
