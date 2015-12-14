#ifndef DUNE_FEM_DG_STOKES_MODEL_HH
#define DUNE_FEM_DG_STOKES_MODEL_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/models/defaultmodel.hh>

namespace Dune
{
namespace Fem
{

  /**********************************************
   * Analytical model                           *
   *********************************************/
  /**
   * \brief Traits class for StokesModel
   */
  template <class GridPartImp, class ProblemImp>
  class StokesModelTraits
    : public DefaultModelTraits< GridPartImp, ProblemImp >
  {
    typedef DefaultModelTraits< GridPartImp, ProblemImp >              BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                       GradientType;

    typedef std::tuple <>                                              ModelParameter;
  };

  /**
   * \brief describes the analytical model
   *
   * This is an description class for the problem
   * \f{eqnarray*}{      V + \nabla a(U)             & = & 0 \\
   * \partial_t U + \nabla\cdot (F(U)+A(U,V)) + S(U) & = & 0 \\
   *                               U                 & = & g_D \\
   *                        \nabla U \cdot n         & = & g_N \f}
   *
   * where each class methods describes an analytical function.
   * <ul>
   * <li> \f$F\f$:   advection() </li>
   * <li> \f$a\f$:   jacobian() </li>
   * <li> \f$A\f$:   diffusion() </li>
   * <li> \f$g_D\f$: boundaryValue() </li>
   * <li> \f$g_N\f$: boundaryFlux() </li>
   * <li> \f$S\f$:   stiffSource()/nonStiffSource() </li>
   * </ul>
   *
   * \attention \f$F(U)\f$ and \f$A(U,V)\f$ are matrix valued, and therefore the
   * divergence is defined as
   *
   * \f[ \Delta M = \nabla \cdot (\nabla \cdot (M_{i\cdot})^t)_{i\in 1\dots n} \f]
   *
   * for a matrix \f$M\in \mathbf{M}^{n\times m}\f$.
   *
   * \param GridPartImp GridPart for extraction of dimension
   * \param ProblemImp Class describing the initial(t=0) and exact solution
   */


  ////////////////////////////////////////////////////////
  //
  //  Analytical model for the Heat Equation
  //      dx(u) + div(uV) - epsilon*lap(u)) = 0
  //  where V is constant vector
  //
  ////////////////////////////////////////////////////////
  template <class GridPartImp, class ProblemImp, bool rightHandSideModel = false >
  class StokesModel : public DefaultModel< StokesModelTraits< GridPartImp, ProblemImp > >
  {
  public:
    typedef StokesModelTraits< GridPartImp, ProblemImp> Traits;

    enum { rhs = 0 };
    typedef std::integral_constant< int, rhs > rhsVar;
    typedef std::tuple < rhsVar > ModelParameter;


    typedef typename Traits::ProblemType                ProblemType;
    typedef typename Traits::GridPartType               GridPartType;
    typedef typename Traits::GridType                   GridType;

    static const int dimDomain = Traits::dimDomain;
    static const int dimRange = Traits::dimRange;
    typedef typename Traits::DomainFieldType            DomainFieldType;
    typedef typename Traits::RangeFieldType             RangeFieldType;
    typedef typename Traits::DomainType                 DomainType;
    typedef typename Traits::RangeType                  RangeType;
    typedef typename Traits::GradientType               GradientType;
    typedef typename Traits::FluxRangeType              FluxRangeType;
    typedef typename Traits::DiffusionRangeType         DiffusionRangeType;
    typedef typename Traits::FaceDomainType             FaceDomainType;
    typedef typename Traits::JacobianRangeType          JacobianRangeType;

    typedef typename Traits::DiffusionMatrixType        DiffusionMatrixType ;

    typedef typename Traits::EntityType                 EntityType;
    typedef typename Traits::IntersectionType           IntersectionType;

    static const bool hasAdvection = true;
    static const bool hasDiffusion = true;

  public:
    /**
     * \brief Constructor
     *
     * initializes model parameter
     *
     * \param problem Class describing the initial(t=0) and exact solution
     */
    StokesModel(const ProblemType& problem)
      : problem_(problem),
        theta_( problem.theta() )
    {}

    inline bool hasFlux() const { return true ; }

    inline bool hasStiffSource() const { return true ; }
    inline bool hasNonStiffSource() const { return false ; }

    struct ComputeRHS
    {
      typedef rhsVar     VarId;
      typedef RangeType  ReturnType;

      template <class LocalEvaluation>
      RangeType operator() (const LocalEvaluation& local) const
      {
        return RangeType( 0 );
      }
    };


    template <class LocalEvaluation>
    inline double stiffSource( const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& du,
                               RangeType & s) const
    {
      if( ! rightHandSideModel )
      {
        s  = u ;
        s /= theta_;
      }
      return 0; //step 2, RhsLaplace ()
    }

    template <class LocalEvaluation>
    inline double nonStiffSource( const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& du,
                                  RangeType & s) const
    {
      s = 0;
      return 0;
    }


    /**
     * \brief advection term \f$F\f$
     *
     * \param local local evaluation
     * \param u \f$U\f$
     * \param jacu \f$\nabla U\f$
     * \param f \f$F(U)\f$
     */
    template <class LocalEvaluation>
    inline void advection (const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& jacu,
                           FluxRangeType & f) const
    {
      f = 0 ;
    }

  protected:

    /**
     * \brief velocity calculation, is called by advection()
     */
    template <class LocalEvaluation>
    inline DomainType velocity(const LocalEvaluation& local, const RangeType& u ) const
    {
      return local.evaluate( ComputeRHS(), local, u);
    }

  public:

    /**
     * \brief diffusion term \f$a\f$
     */
    template <class LocalEvaluation>
    inline void jacobian(const LocalEvaluation& local,
                          const RangeType& u,
                          DiffusionRangeType& a) const
    {
      a = 0;

      assert( a.rows == dimRange * dimDomain );
      assert( a.cols == dimDomain );

      for (int r=0;r<dimRange;r++)
        for (int d=0;d<dimDomain;d++)
          a[dimDomain*r+d][d] = u[r];
    }


    template <class LocalEvaluation>
    inline void eigenValues(const LocalEvaluation& local,
                            const RangeType& u,
                            RangeType& maxValue) const
    {
      DiffusionMatrixType K ;
      DomainType xgl = local.entity().geometry().global( local.point() );
      problem_.K( xgl, K );

      DomainType values ;
      // calculate eigenvalues
      FMatrixHelp :: eigenValues( K, values );

      maxValue = SQR(values[ dimDomain -1 ]) /values[0];
      return ;

      // take max eigenvalue
      maxValue = values.infinity_norm();
    }

    inline double lambdaK( const DiffusionMatrixType& K ) const
    {
      DomainType values ;
      // calculate eigenvalues
      FMatrixHelp :: eigenValues( K, values );

      // value[ 0 ] is smallest ev
      return SQR(values[ dimDomain -1 ]) / values[ 0 ];
    }

    template <class LocalEvaluation>
    inline double penaltyFactor( const LocalEvaluation& left,
                                 const LocalEvaluation& right,
                                 const RangeType& uLeft,
                                 const RangeType& uRight ) const
    {
      DiffusionMatrixType K ;
      double betaK, betaInside, betaOutside ;
      if( problem_.constantK() )
      {
        DiffusionMatrixType Kinside ;
        DiffusionMatrixType Koutside;

        const DomainType xglIn = left.entity().geometry().center();
        problem_.K( xglIn , Kinside );
        const DomainType xglOut = right.entity().geometry().center();
        problem_.K( xglOut , Koutside );

        K = Kinside ;
        K += Koutside ;
        K *= 0.5 ;

        betaInside  = lambdaK( Kinside );
        betaOutside = lambdaK( Koutside );

        betaK = lambdaK( K );
      }
      else
      {
        const DomainType xgl = left.entity().geometry().global( left.point() );
        problem_.K( xgl , K );

        betaK = lambdaK( K );
        betaInside = betaOutside = betaK;
      }

      const double jump = std::tanh( std::abs( betaInside - betaOutside ) );

      // only for small values of betS apply betS in smooth cases
      const double betaN = std :: min( betaK , 1.0 );

      // betS becomes 1 if the eigen values of both matrices are the same
      betaK = betaK * jump + (1.0 - jump) * betaN;

      return betaK ;
    }

    inline double penaltyBoundary( const EntityType& inside,
                                   const double time,
                                   const DomainType& xInside,
                                   const RangeType& uLeft ) const
    {
      return penaltyFactor( inside, inside, time, xInside, uLeft, uLeft );
    }

    /**
     * \brief diffusion term \f$A\f$
     */
    template <class LocalEvaluation>
    inline void diffusion(const LocalEvaluation& local,
                          const RangeType& u,
                          const JacobianRangeType& jac,
                          FluxRangeType& A) const
    {
      // for constant K evalute at center (see Problem 4)
      const DomainType xgl = ( problem_.constantK() ) ?
        local.entity().geometry().center () : local.entity().geometry().global(local.point())  ;

      DiffusionMatrixType K ;

      // fill diffusion matrix
      problem_.K( xgl, K );

      // apply diffusion
      for( int r =0; r<dimRange; ++r )
        K.mv( jac[ r ] , A[ r ] );
    }

    template <class LocalEvaluation>
    inline double diffusionTimeStep (const LocalEvaluation& local,
                                     const double circumEstimate,
                                     const RangeType& u) const
    {
      return 0;
    }

  public:
    /**
     * \brief checks for existence of dirichlet boundary values
     */
    template <class LocalEvaluation>
    inline bool hasBoundaryValue (const LocalEvaluation& local) const
    {
      return true;
    }

    /**
     * \brief dirichlet boundary values
     */
    template <class LocalEvaluation>
    inline void boundaryValue (const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight) const
    {
      DomainType xgl = local.entity().geometry().global( local.point() );
      problem_.g(xgl, uRight);
    }

    /**
     * \brief diffusion boundary flux
     */
    template <class LocalEvaluation>
    inline double diffusionBoundaryFlux (const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const JacobianRangeType& gradLeft,
                                         RangeType& gLeft) const
    {
      return 0.0;
    }

    const ProblemType& problem () const { return problem_; }

  private:
    template <class T>
    T SQR( const T& a ) const
    {
      return (a * a);
    }
   protected:
    const ProblemType& problem_;
    const double theta_;
  };

}
}
#endif
