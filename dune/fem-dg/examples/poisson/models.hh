#ifndef DUNE_FEM_DG_POISSON_MODEL_HH
#define DUNE_FEM_DG_POISSON_MODEL_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/boundaryidprovider.hh>

#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>
#include <dune/fem-dg/misc/error/l2eocerror.hh>

namespace Dune
{
namespace Fem
{
namespace Poisson
{

  /**********************************************
   * Analytical model                           *
   *********************************************/
  /**
   * \brief Traits class for PoissonModel
   */
  template <class GridImp,
            class ProblemImp>
  class PoissonModelTraits
    : public DefaultModelTraits< GridImp, ProblemImp >
  {
    typedef DefaultModelTraits< GridImp, ProblemImp >              BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                       GradientType;

    static const int modelParameterSize = 0;
  };

  /**
   * \brief Describes the analytical model.
   *
   * \ingroup AnalyticalModels
   *
   * The main goal of this class is provide all the analytical data of
   * the partial differential equation.
   *
   * This is an description class for the problem
   * \f{eqnarray*}{ v + \nabla a(u)                 & = & 0 \\
   * \partial_t u + \nabla \cdot (F(u)+A(u,v)) +S_1(u) + S_2(u) & = & 0 \\
   *                                u                & = & g_D \\
   *                               \nabla u \cdot n  & = & g_N \f}
   *
   * where each class methods describes an analytical function.
   *
   * | Method                     | formular                                  |
   * | -------------------------- | ----------------------------------------- |
   * | stiffSource()              | \f$ S_1 \f$                               |
   * | nonStiffSource()           | \f$ S_2 \f$                               |
   * | boundaryValue()            | \f$ g_D\f$                                |
   * | hasBoundaryValue()         | true if \f$ x \in \Gamma_D \f$            |
   * | diffusion()                | \f$ A \f$                                 |
   * | advection()                | \f$ F \f$                                 |
   * | jacobian()                 | \f$ a \f$                                 |
   * | boundaryFlux()             | \f$ g_N \f$                               |
   * | diffusionBoundaryFlux()    | \f$ ??? \f$                               |
   * | hasFlux()                  | true if \f$ F\neq 0\f$, false otherwise   |
   * | hasStiffSource()           | true if \f$ S_1\neq 0\f$, false otherwise |
   * | hasNonStiffSource()        | true if \f$ S_2\neq 0\f$, false otherwise |
   *
   * \attention \f$F(u)\f$ and \f$A(u,v)\f$ are matrix valued, and therefore the
   * divergence is defined as
   * \f[ \Delta M = \nabla \cdot (\nabla \cdot (M_{i\cdot})^t)_{i\in 1\dots n} \f]
   * for a matrix \f$M\in \mathbf{M}^{n\times m}\f$.
   *
   * \param Grid Grid for extraction of dimension
   * \param ProblemType Class describing the initial(t=0) and exact solution
   */
  template <class GridImp, class ProblemImp>
  class Model : public DefaultModel< PoissonModelTraits< GridImp, ProblemImp > >
  {
    typedef DefaultModel< PoissonModelTraits< GridImp, ProblemImp > > BaseType;
  public:
    typedef GridImp                                      GridType;
    typedef BoundaryIdProvider< GridType >               BoundaryIdProviderType;
    typedef PoissonModelTraits< GridImp, ProblemImp >    Traits;
    typedef typename Traits::ProblemType                 ProblemType;

    static const int dimDomain = Traits::dimDomain;
    static const int dimRange = Traits::dimRange;

    typedef typename Traits::DomainFieldType             DomainFieldType;
    typedef typename Traits::RangeFieldType              RangeFieldType;
    typedef typename Traits::DomainType                  DomainType;
    typedef typename Traits::RangeType                   RangeType;
    typedef typename Traits::GradientType                GradientType;
    typedef typename Traits::FluxRangeType               FluxRangeType;
    typedef typename Traits::DiffusionRangeType          DiffusionRangeType;
    typedef typename Traits::FaceDomainType              FaceDomainType;
    typedef typename Traits::JacobianRangeType           JacobianRangeType;

    typedef typename ProblemType::DiffusionMatrixType    DiffusionMatrixType ;

    static const bool hasDiffusion = true;
    static const int ConstantVelocity = false;

    using BaseType::stiffSource;
    using BaseType::nonStiffSource;

    /**
     * \brief Constructor
     *
     * initializes model parameter
     *
     * \param problem Class describing the initial(t=0) and exact solution
     */
    Model (const ProblemType& problem)
      : problem_(problem), K_( 0 )
    {
      if( problem_.constantK() )
      {
        DomainType xgl( 0 );
        problem_.K( xgl, K_ );
      }
    }

    /**
     * \copydoc DefaultModel::hasFlux()
     */
    inline bool hasFlux () const { return true ; }

    /**
     * \copydoc DefaultModel::hasStiffSource()
     */
    inline bool hasStiffSource () const { return true ; }

    /**
     * \copydoc DefaultModel::hasNonStiffSource()
     */
    inline bool hasNonStiffSource () const { return false ; }

    /**
     * \copydoc DefaultModel::stiffSource(LocalEvaluation,RangeType,JacobianRangeType,RangeType)
     */
    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& du,
                               RangeType& s) const
    {
      DomainType xgl = local.entity().geometry().global( local.position() );
      problem_.f( xgl, s );
      return 0.0;
    }

    /**
     * \copydoc DefaultModel::nonStiffSource(const LocalEvaluation& local,const RangeType& u,const JacobianRangeType& jac,RangeType & s))
     */
    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& du,
                                  RangeType & s) const
    {
      DomainType xgl = local.entity().geometry().global( local.position() );
      problem_.f( xgl, s );
      return 0.0;
    }


    /**
     * \brief advection term \f$F\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  jac evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] f the result \f$F(u)\f$
     */
    template <class LocalEvaluation>
    inline void advection (const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& jac,
                           FluxRangeType& f) const
    {
      const DomainType v = velocity( local, u );
      //f = uv;
      for( int r=0; r<dimRange; ++r )
        for( int d=0; d<dimDomain; ++d )
          f[r][d] = v[ d ] * u[ r ];
    }

    /**
     * \copydoc DefaultModel::jacobian()
     */
    template <class LocalEvaluation>
    inline void jacobian (const LocalEvaluation& local,
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

    /**
     * \brief returns the maximum eigen value of \f$ K \f$.
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] the maximal eigen value
     */
    template <class LocalEvaluation>
    inline void eigenValues (const LocalEvaluation& local,
                             const RangeType& u,
                             RangeType& maxValue) const
    {
      DiffusionMatrixType K;
      DomainType xgl = local.entity().geometry().global( local.position() );
      problem_.K( xgl, K );

      DomainType values ;
      // calculate eigenvalues
      Dune::FMatrixHelp :: eigenValues( K, values );

      maxValue = SQR(values[ dimDomain -1 ]) /values[0];
      return ;

      // take max eigenvalue
      maxValue = values.infinity_norm();
    }

    /**
     * \brief diffusion term \f$A\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  jac evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] A The result f$ A(u) \f$
     */
    template <class LocalEvaluation>
    inline void diffusion (const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& jac,
                           FluxRangeType& A) const
    {
      if( problem_.constantK() )
      {
        // apply diffusion
        for( int r =0; r<dimRange; ++r )
          K_.mv( jac[ r ] , A[ r ] );
      }
      else
      {
        // for constant K evalute at center (see Problem 4)
        const DomainType xgl = local.entity().geometry().global(local.position());

        DiffusionMatrixType K ;

        // fill diffusion matrix
        problem_.K( xgl, K );
        // apply diffusion
        for( int r =0; r<dimRange; ++r )
          K.mv( jac[ r ] , A[ r ] );
      }
    }

    /**
     * \brief returns the maximal time step size which is given through
     * the diffusion term.
     *
     * \param[in]  local local evaluation
     * \param[in]  circumEstimate estimation of the circum
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     */
    template <class LocalEvaluation>
    inline double diffusionTimeStep (const LocalEvaluation& local,
                                     const RangeType& u) const
    {
      return 0;
    }

    /**
     * \brief checks for existence of dirichlet boundary values
     *
     * \param[in] local local evaluation
     */
    template <class LocalEvaluation>
    inline bool hasBoundaryValue (const LocalEvaluation& local ) const
    {
      return true;
    }

    /**
     * \brief returns the Dirichlet boundary values
     *
     * \param[in]  local local evaluation
     * \param[in]  uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[out] uRight the Dirichlet boundary value
     */
    template <class LocalEvaluation>
    inline void boundaryValue (const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight) const
    {
      if ( BoundaryIdProviderType::boundaryId( local.intersection() ) == 99) // Dirichlet zero boundary conditions
      {
        uRight = 0;
      }
      else
      {
        DomainType xgl = local.entity().geometry().global( local.position() );
        problem_.g(xgl, uRight);
      }
    }

    /**
     * \brief returns the Diffusion boundary flux
     *
     * \param[in]  local local evaluation
     * \param[in]  uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[in]  gradLeft evaluation of the gradient of the local function, i.e. \f$\nabla u_{E^+}( \hat{x} )\f$
     * \param[out] gLeft the Dirichlet boundary value for the diffusion part
     */
    template <class LocalEvaluation>
    inline double diffusionBoundaryFlux (const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const JacobianRangeType& gradLeft,
                                         RangeType& gLeft) const
    {
      return 0.0;
    }

    /**
     * \brief returns the problem.
     */
    const ProblemType& problem () const { return problem_; }

    /**
     * \brief velocity calculation
     *
     * \note This method is internally called by by advection() and
     * externally called by some numerical fluxes (i.e. UpwindFlux etc.).
     *
     * \param[in] local local evaluation
     */
    template <class LocalEvaluation>
    DomainType velocity (const LocalEvaluation& local,
                         const RangeType& uLeft ) const
    {
      return DomainType( problem().constantAdvection() );
    }

    template <class LocalEvaluation>
    inline double penaltyFactor (const LocalEvaluation& left,
                                 const LocalEvaluation& right,
                                 const RangeType& uLeft,
                                 const RangeType& uRight) const
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
        const DomainType xgl = left.entity().geometry().global( left.position() );
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
    template< class DiscreteFunction, class SigmaDiscreteFunction >
    void eocErrors( const DiscreteFunction& df, const SigmaDiscreteFunction& sigma ) const
    {
      eocErrors( df );
      //TODO: What is the analytical solution for sigma?
      //EOCErrorList::setErrors<H1EOCError>( *this, sigma );
    }

    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
      //default version!
      EOCErrorList::setErrors<L2EOCError>( *this, df );
      EOCErrorList::setErrors<DGEOCError>( *this, df );
    }


  private:
    template <class T>
    T SQR( const T& a ) const
    {
      return (a * a);
    }

    inline double lambdaK( const DiffusionMatrixType& K ) const
    {
      DomainType values ;
      // calculate eigenvalues
      Dune::FMatrixHelp :: eigenValues( K, values );

      // value[ 0 ] is smallest ev
      return SQR(values[ dimDomain -1 ]) / values[ 0 ];
    }

  protected:
    const ProblemType& problem_;
    DiffusionMatrixType K_ ;
  };

}
}
}
#endif
