#ifndef DUNE_FEM_DG_STOKES_MODEL_HH
#define DUNE_FEM_DG_STOKES_MODEL_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/boundaryidprovider.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/dgeocerror.hh>

namespace Dune
{
namespace Fem
{
namespace Stokes
{

  /**********************************************
   * Analytical model                           *
   *********************************************/
  /**
   * \brief Traits class for StokesModel
   */
  template <class GridImp, class ProblemImp>
  class StokesPoissonModelTraits
    : public DefaultModelTraits< GridImp, ProblemImp >
  {
    typedef DefaultModelTraits< GridImp, ProblemImp >       BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                       GradientType;

    static const int modelParameterSize = 0;
  };

  template <class GridImp, class ProblemImp>
  class PoissonModel
    : public DefaultModel< StokesPoissonModelTraits< GridImp, ProblemImp > >
  {
  public:
    typedef ProblemImp                                      ProblemType;
    typedef StokesPoissonModelTraits< GridImp, ProblemType >
                                                            Traits;
    typedef DefaultModel<Traits>                            BaseType;


    typedef typename Traits::GridType                       GridType;
    typedef Dune::Fem::BoundaryIdProvider< GridType >       BoundaryIdProviderType;

    static const int dimDomain = Traits::dimDomain;
    static const int dimRange  = Traits::dimRange;
    typedef typename Traits::DomainFieldType                DomainFieldType;
    typedef typename Traits::RangeFieldType                 RangeFieldType;
    typedef typename Traits::DomainType                     DomainType;
    typedef typename Traits::RangeType                      RangeType;
    typedef typename Traits::GradientType                   GradientType;
    typedef typename Traits::FluxRangeType                  FluxRangeType;
    typedef typename Traits::DiffusionRangeType             DiffusionRangeType;
    typedef typename Traits::FaceDomainType                 FaceDomainType;
    typedef typename Traits::JacobianRangeType              JacobianRangeType;

    typedef typename Traits::DiffusionMatrixType            DiffusionMatrixType ;

    static const bool hasDiffusion = true;
    static const int ConstantVelocity = false;

    using BaseType::time;

    /**
     * \brief Constructor
     *
     * initializes model parameter
     *
     * \param problem Class describing the initial(t=0) and exact solution
     */
    PoissonModel(const ProblemType& problem)
      : BaseType(0.0),
        problem_(problem)
    {}

    inline bool hasFlux() const { return true ; }

    inline bool hasStiffSource() const { return true ; }
    inline bool hasNonStiffSource() const { return false ; }

    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& du,
                               RangeType & s) const
    {
      const DomainType x = local.entity().geometry().global( local.position() );
      // right hand side
      problem_.f( x, s );

      RangeType mass( u );
      mass *= problem_.gamma() ;
      s += mass;
      return 0.0;
    }

    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& du,
                                  RangeType & s) const
    {
      s = 0;
      return 0.0;
    }

    /**
     * \brief advection term \f$F\f$
     *
     * \param local local evaluation
     * \param u \f$U\f$
     * \param jacu \f$\nabla U\f$
     * \param f \f$f(U)\f$
     */
    template <class LocalEvaluation>
    inline void advection (const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& jacu,
                           FluxRangeType & f) const
    {
      // evaluate velocity V
      const DomainType v = velocity( local, u );

      //f = uV;
      for( int r=0; r<dimRange; ++r )
        for( int d=0; d<dimDomain; ++d )
          f[r][d] = v[ d ] * u[ r ];
    }

    /**
     * \brief velocity calculation, is called by advection()
     */
    template <class LocalEvaluation>
    inline DomainType velocity ( const LocalEvaluation& local,
                                 const RangeType& u ) const
    {
      DomainType v;
      problem_.velocity( local.entity().geometry().global( local.position() ), v );
      return v;
    }

    /**
     * \brief diffusion term \f$a\f$
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

    template <class LocalEvaluation>
    inline void eigenValues(const LocalEvaluation& local,
                            const RangeType& u,
                            RangeType& maxValue) const
    {
      DiffusionMatrixType K ;
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
     */
    template <class LocalEvaluation>
    inline void diffusion(const LocalEvaluation& local,
                          const RangeType& u,
                          const JacobianRangeType& jac,
                          FluxRangeType& A) const
    {
      // for constant K evalute at center (see Problem 4)
      const DomainType xgl = ( problem_.constantK() ) ?
        local.entity().geometry().center () : local.entity().geometry().global(local.position())  ;

      DiffusionMatrixType K ;

      // fill diffusion matrix
      problem_.K( xgl, K );
      // apply diffusion
      for( int r =0; r<dimRange; ++r )
        K.mv( jac[ r ] , A[ r ] );
    }

    template <class LocalEvaluation>
    inline double diffusionTimeStep( const LocalEvaluation& local,
                                     const RangeType& u ) const
    {
      return 0;
    }

    /**
     * \brief checks for existence of dirichlet boundary values
     */
    template <class LocalEvaluation>
    inline bool hasBoundaryValue(const LocalEvaluation& local ) const
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
     * \brief diffusion boundary flux
     */
    template <class LocalEvaluation>
    inline double diffusionBoundaryFlux( const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const JacobianRangeType& gradLeft,
                                         RangeType& gLeft ) const
    {
      return 0.0;
    }

    const ProblemType& problem () const { return problem_; }

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
  };



  /**
   * \brief describes the analytical model
   *
   * \ingroup AnalyticalModels
   *
   * Analytical Model for the stokes problem.
   *
   * This is an description class for the problem
   * \f{eqnarray*}{ V + \nabla a(U)      & = & 0 \\
   * \partial_t U + \nabla (F(U)+A(U,V)) & = & 0 \\
   *                          U          & = & g_D \\
   *                   \nabla U \cdot n  & = & g_N \f}
   *
   * where each class methods describes an analytical function.
   * <ul>
   * <li> \f$F\f$:   advection() </li>
   * <li> \f$a\f$:   jacobian() </li>
   * <li> \f$A\f$:   diffusion() </li>
   * <li> \f$g_D\f$  boundaryValue() </li>
   * <li> \f$g_N\f$  boundaryFlux() </li>
   * </ul>
   *
   * \attention \f$F(U)\f$ and \f$A(U,V)\f$ are matrix valued, and therefore the
   * divergence is defined as
   *
   * \f[ \Delta M = \nabla \cdot (\nabla \cdot (M_{i\cdot})^t)_{i\in 1\dots n} \f]
   *
   * for a matrix \f$M\in \mathbf{M}^{n\times m}\f$.
   *
   * \param GridImp Grid for extraction of dimension
   * \param ProblemImp Class describing the initial(t=0) and exact solution
   */
  template <class GridImp, class ProblemImp>
  class StokesModel
  : public PoissonModel< GridImp, typename ProblemImp::PoissonProblemType >
  {
    typedef PoissonModel< GridImp, typename ProblemImp::PoissonProblemType > BaseType;
  public:
    typedef ProblemImp                                      ProblemType;
    typedef StokesPoissonModelTraits< GridImp, typename ProblemType::PoissonProblemType >
                                                            Traits;

    typedef typename Traits::GridType                       GridType;
    static const int dimDomain = Traits::dimDomain;
    static const int dimRange  = Traits::dimRange;
    typedef typename Traits::DomainFieldType                DomainFieldType;
    typedef typename Traits::RangeFieldType                 RangeFieldType;
    typedef typename Traits::DomainType                     DomainType;
    typedef typename Traits::RangeType                      RangeType;
    typedef typename Traits::GradientType                   GradientType;
    typedef typename Traits::FluxRangeType                  FluxRangeType;
    typedef typename Traits::DiffusionRangeType             DiffusionRangeType;
    typedef typename Traits::FaceDomainType                 FaceDomainType;
    typedef typename Traits::JacobianRangeType              JacobianRangeType;

    typedef typename Traits::DiffusionMatrixType            DiffusionMatrixType ;

    static const bool hasDiffusion = true;
    static const int ConstantVelocity = false;

    /**
     * \brief Constructor
     *
     * initializes model parameter
     *
     * \param problem Class describing the initial(t=0) and exact solution
     */
    StokesModel(const ProblemType& problem)
      : BaseType( problem.template get<0>() ),
        problem_( problem )
    {}

    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
      //default version, pressure!
      EOCErrorList::setErrors<L2EOCError>( *this, df );
    }


    const ProblemType& problem () const { return problem_; }

   protected:
    const ProblemType& problem_;
  };

}
}
}
#endif
