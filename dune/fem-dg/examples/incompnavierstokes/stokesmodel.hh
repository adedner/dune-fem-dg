#ifndef DUNE_FEM_DG_STOKES_MODEL_HH
#define DUNE_FEM_DG_STOKES_MODEL_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/io/parameter.hh>

#include "problems/operatorsplitting.hh"
#include <dune/fem-dg/models/defaultmodel.hh>

#ifndef DOXYGEN

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
  template <class GridPartImp, class ProblemImp, class ActivationImp, int maxModelParameterSize >
  class StokesPoissonModelTraits
    : public DefaultModelTraits< GridPartImp, ProblemImp, ActivationImp, maxModelParameterSize >
  {
    typedef DefaultModelTraits< GridPartImp, ProblemImp, ActivationImp, maxModelParameterSize >  BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                       GradientType;
  };


  template <class GridPartImp, class ProblemImp, class ActivationImp = std::tuple<> >
  class PoissonModel
  : public DefaultModel< StokesPoissonModelTraits< GridPartImp, ProblemImp, ActivationImp, 1 > >
  {
  public:

    typedef ProblemImp                                      ProblemType;
    typedef StokesPoissonModelTraits< GridPartImp, ProblemType, ActivationImp, 1 >
                                                            Traits;

    static const int velo = Traits::template IdGenerator<0>::id;
    //dummy atm
    static const int rhs = -1;
    typedef std::integral_constant< int, velo >             velocityVar;

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

    // one function space here, with dimension of range type 'dimDomain'
    typedef typename Traits::template ParameterSpaces<dimDomain>      ParameterSpacesType;

    typedef FractionalStepThetaScheme<0, rhs!=-1> SplitType;

    static const bool hasDiffusion = SplitType::hasDiffusion;
    static const bool hasAdvection = SplitType::hasAdvection;
    static const int ConstantVelocity = false;
  public:
    /**
     * \brief Constructor
     *
     * initializes model parameter
     *
     * \param problem Class describing the initial(t=0) and exact solution
     */
    PoissonModel(const ProblemType& problem)
      : problem_(problem)
    {}

    inline bool hasFlux() const { return true; }

    inline bool hasStiffSource() const { return true; }
    inline bool hasNonStiffSource() const { return false; }

    struct ComputeVelocity
    {
      typedef velocityVar VarId;

      template <class LocalEvaluation>
      const RangeType& operator() (const LocalEvaluation& local, const RangeType& u ) const
      {
        return u;
      }
    };

    template <class LocalEvaluation>
    inline double stiffSource( const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& du,
                               RangeType & s) const
    {
      // mass part
      RangeType mass( u );
      mass /= (problem_.deltaT() * SplitType::mass());
      s = mass;

      assert( mass == mass );

      if( SplitType::hasSource )
      {
        // right hand side
        const DomainType x = local.entity().geometry().global( local.position() );
        problem_.f( x, s );
        s *= SplitType::source();
        //s  = u ;
        //s /= theta_;
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
    inline void advection(const LocalEvaluation& local,
                          const RangeType& u,
                          const JacobianRangeType& jacu,
                          FluxRangeType & f) const
    {
      if( SplitType::hasAdvection )
      {
        DomainType v = velocity( local, u );
        v *= SplitType::advection();

        // f = uV;
        for( int r=0; r<dimRange; ++r )
          for( int d=0; d<dimDomain; ++d )
            f[r][d] = v[ d ] * u[ r ];
      }
      else
        f = 0;
    }

    /**
     * \brief velocity calculation, is called by advection()
     */
    template <class LocalEvaluation>
    inline DomainType velocity(const LocalEvaluation& local, const RangeType& u ) const
    {
      return local.values( ComputeVelocity(), local, u);
    }

    /**
     * \brief diffusion term \f$a\f$
     */
    template <class LocalEvaluation>
    inline void jacobian(const LocalEvaluation& local,
                          const RangeType& u,
                          DiffusionRangeType& a) const
    {
      abort();
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
      DiffusionMatrixType K;
      DomainType xgl = local.entity().geometry().global( local.position() );
      problem_.K( xgl, K );

      DomainType values ;
      // calculate eigenvalues
      FMatrixHelp :: eigenValues( K, values );

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
      if( SplitType::hasDiffusion )
      {
        // for constant K evalute at center (see Problem 4)
        const DomainType xgl = ( problem_.constantK() ) ?
          local.entity().geometry().center () : local.entity().geometry().global(local.position())  ;

        DiffusionMatrixType K ;

        // fill diffusion matrix
        problem_.K( xgl, K );

        // scale with mu
        K *= SplitType::diffusion();

        // apply diffusion
        for( int r =0; r<dimRange; ++r )
        {
          K.mv( jac[ r ] , A[ r ] );
        }
      }
      else
        A = 0.0;
    }

    template <class LocalEvaluation>
    inline double diffusionTimeStep (const LocalEvaluation& local,
                                     const double circumEstimate,
                                     const RangeType& u) const
    {
      return 0;
    }

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
      DomainType xgl = local.entity().geometry().global( local.position() );
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

  private:

    inline double lambdaK( const DiffusionMatrixType& K ) const
    {
      DomainType values ;
      // calculate eigenvalues
      FMatrixHelp :: eigenValues( K, values );

      // value[ 0 ] is smallest ev
      return SQR(values[ dimDomain -1 ]) / values[ 0 ];
    }

    template <class T>
    T SQR( const T& a ) const
    {
      return (a * a);
    }
  protected:
    const ProblemType& problem_;
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
  template <class GridPartImp, class ProblemImp, class ActivationImp = std::tuple<> >
  class StokesModel
  : public PoissonModel< GridPartImp, typename ProblemImp::PoissonProblemType, ActivationImp >
  {
    typedef PoissonModel< GridPartImp, typename ProblemImp::PoissonProblemType, ActivationImp > BaseType;
  public:
    typedef ProblemImp                                      ProblemType;
    typedef StokesPoissonModelTraits< GridPartImp, typename ProblemType::PoissonProblemType, ActivationImp, 1 >
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
    static const bool hasAdvection = true;
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


    const ProblemType& problem () const { return problem_; }

   protected:
    const ProblemType& problem_;
  };
#endif

}
}
#endif
