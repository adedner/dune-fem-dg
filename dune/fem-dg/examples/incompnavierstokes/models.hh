#ifndef HEAT_MODELS_HH
#define HEAT_MODELS_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/io/parameter.hh>


#include <dune/fem-dg/operator/limiter/limitpass.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/pass/dgpass.hh>

#include "problems/operatorsplitting.hh"

/**********************************************
 * Analytical model                           *
 *********************************************/

namespace Dune
{
namespace Fem
{

  /**
   * \brief Traits class for NavierStokesModel
   */
  template <class GridPartImp, class ProblemImp, class ActivationImp, int maxModelParameterSize  >
  class NavierStokesModelTraits
    : public DefaultModelTraits< GridPartImp, ProblemImp, ActivationImp, maxModelParameterSize  >
  {
    typedef DefaultModelTraits< GridPartImp, ProblemImp, ActivationImp, maxModelParameterSize > BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                       GradientType;

    //typedef Fem::MinModLimiter< FieldType > LimiterFunctionType ;
    //typedef SuperBeeLimiter< FieldType > LimiterFunctionType ;
    //typedef VanLeerLimiter< FieldType > LimiterFunctionType ;
  };

  /**
   * \brief describes the analytical model
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
   * \param GridPart GridPart for extraction of dimension
   * \param ProblemType Class describing the initial(t=0) and exact solution
   */


  ////////////////////////////////////////////////////////
  //
  //  Analytical model for the Heat Equation
  //      dt(u) + div(uV) - epsilon*lap(u)) = 0
  //  where V is constant vector
  //
  ////////////////////////////////////////////////////////
  template <class GridPartImp, class ProblemImp, class ActivationImp = std::tuple<> >
  class NavierStokesModel :
    public DefaultModel < NavierStokesModelTraits< GridPartImp, ProblemImp, ActivationImp, 2 > >
  {
  public:

    typedef NavierStokesModelTraits< GridPartImp, ProblemImp, ActivationImp, 2 > Traits;
    typedef DefaultModel< Traits >                         BaseType;

    static const int velo = Traits::template IdGenerator<0>::id;
    static const int rhs = Traits::template IdGenerator<1>::id;

    typedef std::integral_constant< int, velo >           velocityVar;
    typedef std::integral_constant< int, rhs  >           rhsVar;

    typedef typename Traits::ProblemType                  ProblemType;
    static_assert( ProblemType::dimRange == Traits::dimRange, "dimRange of Problem and Model does not fit.");

    static const int ConstantVelocity = ProblemType::ConstantVelocity;
    typedef typename Traits::GridType                     GridType;
    static const int dimDomain = Traits::dimDomain;
    static const int dimRange  = Traits::dimRange;
    typedef typename Traits::DomainType                   DomainType;
    typedef typename Traits::RangeType                    RangeType;
    typedef typename Traits::GradientType                 GradientType;
    typedef typename Traits::FluxRangeType                FluxRangeType;
    typedef typename Traits::DiffusionRangeType           DiffusionRangeType;
    typedef typename Traits::DiffusionMatrixType          DiffusionMatrixType;
    typedef typename Traits::FaceDomainType               FaceDomainType;
    typedef typename Traits::JacobianRangeType            JacobianRangeType;


    typedef FractionalStepThetaScheme<1,rhs!=-1> SplitType;

    static const bool hasDiffusion = SplitType::hasDiffusion;
    static const bool hasAdvection = SplitType::hasAdvection;

    // two function spaces here...
    typedef typename Traits::template ParameterSpaces<dimDomain,dimRange>      ParameterSpacesType;

    NavierStokesModel(const NavierStokesModel& other);
    const NavierStokesModel &operator=(const NavierStokesModel &other);
  public:
    using BaseType::time;

    /**
     * \brief Constructor
     *
     * initializes model parameter
     *
     * \param problem Class describing the initial(t=0) and exact solution
     */
    NavierStokesModel(const ProblemType& problem)
    : BaseType( problem.startTime() ),
      problem_(problem),
      epsilon_(problem.epsilon()),
      tstepEps_( problem.beta()*problem.mu() ),
      theta_( problem.theta() )
    {}

    inline const ProblemType& problem() const { return problem_; }

    inline bool hasFlux() const { return true; }
    inline bool hasStiffSource() const { return true; }
    inline bool hasNonStiffSource() const { return false; }

    template <class LocalEvaluation>
    inline double nonStiffSource( const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& du,
                                  RangeType & s) const
    {
      s = 0;
      return 0;
    }

    struct ComputeRHS
    {
      typedef rhsVar     VarId;

      //analytical rhs, if not specified in extra parameter tuple
      template <class LocalEvaluation>
      RangeType operator() (const LocalEvaluation& local) const
      {
        return RangeType( 0 );
      }
    };

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
                               const JacobianRangeType& jac,
                               RangeType & s) const
    {
      // no mass part

      if( SplitType::hasSource )
      {
        // right hand side
        problem_.stiffSource( local.position(), time(), u, s );
        s *= SplitType::source();
      }

      // right hand side f
      if( rhs != -1 )
      {
        // + \alpha \mu \Delta u^n+\theta - \nabla p
        //step 2: dgOperator (u*,rhs)  -> rhs_
        s += local.values( ComputeRHS(), local ) ;
      }
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
      FluxRangeType A(0);
      maxValue = RangeType( problem_.diffusion( u,A ) );
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
        // copy v to A
        A = jac;

        // apply diffusion coefficient
        A *= SplitType::diffusion() * problem_.mu();// diffusion( u, A );//*(1.+d);
      }
      else
        A = 0.0;

    }

    template <class LocalEvaluation>
    inline double diffusionTimeStep(const LocalEvaluation& local,
                                    const RangeType& u ) const
    {
      return tstepEps_;
    }

    /** \brief convert conservative to primitive variables
     *  \note This method is used only for additional output to hdd
     */
    inline void conservativeToPrimitive( const double time,
                                         const DomainType& xgl,
                                         const RangeType& cons,
                                         RangeType& prim,
                                         const bool ) const
    {
      prim = cons;
    }

  public:
    /**
     * \brief checks for existence of dirichlet boundary values
     */
    template <class LocalEvaluation>
    inline bool hasBoundaryValue(const LocalEvaluation& local ) const
    {
      return true;
    }

    /**
     * \brief neuman boundary values \f$g_N\f$ for pass1
     */
    template <class LocalEvaluation>
    inline double boundaryFlux(const LocalEvaluation& local,
                               const RangeType& uLeft,
                               const JacobianRangeType& jacLeft,
                               RangeType& gLeft) const
    {
      gLeft = 0.;
      return 0.;
    }

    /** \brief boundary flux for the diffusion part
     */
    template <class LocalEvaluation, class JacobianRangeImp>
    inline double diffusionBoundaryFlux( const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const JacobianRangeImp& jacLeft,
                                         RangeType& gLeft ) const
    {
      std::cerr <<"diffusionBoundaryFlux shouldn't be used for this testcase" <<std::endl;
      abort();
    }

    /**
     * \brief dirichlet boundary values
     */
    template <class LocalEvaluation>
    inline  void boundaryValue(const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight) const
    {
  #ifdef TESTOPERATOR
    uRight = 0;
    return;
  #endif
      DomainType xgl = local.intersection().geometry().global( local.localPosition() );
      problem_.evaluate(xgl, time(), uRight);
    }


    // here x is in global coordinates
    template <class LocalEvaluation>
    inline void maxWaveSpeed( const LocalEvaluation& local,
                          const DomainType& normal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      const DomainType& v = velocity( local, u );
      advspeed   = std::abs( v * normal );
      totalspeed = advspeed;
    }

   protected:
    const ProblemType& problem_;
    const double epsilon_;
    const double tstepEps_;
    double theta_;
  };

}
}
#endif
