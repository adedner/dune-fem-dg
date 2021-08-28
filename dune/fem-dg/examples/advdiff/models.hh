#ifndef HEAT_MODELS_HH
#define HEAT_MODELS_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/io/parameter.hh>


#include <dune/fem-dg/operator/limiter/limitpass.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/pass/dgpass.hh>


/**********************************************
 * Analytical model                           *
 *********************************************/

namespace Dune
{
namespace Fem
{
  template <class GridImp, class ProblemImp, class ActivationImp, int maxModelParameterSize >
  class HeatEqnModelTraits
    : public DefaultModelTraits< GridImp, ProblemImp, ActivationImp, maxModelParameterSize >
  {
    typedef DefaultModelTraits< GridImp, ProblemImp, ActivationImp, maxModelParameterSize >     BaseType;
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
   * \ingroup AnalyticalModels
   *
   * This is an description class for the problem
   * \f{eqnarray*}{ V + \nabla a(U)      & = & 0 \\
   * \partial_t U + \nabla \cdot(F(U)+A(U,V)) & = & 0 \\
   *                               U          & = & g_D \\
   *                        \nabla U \cdot n  & = & g_N \f}
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
   * \param Grid Grid for extraction of dimension
   * \param ProblemType Class describing the initial(t=0) and exact solution
   */
  template <class GridImp, class ProblemImp, class ActivationImp = std::tuple<> >
  class HeatEqnModel
  : public DefaultModel < HeatEqnModelTraits< GridImp, ProblemImp, ActivationImp, 1 > >
  {
  public:
    typedef HeatEqnModelTraits< GridImp, ProblemImp, ActivationImp, 1 > Traits;
    typedef DefaultModel< Traits >                                    BaseType;

    static const int velo = Traits::template IdGenerator<0>::id;

    typedef std::integral_constant< int, velo >                       velocityVar;

    typedef typename Traits::ProblemType                              ProblemType;

    static const int ConstantVelocity = ProblemType::ConstantVelocity;
    typedef typename Traits::GridType                                 GridType;
    static const int dimDomain = Traits::dimDomain;
    static const int dimRange  = Traits::dimRange;
    typedef typename Traits::DomainType                               DomainType;
    typedef typename Traits::RangeType                                RangeType;
    typedef typename Traits::GradientType                             GradientType;
    typedef typename Traits::FluxRangeType                            FluxRangeType;
    typedef typename Traits::DiffusionRangeType                       DiffusionRangeType;
    typedef typename Traits::DiffusionMatrixType                      DiffusionMatrixType;
    typedef typename Traits::FaceDomainType                           FaceDomainType;
    typedef typename Traits::JacobianRangeType                        JacobianRangeType;

    // for heat equations advection is disabled
    static const bool hasAdvection = true;
    static const bool hasDiffusion = true;


    // one function space here, with dimension of range type 'dimDomain'
    typedef typename Traits::template ParameterSpaces<dimDomain>      ParameterSpacesType;

    using BaseType::time_;
    using BaseType::time;

  public:
    /**
     * \brief Constructor
     *
     * initializes model parameter
     *
     * \param problem Class describing the initial(t=0) and exact solution
     */
    HeatEqnModel(const ProblemType& problem)
    : BaseType( problem.startTime() ),
      problem_(problem),
      epsilon_(problem.epsilon()),
      tstepEps_( getTStepEps() )
    {}

    inline const ProblemType& problem() const { return problem_; }

    inline bool hasMass() const { return problem_.hasMass(); }
    inline bool hasFlux() const { return true; }
    inline bool hasStiffSource() const { return problem_.hasStiffSource(); }
    inline bool hasNonStiffSource() const { return problem_.hasNonStiffSource(); }

    template <class LocalEvaluation>
    inline double nonStiffSource( const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& du,
                                  RangeType & s) const
    {
      DomainType xgl = local.entity().geometry().global( local.position() );
      return problem_.nonStiffSource( xgl, time_, u, s );
    }


    template <class LocalEvaluation>
    inline double stiffSource( const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& jac,
                               RangeType & s) const
    {
      DomainType xgl = local.entity().geometry().global( local.position() );
      return problem_.stiffSource( xgl, time_, u, s );
    }

    template <class LocalEvaluation>
    inline void mass( const LocalEvaluation& local,
                      const RangeType& u,
                      RangeType& diag ) const
    {
      problem_.mass( local.entity().geometry().global( local.position() ), time_, u, diag );
    }

  private:
    struct GetVelocity
    {
      typedef velocityVar VarId;

      template <class LocalEvaluation>
      DomainType operator() (const LocalEvaluation& local, double time, const ProblemType& problem ) const
      {
        DomainType v;
        problem.velocity( local.entity().geometry().global( local.position() ), time, v);
        return v;
      }
    };
  public:

    /**
     * \brief advection term \f$F\f$
     *
     * \param local local evaluation
     * \param u \f$U\f$
     * \param jacu \f$\nabla U\f$
     * \param f \f$f(U)\f$
     */
    template <class LocalEvaluation>
    inline void advection(const LocalEvaluation& local,
                          const RangeType& u,
                          const JacobianRangeType& jacu,
                          FluxRangeType & f) const
    {
      const DomainType& v = velocity( local, u );

      // f = uV;
      for( int r=0; r<dimRange; ++r )
        for( int d=0; d<dimDomain; ++d )
          f[r][d] = v[ d ] * u[ r ];
    }

    /**
     * \brief velocity calculation, is called by advection()
     */
    template <class LocalEvaluation>
    inline DomainType velocity(const LocalEvaluation& local,
                               const RangeType& u ) const
    {
      return local.values( GetVelocity(), local, time_, problem_ );
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
      // copy v to A
      A = jac;

      // apply diffusion coefficient
      //double d =  0.; //(1.-en.geometry().global(x)[0])*en.geometry().global(x)[0]+
      //                //(1.-en.geometry().global(x)[1])*en.geometry().global(x)[1];
      A *= problem_.diffusion( u, A );//*(1.+d);
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
    template <class LocalEvaluation, class JacobianRangeImp>
    inline double boundaryFlux (const LocalEvaluation& local,
                                const RangeType& uLeft,
                                const JacobianRangeImp& jacLeft,
                                RangeType& gLeft) const
    {
      gLeft = 0;
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
    inline void boundaryValue(const LocalEvaluation& local,
                              const RangeType& uLeft,
                              RangeType& uRight) const
    {
  #ifdef TESTOPERATOR
    uRight = 0;
    return;
  #endif
      DomainType xgl = local.intersection().geometry().global( local.localPosition() );
      problem_.evaluate(xgl, time_, uRight);
    }


    // here x is in global coordinates
    template <class LocalEvaluation>
    inline void maxWaveSpeed( const LocalEvaluation& local,
                          const DomainType& normal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      const DomainType& v = velocity( local );
      advspeed   = std::abs( v * normal );
      totalspeed = advspeed;
    }

    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
      //default version!
      EOCErrorList::setErrors<L2EOCError>( *this, df );
    }

   protected:
    double getTStepEps() const
    {
      // if diffusionTimeStep is set to non-zero in the parameterfile, the
      // deltaT in the timeprovider is updated according to the diffusion
      // parameter epsilon.
      bool diff_tstep;
      Fem::Parameter::get("femdg.stepper.diffusiontimestep", diff_tstep);
      return diff_tstep ? epsilon_ : 0;
    }

   protected:
    const ProblemType& problem_;
    const double epsilon_;
    const double tstepEps_;
  };

}
}
#endif
