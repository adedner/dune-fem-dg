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

/**
 * @brief Traits class for NavierStokesModel
 */
template <class GridPart, int dimR>
class NavierStokesModelTraits
  : public Dune::Fem::FunctionSpace< typename GridPart::GridType::ctype,
                                     double,
                                     GridPart::GridType::dimensionworld,
                                     dimR >
{
  typedef Dune::Fem::FunctionSpace< typename GridPart::GridType::ctype,
                                    double,
                                    GridPart::GridType::dimensionworld,
                                    dimR >  BaseType;
public:
  enum { velo = 0, press = 1, blabla = 2 };
  typedef std::integral_constant< int, velo   > velocityVar;
  typedef std::integral_constant< int, press  > pressure;
  typedef std::integral_constant< int, blabla > blablabla;
  typedef std::tuple < velocityVar, pressure, blablabla > ModelParameter;

  typedef GridPart                                                      GridPartType;
  typedef typename GridPartType :: GridType                             GridType;

  typedef typename BaseType::RangeFieldType                             RangeFieldType;
  typedef typename BaseType::DomainFieldType                            DomainFieldType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };
  static const int dimGradRange = dimRange * dimDomain ;

  // Definition of domain and range types
  typedef Dune::FieldVector< DomainFieldType, dimDomain-1 >             FaceDomainType;
  typedef Dune::FieldVector< RangeFieldType, dimGradRange >             GradientType;
  // ATTENTION: These are matrices (c.f. NavierStokesModel)
  typedef typename BaseType :: JacobianRangeType                        FluxRangeType;
  typedef Dune::FieldMatrix< RangeFieldType, dimGradRange, dimDomain >  DiffusionRangeType;
  typedef Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain >     DiffusionMatrixType;
  typedef typename GridType :: template Codim< 0 > :: Entity            EntityType;
  typedef typename GridPartType :: IntersectionIteratorType             IntersectionIterator;
  typedef typename IntersectionIterator :: Intersection                 IntersectionType;

  //typedef Dune::Fem::MinModLimiter< FieldType > LimiterFunctionType ;
  //typedef SuperBeeLimiter< FieldType > LimiterFunctionType ;
  //typedef VanLeerLimiter< FieldType > LimiterFunctionType ;
};

/**
 * @brief describes the analytical model
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
 * <li> \f$a\f$:   diffusion1() </li>
 * <li> \f$A\f$:   diffusion2() </li>
 * <li> \f$g_D\f$  boundaryValue() </li>
 * <li> \f$g_N\f$  boundaryFlux1(), boundaryFlux2() </li>
 * </ul>
 *
 * \attention \f$F(U)\f$ and \f$A(U,V)\f$ are matrix valued, and therefore the
 * divergence is defined as
 *
 * \f[ \Delta M = \nabla \cdot (\nabla \cdot (M_{i\cdot})^t)_{i\in 1\dots n} \f]
 *
 * for a matrix \f$M\in \mathbf{M}^{n\times m}\f$.
 *
 * @param GridPart GridPart for extraction of dimension
 * @param ProblemType Class describing the initial(t=0) and exact solution
 */


////////////////////////////////////////////////////////
//
//  Analytical model for the Heat Equation
//      dt(u) + div(uV) - epsilon*lap(u)) = 0
//  where V is constant vector
//
////////////////////////////////////////////////////////
template <class GridPartType, class ProblemImp>
class NavierStokesModel :
  public DefaultModel < NavierStokesModelTraits< GridPartType,ProblemImp::dimRange> >
{
public:
  enum { velo = 0, laplaceU = 1, gradP = 2 };
  typedef std::integral_constant< int, velo      > velocityVar;
  typedef std::integral_constant< int, laplaceU  > laplaceUVar;
  typedef std::integral_constant< int, gradP     > gradientPressureVar;
  typedef std::tuple < velocityVar, laplaceUVar, gradientPressureVar > ModelParameter;

  //typedef Dune::Fem::Selector< velo >  ModelParameterSelectorType;
  //typedef std::tuple< VelocityType* >  ModelParameterTypes;
  //typedef Dune::Fem::Selector< >  ModelParameterSelectorType;
  //typedef std::tuple< >  ModelParameterTypes;

  // for heat equations advection is disabled
  static const bool hasAdvection = true ;
  static const bool hasDiffusion = true ;

  typedef ProblemImp ProblemType ;

  static const int ConstantVelocity = ProblemType :: ConstantVelocity;
  typedef typename GridPartType :: GridType                        GridType;
  typedef NavierStokesModelTraits< GridPartType,ProblemImp::dimRange >  Traits;
  static const int dimDomain = Traits :: dimDomain ;
  static const int dimRange  = Traits :: dimRange ;
  typedef typename Traits :: DomainType                          DomainType;
  typedef typename Traits :: RangeType                           RangeType;
  typedef typename Traits :: GradientType                        GradientType;
  typedef typename Traits :: FluxRangeType                       FluxRangeType;
  typedef typename Traits :: DiffusionRangeType                  DiffusionRangeType;
  typedef typename Traits :: DiffusionMatrixType                 DiffusionMatrixType;
  typedef typename Traits :: FaceDomainType                      FaceDomainType;
  typedef typename Traits :: JacobianRangeType                   JacobianRangeType;

  typedef typename Traits :: EntityType                          EntityType;
  typedef typename Traits :: IntersectionType                    IntersectionType;

  NavierStokesModel(const NavierStokesModel& otehr);
  const NavierStokesModel &operator=(const NavierStokesModel &other);
public:
  /**
   * @brief Constructor
   *
   * initializes model parameter
   *
   * @param problem Class describing the initial(t=0) and exact solution
   */
  NavierStokesModel(const ProblemType& problem, const bool rightHandSideModel ) :
    problem_(problem),
    epsilon_(problem.epsilon()),
    tstepEps_( getTStepEps() ),
    rightHandSideModel_( rightHandSideModel )
  {}

  inline const ProblemType& problem() const { return problem_; }

  inline bool hasFlux() const { return true ; }
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

  struct ComputeLaplaceU
  {
    typedef laplaceUVar VarId;
    typedef RangeType  ReturnType;

    template <class LocalEvaluation>
    RangeType operator() (const LocalEvaluation& local) const
    {
      return RangeType( 0 );
    }
  };

  struct ComputeGradientPressure
  {
    typedef gradientPressureVar VarId;
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
                             const JacobianRangeType& jac,
                             RangeType & s) const
  {
    DomainType xgl = local.entity().geometry().global( local.point() );
    // right hand side f
    problem_.stiffSource( xgl, local.time(), u, s );

    if( ! rightHandSideModel_ )
    {
      // + \alpha \mu \Delta u^n+\theta
      RangeType laplaceU ( local.evaluate( ComputeLaplaceU, local ) );
      laplaceU *= problem_.alphaMu();
      s += laplaceU;

      // + \nable p^n+\theta
      s += local.evaluate( ComputeGradientPressure, local );
    }
    else
    {
      // + u^n * \theta \Delta t
      RangeType uS ( u );
      uS *= theta_;
      s += uS;
    }
  }

  struct ComputeVelocity
  {
    typedef velocityVar VarId;
    typedef DomainType  ReturnType;

    template <class LocalEvaluation>
    const RangeType& operator() (const LocalEvaluation& local, const RangeType& u ) const
    {
      return u;
    }
  };

  /**
   * @brief advection term \f$F\f$
   *
   * @param en entity on which to evaluate the advection term
   * @param time current time of TimeProvider
   * @param x coordinate local to entity
   * @param u \f$U\f$
   * @param f \f$f(U)\f$
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
   * @brief velocity calculation, is called by advection()
   */
  template <class LocalEvaluation>
  inline DomainType velocity(const LocalEvaluation& local, const RangeType& u ) const
  {
    return local.evaluate( ComputeVelocity(), local, u);
  }


  /**
   * @brief diffusion term \f$a\f$
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
   * @brief diffusion term \f$A\f$
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
    A *= problem_.betaMu();// diffusion( u, A );//*(1.+d);
  }

  template <class LocalEvaluation>
  inline double diffusionTimeStep(const LocalEvaluation& local,
                                  const double circumEstimate,
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
   * @brief checks for existence of dirichlet boundary values
   */
  template <class LocalEvaluation>
  inline bool hasBoundaryValue(const LocalEvaluation& local ) const
  {
    return true;
  }

  /**
   * @brief neuman boundary values \f$g_N\f$ for pass2
   */
  template <class LocalEvaluation>
  inline double boundaryFlux(const LocalEvaluation& local,
                             const RangeType& uLeft,
                             const GradientType& vLeft,
                             RangeType& gLeft) const
  {
    gLeft = 0.;
    return 0.;
  }

  /**
   * @brief neuman boundary values \f$g_N\f$ for pass1
   */
  template <class LocalEvaluation>
  inline double boundaryFlux(const LocalEvaluation& local,
                             const RangeType& uLeft,
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
   * @brief dirichlet boundary values
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
    DomainType xgl = local.intersection().geometry().global( local.localPoint() );
    problem_.evaluate(xgl, local.time(), uRight);
  }


  // here x is in global coordinates
  template <class LocalEvaluation>
  inline void maxSpeed( const LocalEvaluation& local,
                        const DomainType& normal,
                        const RangeType& u,
                        double& advspeed,
                        double& totalspeed ) const
  {
    const DomainType& v = velocity( local );
    advspeed   = std::abs( v * normal );
    totalspeed = advspeed;
  }

 protected:
  double getTStepEps() const
  {
    // if diffusionTimeStep is set to non-zero in the parameterfile, the
    // deltaT in the timeprovider is updated according to the diffusion
    // parameter epsilon.
    bool diff_tstep;
    Dune::Fem::Parameter::get("femdg.stepper.diffusiontimestep", diff_tstep);
    return diff_tstep ? epsilon_ : 0;
  }

 protected:
  const ProblemType& problem_;
  const double epsilon_;
  const double tstepEps_;
  const bool rightHandSideModel_;
};


#if 0
/**
 * @brief defines the advective flux
 */
template <class ModelType>
class UpwindFlux {
public:
  typedef ModelType Model;
  typedef typename Model::Traits Traits;
  enum { dimRange = Model::dimRange };
  typedef typename Model :: DomainType DomainType;
  typedef typename Model :: RangeType RangeType;
  typedef typename Model :: JacobianRangeType JacobianRangeType;
  typedef typename Model :: FluxRangeType FluxRangeType;
  typedef typename Model :: DiffusionRangeType DiffusionRangeType;
  typedef typename Model :: FaceDomainType  FaceDomainType;
  typedef typename Model :: EntityType  EntityType;
  typedef typename Model :: IntersectionType  IntersectionType;

public:
  /**
   * @brief Constructor
   */
  UpwindFlux(const Model& mod) : model_(mod) {}

  static std::string name () { return "UpwindFlux"; }

  const Model& model() const {return model_;}

  /**
   * @brief evaluates the flux \f$g(u,v)\f$
   *
   * @return maximum wavespeed * normal
   */
  template <class LocalEvaluation>
  inline double numericalFlux( const LocalEvaluation& left,
                               const LocalEvaluation& right,
                               const RangeType& uLeft,
                               const RangeType& uRight,
                               const JacobianRangeType& jacLeft,
                               const JacobianRangeType& jacRight,
                               RangeType& gLeft,
                               RangeType& gRight ) const
  {
    const FaceDomainType& x = left.localPoint();

    // get normal from intersection
    const DomainType normal = left.intersection().integrationOuterNormal(x);

    // get velocity
    const DomainType v = model_.velocity( left );
    const double upwind = normal * v;

    if (upwind>0)
      gLeft = uLeft;
    else
      gLeft = uRight;
    gLeft *= upwind;
    gRight = gLeft;
    return std::abs(upwind);
  }
protected:
  const Model& model_;
};
#endif

#endif
