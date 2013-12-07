#ifndef HEAT_MODELS_HH
#define HEAT_MODELS_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/io/parameter.hh>


#include <dune/fem-dg/operator/limiter/limitpass.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>


/**********************************************
 * Analytical model                           *
 *********************************************/

/**
 * @brief Traits class for HeatEqnModel
 */
template <class GridPart > 
class HeatEqnModelTraits
{
public:
  typedef GridPart                                                   GridPartType;
  typedef typename GridPartType :: GridType                          GridType;
  static const int dimDomain = GridType::dimensionworld;
  static const int dimRange = DIMRANGE;
  static const int dimGradRange = dimRange * dimDomain ;
  // Definition of domain and range types
  typedef Dune::FieldVector< double, dimDomain >                     DomainType;
  typedef Dune::FieldVector< double, dimDomain-1 >                   FaceDomainType;
  typedef Dune::FieldVector< double, dimRange >                      RangeType;
  typedef Dune::FieldVector< double, dimGradRange >                  GradientType;
  // ATTENTION: These are matrices (c.f. HeatEqnModel)
  typedef Dune::FieldMatrix< double, dimRange, dimDomain >           FluxRangeType;
  typedef Dune::FieldMatrix< double, dimRange, dimDomain >           JacobianRangeType;
  typedef Dune::FieldMatrix< double, dimGradRange, dimDomain >       DiffusionRangeType;
  typedef Dune::FieldMatrix< double, dimDomain, dimDomain >          DiffusionMatrixType;
  typedef typename GridType :: template Codim< 0 > :: Entity         EntityType;
  typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
  typedef typename IntersectionIterator :: Intersection              IntersectionType;

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
//      dx(u) + div(uV) - epsilon*lap(u)) = 0
//  where V is constant vector
//
////////////////////////////////////////////////////////
template <class GridPartType,class ProblemImp>
class HeatEqnModel : public DefaultModel < HeatEqnModelTraits< GridPartType > >
{
public:
  typedef ProblemImp ProblemType ;

  static const int ConstantVelocity = ProblemType :: ConstantVelocity;
  typedef typename GridPartType :: GridType                          GridType;
  typedef HeatEqnModelTraits< GridPartType >        Traits;
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

  typedef typename Traits :: EntityType                       EntityType;
  typedef typename Traits :: IntersectionType                 IntersectionType;

  HeatEqnModel(const HeatEqnModel& otehr);
  const HeatEqnModel &operator=(const HeatEqnModel &other);
public:
  /**
   * @brief Constructor
   *
   * initializes model parameter
   *
   * @param problem Class describing the initial(t=0) and exact solution
   */
  HeatEqnModel(const ProblemType& problem) :
    problem_(problem),
    epsilon_(problem.epsilon()),
    tstepEps_( getTStepEps() ),
    velocity_( getVelocity() )
  {}

  inline bool hasFlux() const { return true ; }
  inline const ProblemType& problem() const { return problem_; }
  inline bool hasStiffSource() const { return problem_.hasStiffSource(); }
  inline bool hasNonStiffSource() const { return problem_.hasNonStiffSource(); }

  inline double nonStiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const GradientType& du,
                        RangeType & s) const
  {
    //FieldMatrixConverter< GradientType, JacobianRangeType> jac( du );
    return nonStiffSource( en, time, x, u, s );
  }

  inline double nonStiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const JacobianRangeType& jac,
                        RangeType & s) const
  {
    return nonStiffSource( en, time, x, u, s );
  }

  inline double nonStiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        RangeType & s) const
  {
    DomainType xgl = en.geometry().global( x );
    return problem_.nonStiffSource( xgl, time, u, s );
  }

  inline double stiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const GradientType& du,
                        RangeType & s) const
  {
    //FieldMatrixConverter< GradientType, JacobianRangeType> jac( du );
    return stiffSource( en, time, x, u, s );
  }


  inline double stiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const JacobianRangeType& jac,
                        RangeType & s) const
  {
    return stiffSource( en, time, x, u, s );
  }

  inline double stiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        RangeType & s) const
  {
    DomainType xgl = en.geometry().global( x );
    return problem_.stiffSource( xgl, time, u, s );
  }

  /**
   * @brief advection term \f$F\f$
   *
   * @param en entity on which to evaluate the advection term
   * @param time current time of TimeProvider
   * @param x coordinate local to entity
   * @param u \f$U\f$
   * @param f \f$f(U)\f$
   */
  inline  void advection(const EntityType& en,
                         const double time,
                         const DomainType& x,
                         const RangeType& u,
                         FluxRangeType & f) const
  {
    // evaluate velocity V
    DomainType v;
    velocity( en, time, x, u, v );

    //f = uV;
    for( int r=0; r<dimRange; ++r )
      for( int d=0; d<dimDomain; ++d )
        f[r][d] = v[ d ] * u[ r ];
  }

  /**
   * @brief velocity calculation, is called by advection()
   */
  inline  void velocity(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u, 
                        DomainType& v) const
  {
    velocity( en.geometry().global(x), time, v );
  }

  /*
   * @brief velocity calculation, is called by advection()
   */
  inline  void velocity(const DomainType& xGlobal,
                        const double time,
                        DomainType& v) const
  {
    problem_.velocity(xGlobal, time, v);
  }

  /**
   * @brief diffusion term \f$a\f$
   */
  inline void jacobian(const EntityType& en,
                        const double time,
                        const DomainType& x,
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

  inline void eigenValues(const EntityType& en,
                          const double time,
                          const DomainType& x,
                          const RangeType& u,
                          RangeType& maxValue) const
  {
    FluxRangeType A(0);
    maxValue = RangeType( problem_.diffusion( u,A ) );
    return;
    abort();
    /*
    DiffusionMatrixType K(0) ;
    DomainType values ;
    for( int r = 0; r<dimRange; ++r )
    {
      // setup matrix 
      for(int i=0; i<dimDomain; ++i) 
        K[i][i] = problem_.diffusion( u,A );

      // calculate eigenvalues 
      Dune::FMatrixHelp :: eigenValues( K, values );
      // take max eigenvalue 
      maxValue[ r ] = values.infinity_norm();
    }
    */
  }

  /**
   * @brief diffusion term \f$A\f$
   */
  template <class JacobianType>
  inline void diffusion(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const JacobianType& jac,
                        FluxRangeType& A) const
  {
    // copy v to A 
    A = jac;

    // apply diffusion coefficient 
    //double d =  0.; //(1.-en.geometry().global(x)[0])*en.geometry().global(x)[0]+
    //                //(1.-en.geometry().global(x)[1])*en.geometry().global(x)[1];
    A *= problem_.diffusion( u, A );//*(1.+d);
  }

  inline void diffusion(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const GradientType& vecJac,
                        FluxRangeType& A) const
  {
    Dune::Fem::FieldMatrixConverter< GradientType, FluxRangeType> jac( vecJac );
    diffusion( en, time, x, u, jac, A );
  }

  inline double diffusionTimeStep( const IntersectionType &it,
                                   const double enVolume,
                                   const double circumEstimate,
                                   const double time,
                                   const FaceDomainType& x,
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
  inline bool hasBoundaryValue(const IntersectionType& it,
                               const double time,
                               const FaceDomainType& x) const
  {
    return true;
  }

  /**
   * @brief neuman boundary values \f$g_N\f$ for pass2
   */
  inline double boundaryFlux(const IntersectionType& it,
                             const double time,
                             const FaceDomainType& x,
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
  inline double boundaryFlux(const IntersectionType& it,
                             const double time,
                             const FaceDomainType& x,
                             const RangeType& uLeft,
                             RangeType& gLeft) const
  {
    gLeft = 0.;
    return 0.;
  }

  inline double diffusionBoundaryFlux( const IntersectionType& it,
                                       const double time,
                                       const FaceDomainType& x,
                                       const RangeType& uLeft,
                                       const GradientType& gradLeft,
                                       RangeType& gLeft ) const  
  {
    Dune::Fem::FieldMatrixConverter< GradientType, JacobianRangeType> jacLeft( gradLeft );
    return diffusionBoundaryFlux( it, time, x, uLeft, jacLeft, gLeft );
  }

  /** \brief boundary flux for the diffusion part
   */
  template <class JacobianRangeImp>
  inline double diffusionBoundaryFlux( const IntersectionType& it,
                                       const double time,
                                       const FaceDomainType& x,
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
  inline  void boundaryValue(const IntersectionType& it,
                             const double time,
                             const FaceDomainType& x,
                             const RangeType& uLeft,
                             RangeType& uRight) const
  {
#ifdef TESTOPERATOR
  uRight = 0;
  return;
#endif
    DomainType xgl = it.geometry().global( x );
    problem_.evaluate(xgl, time, uRight);
  }


  // here x is in global coordinates
  inline void maxSpeed( const EntityType& entity,
                        const double time,
                        const DomainType& xGlobal,
                        const DomainType& normal,
                        const RangeType& u,
                        double& advspeed,
                        double& totalspeed ) const
  {
    DomainType v;
    problem_.velocity( xGlobal, time, v );
    advspeed   = v * normal ;
    totalspeed = advspeed;
  }

 protected:
  DomainType getVelocity() const 
  {
    DomainType velocity(0);
    if(ConstantVelocity) {
      problem_.velocity(velocity,velocity);
    }
    return velocity;  
  }
  double getTStepEps() const 
  {
    // if diffusionTimeStep is set to non-zero in the parameterfile, the
    // deltaT in the timeprovider is updated according to the diffusion
    // parameter epsilon.
    bool diff_tstep;
    Dune::Fem::Parameter::get("femhowto.diffusionTimeStep", diff_tstep);
    return diff_tstep ? epsilon_ : 0;
  }

 protected:
  const ProblemType& problem_;
  const double epsilon_;
  const double tstepEps_;
 public:
  const DomainType velocity_;

};


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
  typedef typename Model :: FluxRangeType FluxRangeType;
  typedef typename Model :: DiffusionRangeType DiffusionRangeType;
  typedef typename Model :: FaceDomainType  FaceDomainType;
  typedef typename Model :: EntityType  EntityType;
  typedef typename Model :: IntersectionType  IntersectionType;
protected:
  template <class Model, bool constVelo>
    struct Velocity
    {
      /**
       * @brief computes and returns the wind direction
       */
      static inline double upwind(const Model& model,
                                  const IntersectionType& it,
                                  const double time,
                                  const FaceDomainType& x,
                                  const RangeType& uLeft)
      {
        const DomainType normal = it.integrationOuterNormal(x);
        DomainType velocity;
        model.velocity(*it.inside(),time,
                       it.geometryInInside().global(x),
                       uLeft,velocity);
        return normal*velocity;
      }
    };

  template <class Model>
    struct Velocity<Model,true>
    {
      /**
       * @brief computes and returns the wind direction for models with
       * constant velocity
       */
      static inline double upwind( const Model& model,
                                   const IntersectionType& it,
                                   const double time,
                                   const FaceDomainType& x,
                                   const RangeType& uLeft )
      {
        const DomainType normal = it.integrationOuterNormal(x);
        return normal * model.velocity_;
      }
    };

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
  template <class QuadratureImp>
  inline double numericalFlux( const IntersectionType& it,
                               const EntityType& inside,
                               const EntityType& outside,
                               const double time,
                               const QuadratureImp& faceQuadInner,
                               const QuadratureImp& faceQuadOuter,
                               const int quadPoint,
                               const RangeType& uLeft,
                               const RangeType& uRight,
                               RangeType& gLeft,
                               RangeType& gRight ) const
  {
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
    const double upwind = Velocity<Model,Model::ConstantVelocity>::
      upwind(model_,it,time,x,uLeft);

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