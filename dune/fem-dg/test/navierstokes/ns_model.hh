#ifndef NS_MODEL_HH
#define NS_MODEL_HH

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include "thermodynamics.hh"
#include "ns_model_spec.hh"


namespace Dune {

//////////////////////////////////////////////////////
//
//                 NAVIER-STOKES EQNS
//
//////////////////////////////////////////////////////
template< class GridPart, class Problem >
class NSModelTraits
  : public Fem :: FunctionSpace< double, double,
                                 GridPart::GridType::dimensionworld,
                                 GridPart::GridType::dimensionworld + 2 >
{
  typedef Fem :: FunctionSpace< double, double,
                                 GridPart::GridType::dimensionworld,
                                 GridPart::GridType::dimensionworld + 2 > BaseType;

public:
  typedef Problem  ProblemType;
  typedef GridPart GridPartType;
  typedef typename GridPart :: GridType                     GridType;

  enum { dimDomain = GridType :: dimensionworld };
  enum { dimRange = dimDomain + 2 };
  enum { dimGradRange = dimRange * dimDomain };
  enum { dimGradient = dimDomain + 1 };

  typedef BaseType FunctionSpaceType ;

  typedef typename FunctionSpaceType :: DomainFieldType    DomainFieldType ;
  typedef FieldVector< DomainFieldType, FunctionSpaceType::dimDomain - 1 >    FaceDomainType;

  typedef FieldVector< double, dimGradRange >               GradientType;
  typedef typename FunctionSpaceType :: JacobianRangeType   FluxRangeType;

  typedef FieldVector< double, dimGradRange >               GradientRangeType;
  typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;

  typedef typename GridPart :: IntersectionIteratorType     IntersectionIterator;
  typedef typename IntersectionIterator::Intersection       Intersection;
  typedef Intersection       IntersectionType;
  typedef typename GridPartType :: template Codim<0> :: EntityType  EntityType;

  typedef Thermodynamics< dimDomain >                       ThermodynamicsType;

  typedef MinModLimiter< DomainFieldType > LimiterFunctionType;
};



template< class GridPartType , class ProblemImp >
class NSModel : public DefaultModel < NSModelTraits< GridPartType, ProblemImp > >
{
  public:
  typedef ProblemImp                                        ProblemType;
  typedef typename GridPartType::GridType                   GridType;
  typedef NSModelTraits< GridPartType, ProblemType >        Traits;

  typedef NSFlux< Traits >  FluxType ;

  enum { dimDomain = Traits :: dimDomain };
  enum { dimRange  = Traits :: dimRange  };
  enum { dimGradRange = Traits::dimGradRange } ;

  typedef typename Traits :: EntityType                     EntityType;
  typedef typename Traits :: IntersectionIterator           IntersectionIterator;
  typedef typename Traits :: Intersection                   IntersectionType;
  typedef typename Traits :: FaceDomainType                 FaceDomainType;

  typedef typename Traits :: DomainType                     DomainType;
  typedef typename Traits :: RangeType                      RangeType;
  typedef typename Traits :: GradientRangeType              GradientRangeType;
  typedef typename Traits :: JacobianRangeType              JacobianRangeType;
  typedef typename Traits :: JacobianFluxRangeType          JacobianFluxRangeType;
  typedef typename Traits :: ThermodynamicsType             ThermodynamicsType;

  // for heat equations advection is disabled
  static const bool hasAdvection = true ;
  static const bool hasDiffusion = true ;

 public:
  NSModel( const ProblemType& problem )
    : thermodynamics_( problem.thermodynamics() )
    , problem_( problem )
    , nsFlux_( problem )
    , alpha_( std::pow( problem.gamma(), 1.5 ) * (problem.Re_inv() * problem.Pr_inv()) )
  {
  }

  double gamma() const { return problem_.gamma() ; }

  inline bool hasStiffSource() const { return problem_.hasStiffSource(); }
  inline bool hasNonStiffSource() const { return problem_.hasNonStiffSource(); }
  inline bool hasFlux() const { return true ; }

  template <class LocalEvaluation>
  inline double stiffSource( const LocalEvaluation& local,
                             const RangeType& u,
                             const JacobianRangeType& jac,
                             RangeType& s ) const
  {
    // some special RHS for testcases/NSWaves
    const DomainType& xgl = local.entity().geometry().global( local.point() );
    return problem_.stiffSource( local.time(), xgl, u, s );
  }

  template <class LocalEvaluation>
  inline double nonStiffSource( const LocalEvaluation& local,
                                const RangeType& u,
                                const JacobianRangeType& jac,
                                RangeType& s ) const
  {
    // some special RHS for testcases/NSWaves
    const DomainType& xgl = local.entity().geometry().global( local.point() );
    return problem_.nonStiffSource( local.time(), xgl, u, s );
  }

  template <class LocalEvaluation>
  inline void advection( const LocalEvaluation& local,
                         const RangeType& u,
                         JacobianRangeType& f ) const
  {
    nsFlux_.analyticalFlux( u, f );
  }

  /////////////////////////////////////////////////////////////////
  // Limiter section
  ////////////////////////////////////////////////////////////////
  inline void velocity( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        DomainType& velocity) const
  {
    for(int i=0; i<dimDomain; ++i)
    {
      // we store \rho u but do not need to divide by \rho here since only
      // sign is needed.
      velocity[i] = u[i+1];
    }
  }

  // adjust average value if necessary
  // (e.g. transform from conservative to primitive variables )
  void adjustAverageValue( const EntityType& entity,
                           const DomainType& xLocal,
                           RangeType& u ) const
  {
    // nothing to be done here for this test case
  }

  // calculate jump between left and right value
  inline void jump(const IntersectionType& it,
                   const double time,
                   const FaceDomainType& x,
                   const RangeType& uLeft,
                   const RangeType& uRight,
                   RangeType& jump) const
  {
    const double pl = problem_.pressure( uLeft );
    const double pr = problem_.pressure( uRight );

    // take pressure as shock detection values
    jump  = (pl-pr)/(0.5*(pl+pr));
  }

  // calculate jump between left and right value
  inline void adaptationIndicator(
                   const IntersectionType& it,
                   const double time,
                   const FaceDomainType& x,
                   const RangeType& uLeft,
                   const RangeType& uRight,
                   RangeType& indicator) const
  {
    jump( it, time, x, uLeft, uRight, indicator );
  }


  //! \brief returns true if physical check does something useful */
  bool hasPhysical () const { return true; }

  template <class LocalDomainType>
  bool physical( const EntityType& en,
                 const LocalDomainType& x,
                 const RangeType& u) const
  {
    if (u[0]<1e-8)
      return false;
    else
    {
      const double p = problem_.pressure( u );
      return p > 1e-8;
    }


    return u[ 0 ] > 1e-8 ;
  }


  template <class LocalEvaluation>
  inline double diffusionTimeStep( const LocalEvaluation& local,
                                   const double circumEstimate,
                                   const RangeType& u ) const
  {
    // look at Ch. Merkle Diplom thesis, pg. 38
    // get value of mu at temperature T

    // get mu
    const double mu = problem_.mu( u );

    // ksi = 0.25
    return mu * circumEstimate * alpha_ / (0.25 * u[0] * local.volume());
  }


  //! return analyticalFlux for 1st pass
  template <class LocalEvaluation>
  inline void jacobian( const LocalEvaluation& local,
                        const RangeType& u,
                        JacobianFluxRangeType& a ) const
  {
    nsFlux_.jacobian( u, a );
  }


  template <class LocalEvaluation>
  inline bool hasBoundaryValue( const LocalEvaluation& local ) const
  {
    return true;
  }


  template <class LocalEvaluation>
  inline double boundaryFlux( const LocalEvaluation& local,
                              const RangeType& uLeft,
                              RangeType& gLeft ) const
  {
    abort();
    gLeft = 0;
    return 0.0;
  }

  template <class LocalEvaluation>
  inline double boundaryFlux( const LocalEvaluation& local,
                              const RangeType& uLeft,
                              const GradientRangeType& duLeft,
                              RangeType& gLeft ) const
  {
    abort();
    DomainType xgl= local.intersection().intersectionGlobal().global( local.localPoint() );
    const DomainType normal = local.intersection().integrationOuterNormal( local.localPoint() );
    double p;
    double T;
    pressAndTemp( uLeft, p, T );
    gLeft = 0;

    // bnd. cond. from euler part
    for (int i=0;i<dimDomain; ++i)
      gLeft[i+1] = normal[i]*p;

    return 0.0;
  }

  /** \brief boundary flux for the diffusion part
   */
  template <class LocalEvaluation>
  inline double diffusionBoundaryFlux( const LocalEvaluation& local,
                                       const RangeType& uLeft,
                                       const JacobianRangeType& jacLeft,
                                       RangeType& gLeft ) const
  {
    std::cerr <<"diffusionBoundaryFlux shouldn't be used for this testcase" <<std::endl;
    abort();
  }


  template <class LocalEvaluation>
  inline void boundaryValue( const LocalEvaluation& local,
                             const RangeType& uLeft,
                             RangeType& uRight ) const
  {
    const DomainType xgl = local.intersection().geometry().global( local.localPoint() );
    problem_.evaluate( local.time(), xgl, uRight );
  }

  // here x is in global coordinates
  template <class LocalEvaluation>
  inline void maxSpeed( const LocalEvaluation& local,
                        const DomainType& normal,
                        const RangeType& u,
                        double& advspeed,
                        double& totalspeed ) const
  {
    advspeed   = nsFlux_.maxSpeed( normal , u );
    totalspeed = advspeed;
  }


  template <class LocalEvaluation>
  void diffusion( const LocalEvaluation& local,
                  const RangeType& u,
                  const JacobianRangeType& jac,
                  JacobianRangeType& diff ) const
  {
    nsFlux_.diffusion( u, jac, diff );
  }

  inline void pressAndTemp( const RangeType& u, double& p, double& T ) const
  {
    thermodynamics_.pressAndTempEnergyForm( u, p, T );
  }

  inline void conservativeToPrimitive( const double time,
                                       const DomainType& xgl,
                                       const RangeType& cons,
                                       RangeType& result,
                                       bool ) const
  {
    problem_.evaluate( time, xgl, result );
    //thermodynamics_.conservativeToPrimitiveEnergyForm( cons, result );
  }


  inline const ProblemType& problem() const { return problem_; }

 protected:
  const ThermodynamicsType& thermodynamics_;
  const ProblemType& problem_;
  const FluxType  nsFlux_;
  const double alpha_;
};

} // end namespace Dune

#endif // file definition
