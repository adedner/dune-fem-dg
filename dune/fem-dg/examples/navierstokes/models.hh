#ifndef NS_MODEL_HH
#define NS_MODEL_HH

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include "thermodynamics.hh"
#include "navierstokesflux.hh"

#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/l1eocerror.hh>

namespace Dune
{
namespace Fem
{


  template <class GridImp, class ProblemImp >
  class NSModelTraits
    : public DefaultModelTraits< GridImp, ProblemImp >
  {
    typedef DefaultModelTraits< GridImp, ProblemImp >           BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::RangeFieldType, BaseType::dimGradRange >
                                                                    GradientType;

    typedef MinModLimiter< typename BaseType::RangeFieldType >     LimiterFunctionType;
    //typedef Fem::MinModLimiter< FieldType > LimiterFunctionType ;
    //typedef SuperBeeLimiter< FieldType > LimiterFunctionType ;
    //typedef VanLeerLimiter< FieldType > LimiterFunctionType ;

    typedef Thermodynamics< BaseType::dimDomain, typename BaseType::RangeFieldType >  ThermodynamicsType;
  };


  /**
   *  \brief Analytical Model for the compressible Navier-Stokes equations
   *
   *  \ingroup AnalyticalModels
   */
  template< class GridImp, class ProblemImp >
  class NSModel
    : public DefaultModel < NSModelTraits< GridImp, ProblemImp > >
  {
    public:
    typedef NSModelTraits< GridImp, ProblemImp >       Traits;
    typedef DefaultModel< Traits >                     BaseType;

    typedef GridImp                                    GridType;
    typedef typename Traits::ProblemType               ProblemType;

    typedef NSFlux< Traits >                           FluxType;

    enum { dimDomain = Traits::dimDomain };
    enum { dimRange  = Traits::dimRange  };
    enum { dimGradRange = Traits::dimGradRange };

    typedef typename Traits::DomainType                DomainType;
    typedef typename Traits::FaceDomainType            FaceDomainType;
    typedef typename Traits::RangeType                 RangeType;
    typedef typename Traits::RangeFieldType            RangeFieldType;
    typedef typename Traits::GradientType              GradientType;
    typedef typename Traits::JacobianRangeType         JacobianRangeType;
    typedef typename Traits::DiffusionRangeType        DiffusionRangeType;

    typedef typename Traits::ThermodynamicsType        ThermodynamicsType;

    // for heat equations advection is disabled
    static const bool hasAdvection = true;
    static const bool hasDiffusion = true;

    using BaseType::time;

   public:
    NSModel( const ProblemType& problem )
      : BaseType( problem.startTime() )
      , thermodynamics_( problem.thermodynamics() )
      , problem_( problem )
      , nsFlux_( problem )
      , alpha_( std::pow( problem.gamma(), 1.5 ) * (problem.Re_inv() * problem.Pr_inv()) )
    {
    }

    RangeFieldType gamma() const { return problem_.gamma() ; }

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
      const DomainType& xgl = local.entity().geometry().global( local.position() );
      return problem_.stiffSource( time(), xgl, u, s );
    }

    template <class LocalEvaluation>
    inline double nonStiffSource( const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& jac,
                                  RangeType& s ) const
    {
      // some special RHS for testcases/NSWaves
      const DomainType& xgl = local.entity().geometry().global( local.position() );
      return problem_.nonStiffSource( time(), xgl, u, s );
    }

    template <class LocalEvaluation>
    inline void advection( const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& jacu,
                           JacobianRangeType& f ) const
    {
      nsFlux_.analyticalFlux( u, f );
    }

    /////////////////////////////////////////////////////////////////
    // Limiter section
    ////////////////////////////////////////////////////////////////
    template< class Entity >
    inline void velocity( const Entity& en,
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
    template< class Entity >
    void adjustAverageValue( const Entity& entity,
                             const DomainType& xLocal,
                             RangeType& u ) const
    {
      // nothing to be done here for this test case
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void jump(const Intersection& it,
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
    template< class Intersection >
    inline void adaptationIndicator(
                     const Intersection& it,
                     const FaceDomainType& x,
                     const RangeType& uLeft,
                     const RangeType& uRight,
                     RangeType& indicator) const
    {
      jump( it, x, uLeft, uRight, indicator );
    }


    //! \brief returns true if physical check does something useful */
    bool hasPhysical () const { return true; }

    template <class Entity, class LocalDomainType>
    bool physical( const Entity& en,
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
                                     const RangeType& u ) const
    {
      // look at Ch. Merkle Diplom thesis, pg. 38
      // get value of mu at temperature T

      // get mu
      const RangeFieldType mu = problem_.mu( u );

      // ksi = 0.25
      return mu * alpha_ / (0.25 * u[0]);
    }


    //! return analyticalFlux for 1st pass
    template <class LocalEvaluation>
    inline void jacobian( const LocalEvaluation& local,
                          const RangeType& u,
                          DiffusionRangeType& a ) const
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
                                const JacobianRangeType&,
                                RangeType& gLeft ) const
    {
      abort();
      gLeft = 0;
      return 0.0;
    }

    template <class LocalEvaluation>
    inline double boundaryFlux( const LocalEvaluation& local,
                                const RangeType& uLeft,
                                const GradientType& duLeft,
                                RangeType& gLeft ) const
    {
      abort();
      DomainType xgl= local.intersection().intersectionGlobal().global( local.localPosition() );
      const DomainType normal = local.intersection().integrationOuterNormal( local.localPosition() );
      RangeFieldType p;
      RangeFieldType T;
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
      const DomainType xgl = local.intersection().geometry().global( local.localPosition() );
      problem_.evaluate( time(), xgl, uRight );
    }

    // here x is in global coordinates
    template <class LocalEvaluation>
    inline void maxWaveSpeed( const LocalEvaluation& local,
                          const DomainType& normal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      advspeed   = nsFlux_.maxWaveSpeed( normal , u );
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

    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
      //default version!
      EOCErrorList::setErrors<L2EOCError>( *this, df );
    }

    inline const ProblemType& problem() const { return problem_; }

   protected:
    const ThermodynamicsType& thermodynamics_;
    const ProblemType& problem_;
    const FluxType  nsFlux_;
    const RangeFieldType alpha_;
  };

} // end namespace
} // end namespace

#endif // file definition
