#ifndef DUNE_EULERMODEL_HH
#define DUNE_EULERMODEL_HH

// system includes
#include <config.h>
#include <cmath>
#include <type_traits>

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/space/common/functionspace.hh>

#include "../navierstokes/thermodynamics.hh"

#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/rotator.hh>
#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include <dune/fem-dg/operator/fluxes/analyticaleulerflux.hh>

namespace Dune
{
namespace Fem
{
  template< class GridPartImp, class ProblemImp >
  class EulerModelTraits
    : public DefaultModelTraits< GridPartImp, ProblemImp >
  {
    typedef DefaultModelTraits< GridPartImp, ProblemImp >           BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                    GradientType;

    typedef std::tuple<>                                            ModelParameter;

    typedef MinModLimiter< typename BaseType::DomainFieldType >     LimiterFunctionType ;
    //typedef SuperBeeLimiter< typename BaseType::DomainFieldType > LimiterFunctionType ;
    //typedef VanLeerLimiter< typename BaseType::DomainFieldType >  LimiterFunctionType ;

    typedef Thermodynamics< BaseType::dimDomain >                   ThermodynamicsType;
  };


  /**
   * \brief Euler equations for dry atmosphere
   *
   * \ingroup AnalyticalModels
   */
  template< class GridPartType, class ProblemImp >
  class EulerModel :
    public DefaultModel< EulerModelTraits< GridPartType, ProblemImp > >
  {
  public:
    typedef EulerModelTraits< GridPartType, ProblemImp > Traits;
    typedef typename Traits::ProblemType                 ProblemType;

    enum { dimDomain = Traits::dimDomain };
    enum { dimRange = Traits::dimRange };

    typedef typename Traits::GridType                    GridType;
    typedef typename Traits::EntityType                  EntityType;
    typedef typename Traits::IntersectionType            IntersectionType;
    typedef typename Traits::FaceDomainType              FaceDomainType;

    typedef typename Traits::RangeType                   RangeType;
    typedef typename Traits::RangeFieldType              RangeFieldType ;
    typedef typename Traits::DomainType                  DomainType;
    typedef typename Traits::FluxRangeType               FluxRangeType;
    typedef typename Traits::GradientType                GradientType;
    typedef typename Traits::JacobianRangeType           JacobianRangeType;
    typedef typename Traits::DiffusionRangeType          DiffusionRangeType;

    typedef typename Traits::ThermodynamicsType          ThermodynamicsType;

    // for Euler equations diffusion is disabled
    static const bool hasAdvection = true;
    static const bool hasDiffusion = false;

   public:
    EulerModel( const ProblemType& problem )
      : gamma_( problem.gamma() )
      , problem_( problem )
      , fieldRotator_( 1 ) // insert number of fist velocity component
    {}

    double gamma () const { return gamma_; }

    inline bool hasStiffSource() const { return false; }
    inline bool hasNonStiffSource() const { return false; }
    inline bool hasFlux() const { return true ; }

    template <class LocalEvaluation>
    inline double stiffSource( const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& du,
                               RangeType & s) const
    {
      s = 0;
      return 0;
    }


    template< class LocalEvaluation >
    inline double nonStiffSource( const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& jac,
                                  RangeType& s) const
    {
      s = 0;
      return 0;
    }

    inline double pressure( const RangeType& u ) const
    {
      return EulerAnalyticalFlux< dimDomain >().pressure( gamma_ , u );
    }

    inline void conservativeToPrimitive( const double time,
                                         const DomainType& xgl,
                                         const RangeType& cons,
                                         RangeType& prim,
                                         const bool ) const
    {
      problem_.evaluate( xgl, time, prim );
      //thermodynamics_.conservativeToPrimitiveEnergyForm( cons, prim );
    }

    template <class LocalEvaluation>
    inline void advection( const LocalEvaluation& local,
                           const RangeType& u,
                           const FluxRangeType& jacu,
                           FluxRangeType& f ) const
    {
      EulerAnalyticalFlux<dimDomain>().analyticalFlux( gamma_ , u , f );
    }

    template <class LocalEvaluation>
    inline void eigenValues(const LocalEvaluation& local,
                            const RangeType& u,
                            RangeType& maxValue) const
    {
      std::cerr <<"eigenValues for problems/euler not implemented\n";
      abort();
    }

    template <class LocalEvaluation>
    inline double diffusionTimeStep( const LocalEvaluation& local,
                                     const double circumEstimate,
                                     const RangeType& u ) const
    {
      return 0;
    }

    // is not used
    template <class LocalEvaluation>
    inline  void jacobian( const LocalEvaluation& local,
                           const RangeType& u,
                           const FluxRangeType& du,
                           RangeType& A ) const
    {
      EulerAnalyticalFlux<dimDomain>().jacobian( gamma_ , u , du , A );
    }


    enum { Inflow = 1, Outflow = 2, Reflection = 3 , Slip = 4 };
    enum { MaxBnd = Slip };

    template <class LocalEvaluation>
    inline bool hasBoundaryValue( const LocalEvaluation& local ) const
    {
      const int bndId = problem_.boundaryId( local.intersection().boundaryId() );
      // on slip boundary we use boundaryFlux
      return bndId != Slip;
    }

    // return iRight for insertion into the numerical flux
    template <class LocalEvaluation>
    inline void boundaryValue( const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight ) const
    {
      // Neumann boundary condition
      //uRight = uLeft;
      // 5 and 6 is also Reflection
      //const int bndId = (it.boundaryId() > MaxBnd) ? MaxBnd : it.boundaryId();
      const int bndId = problem_.boundaryId( local.intersection().boundaryId() );

      assert( bndId > 0 );
      if( bndId == Inflow )
      {
        const DomainType xgl = local.intersection().geometry().global( local.localPosition() );
        problem_.evaluate(xgl, local.time(), uRight);
        return ;
      }
      else if ( bndId == Outflow )
      {
        uRight = uLeft;
        return ;
      }
      else if ( bndId == Reflection )
      {
        uRight = uLeft;
        const DomainType unitNormal = local.intersection().unitOuterNormal( local.localPosition() );
        fieldRotator_.rotateForth( uRight , unitNormal );
        // Specific for euler: opposite sign for first component of momentum
        uRight[1] = -uRight[1];
        fieldRotator_.rotateBack( uRight, unitNormal );
        return ;
      }
      /*
      else if ( bndId == Slip )
      {
        RangeType tmp(uLeft);
        const DomainType unitNormal = it.unitOuterNormal(x);
        this->rot_.rotateForth( tmp , unitNormal );
        tmp[1] = 0.;
        tmp[2] /= tmp[0];
        tmp[3] = pressure(uLeft);
        prim2cons(tmp,uRight);
        this->rot_.rotateBack( uRight, unitNormal );
        return ;
      }
      */
      else
      {
        uRight = uLeft;
        return ;
        assert( false );
        abort();
      }
    }

    // boundary condition here is slip boundary cond. <u,n>=0
    // gLeft= p*[0 n(global(x)) 0]
    template <class LocalEvaluation>
    inline double boundaryFlux( const LocalEvaluation& local,
                                const RangeType& uLeft,
                                const JacobianRangeType&,
                                RangeType& gLeft ) const
    {
      // Slip boundary condition
      const DomainType normal = local.intersection().integrationOuterNormal( local.localPosition() );

      const double p = EulerAnalyticalFlux< dimDomain >().pressure( gamma_ , uLeft );
      gLeft = 0;
      for ( int i = 0 ; i < dimDomain ; ++i )
        gLeft[i+1] = normal[i] * p;
      return 0.;
    }

    template <class LocalEvaluation>
    void diffusion( const LocalEvaluation& local,
                    const RangeType& u,
                    const JacobianRangeType& v,
                    JacobianRangeType& diff ) const
    {
    }


    /** \brief boundary flux for the diffusion part
     */
    template <class LocalEvaluation>
    inline double diffusionBoundaryFlux( const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const JacobianRangeType& jacLeft,
                                         RangeType& gLeft ) const
    {
      return 0;
    }

    // here x is in global coordinates
    template <class LocalEvaluation>
    inline void maxSpeed( const LocalEvaluation& local,
                          const DomainType& normal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      advspeed = EulerAnalyticalFlux< dimDomain >().maxSpeed( gamma_ , normal , u );
      totalspeed = advspeed;
    }

    inline const ProblemType& problem() const
    {
      return problem_;
    }

    const double gamma_;

    /**
     * \brief velocity calculation
     */
    template <class LocalEvaluation>
    inline DomainType velocity(const LocalEvaluation& local) const
    {
      std::integral_constant< int, 2 > uVar;
      DomainType v;
      for(int i=0; i<dimDomain; ++i)
        v[i] = local.values()[ uVar ][i];
      v *= EulerAnalyticalFlux<dimDomain>().rhoeps(local.values()[uVar]);
      return v;
    }


    /////////////////////////////////////////////////////////////////
    // Limiter section
    ////////////////////////////////////////////////////////////////
    inline void velocity (const EntityType& en,
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

    // we have physical check for this model
    bool hasPhysical() const
    {
      return true;
    }

    // calculate jump between left and right value
    inline bool physical(const EntityType& entity,
                         const DomainType& xGlobal,
                         const RangeType& u) const
    {
      if (u[0]<1e-8)
        return false;
      else
      {
        //std::cout << EulerAnalyticalFlux<dimDomain>().rhoeps(u) << std::endl;
        return (EulerAnalyticalFlux<dimDomain>().rhoeps(u) > 1e-8);
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
      // take pressure as shock detection values
      const RangeFieldType pl = pressure( uLeft );
      const RangeFieldType pr = pressure( uRight );
      jump  = (pl-pr)/(0.5*(pl+pr));
    }

    // calculate jump between left and right value
    inline void adaptationIndicator (const IntersectionType& it,
                                     const double time,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const RangeType& uRight,
                                     RangeType& indicator) const
    {
      // take density as shock detection values
      indicator = (uLeft[0] - uRight[0])/(0.5 * (uLeft[0]+uRight[0]));

      const DomainType unitNormal = it.unitOuterNormal(x);
      RangeType ul(uLeft),ur(uRight);
      fieldRotator_.rotateForth( ul , unitNormal );
      fieldRotator_.rotateForth( ur , unitNormal );
      for (int i=1; i<dimDomain+1; ++i)
      {
        indicator += ((ul[i+1]/ul[0] - ur[i+1]/ur[0])/problem_.V());
      }
    }

   protected:
    const ThermodynamicsType thermodynamics_;
    const ProblemType& problem_;
    FieldRotator< DomainType, RangeType > fieldRotator_;
  };

}

}

#endif
