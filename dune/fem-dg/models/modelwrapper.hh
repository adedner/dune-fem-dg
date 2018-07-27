#ifndef DUNE_MODELWRAPPER_HH
#define DUNE_MODELWRAPPER_HH

// system includes
#include <config.h>
#include <cmath>
#include <type_traits>

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/boundaryidprovider.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/l1eocerror.hh>

namespace Dune
{
namespace Fem
{
  namespace detail {

    template <class ModelImp>
    class ProblemWrapper :
      public Dune::Fem::EvolutionProblemInterface< typename ModelImp :: RFunctionSpaceType, false >
    {
      typedef Dune::Fem::EvolutionProblemInterface< typename ModelImp :: RFunctionSpaceType, false > BaseType;

    public:
      using BaseType :: evaluate ;
      using BaseType :: fixedTimeFunction;

      typedef Dune::Fem::Parameter ParameterType;

      typedef typename BaseType :: FunctionSpaceType  FunctionSpaceType;

      enum { dimDomain = FunctionSpaceType::dimDomain };
      typedef typename FunctionSpaceType :: RangeType  RangeType ;
      typedef typename FunctionSpaceType :: DomainType DomainType ;

      typedef typename FunctionSpaceType :: RangeFieldType  RangeFieldType ;
      typedef RangeFieldType FieldType ;

      ProblemWrapper() {}

      void init () {}

      virtual double endTime () const { return 0.1; }

      void bg ( const DomainType&, RangeType& ) const {}

      //! methods for gradient based indicator
      bool twoIndicators() const { return false ; }

      //! methods for gradient based indicator
      double indicator1( const DomainType& xgl, const RangeType& u ) const
      {
        // use density as indicator
        return u[ 0 ];
      }

      virtual int boundaryId ( const int id ) const
      {
        return 1;
      }

      void evaluate(const DomainType& x, const double time, RangeType& res) const
      {
        // TODO: extract from ModelImp
        res = 0 ;
      }
    };

  } // end namespace detail


  template <class ModelImp>
  struct ModelImplementationWrapper
    : public ModelImp,
      public detail::ProblemWrapper< ModelImp >
  {
    typedef typename ModelImp :: RFunctionSpaceType   FunctionSpaceType ;
    ModelImplementationWrapper() : ModelImp() {}

    using ModelImp :: init;
  };

  template< class GridImp, class ProblemImp >
  class ModelWrapperTraits
    : public DefaultModelTraits< GridImp, ProblemImp >
  {
    typedef DefaultModelTraits< GridImp, ProblemImp >           BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                    GradientType;
    static const int modelParameterSize = 0;

    typedef MinModLimiter< typename BaseType::DomainFieldType >     LimiterFunctionType ;
    //typedef SuperBeeLimiter< typename BaseType::DomainFieldType > LimiterFunctionType ;
    //typedef VanLeerLimiter< typename BaseType::DomainFieldType >  LimiterFunctionType ;
  };


  /**
   * \brief Euler equations for dry atmosphere
   *
   * \ingroup AnalyticalModels
   */
  template< class GridImp, class ProblemImp >
  class ModelWrapper :
    public DefaultModel< ModelWrapperTraits< GridImp, ProblemImp > >
  {
  public:
    typedef GridImp                                      GridType;
    typedef ModelWrapperTraits< GridType, ProblemImp >     Traits;
    typedef DefaultModel< Traits >                       BaseType;
    typedef typename Traits::ProblemType                 ProblemType;

    enum { dimDomain = Traits::dimDomain };
    enum { dimRange = Traits::dimRange };

    typedef typename Traits::FaceDomainType              FaceDomainType;

    typedef typename Traits::RangeType                   RangeType;
    typedef typename Traits::RangeFieldType              RangeFieldType ;
    typedef typename Traits::DomainType                  DomainType;
    typedef typename Traits::FluxRangeType               FluxRangeType;
    typedef typename Traits::GradientType                GradientType;
    typedef typename Traits::JacobianRangeType           JacobianRangeType;
    typedef typename Traits::DiffusionRangeType          DiffusionRangeType;

    typedef Dune::Fem::BoundaryIdProvider < GridType >   BoundaryIdProviderType;

    typedef Dune::FieldVector< int, 2 > ModifiedRangeType;

    // for Euler equations diffusion is disabled
    static const bool hasAdvection = true;
    static const bool hasDiffusion = false;

    using BaseType::time;
    using BaseType::time_;

   public:
#if 0
    ModelWrapper()
      : impl()
    {
    }
#endif

    ModelWrapper( const ProblemType& problem )
      : impl_( problem )
    {
      modified_[ 0 ] = 0;
      modified_[ 1 ] = dimRange-1;
    }

    template <class Entity>
    void setEntity( const Entity& entity ) const
    {
      impl_.init( entity );
    }

    double gamma () const { return 1.4; }

    inline bool hasStiffSource() const { return false; }
    inline bool hasNonStiffSource() const { return false; }
    inline bool hasFlux() const { return true ; }

    const ModifiedRangeType& modifiedRange() const { return modified_; }

    void obtainBounds( RangeType& globalMin, RangeType& globalMax) const
    {
      globalMin = 0;
      globalMax = std::numeric_limits<double>::max();
    }

    bool isConstant( const RangeType& min, const RangeType& max ) const
    {
      return (min - max).infinity_norm() < 1e-10;
    }


    template <class LocalEvaluation>
    inline double stiffSource( const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& du,
                               RangeType & s) const
    {
      impl_.source( local.quadraturePoint(), u, du, s );
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

    inline void conservativeToPrimitive( const double time,
                                         const DomainType& xgl,
                                         const RangeType& cons,
                                         RangeType& prim,
                                         const bool ) const
    {
      // impl_.evaluate( xgl, time, prim );
    }

    template <class LocalEvaluation>
    inline void advection( const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& du,
                           JacobianRangeType& f ) const
    {
      impl_.diffusiveFlux( local.quadraturePoint(), u, du, f);
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
      // TODO: u != ubar and du != dubar
      impl_.linDiffusiveFlux( u, du, local.quadraturePoint(), u, du, A);
    }

    template <class LocalEvaluation>
    inline bool hasBoundaryValue( const LocalEvaluation& local ) const
    {
      Dune::FieldVector< int, dimRange > bndIds;
      return impl_.isDirichletIntersection( local.intersection(), bndIds );
    }

    // return iRight for insertion into the numerical flux
    template <class LocalEvaluation>
    inline void boundaryValue( const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight ) const
    {
      Dune::FieldVector< int, dimRange > bndIds;
#ifndef NDEBUG
      const bool isDirichlet =
#endif
      impl_.isDirichletIntersection( local.intersection(), bndIds );
      assert( isDirichlet );
      impl_.dirichlet( bndIds[ 0 ], local.quadraturePoint(), uRight );
    }

    // boundary condition here is slip boundary cond. <u,n>=0
    // gLeft= p*[0 n(global(x)) 0]
    template <class LocalEvaluation>
    inline double boundaryFlux( const LocalEvaluation& local,
                                const RangeType& uLeft,
                                const JacobianRangeType&,
                                RangeType& gLeft ) const
    {
      /*
      // Slip boundary condition
      const DomainType normal = local.intersection().integrationOuterNormal( local.localPosition() );

      const double p = eulerFlux_.pressure( gamma_ , uLeft );
      gLeft = 0;
      for ( int i = 0 ; i < dimDomain ; ++i )
        gLeft[i+1] = normal[i] * p;
        */

      gLeft = 0;
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
      // TODO
      std::abort();
      // advspeed = eulerFlux_.maxSpeed( gamma_ , normal , u );
      totalspeed = advspeed;
    }

    inline const ProblemType& problem() const
    {
      return impl_;
    }

    /////////////////////////////////////////////////////////////////
    // Limiter section
    ////////////////////////////////////////////////////////////////
    template< class Entity >
    inline void velocity (const Entity& en,
                          const DomainType& x,
                          const RangeType& u,
                          DomainType& velocity) const
    {
      for(int i=0; i<dimDomain; ++i)
      {
        // U = (rho, rho v_0,...,rho v_(d-1), e )
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
    template< class Entity >
    inline bool physical(const Entity& entity,
                         const DomainType& xGlobal,
                         const RangeType& u) const
    {
      if (u[0]<1e-8)
        return false;
      else
      {
        //std::cout << eulerFlux_.rhoeps(u) << std::endl;
        // return (eulerFlux_.rhoeps(u) > 1e-8);
        return true ;
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
      // take pressure as shock detection values
      //const RangeFieldType pl = pressure( uLeft );
      //const RangeFieldType pr = pressure( uRight );
      //jump  = (pl-pr)/(0.5*(pl+pr));
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void adaptationIndicator (const Intersection& it,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const RangeType& uRight,
                                     RangeType& indicator) const
    {
    }

    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
      EOCErrorList::setErrors<L2EOCError>( *this, df );
      EOCErrorList::setErrors<L1EOCError>( *this, df );
    }

  protected:
    ProblemType impl_;
    ModifiedRangeType modified_;
  };

}

}

#endif