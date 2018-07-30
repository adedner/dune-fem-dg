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

//#include <dune/fem-dg/operator/fluxes/analyticaleulerflux.hh>
#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem-dg/misc/error/l1eocerror.hh>

#if HAVE_DUNE_FEMPY
#include <dune/fem/schemes/diffusionmodel.hh>
#endif

namespace Dune
{
namespace Fem
{
  namespace detail {

    template <class ModelImp>
    class EmptyProblem :
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

      EmptyProblem() {}

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
        res = 0;
        if( x[0] < 0.5 )
        {
          // default is sod's rp
          res[ 0 ] = 1.0;
          res[ dimDomain+1 ] = 1.0;
        }
        else
        {
          res[ 0 ] = 0.125;
          res[ dimDomain + 1 ] = 0.1;
        }
      }
    };

  } // end namespace detail

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
  template< class GridImp,
            class AdvectionModelImp,
            class DiffusionModelImp,
            class AdditionalImp >
  class ModelWrapper :
    public DefaultModel< ModelWrapperTraits< GridImp, detail::EmptyProblem< AdvectionModelImp > > >
  {
  public:
    typedef GridImp                                      GridType;
    typedef ModelWrapperTraits< GridType, detail::EmptyProblem< AdvectionModelImp > > Traits;
    typedef DefaultModel< Traits >                       BaseType;
    typedef typename Traits::ProblemType                 ProblemType;
    typedef AdvectionModelImp                            AdvectionModelType ;
    typedef DiffusionModelImp                            DiffusionModelType ;

    typedef AdditionalImp                                AdditionalType;

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
    static const bool hasAdvection = AdditionalType::hasAdvection;
    static const bool hasDiffusion = AdditionalType::hasDiffusion;

    using BaseType :: time;
    using BaseType :: hasMass;

    void setTime (double t)
    {
      BaseType::setTime(t);
    }

    ModelWrapper( const AdvectionModelType& advModel, const DiffusionModelType& diffModel )
      : advection_( advModel ),
        diffusion_( diffModel ),
        //eulerFlux_(),
        problem_()
    {
      modified_[ 0 ] = 0;
      modified_[ 1 ] = dimRange-1;
    }

    template <class Entity>
    void setEntity( const Entity& entity ) const
    {
      if( hasAdvection )
        advection_.init( entity );
      if( hasDiffusion )
        diffusion_.init( entity );
    }

    inline bool hasStiffSource() const { return AdditionalType::hasStiffSource; }
    inline bool hasNonStiffSource() const { return AdditionalType::hasNonStiffSource; }
    inline bool hasFlux() const { return AdditionalType::hasFlux; }

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
      assert( hasDiffusion );
      diffusion_.source( local.quadraturePoint(), u, du, s );
      return 0;
    }


    template< class LocalEvaluation >
    inline double nonStiffSource( const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& du,
                                  RangeType& s) const
    {
      assert( hasAdvection );
      advection_.source( local.quadraturePoint(), u, du, s );
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
      assert( hasAdvection );
      advection_.diffusiveFlux( local.quadraturePoint(), u, du, f);

      //JacobianRangeType fTest;
      //eulerFlux_.analyticalFlux( gamma() , u , fTest );

      /*
      if( ( f - fTest ).infinity_norm() > 1e-10 )
      {
        std::cout << f << " diffFlux" << std::endl;
        std::cout << fTest << " eulerFlux" << std::endl;
      }
      */
    }

    template <class LocalEvaluation>
    inline void eigenValues(const LocalEvaluation& local,
                            const RangeType& u,
                            RangeType& maxValue) const
    {
      std::cerr <<"eigenValues for problems/euler not implemented\n";
      std::abort();
    }

    template <class LocalEvaluation>
    inline double diffusionTimeStep( const LocalEvaluation& local,
                                     const double circumEstimate,
                                     const RangeType& u ) const
    {
      // TODO: implement using diffusion model - something on Additional?
      return 0;
    }

    // is not used
    template <class LocalEvaluation>
    inline  void jacobian( const LocalEvaluation& local,
                           const RangeType& u,
                           const FluxRangeType& du,
                           RangeType& A ) const
    {
      assert( hasAdvection );
      // TODO: u != ubar and du != dubar
      advection_.linDiffusiveFlux( u, du, local.quadraturePoint(), u, du, A);
    }

    template <class LocalEvaluation>
    int getBoundaryId( const LocalEvaluation& local ) const
    {
      return BoundaryIdProviderType::boundaryId( local.intersection() );
    }

    template <class LocalEvaluation>
    inline bool hasBoundaryValue( const LocalEvaluation& local ) const
    {
      RangeType u;
      int id = getBoundaryId( local );
      return AdditionalType::boundaryValue(id, time(), local.entity(), local.quadraturePoint(), u, u);
    }

    // return uRight for insertion into the numerical flux
    template <class LocalEvaluation>
    inline void boundaryValue( const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight ) const
    {
      int id = getBoundaryId( local );
#ifndef NDEBUG
      const bool isDirichlet =
#endif
      AdditionalType::boundaryValue(id, time(), local.entity(), local.quadraturePoint(), uLeft, uRight);
      assert( isDirichlet );
    }

    // boundary condition here is slip boundary cond. <u,n>=0
    // gLeft= p*[0 n(global(x)) 0]
    template <class LocalEvaluation>
    inline double boundaryFlux( const LocalEvaluation& local,
                                const RangeType& uLeft,
                                const JacobianRangeType&,
                                RangeType& gLeft ) const
    {
      const DomainType normal = local.intersection().integrationOuterNormal( local.localPosition() );
      int id = getBoundaryId( local );
#ifndef NDEBUG
      const bool isFluxBnd =
#endif
      AdditionalType::boundaryFlux(id, time(), local.entity(), local.quadraturePoint(), normal, uLeft, gLeft);
#if 0
      std::cout << "BoundaryFlux with id=" << id << " with uLeft=" << uLeft
                << " results in flux=" << gLeft
                << " --- " << isFluxBnd
                << std::endl;
#endif
      assert( isFluxBnd );
      return 0; // QUESTION: do something better here? Yes, return time step restriction if possible
    }

    template <class LocalEvaluation>
    void diffusion( const LocalEvaluation& local,
                    const RangeType& u,
                    const JacobianRangeType& du,
                    JacobianRangeType& diff ) const
    {
      assert( hasDiffusion );
      diffusion_.diffusiveFlux( local.quadraturePoint(), u, du, diff);
    }


    /** \brief boundary flux for the diffusion part
     */
    template <class LocalEvaluation>
    inline double diffusionBoundaryFlux( const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const JacobianRangeType& jacLeft,
                                         RangeType& gLeft ) const
    {
      // TODO: need to add something to 'Addtional'?
      assert( hasDiffusion );
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
      // TODO: add a max speed for the diffusion time step control
      assert( hasAdvection );
      double len = normal.two_norm();
      DomainType unitNormal(normal);
      unitNormal /= len;
      advspeed = AdditionalType::maxSpeed( time(), local.entity(), local.quadraturePoint(), unitNormal, u );
      advspeed *= len;
      totalspeed = advspeed;
    }

    inline const ProblemType& problem() const
    {
      return problem_;
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
      velocity = AdditionalType :: velocity( time(), en, x, u );
    }

    // we have physical check for this model
    bool hasPhysical() const
    {
      return true;
    }

    // calculate jump between left and right value
    template< class Entity >
    inline bool physical(const Entity& entity,
                         const DomainType& x,
                         const RangeType& u) const
    {
      return AdditionalType :: physical( entity, x, u ) > 0;
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
      jump = AdditionalType :: jump( it, x, uLeft, uRight );
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void adaptationIndicator (const Intersection& it,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const RangeType& uRight,
                                     RangeType& indicator) const
    {
      indicator = AdditionalType :: jump( it, x, uLeft, uRight );
    }

    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
      EOCErrorList::setErrors<L2EOCError>( *this, df );
      EOCErrorList::setErrors<L1EOCError>( *this, df );
    }

  protected:
    const AdvectionModelType& advection_;
    const DiffusionModelType& diffusion_;

    //EulerAnalyticalFlux<dimDomain, RangeFieldType > eulerFlux_;
    ProblemType problem_;
    ModifiedRangeType modified_;
  };

}

}

#endif
