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
#include <dune/fem/schemes/diffusionmodel.hh>

// fem-dg includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>
#include <dune/fem-dg/operator/limiter/limiterutility.hh>
#include <dune/fem-dg/examples/euler/problems.hh>

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
        std::abort();
      }

      double gamma () const { return 1.4; }
    };

  } // end namespace detail

  template< class GridImp, class ProblemImp, class LimiterFunction = MinModLimiter< typename GridImp::ctype> >
  class ModelWrapperTraits
    : public DefaultModelTraits< GridImp, ProblemImp >
  {
    typedef DefaultModelTraits< GridImp, ProblemImp >           BaseType;
  public:
    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                    GradientType;
    static const int modelParameterSize = 0;

    typedef LimiterFunction LimiterFunctionType;
    // possible options are
    //typedef MinModLimiter< typename BaseType::DomainFieldType >     LimiterFunctionType ;
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
            class AdditionalImp,
            class LimiterFunction,
            class Problem = detail::EmptyProblem< AdvectionModelImp > >
  class ModelWrapper :
    public DefaultModel< ModelWrapperTraits< GridImp, Problem, LimiterFunction > >
  {
  public:
    typedef GridImp                                      GridType;
    typedef ModelWrapperTraits< GridType, Problem, LimiterFunction > Traits;
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

    static const int limitedDimRange = AdditionalType :: limitedDimRange ;
    typedef Dune::FieldVector< int, limitedDimRange > LimitedRangeType;

    // for Euler equations diffusion is disabled
    static const bool hasAdvection = AdditionalType::hasAdvection;
    static const bool hasDiffusion = AdditionalType::hasDiffusion;

    using BaseType :: time;
    using BaseType :: hasMass;
    using BaseType :: velocity ;

    ModelWrapper( const AdvectionModelType& advModel, const DiffusionModelType& diffModel )
      : advection_( advModel ),
        diffusion_( diffModel ),
        problem_(),
        limitedRange_()
    {
      // by default this should be the identity
      for( int i=0; i<limitedDimRange; ++i )
        limitedRange_[ i ] = i;

      // if method has been filled then modified will be set differently
      // AdditionalType :: limitedRange( limitedRange_ );
      advection_.limitedRange( limitedRange_ );
    }

    void setTime (const double t)
    {
      BaseType::setTime(t);
      // update model times (only if time method is available on these models)
      //! TODO problem without virtualization advection_.setTime(t);
      //! TODO problem without virtualization diffusion_.setTime(t);
      ::detail::CallSetTime< AdvectionModelType,
                           ::detail::CheckTimeMethod< AdvectionModelType >::value >
        ::setTime( const_cast< AdvectionModelType& > (advection_), t );
      ::detail::CallSetTime< DiffusionModelType,
                           ::detail::CheckTimeMethod< DiffusionModelType >::value >
        ::setTime( const_cast< DiffusionModelType& > (diffusion_), t );
    }

    double gamma () const { return problem_.gamma(); }

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

    const LimitedRangeType& limitedRange() const { return limitedRange_; }

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

    template <class LocalEvaluation>
    inline void advection( const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& du,
                           JacobianRangeType& result ) const
    {
      assert( hasAdvection );
      advection_.diffusiveFlux( local.quadraturePoint(), u, du, result );
    }

    template <class LocalEvaluation>
    inline void eigenValues(const LocalEvaluation& local,
                            const RangeType& u,
                            RangeType& maxValue) const
    {
      std::cerr <<"eigenValues for problems/euler not implemented\n";
      std::abort();
    }

    template <class LocalEvaluation, class T>
    inline double diffusionTimeStep( const LocalEvaluation& local,
                                     const T& circumEstimate,
                                     const RangeType& u ) const
    {
      // return AdditionalType :: diffusionTimeStep( local.entity(), local.quadraturePoint(), circumEstimate, u );
      return diffusion_.diffusionTimeStep( local.entity(), local.quadraturePoint(), circumEstimate, u );
    }

    // is not used
    template <class LocalEvaluation>
    inline  void jacobian( const LocalEvaluation& local,
                           const RangeType& u,
                           const FluxRangeType& du,
                           RangeType& A ) const
    {
      assert( hasAdvection );
      assert( 0 ); // if this is used we have to check if this is correct
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
      // note: the local coordinate used here is not important since we
      // only check if the intersection as a whole is part of the dirichlet
      // boundary - so the center is used
      // bool bndVal = AdditionalType::boundaryValue(id, time(), local.entity(), local.intersection().geometryInInside().center(), u, u);
      bool bndVal = advection_.boundaryValue(id, time(), local.entity(), local.intersection().geometryInInside().center(), u, u);

      return bndVal;
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
      // AdditionalType::boundaryValue(id, time(), local.entity(), local.quadraturePoint(), uLeft, uRight);
      advection_.boundaryValue(id, time(), local.entity(), local.quadraturePoint(), uLeft, uRight);
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
      DomainType normal = local.intersection().integrationOuterNormal( local.localPosition() );
      double len = normal.two_norm();
      normal /= len;
      int id = getBoundaryId( local );
#ifndef NDEBUG
      const bool isFluxBnd =
#endif
      // AdditionalType::boundaryFlux(id, time(), local.entity(), local.quadraturePoint(), normal, uLeft, gLeft);
      advection_.boundaryFlux(id, time(), local.entity(), local.quadraturePoint(), normal, uLeft, gLeft);
      gLeft *= len;
      assert( isFluxBnd );
      return 0; // TODO: do something better here
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
      DomainType normal = local.intersection().integrationOuterNormal( local.localPosition() );
      double len = normal.two_norm();
      normal /= len;
      int id = getBoundaryId( local );
#ifndef NDEBUG
      const bool isFluxBnd =
#endif
      // AdditionalType::diffusionBoundaryFlux(id, time(), local.entity(), local.quadraturePoint(), normal, uLeft, jacLeft, gLeft);
      diffusion_.diffusionBoundaryFlux(id, time(), local.entity(), local.quadraturePoint(), normal, uLeft, jacLeft, gLeft);
      gLeft *= len;
      assert( isFluxBnd );
      return 0; // QUESTION: do something better here? Yes, return time step restriction if possible
    }

    template <class LocalEvaluation>
    inline void maxSpeed( const LocalEvaluation& local,
                          const DomainType& unitNormal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      // TODO: add a max speed for the diffusion time step control
      // this needs to be added in diffusionTimeStep
      assert( hasAdvection );
      // advspeed = AdditionalType::maxSpeed( time(), local.entity(), local.quadraturePoint(), unitNormal, u );
      advection_.maxSpeed( time(), local.entity(), local.quadraturePoint(), unitNormal, u );
      totalspeed = advspeed;
    }

    inline const ProblemType& problem() const
    {
      return problem_;
    }

    // velocity method for Upwind flux
    template< class LocalEvaluation >
    inline DomainType velocity (const LocalEvaluation& local,
                                RangeType& u) const
    {
      // return AdditionalType :: velocity( time(), local.entity(), local.quadraturePoint(), u );
      return advection_.velocity( time(), local.entity(), local.quadraturePoint(), u );
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
      // velocity = AdditionalType :: velocity( time(), en, x, u );
      velocity = advection_.velocity( time(), en, x, u );
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
      // return AdditionalType :: physical( entity, x, u ) > 0;
      return advection_.physical( entity, x, u ) > 0;
    }

    // adjust average value if necessary
    // (e.g. transform from conservative to primitive variables )
    template< class Entity >
    void adjustAverageValue( const Entity& entity,
                             const DomainType& xLocal,
                             RangeType& u ) const
    {
      // nothing to be done here for this test case
      // AdditionalType :: adjustAverageValue( entity, xLocal, u );
      advection_.adjustAverageValue( entity, xLocal, u );
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void jump(const Intersection& it,
                     const FaceDomainType& x,
                     const RangeType& uLeft,
                     const RangeType& uRight,
                     RangeType& jump) const
    {
      // jump = AdditionalType :: jump( it, x, uLeft, uRight );
      jump = advection_.jump( it, x, uLeft, uRight );
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void adaptationIndicator (const Intersection& it,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const RangeType& uRight,
                                     RangeType& indicator) const
    {
      // indicator = AdditionalType :: jump( it, x, uLeft, uRight );
      indicator = advection_.jump( it, x, uLeft, uRight );
    }

    // misc for eoc calculation
    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
    }

  protected:
    const AdvectionModelType& advection_;
    const DiffusionModelType& diffusion_;

    ProblemType problem_;
    LimitedRangeType limitedRange_;
  };

}

}

#endif
