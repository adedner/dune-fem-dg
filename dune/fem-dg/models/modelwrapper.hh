#ifndef DUNE_MODELWRAPPER_HH
#define DUNE_MODELWRAPPER_HH

// system includes
#include <config.h>
#include <cmath>
#include <type_traits>

// DUNE includes
#include <dune/common/version.hh>
#include <dune/common/fmatrixev.hh>

#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/boundaryidprovider.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/fem/schemes/conservationlawmodel.hh>

// fem-dg includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>
#include <dune/fem-dg/operator/limiter/limiterutility.hh>
//#include <dune/fem-dg/examples/euler/problems.hh>

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

      virtual double endTime () const
      {
        return ParameterType::getValue< double >( "femdg.stepper.endtime" );
      }

      void bg ( const DomainType&, RangeType& ) const
      {
        DUNE_THROW(InvalidStateException,"EmptyProblem::endTime: this method should not be called");
      }

      //! methods for gradient based indicator
      double indicator1( const DomainType& xgl, const RangeType& u ) const
      {
        //DUNE_THROW(InvalidStateException,"EmptyProblem::endTime: this method should not be called");
        // use density as indicator
        return u[ 0 ];
      }

      virtual int boundaryId ( const int id ) const
      {
        DUNE_THROW(InvalidStateException,"EmptyProblem::endTime: this method should not be called");
        return -1;
      }

      void evaluate(const DomainType& x, const double time, RangeType& res) const
      {
        DUNE_THROW(InvalidStateException,"EmptyProblem::endTime: this method should not be called");
      }

    private:
      template <int, bool>
      struct Gamma
      {
        static double value()
        {
          DUNE_THROW(NotImplemented,"Python Model does not implement a gamma() method, please add that!");
          return 0.0;
        }
      };

      template <int d>
      struct Gamma<d, true>
      {
        static double value() { return ModelImp::gamma; }
      };

    public:
      double gamma() const
      {
        return Gamma< 0, ModelImp::hasGamma >::value();
      }
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
    //typedef NoLimiter< typename BaseType::DomainFieldType >  LimiterFunctionType ;

    // check whether LimiterFunction is NoLimiter, in that case scalingLimiter
    // might have been chosen
    static const bool scalingLimiter =
      std::is_same< LimiterFunctionType, NoLimiter< typename LimiterFunctionType::FieldType > > :: value ;
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

    typedef typename AdvectionModelImp::GridPartType     GridPartType;
    typedef typename Traits::ProblemType                 ProblemType;
    typedef AdvectionModelImp                            AdvectionModelType ;
    typedef DiffusionModelImp                            DiffusionModelType ;

    typedef AdditionalImp                                AdditionalType;

    enum { dimDomain = Traits::dimDomain };
    enum { dimRange  = Traits::dimRange  };

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

    static const bool scalingLimiter = Traits::scalingLimiter;

    // for Euler equations diffusion is disabled
    static const bool hasAdvection = AdditionalType::hasAdvection;
    static const bool hasDiffusion = AdditionalType::hasDiffusion;

    using BaseType :: time;
    using BaseType :: hasMass;
    using BaseType :: velocity ;

    // default constructor called by DGOperator
    ModelWrapper( const AdvectionModelType& advModel, const DiffusionModelType& diffModel, const ProblemType& problem )
      : advection_( advModel ),
        diffusion_( diffModel ),
        problem_( problem ),
        limitedRange_()
    {
      // by default this should be the identity
      for( int i=0; i<limitedDimRange; ++i )
        limitedRange_[ i ] = i;

      // if method has been filled then modified will be set differently
      advection().limitedRange( limitedRange_ );
    }

#ifdef EULER_WRAPPER_TEST
    ModelWrapper( const ProblemType& problem )
      : ModelWrapper( *(new AdvectionModelType()), *(new DiffusionModelType()), problem )
    {}
#endif

    void setTime (const double t)
    {
      BaseType::setTime(t);
      // update model times (only if time method is available on these models)
      //! TODO problem without virtualization advection_.setTime(t);
      //! TODO problem without virtualization diffusion().setTime(t);
      ::detail::CallSetTime< AdvectionModelType,
                             ::detail::CheckTimeMethod< AdvectionModelType >::value >
        ::setTime( const_cast< AdvectionModelType& > (advection()), t );
      ::detail::CallSetTime< DiffusionModelType,
                             ::detail::CheckTimeMethod< DiffusionModelType >::value >
        ::setTime( const_cast< DiffusionModelType& > (diffusion()), t );
    }

    double gamma () const { return problem_.gamma(); }

    template <class Entity>
    void setEntity( const Entity& entity ) const
    {
      if constexpr ( hasAdvection )
        advection().init( entity );
      if constexpr ( hasDiffusion )
        diffusion().init( entity );
    }

    inline bool hasStiffSource() const { return AdditionalType::hasStiffSource; }
    inline bool hasNonStiffSource() const { return AdditionalType::hasNonStiffSource; }
    inline bool hasFlux() const { return AdditionalType::hasFlux; }

    const LimitedRangeType& limitedRange() const { return limitedRange_; }

    void obtainBounds( RangeType& globalMin, RangeType& globalMax) const
    {
      assert( hasAdvection );
      advection().obtainBounds(globalMin, globalMax);
    }

    bool isConstant( const RangeType& min, const RangeType& max ) const
    {
      return (min - max).infinity_norm() < 1e-10;
    }

    template <class State> // double[dim+2]
    void pressureTemperature( const State& u, double& p, double& T ) const
    {
      if constexpr ( hasAdvection )
      {
        const auto pT = advection().pressureTemperature( u );
        p = pT[0];
        T = pT[1];
        return ;
      }

      if constexpr ( hasDiffusion )
      {
        const auto pT = diffusion().pressureTemperature( u );
        p = pT[0];
        T = pT[1];
        return ;
      }

      p = 0.0;
      T = 0.0;
    }

    template <class LocalEvaluation>
    inline double stiffSource( const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& du,
                               RangeType & s) const
    {
      assert( hasDiffusion );
      diffusion().source( local.quadraturePoint(), u, du, s );
      return 0;
    }


    template< class LocalEvaluation >
    inline double nonStiffSource( const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& du,
                                  RangeType& s) const
    {
      assert( hasAdvection );
      advection().source( local.quadraturePoint(), u, du, s );
      return 0;
    }

    template <class LocalEvaluation>
    inline void advection( const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianRangeType& du,
                           JacobianRangeType& result ) const
    {
      assert( hasAdvection );
      advection().flux( local.quadraturePoint(), u, du, result );
    }

    template <class LocalEvaluation>
    inline void eigenValues(const LocalEvaluation& local,
                            const RangeType& u,
                            RangeType& maxValue) const
    {
      if constexpr ( hasDiffusion )
      {
        maxValue = diffusion().diffusionTimeStep( local.entity(), local.quadraturePoint(), u );
      }
    }

    template <class LocalEvaluation>
    inline double diffusionTimeStep( const LocalEvaluation& local,
                                     const RangeType& u ) const
    {
      return diffusion().diffusionTimeStep( local.entity(), local.quadraturePoint(), u );
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
      advection().linDiffusiveFlux( u, du, local.quadraturePoint(), u, du, A);
    }

    template <class LocalEvaluation>
    int getBoundaryId( const LocalEvaluation& local ) const
    {
      return BoundaryIdProviderType::boundaryId( local.intersection() );
    }

    template <class LocalEvaluation>
    inline bool hasBoundaryValue( const LocalEvaluation& local ) const
    {
      RangeType u; // fake return variable
      const int id = getBoundaryId( local );
      // the following fails since is is called with
      return advection().hasBoundaryValue(id, time(), local.entity(), local.position(), u, u);
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
      advection().boundaryValue(id, time(), local.entity(), local.quadraturePoint(),
                                    local.intersection().unitOuterNormal( local.localPosition() ),
                                    uLeft, uRight);
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
      advection().boundaryFlux(id, time(), local.entity(), local.quadraturePoint(), normal, uLeft, gLeft);
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
      // TODO: This should be done outside when calling setEntity,
      // however, this seems not to work and will trigger the assertion
      diffusion().init( local.entity() );
      //assert( &(diffusion().entity()) == &(local.entity()) );
      diffusion().flux( local.quadraturePoint(), u, du, diff);
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
      diffusion().diffusionBoundaryFlux(id, time(), local.entity(), local.quadraturePoint(), normal, uLeft, jacLeft, gLeft);
      gLeft *= len;
      assert( isFluxBnd );
      return 0; // QUESTION: do something better here? Yes, return time step restriction if possible
    }

    template <class LocalEvaluation>
    inline void maxWaveSpeed( const LocalEvaluation& local,
                          const DomainType& unitNormal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      // TODO: add a max speed for the diffusion time step control
      // this needs to be added in diffusionTimeStep
      assert( hasAdvection );
      advspeed = advection().maxWaveSpeed( time(), local.entity(), local.quadraturePoint(), unitNormal, u );
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
      return advection().velocity( time(), local.entity(), local.quadraturePoint(), u );
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
      velocity = advection().velocity( time(), en, x, u );
    }

    // we have physical check for this model
    constexpr bool hasPhysical() const
    {
      return AdvectionModelImp::hasPhysical;
    }

    // calculate jump between left and right value
    template< class Entity >
    inline bool physical(const Entity& entity,
                         const DomainType& x,
                         const RangeType& u) const
    {
      if constexpr ( AdvectionModelImp::hasPhysical )
      {
        return advection().physical( entity, x, u ) > 0;
      }
      else
        return true;
    }

    // adjust average value if necessary
    // (e.g. transform from conservative to primitive variables )
    template< class Entity >
    void adjustAverageValue( const Entity& entity,
                             const DomainType& xLocal,
                             RangeType& u ) const
    {
      // nothing to be done here for this test case
      advection().adjustAverageValue( entity, xLocal, u );
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void jump(const Intersection& it,
                     const FaceDomainType& x,
                     const RangeType& uLeft,
                     const RangeType& uRight,
                     RangeType& jump) const
    {
      jump = advection().jump( it, x, uLeft, uRight );
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void adaptationIndicator (const Intersection& it,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const RangeType& uRight,
                                     RangeType& indicator) const
    {
      indicator = advection().jump( it, x, uLeft, uRight );
    }

    // misc for eoc calculation
    template< class DiscreteFunction >
    void eocErrors( const DiscreteFunction& df ) const
    {
    }

  protected:
    const AdvectionModelType& advection () const { return advection_; }
    const DiffusionModelType& diffusion () const { return diffusion_; }

    // store a copy of the models here for thread parallel runs
    // where we need class variables to be thread safe
    AdvectionModelType advection_;
    DiffusionModelType diffusion_;

    const ProblemType& problem_;
    LimitedRangeType limitedRange_;
  };

  template< class GridImp,
            class AdvectionModelImp,
            class AdditionalImp,
            class LimiterFunction >
  class AdvectionModelWrapper
    : public ModelWrapper< GridImp, AdvectionModelImp, AdvectionModelImp, AdditionalImp, LimiterFunction >
  {
    typedef ModelWrapper< GridImp, AdvectionModelImp, AdvectionModelImp, AdditionalImp, LimiterFunction >  BaseType;

  public:
    typedef typename BaseType :: ProblemType          ProblemType;
    typedef typename BaseType :: AdvectionModelType   AdvectionModelType;

    // default constructor called by DGOperator
    AdvectionModelWrapper( const AdvectionModelType& advModel )
      : BaseType( advModel, advModel, *(new ProblemType() )),
        problemPtr_( &this->problem_ )
    {}

    // copy constructor
    AdvectionModelWrapper( const AdvectionModelWrapper& other )
      : AdvectionModelWrapper( other.advection() )
    {}

  protected:
    std::unique_ptr< const ProblemType > problemPtr_;
  };

} // end namespace Fem
} // end namespace Dune

#endif
