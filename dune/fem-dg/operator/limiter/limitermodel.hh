#ifndef DUNE_FEMDG_LIMITERDEFAULTMODEL_HH
#define DUNE_FEMDG_LIMITERDEFAULTMODEL_HH

// system includes
#include <config.h>
#include <cmath>
#include <type_traits>
#include <limits>

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/boundaryidprovider.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/operator/limiter/limitpass.hh>

namespace Dune
{
namespace Fem
{
  template <class FS>
  struct PhonyProblem
  {
    typedef FS FunctionSpaceType;
  };

  template< class GridImp, class FunctionSpace >
  class LimiterDefaultModelTraits
    : public DefaultModelTraits< GridImp, PhonyProblem< FunctionSpace > >
  {
    typedef DefaultModelTraits< GridImp, PhonyProblem< FunctionSpace > >           BaseType;
  public:
    typedef PhonyProblem< FunctionSpace >  ProblemType;

    typedef Dune::FieldVector< typename BaseType::DomainFieldType, BaseType::dimGradRange >
                                                                    GradientType;
    static const int modelParameterSize = 0;

    typedef NoLimiter< typename BaseType::DomainFieldType >    LimiterFunctionType;
    //typedef MinModLimiter< typename BaseType::DomainFieldType >     LimiterFunctionType ;
    //typedef SuperBeeLimiter< typename BaseType::DomainFieldType > LimiterFunctionType ;
    //typedef VanLeerLimiter< typename BaseType::DomainFieldType >  LimiterFunctionType ;

  };


  /**
   * \brief Euler equations for dry atmosphere
   *
   * \ingroup AnalyticalModels
   */
  template< class GridImp, class FunctionSpace >
  class LimiterDefaultModel :
    public DefaultModel< LimiterDefaultModelTraits< GridImp, FunctionSpace > >
  {
  public:
    typedef LimiterDefaultModelTraits< GridImp, FunctionSpace >  Traits;
    typedef DefaultModel< Traits >                               BaseType;
    typedef typename Traits::ProblemType                         ProblemType;

    enum { dimDomain = Traits::dimDomain };
    enum { dimRange = Traits::dimRange };

    typedef typename Traits::GridType                    GridType;
    typedef typename Traits::FaceDomainType              FaceDomainType;

    typedef typename Traits::RangeType                   RangeType;
    typedef typename Traits::RangeFieldType              RangeFieldType ;
    typedef typename Traits::DomainType                  DomainType;
    typedef typename Traits::FluxRangeType               FluxRangeType;
    typedef typename Traits::GradientType                GradientType;
    typedef typename Traits::JacobianRangeType           JacobianRangeType;
    typedef typename Traits::DiffusionRangeType          DiffusionRangeType;

    typedef Dune::Fem::BoundaryIdProvider < GridType >   BoundaryIdProviderType;

    typedef Dune::FieldVector< int, dimRange-1 > LimitedRangeType;

    // for Euler equations diffusion is disabled
    static const bool hasAdvection = true;
    static const bool hasDiffusion = false;

    using BaseType::time;
    using BaseType::time_;

    // saturation component in RangeType vector, saturation in the last component
    static const int sat = dimRange - 1;

  protected:
    RangeType lower_;
    RangeType upper_;
    LimitedRangeType limitedRange_;

   public:
    LimiterDefaultModel( const double lower, const double upper,
                         LimitedRangeType mod = LimitedRangeType( sat ) )
      : lower_( lower ),
        upper_( std::numeric_limits< double >::max() ),
        limitedRange_( mod )
    {
      lower_[ sat ] = lower;
      upper_[ sat ] = upper;

      //for( int d=0; d<dimRange; ++d )
      //  modified_[ d ] = d;
    }

    inline bool hasStiffSource() const { return false; }
    inline bool hasNonStiffSource() const { return false; }
    inline bool hasFlux() const { return true ; }
    // we only need physical check
    inline bool calculateIndicator() const { return false ; }

    void obtainBounds( RangeType& globalMin, RangeType& globalMax) const
    {
      globalMin = lower_;
      globalMax = upper_;
    }

    // return set with components to be modified by limiter
    const LimitedRangeType& limitedRange() const { return limitedRange_; }

    /*
    // return true if solution is constant, i.e. min == max
    bool isConstant( const RangeType& minVal, const RangeType& maxVal ) const
    {
      return false ;//(std::abs( maxVal[ sat ] - minVal[ sat ] ) / (upper_[ sat ] - lower_[ sat ]))  < 1e-10;
    }
    */

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
      return 0.0;
    }

    inline void conservativeToPrimitive( const double time,
                                         const DomainType& xgl,
                                         const RangeType& cons,
                                         RangeType& prim,
                                         const bool ) const
    {
    }

    template <class LocalEvaluation>
    inline void advection( const LocalEvaluation& local,
                           const RangeType& u,
                           const FluxRangeType& jacu,
                           FluxRangeType& f ) const
    {
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
    }

    template <class LocalEvaluation>
    int getBoundaryId( const LocalEvaluation& local ) const
    {
      return BoundaryIdProviderType::boundaryId( local.intersection() );
    }

    template <class LocalEvaluation>
    inline bool hasBoundaryValue( const LocalEvaluation& local ) const
    {
      return false;
    }

    // return iRight for insertion into the numerical flux
    template <class LocalEvaluation>
    inline void boundaryValue( const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight ) const
    {
      uRight = uLeft;
    }

    // boundary condition here is slip boundary cond. <u,n>=0
    // gLeft= p*[0 n(global(x)) 0]
    template <class LocalEvaluation>
    inline double boundaryFlux( const LocalEvaluation& local,
                                const RangeType& uLeft,
                                const JacobianRangeType&,
                                RangeType& gLeft ) const
    {
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
    inline void maxWaveSpeed( const LocalEvaluation& local,
                          const DomainType& normal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      totalspeed = advspeed = 0;
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
      velocity = 0;
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
      //return (u[ 0 ] >= lower_[ 0 ]) && (u[ sat ] >= lower_[ sat ]) && (u[ sat ] <= upper_[ sat ]);
      return (u[ sat ]>= lower_[ sat ]) && (u[ sat ] <= upper_[ sat ]);
    }

    // adjust average value if necessary
    // (e.g. transform from conservative to primitive variables )
    template< class Entity >
    bool adjustAverageValue( const Entity& entity,
                             const DomainType& xLocal,
                             RangeType& u ) const
    {

      if( u[ sat ] < lower_[ sat])
      {
        u[ sat ] = lower_[ sat ] + 1e-14 ;
        return false ;
      }

      // nothing to be done here for this test case
      return true;
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void jump(const Intersection& it,
                     const FaceDomainType& x,
                     const RangeType& uLeft,
                     const RangeType& uRight,
                     RangeType& jump) const
    {
      assert( calculateIndicator() );
      jump = 0;
      //jump[ sat ] = (uLeft[ sat ] - uRight[ sat ])/(0.5*(uLeft[ sat ] + uRight[ sat ]));
    }

    // calculate jump between left and right value
    template< class Intersection >
    inline void adaptationIndicator (const Intersection& it,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const RangeType& uRight,
                                     RangeType& indicator) const
    {
      jump( it, x, uLeft, uRight, indicator );
    }
  };

} // end namespace Fem

} // end namespace Dune

#endif
