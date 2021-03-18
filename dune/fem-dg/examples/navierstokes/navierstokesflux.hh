#ifndef NS_MODEL_SPEC_HH
#define NS_MODEL_SPEC_HH

#include <dune/common/fvector.hh>

#include <dune/fem-dg/operator/fluxes/analyticaleulerflux.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

namespace Dune
{
namespace Fem
{

  template< class Traits >
  class NSFlux
  {
    enum { dimDomain = Traits :: dimDomain };

    enum { e = dimDomain + 1 };
    enum { e_ = e };
    enum { dimRange = Traits :: dimRange };
    enum { dimGradRange = dimRange * dimDomain };

  public:
    typedef typename Traits::ProblemType                      ProblemType;
    typedef typename Traits::DomainType                       DomainType;
    typedef typename Traits::DomainFieldType                  DomainFieldType;
    typedef typename Traits::RangeFieldType                   RangeFieldType;
    typedef typename Traits::RangeType                        RangeType;
    typedef Dune::FieldVector< RangeFieldType, dimGradRange > GradientRangeType;
    typedef typename Traits::JacobianRangeType                JacobianRangeType;

    typedef Dune::FieldMatrix< RangeFieldType, dimGradRange, dimDomain >
                                                              DiffusionRangeType;
    typedef Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType >
                                                              ConvertedJacobianRangeType;
    typedef Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain >
                                                              VelocityGradientType;

    NSFlux( const ProblemType& problem )
      : eulerFlux_()
      , problem_( problem )
      , gamma_( problem.gamma() )
      , R_d_inv_( problem.thermodynamics().R_d_inv() )
      , Re_inv_( problem.Re_inv() )
      , c_v_inv_( problem.thermodynamics().c_vd_inv() )
      , c_v_( problem.thermodynamics().c_vd() )
    {}

    inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const
    {
      eulerFlux_.analyticalFlux( gamma_ , u , f );
    }

    inline void jacobian( const RangeType& u, DiffusionRangeType& du ) const
    {
      du = 0;

      assert( int(du.rows) == int(dimRange * dimDomain) );
      assert( int(du.cols) == int(dimDomain) );

      for (int r=0;r<dimRange;r++)
        for (int d=0;d<dimDomain;d++)
          du[dimDomain*r+d][d] = u[r];
    }

    inline double maxWaveSpeed( const DomainType& n, const RangeType& u ) const
    {
      return eulerFlux_.maxWaveSpeed( gamma_, n, u );
    }

    inline void diffusion( const RangeType& u,
                           const GradientRangeType& du,
                           JacobianRangeType& f ) const;

    template <class JacobianRangeImp>
    inline void diffusion( const RangeType& u,
                           const JacobianRangeImp& du,
                           JacobianRangeType& f ) const;

    inline double mu( const double T ) const { return problem_.mu(T); }
    inline double lambda( const double T ) const { return problem_.lambda(T); }
    inline double k( const double T ) const { return problem_.k(T); }

  protected:
    const EulerAnalyticalFlux<dimDomain, RangeFieldType > eulerFlux_;
    const ProblemType& problem_;
    const double gamma_;
    const double R_d_inv_;
    const double Re_inv_;
    const double c_v_inv_;
    const double c_v_;
  };


  /////////////////////////////////////////////
  // Implementation of flux functions
  /////////////////////////////////////////////
  template< class Traits >
  template< class JacobianRangeImp >
  void NSFlux< Traits >
  :: diffusion( const RangeType& u,
                const JacobianRangeImp& du,
                JacobianRangeType& diff ) const
  {
    //std::cout << du << " du " << std::endl;
    assert( u[0] > 1e-10 );
    const double rho_inv = 1. / u[0];

    // get velocity field
    DomainType v;
    for( int i = 0; i < dimDomain; ++i )
      v[i] = u[ i+1 ] * rho_inv;

    // | v |^2
    const double vTwoNorm = v.two_norm2();

    // get all partial derivatives of all velocities
    VelocityGradientType dVel;
    for( int j = 0; j < dimDomain; ++j ) // v components
    {
      for( int i = 0; i < dimDomain; ++i ) // space derivatives
      {
        // substract d_x rho * v from the derivative of the conservative variables
        dVel[ j ][ i ] = rho_inv * (du[ j+1 ][ i ] - v[ j ]*du[ 0 ][ i ]);
      }
    }

    DomainType vGradVel;
    for( int i = 0; i < dimDomain; ++i )
    {
      vGradVel[ i ] = 0.;
      for( int j = 0; j < dimDomain; ++j )
        vGradVel[ i ] += v[ j ] * dVel[ j ][ i ];
    }

    // get the absolute temperature
    const double T = problem_.temperature( u );

    // get all derivatives of the absolute temperature
    DomainType dTemp;
    for( int i = 0; i < dimDomain; ++i )
    {
      dTemp[ i ]  = du[ e_ ][ i ] - T * c_v_ * du[ 0 ][ i ] - 0.5*du[ 0 ][ i ]*vTwoNorm - u[ 0 ] * vGradVel[ i ];
      dTemp[ i ] *= rho_inv * c_v_inv_;
    }

    // get temperature dependend viscosity coefficients
    const double muLoc = mu( T );
    const double kLoc = k( T );

    // get divergence of the velocity field
    double divVel = 0;
    for( int i = 0; i < dimDomain; ++i )
      divVel += dVel[ i ][ i ];

    // apply lambda to divergence
    divVel *= lambda( T );

    for( int i = 0; i < dimDomain; ++i )
    {
      // assemble the diffusion part for the density equation
      // this equation is purely hyperbolic, so no diffusion
      diff[ 0 ][ i ] = 0;

      // assemble the diffusion part for the momentum equation
      // the stesstensor tau
      for( int j = 0; j < dimDomain; ++j )
      {
        diff[ j+1 ][ i ] = muLoc*(dVel[ j ][ i ] + dVel[ i ][ j ]);
      }

      // add   lambda * div v * 1I
      diff[ i+1 ][ i ] += divVel;
    }

    for( int i = 0; i < dimDomain; ++i )
    {
      // assemble the diffusion part for the energy (internal+kinetic) equation
      for( int j = 0; j < dimDomain; ++j )
        diff[e_][ i ] = v[ j ] * diff[ i+1 ][ j ];

      //diff[e_][i] += problem_.Pr_inv() * kLoc * dTemp[i];
      diff[e_][i] += kLoc * dTemp[ i ];
    }

    // scale diffusion part with Reynold's number
    diff *= Re_inv_;
  }

} // end namespace
} // end namespace
#endif
