#ifndef NSWAVES_HH
#define NSWAVES_HH

#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>


// local includes
#include <dune/fem-dg/examples/navierstokes/thermodynamics.hh>
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


/***********************************************************
 *
 * 2d problem
 * PhD thesis Gregor Gassner, pg. 99
 *   Analytical solution for the Navier-Stokes equations
 *
 **********************************************************/


namespace Dune
{
namespace Fem
{
  /**
   * \brief Compressible Navier-Stokes problem.
   *
   * \ingroup NavierStokesProblems
   */

  template <class GridType>
  class NSWaves : public EvolutionProblemInterface<
                    Dune::Fem::FunctionSpace< typename GridType::ctype,
                                              typename GridType::ctype,
                                              GridType::dimension, GridType::dimension + 2 >,
                    true >,
                  public Thermodynamics< GridType::dimensionworld, typename GridType::ctype>
  {
    NSWaves( const NSWaves& );
  public:
    typedef Fem::FunctionSpace<typename GridType::ctype,
                               typename GridType::ctype,
                               GridType::dimensionworld,
                               GridType::dimensionworld + 2 > FunctionSpaceType ;

    typedef Fem::Parameter  ParameterType ;

    enum { dimension = GridType::dimensionworld };
    enum { energyId = dimension + 1 };
    typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
    typedef typename FunctionSpaceType :: DomainType        DomainType;
    typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
    typedef typename FunctionSpaceType :: RangeType         RangeType;
    typedef Thermodynamics< dimension, RangeFieldType >     ThermodynamicsType;

    NSWaves() : ThermodynamicsType(),
      myName_( "NSWaves" ),
      omegaGNS_( ParameterType::template getValue< double >( "omegaGNS" ) ),
      kGNS_( ParameterType::template getValue< double >( "kGNS" ) ),
      gammaGNS_( ParameterType::template getValue< double >( "gammaGNS" ) ),
      endTime_ ( ParameterType::template getValue<double>( "femdg.stepper.endtime" )),
      mu_( ParameterType::template getValue< double >( "mu" )),
      k_ ( c_pd() * Pr_inv() * mu_),
      A_( init( true ) ),
      B_( init( false ) )
    {
    }


    // initialize A and B
    RangeFieldType init(const bool returnA ) const ;

    // print info
    void printInitInfo() const;

    // source implementations
    inline bool hasStiffSource() const { return false; }
    inline bool hasNonStiffSource() const  { return true; }
    inline double stiffSource( const double t, const DomainType& x, const RangeType& u, RangeType& res ) const;
    inline double nonStiffSource( const double t, const DomainType& x, const RangeType& u, RangeType& res ) const;

    // this is the initial data
    inline void evaluate( const DomainType& arg , RangeType& res ) const
    {
      evaluate( 0., arg, res );
    }

    // evaluate function
    inline void evaluate( const double t, const DomainType& x, RangeType& res ) const;

    // cloned method
    inline void evaluate( const DomainType& x, const double t, RangeType& res ) const
    {
      evaluate( t, x, res );
    }

    //template< class RangeImp >
    RangeFieldType pressure( const RangeType& u ) const
    {
      return thermodynamics().pressureEnergyForm( u );
    }

    RangeFieldType temperature( const RangeType& u ) const
    {
      return thermodynamics().temperatureEnergyForm( u, pressure( u ) );
    }

    // pressure and temperature
    template< class RangeImp >
    inline void pressAndTemp( const RangeImp& u, RangeFieldType& p, RangeFieldType& T ) const;


    const ThermodynamicsType& thermodynamics() const { return *this; }
    using ThermodynamicsType :: Re ;
    using ThermodynamicsType :: Re_inv;
    using ThermodynamicsType :: Pr;
    using ThermodynamicsType :: Pr_inv;
    using ThermodynamicsType :: c_pd;
    using ThermodynamicsType :: c_pd_inv;
    using ThermodynamicsType :: c_vd;
    using ThermodynamicsType :: c_vd_inv;
    using ThermodynamicsType :: gamma;
    using ThermodynamicsType :: g;
    using ThermodynamicsType :: R_d_inv;
    void printmyInfo( std::string filename ) const {}
    inline double endtime() const { return endTime_; }
    inline std::string myName() const { return myName_; }
    void paraview_conv2prim() const {}
    std::string description() const;

    inline RangeFieldType mu( const RangeType& ) const { return mu_; }
    inline RangeFieldType mu( const RangeFieldType T ) const { return mu_; }
    inline RangeFieldType lambda( const RangeFieldType T ) const { return -2./3.*mu(T); }
    inline RangeFieldType k( const RangeFieldType T ) const { return c_pd() *mu(T) * Pr_inv(); }

  protected:
    const std::string myName_;
    const RangeFieldType omegaGNS_;
    const RangeFieldType kGNS_;
    const RangeFieldType gammaGNS_;
    const double endTime_;
    const RangeFieldType mu_;
    const RangeFieldType k_;
    const RangeFieldType A_;
    const RangeFieldType B_;
  };


  template <class GridType>
  inline typename NSWaves<GridType> :: RangeFieldType
  NSWaves<GridType>:: init(const bool returnA ) const
  {
    if( dimension == 1 )
    {
      if( returnA ) // A
        return (-omegaGNS_ + kGNS_*(3.5*gamma()-2.5));
      else // B
        return ( kGNS_*(0.5+3.5*gamma()) - 4.*omegaGNS_ );
    }
    if( dimension == 2 )
    {
      if( returnA ) // A
        return (-omegaGNS_ + kGNS_*(3.*gamma() - 1.));
      else // B
        return ( kGNS_*(2.+6.*gamma()) - 4.*omegaGNS_ );
    }
    else if( dimension == 3 )
    {
      if( returnA ) // A
        return (-omegaGNS_ + 0.5*kGNS_*(5.*gamma() + 1.));
      else // B
       return (kGNS_*(4.5+7.5*gamma()) - 4.*omegaGNS_);
    }

    abort();
    return RangeFieldType(0);
  }



  template <class GridType>
  inline void NSWaves<GridType>
  :: printInitInfo() const
  {}


  template <class GridType>
  inline double NSWaves<GridType>
  :: stiffSource( const double t, const DomainType& x, const RangeType& u, RangeType& res ) const
  {
    res = 0.;
    return 0.;
  }


  template <class GridType>
  inline double NSWaves<GridType>
  :: nonStiffSource( const double t, const DomainType& x, const RangeType& u, RangeType& res ) const
  {
    /*
    res = 0;
    double sumX = 0;
    for( int i=0; i< dimension; ++i ) sumX += x[i];

    double Frequency=1. ;
    double Amplitude=0.1 ;
    double Pi = M_PI;
    double Omega=Pi*Frequency ;
    double a=1.0*2.*Pi;
    double Kappa = gamma();
    double KappaM1 = gamma() - 1.0;
    double mu0 = mu_ ;

    res[ 0 ] = (-a+3*Omega)*std::cos(Omega*sumX-a*t);
    res[ 1 ] = (-a+Omega/(2.)*(1.+Kappa*5.))*std::cos(Omega*sumX-a*t)
                + Amplitude*Omega*KappaM1*std::sin(2.*(Omega*sumX-a*t));
    for(int i=2; i<energyId; ++i )
      res[ i ] = res[ 1 ];

    res[ energyId ]  = 0.5*((9.+Kappa*15.)*Omega-8.*a)*std::cos(Omega*sumX-a*t)
                +Amplitude*(3.*Omega*Kappa-a)*std::sin(2.*(Omega*sumX-a*t))
                +3.*mu0*Kappa*(Omega*Omega)*Pr_inv()*std::sin(Omega*sumX-a*t);

    res *= Amplitude;
    */

    DomainFieldType sumX = 0;
    for( int i=0; i< dimension; ++i ) sumX += x[i];
    const DomainFieldType beta = kGNS_ * sumX - omegaGNS_*t;
    const DomainFieldType cosBeta = std::cos( beta );
    const DomainFieldType sinBeta = std::sin( beta );
    const DomainFieldType sin2BetaGamma = std::sin( 2.*beta ) * gammaGNS_;
    const DomainFieldType sinGammaKappa = sin2BetaGamma * kGNS_ * (gamma() - 1.);

    res[0] = gammaGNS_*( cosBeta * ( dimension * kGNS_ - omegaGNS_) );
    res[1] = gammaGNS_*( cosBeta * A_ + sinGammaKappa );

    // set other velocity components
    res[2] = res[ 1 ];
    res[ energyId - 1 ] = res[ 1 ];

    res[ energyId ] = gammaGNS_*( sin2BetaGamma*( dimension * gamma()*kGNS_ - omegaGNS_ ) );
    res[ energyId ] += gammaGNS_*( cosBeta *B_ );
    res[ energyId ] += dimension * gammaGNS_ * kGNS_ * kGNS_* k_ * c_vd_inv()
                       * Re_inv() * sinBeta;

    // time step restriction
    return 0.0;
  }



  template <class GridType>
  inline void NSWaves<GridType>
  :: evaluate( const double t, const DomainType& x, RangeType& res ) const
  {
#ifdef WBPROBLEM
    // 10-0.7*(10.*x-x*x*2.) , 10.-x*4.
    res = 0.;
    double z = x[dimension-1];
    res[0] = 10.-4.*z;
    double p = thermodynamics_.p0() - g_*(10.*z - 2.*z*z);
    res[energyId] = p/(gamma()-1.);
    for (int i=0;i<dimension;++i)
      res[energyId] += 0.5*res[i+1]*res[i+1]/res[0];
#else
    /*
    res = 0;

    double sumX = 0;
    for( int i=0; i< dimension; ++i ) sumX += x[i];

    double Frequency=1. ;
    double Amplitude=0.1 ;
    double Pi = M_PI;
    double Omega=Pi*Frequency ;
    double a=1.0*2.*Pi;

    res[ 0 ] = 2.+ Amplitude*std::sin(Omega*sumX - a*t);
    for( int i=1; i<=energyId; ++i )
      res[ i ] = res[ 0 ];
    res[ energyId ] *= res[ 0 ];
    */

    DomainFieldType sumX = 0;
    for( int i=0; i< dimension; ++i ) sumX += x[i];
    const DomainFieldType beta = kGNS_ * sumX - omegaGNS_*t;
    const DomainFieldType sinBeta = std::sin( beta );
    const DomainFieldType sinGamma2 = sinBeta * gammaGNS_ + 2.;

    for(int i=0; i<energyId; ++i)
    {
      res[ i ] = sinGamma2; // constant velocity field
    }
    res[ energyId ] = sinGamma2 * sinGamma2;
#endif
  }



  template <class GridType>
  template< class RangeImp >
  inline void NSWaves<GridType>
  :: pressAndTemp( const RangeImp& u, RangeFieldType& p, RangeFieldType& T ) const
  {
    thermodynamics().pressAndTempEnergyForm( u, p, T );
  }


  template <class GridType>
  inline std::string NSWaves<GridType>
  :: description() const
  {
    std::ostringstream stream;

    stream <<"{\\bf Problem:}" <<myName_
           <<", {\\bf $\\mu$:} " <<mu_
           <<", {\\bf End time:} " <<endTime_
           <<", {$\\gamma_{GNS}$:} " <<gammaGNS_
           <<"\n"
           <<", {$\\omega_{GNS}$:} " <<omegaGNS_
           <<", {$k_{GNS}$:} " <<kGNS_;

    std::string returnString = stream.str();

    return returnString;
  }

} // end namespace
} // end namespace Dune
#endif
