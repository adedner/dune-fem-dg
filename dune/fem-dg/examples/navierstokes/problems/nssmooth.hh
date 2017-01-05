#ifndef NS_SMOOTH_SOLUTION_HH
#define NS_SMOOTH_SOLUTION_HH

#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/fem-dg/examples/navierstokes/thermodynamics.hh>
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


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
  class NSSmoothSolution : public EvolutionProblemInterface<
                              Dune::Fem::FunctionSpace< double, double, GridType::dimension, GridType::dimension + 2 >,
                              true >,
                           public Thermodynamics< GridType::dimensionworld >
  {
    NSSmoothSolution( const NSSmoothSolution& );
    public:
    typedef Fem::FunctionSpace<typename GridType::ctype,
                          double, GridType::dimensionworld,
                          GridType::dimensionworld + 2 > FunctionSpaceType ;

    enum{ dimension = GridType::dimensionworld };
    enum { energyId = dimension + 1 };
    typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
    typedef typename FunctionSpaceType :: DomainType        DomainType;
    typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
    typedef typename FunctionSpaceType :: RangeType         RangeType;
    typedef Thermodynamics< dimension >                     ThermodynamicsType;

    typedef Fem :: Parameter ParameterType;

    NSSmoothSolution() : ThermodynamicsType(),
      myName_( "NSSmoothSolution" ),
      omegaGNS_( ParameterType::template getValue< double >( "omegaGNS" ) ),
      kGNS_( ParameterType::template getValue< double >( "kGNS" ) ),
      gammaGNS_( ParameterType::template getValue< double >( "gammaGNS" ) ),
      endTime_ ( ParameterType::template getValue<double>( "femdg.stepper.endtime" )),
      mu_( ParameterType::template getValue< double >( "mu" )),
      //k_ ( c_pd() * Pr_inv() * mu_),
      k_ ( 2.0 ),
      alpha_( 1.0 ),
      omega_( dimension ),
      A_( init( true ) ),
      B_( init( false ) ),
      e0_( dimension ) // e0 > d/2
    {
    }


    // initialize A and B
    double init(const bool returnA ) const ;

    // print info
    void printInitInfo() const {}

    // source implementations
    inline bool hasStiffSource() { return false; }
    inline bool hasNonStiffSource() { return true; }
    inline double stiffSource( const double t, const DomainType& x, const RangeType& u, RangeType& res ) const
    {
      res = 0;
      return 0;
    }

    //! beta = k * sum_i x_i  - omega * t
    const double beta(const double t, const DomainType& x ) const
    {
      double sumX = 0;
      for( int i=0; i< dimension; ++i ) sumX += x[i];
      return k_ * sumX - omega_ * t ;
    }

    inline double nonStiffSource( const double t, const DomainType& x, const RangeType& u, RangeType& res ) const
    {
      const double betaa = beta( t, x );
      const double cosBeta = std::cos( betaa );
      res = cosBeta * alpha_ * ( dimension * k_ - omega_);

      // velocity components
      const double add = (gamma() - 1.0) * ( e0_ - double(dimension) * 0.5 )
                         * ( cosBeta * alpha_ * k_ );

      res[ energyId ] *= e0_ ;

      for( int i=1; i< RangeType :: dimension ; ++i ) res[ i ] += add ;

      // TODO: calculate time step restriction due to source term
      return 0;
    }

    // this is the initial data
    inline void evaluate( const DomainType& arg , RangeType& res ) const
    {
      evaluate( 0., arg, res );
    }

    // evaluate function
    inline void evaluate( const double t, const DomainType& x, RangeType& res ) const
    {
      // sinus waves
      res = std::sin( beta( t, x ) ) * alpha_ + 2.0 ;

      // add factor for energy component
      res[ energyId ] *= e0_ ;

      return ;
    }

    // cloned method
    inline void evaluate( const DomainType& x, const double t, RangeType& res ) const
    {
      evaluate( t, x, res );
    }

    //template< class RangeImp >
    double pressure( const RangeType& u ) const
    {
      return thermodynamics().pressureEnergyForm( u );
    }

    double temperature( const RangeType& u ) const
    {
      return thermodynamics().temperatureEnergyForm( u, pressure( u ) );
    }

    // pressure and temperature
    template< class RangeImp >
    inline void pressAndTemp( const RangeImp& u, double& p, double& T ) const
    {
      thermodynamics().pressAndTempEnergyForm( u, p, T );
    }

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

    inline double mu( const RangeType& ) const { return mu_; }
    inline double mu( const double T ) const { return mu_; }
    inline double lambda( const double T ) const { return -2./3.*mu(T); }
    inline double k( const double T ) const { return c_pd() *mu(T) * Pr_inv(); }

  protected:
    const std::string myName_;
    const double omegaGNS_;
    const double kGNS_;
    const double gammaGNS_;
    const double endTime_;
    const double mu_;
    const double k_;
    const double alpha_;
    const double omega_ ;
    const double A_;
    const double B_;
    const double e0_;
  };


  template <class GridType>
  inline double NSSmoothSolution<GridType>
  :: init(const bool returnA ) const
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
    return 0;
  }


  template <class GridType>
  inline std::string NSSmoothSolution<GridType>
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
