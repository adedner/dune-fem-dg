#ifndef DUNE_THERMODYNAMICS_HH
#define DUNE_THERMODYNAMICS_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>

using namespace Dune;

// Thermodynamics
// --------------

/** \class Thermodynamics
 *  \brief deals with physics used for the atmosphere, in particular
 *         physical constants, equation of state etc.
 *
 *  \tparam dimDomain dimension of the domain
 */
template< int dimDomain >
class Thermodynamics
{
  enum{ energyId = dimDomain+1 };

  public:
  typedef Fem::Parameter  ParameterType ;

  Thermodynamics() :
    Re_( ParameterType::getValue< double >( "Re" ) ),
    Pr_( ParameterType::getValue< double >( "Pr" ) ),
    g_( ParameterType::getValue< double >( "g", 9.81 ) ),
    p0_(ParameterType::getValue< double >( "p0", 100000. )),
    p0Star_( ParameterType::getValue< double >( "p0Star", 610.7 )),
    T0_( ParameterType::getValue< double >( "T0", 273.15 )),
    c_pd_     ( ParameterType::getValue< double >( "c_pd", 1004. )),
    c_vd_     ( ParameterType::getValue< double >( "c_vd", 717. )),
    c_pl_     ( ParameterType::getValue< double >( "c_pl", 4186. )),
    L0_       ( ParameterType::getValue< double >( "L0", 2500000. )),
    relCloud_ ( ParameterType::getValue< double >( "relCloud", 1. )),
    Re_inv_     ( 1. / Re_ ),
    Pr_inv_     ( 1. / Pr_ ),
    c_pd_inv_   ( 1. / c_pd_ ),
    c_vd_inv_   ( 1. / c_vd_ ),
    c_pl_inv_   ( 1. / c_pl_ ),
    R_d_        ( c_pd_ - c_vd_ ),
    R_d_inv_    ( 1. / R_d_ ),
    p0_inv_     ( 1. / p0_ ),
    T0_inv_     ( 1. / T0_ ),
    kappa_      ( R_d_ * c_pd_inv_ ),
    kappa_inv_  ( 1. / kappa_ ),
    gamma_      ( c_pd_ * c_vd_inv_ ),
    gammaM1_    ( gamma_ - 1.0 ),
    gamma_inv_  ( 1. / gamma_ )
  {
    assert( gamma_ > 1. );
    assert( R_d_ > 200. );
  }

  /** \brief calculate the pressure and the temperature assuming the energy form
   *         for conservative variables: \f$[\rho,\rho\boldsymbol{v},\rho e]\f$
   *
   *  Calculate the pressure and temperature from the conservative variables
   *  \f$[\rho,\rho\boldsymbol{v},\rho e]\f$, where \f$ e \f$ is the sum of 
   *  the internal and the kinetic energy.
   *
   *  \param[in] cons Conervative variables
   *
   *  \return pressure in energy form 
   */
  template< class RangeType >
  inline double pressureEnergyForm( const RangeType& cons ) const
  {
    // cons = [rho, rho*v, rho*e]  
    assert( cons[0] > 1e-20 );
    assert( cons[energyId] > 1e-20 ); 

    // kinetic energy
    double kin = 0.;
    for( int i=1; i<=dimDomain; ++i )
      kin += cons[ i ] * cons[ i ];
    kin *= 0.5 / cons[ 0 ];

    return gammaM1_ * ( cons[ energyId ] - kin );
  }

  /** \brief calculate the pressure and the temperature assuming the energy form
   *         for conservative variables: \f$[\rho,\rho\boldsymbol{v},\rho e]\f$
   *
   *  Calculate the pressure and temperature from the conservative variables
   *  \f$[\rho,\rho\boldsymbol{v},\rho e]\f$, where \f$ e \f$ is the sum of 
   *  the internal and the kinetic energy.
   *
   *  \param[in] cons Conervative variables
   *  \param[out] p Pressure
   *  \param[out] T temperature
   *
   *  \tparam RangeType Type of the range value
   */
  template< class RangeType >
  inline double temperatureEnergyForm( const RangeType& cons, 
                                       const double p ) const 
  {
    assert( cons[0] > 1e-20 );
    return R_d_inv_ * p / cons[ 0 ];
  }

  /** \brief calculate the pressure and the temperature assuming the energy form
   *         for conservative variables: \f$[\rho,\rho\boldsymbol{v},\rho e]\f$
   *
   *  Calculate the pressure and temperature from the conservative variables
   *  \f$[\rho,\rho\boldsymbol{v},\rho e]\f$, where \f$ e \f$ is the sum of 
   *  the internal and the kinetic energy.
   *
   *  \param[in] cons Conervative variables
   *
   *  \return pressure in energy form 
   */
  template< class RangeType >
  inline double temperatureEnergyForm( const RangeType& cons ) const
  {
    return temperatureEnergyForm( cons, pressureEnergyForm( cons ) );
  }

  /** \brief calculate the pressure and the temperature assuming the energy form
   *         for conservative variables: \f$[\rho,\rho\boldsymbol{v},\rho e]\f$
   *
   *  \param[in] cons Conervative variables
   *  \param[out] p Pressure
   *  \param[out] T temperature
   *
   *  \tparam RangeType Type of the range value
   */
  template< class RangeType >
  inline double densityThetaForm( const RangeType& prim ) const 
  {
    const double p     = prim[ energyId - 1 ];
    const double theta = prim[ energyId ];

    assert( p > 1e-12 );
    assert( theta > 1e-12 ); 

    const double rho = std::pow( p/p0_ , gamma_inv_ ) * p0_ * R_d_inv_ / theta ;

    assert( rho > 0.0 );  
    return rho;
  }


  template< class RangeType >
  void conservativeToPrimitiveThetaForm( const RangeType& cons, RangeType& prim ) const;

  template< class RangeType >
  void primitiveToConservativeThetaForm( const RangeType& prim, RangeType& cons ) const;

  template< class RangeType >
  void conservativeToPrimitiveEnergyForm( const RangeType& cons, RangeType& prim ) const;

public:
  inline double g() const { return g_; }
  inline double c_pd() const { return c_pd_; }
  inline double c_pd_inv() const { return c_pd_inv_; }
  inline double c_vd() const { return c_vd_; }
  inline double c_vd_inv() const { return c_vd_inv_; }
  inline double c_pl() const { return c_pl_; }
  inline double c_pl_inv() const { return c_pl_inv_; }
  inline double R_d() const { return R_d_; }
  inline double R_d_inv() const { return R_d_inv_; }
  inline double T0() const { return T0_; }
  inline double T0_inv() const { return T0_inv_; }
  inline double L0() const { return L0_; }
  inline double p0Star() const { return p0Star_; }
  inline double p0() const { return p0_; }
  inline double p0_inv() const { return p0_inv_; }
  inline double gamma() const { return gamma_; }
  inline double gamma_inv() const { return gamma_inv_; }
  inline double kappa() const { return kappa_; }
  inline double kappa_inv() const { return kappa_inv_; }

  //! the Reynolds number Re
  inline double Re() const { return Re_; }
  //! the inverse Reynolds number (1/Re)
  inline double Re_inv() const { return Re_inv_; }
  //! the Prandtl number Pr 
  inline double Pr() const { return Pr_; }
  //! the inverse Prandtl number (1/Pr) 
  inline double Pr_inv() const { return Pr_inv_; }

private:
  const double Re_;
  const double Pr_;
  const double g_;
  const double p0_;               // surface pressure
  const double p0Star_;
  const double T0_;               // freezing temperature
  const double c_pd_;             // specific heat capacity of dry air w.r.t. pressure
  const double c_vd_;             // specific heat capacity of water vapour w.r.t. pressure
  const double c_pl_;             // specific heat capacity of liquid water w.r.t. pressure
  const double L0_;               // latent heat of evaporasion at 0 Celsius [J/kg]
  const double relCloud_;

  const double Re_inv_;
  const double Pr_inv_;
  const double c_pd_inv_;
  const double c_vd_inv_;
  const double c_pl_inv_;
  const double R_d_;              // gas constant for dry air
  const double R_d_inv_;
  const double p0_inv_;
  const double T0_inv_;
  const double kappa_;
  const double kappa_inv_;
  const double gamma_;
  const double gammaM1_;
  const double gamma_inv_;
};


/** \brief converts conservative variables in the energy form to primitive ones
 *
 *  Converts conservative variables \f$[\rho\boldsymbol{v},p,\rho e\f$
 *  to primitive ones \f$[\boldsymbol{v},p,\theta]\f$, where \f$e\f$ is the sum
 *  of internal and kinetic energy, and \f$\theta\f$ potential temperature.
 *
 *  \param[in] cons Conservative variables
 *  \param[out] prim Primitive variables
 *
 *  \tparam dimDomain dimension of the domain
 *  \tparam RangeType type of the range value
 */
template< int dimDomain >
template< class RangeType >
void Thermodynamics< dimDomain >
:: conservativeToPrimitiveEnergyForm( const RangeType& cons, RangeType& prim ) const
{
  std::cerr <<"conservativeToPrimitiveEnergyForm not implemented" <<std::endl;
  abort();

  //const double rho_inv = 1./ cons[0];

  //double p, T;
  //pressAndTempEnergyForm( cons, p, T );

  //prim[energyId-1] = p;
  // this is not pot. temp !!!!!!!
  //prim[energyId] = cons[energyId]/cons[0];
}


/** \brief converts conservative variables in the theta form to primitive ones
 *
 *  Converts conservative variables \f$[\rho\boldsymbol{v},p,\rho\theta]\f$
 *  to primitive ones \f$[\boldsymbol{v},p,\theta]\f$, where \f$\theta\f$ is
 *  potential temperature
 *
 *  \param[in] cons Conservative variables
 *  \param[out] prim Primitive variables
 *
 *  \tparam dimDomain dimension of the domain
 *  \tparam RangeType type of the range value
 */
template< int dimDomain >
template< class RangeType >
void Thermodynamics< dimDomain >
:: conservativeToPrimitiveThetaForm( const RangeType& cons, RangeType& prim ) const
{
  assert( cons[0] > 0. );
  assert( cons[energyId] > 0. );

  double p, T;
  pressAndTempThetaForm( cons, p, T );

  for( int i = 0; i < dimDomain; ++i )
    prim[i] = cons[i+1]/cons[0];

  prim[energyId-1] = p;
  prim[energyId] = cons[energyId] / cons[0];
}

template< int dimDomain >
template< class RangeType >
void Thermodynamics< dimDomain >
:: primitiveToConservativeThetaForm( const RangeType& prim, RangeType& cons ) const
{
  // p,theta  --> rho 
  cons[ 0 ] = densityThetaForm( prim );

  // v_i  --> v_i+1 * rho 
  for( int i = 0; i < dimDomain; ++i )
    cons[ i+1 ] = prim[ i ] * cons[0];

  // theta --> theta * rho 
  cons[ energyId ] = prim[ energyId ] * cons[ 0 ];
}

#endif // file define
