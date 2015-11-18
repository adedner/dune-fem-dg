#ifndef DUNE_THERMODYNAMICS_HH
#define DUNE_THERMODYNAMICS_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{
  // Thermodynamics
  // --------------

  /** \class Thermodynamics
   *  \brief deals with physics used for the atmosphere, in particular
   *         physical constants, equation of state etc.
   *
   *  \tparam dimDomain dimension of the domain
   */
  template< int dimDomain, class Field = double >
  class Thermodynamics
  {
    enum{ energyId = dimDomain+1 };

    public:
    typedef Fem::Parameter  ParameterType ;

    Thermodynamics() :
      Re_( ParameterType::template getValue< double >( "Re" ) ),
      Pr_( ParameterType::template getValue< double >( "Pr" ) ),
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
    inline Field pressureEnergyForm( const RangeType& cons ) const
    {
      // cons = [rho, rho*v, rho*e]
      assert( cons[0] > 1e-20 );
      assert( cons[energyId] > 1e-20 );

      // kinetic energy
      Field kin = 0.;
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
    inline Field temperatureEnergyForm( const RangeType& cons,
                                         const Field p ) const
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
    inline Field temperatureEnergyForm( const RangeType& cons ) const
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
    inline Field densityThetaForm( const RangeType& prim ) const
    {
      const Field p     = prim[ energyId - 1 ];
      const Field theta = prim[ energyId ];

      assert( p > 1e-12 );
      assert( theta > 1e-12 );

      const Field rho = std::pow( double(p/p0_) , double(gamma_inv_) ) * p0_ * R_d_inv_ / theta ;

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
    inline Field g() const { return g_; }
    inline Field c_pd() const { return c_pd_; }
    inline Field c_pd_inv() const { return c_pd_inv_; }
    inline Field c_vd() const { return c_vd_; }
    inline Field c_vd_inv() const { return c_vd_inv_; }
    inline Field c_pl() const { return c_pl_; }
    inline Field c_pl_inv() const { return c_pl_inv_; }
    inline Field R_d() const { return R_d_; }
    inline Field R_d_inv() const { return R_d_inv_; }
    inline Field T0() const { return T0_; }
    inline Field T0_inv() const { return T0_inv_; }
    inline Field L0() const { return L0_; }
    inline Field p0Star() const { return p0Star_; }
    inline Field p0() const { return p0_; }
    inline Field p0_inv() const { return p0_inv_; }
    inline Field gamma() const { return gamma_; }
    inline Field gamma_inv() const { return gamma_inv_; }
    inline Field kappa() const { return kappa_; }
    inline Field kappa_inv() const { return kappa_inv_; }

    //! the Reynolds number Re
    inline Field Re() const { return Re_; }
    //! the inverse Reynolds number (1/Re)
    inline Field Re_inv() const { return Re_inv_; }
    //! the Prandtl number Pr
    inline Field Pr() const { return Pr_; }
    //! the inverse Prandtl number (1/Pr)
    inline Field Pr_inv() const { return Pr_inv_; }

  private:
    const Field Re_;
    const Field Pr_;
    const Field g_;
    const Field p0_;               // surface pressure
    const Field p0Star_;
    const Field T0_;               // freezing temperature
    const Field c_pd_;             // specific heat capacity of dry air w.r.t. pressure
    const Field c_vd_;             // specific heat capacity of water vapour w.r.t. pressure
    const Field c_pl_;             // specific heat capacity of liquid water w.r.t. pressure
    const Field L0_;               // latent heat of evaporasion at 0 Celsius [J/kg]
    const Field relCloud_;

    const Field Re_inv_;
    const Field Pr_inv_;
    const Field c_pd_inv_;
    const Field c_vd_inv_;
    const Field c_pl_inv_;
    const Field R_d_;              // gas constant for dry air
    const Field R_d_inv_;
    const Field p0_inv_;
    const Field T0_inv_;
    const Field kappa_;
    const Field kappa_inv_;
    const Field gamma_;
    const Field gammaM1_;
    const Field gamma_inv_;
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
  template< int dimDomain, class Field >
  template< class RangeType >
  void Thermodynamics< dimDomain, Field >
  :: conservativeToPrimitiveEnergyForm( const RangeType& cons, RangeType& prim ) const
  {
    std::cerr <<"conservativeToPrimitiveEnergyForm not implemented" <<std::endl;
    abort();

    //const Field rho_inv = 1./ cons[0];

    //Field p, T;
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
  template< int dimDomain, class Field >
  template< class RangeType >
  void Thermodynamics< dimDomain, Field >
  :: conservativeToPrimitiveThetaForm( const RangeType& cons, RangeType& prim ) const
  {
    assert( cons[0] > 0. );
    assert( cons[energyId] > 0. );

    Field p, T;
    pressAndTempThetaForm( cons, p, T );

    for( int i = 0; i < dimDomain; ++i )
      prim[i] = cons[i+1]/cons[0];

    prim[energyId-1] = p;
    prim[energyId] = cons[energyId] / cons[0];
  }

  template< int dimDomain, class Field >
  template< class RangeType >
  void Thermodynamics< dimDomain, Field >
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

}
}

#endif // file define
