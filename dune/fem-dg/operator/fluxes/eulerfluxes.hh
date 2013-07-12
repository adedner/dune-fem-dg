#ifndef DUNE_FEM_DG_EULERFLUXES_HH
#define DUNE_FEM_DG_EULERFLUXES_HH
  
// system includes
#include <string>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

// dune-grid includes
#if WELLBALANCE
#include <dune/grid/common/genericreferenceelements.hh>
#endif

#include <dune/fem-dg/operator/fluxes/rotator.hh>

///////////////////////////////////////////////////////////
//
//  EulerAnalyticalFlux
//
///////////////////////////////////////////////////////////
template < int dimDomain >
class EulerAnalyticalFlux
{
 public:
  enum { e = dimDomain+1 };
  typedef double FieldType;
  typedef Dune :: FieldVector< FieldType, dimDomain> DomainType;

  template <class RangeType, class FluxRangeType>
  inline void analyticalFlux(const FieldType gamma,
                             const RangeType& u,
                             FluxRangeType& f) const;

  template <class RangeType, class FluxRangeType>
  inline  void jacobian(const FieldType gamma,
                        const RangeType& u,
                        const FluxRangeType& du,
                        RangeType& A) const;

  template <class RangeType>
  inline FieldType rhoeps(const RangeType& u) const;

  template <class RangeType>
  inline FieldType pressure(const FieldType gamma,const RangeType& u) const
  {
    assert(u[0]>1e-10);
    const FieldType rhoe = rhoeps(u);
    assert( rhoe>1e-10 );
    return (gamma-1)*rhoe;
  }

  template <class RangeType>
  inline FieldType maxSpeed(const FieldType gamma, const DomainType& n, const RangeType& u) const;
};

// ***********************
template <>
template <class RangeType> 
inline
EulerAnalyticalFlux<1>::FieldType EulerAnalyticalFlux<1>::rhoeps(const RangeType& u) const 
{
  assert( u[0] >= 1e-10 );
  return u[e]-0.5*(u[1]*u[1])/u[0];
}
template <>
template <class RangeType>
inline
EulerAnalyticalFlux<1>::FieldType EulerAnalyticalFlux<1>::pressure(const FieldType gamma,const RangeType& u) const {
  FieldType rhoe = rhoeps(u);
  assert(rhoe>1e-10);
  return (gamma-1)*rhoe;
}
template <>
template <class RangeType, class FluxRangeType>
inline
void EulerAnalyticalFlux<1>::analyticalFlux(const FieldType gamma,
          const RangeType& u,
          FluxRangeType& f) const {
  assert(u[0]>1e-10);
  const FieldType p = pressure(gamma, u);
  f[0][0] = u[1];
  f[1][0] = u[1]/u[0]*u[1]+p;
  f[e][0] = u[1]/u[0]*(u[e]+p);
}

template <>
template <class RangeType> 
inline
EulerAnalyticalFlux<2>::FieldType EulerAnalyticalFlux<2>::rhoeps(const RangeType& u) const 
{
  assert( u[0] >= 1e-10 );
  return u[e]-0.5*(u[1]*u[1]+u[2]*u[2])/u[0];
}

template <>
template <class RangeType>
inline
EulerAnalyticalFlux<2>::FieldType EulerAnalyticalFlux<2>::pressure(const FieldType gamma,const RangeType& u) const 
{
  const FieldType re = rhoeps(u);
  assert(re>=1e-10);
  return (gamma-1)*re;
}

template <>
template <class RangeType, class FluxRangeType>
inline
void EulerAnalyticalFlux<2>::
analyticalFlux(const FieldType gamma,
               const RangeType& u,
               FluxRangeType& f) const 
{
  assert( u[0] >= 1e-10 );
  const FieldType v[2] = {u[1]/u[0],u[2]/u[0]};
  const FieldType p = pressure( gamma, u );
  const FieldType ue_p = (u[e]+p);
  f[0][0] = u[1];            f[0][1] = u[2];
  f[1][0] = v[0]*u[1]+p;     f[1][1] = v[1]*u[1];
  f[2][0] = v[0]*u[2];       f[2][1] = v[1]*u[2]+p;
  f[e][0] = v[0]*ue_p;       f[e][1] = v[1]*ue_p;
}

template <>
template <class RangeType>
inline
EulerAnalyticalFlux<3>::FieldType EulerAnalyticalFlux<3>::rhoeps(const RangeType& u) const 
{
  assert( u[0] >= 1e-10 );
  return u[e]-0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3])/u[0];
}

template <>
template <class RangeType>
inline
EulerAnalyticalFlux<3>::FieldType EulerAnalyticalFlux<3>::
pressure(const FieldType gamma,const RangeType& u) const 
{
  assert(u[0]>1e-10);
  const FieldType rhoe = rhoeps(u);
  assert( rhoe>1e-10 );
  return (gamma-1)*rhoe;
}

template <>
template <class RangeType, class FluxRangeType>
inline
void EulerAnalyticalFlux<3>::analyticalFlux(const FieldType gamma,
                                            const RangeType& u,
                                            FluxRangeType& f) const 
{
  assert(u[0]>1e-10);
  const FieldType p = pressure(gamma,u);
  // get the velocity 
  const FieldType v[3] = {u[1]/u[0], u[2]/u[0], u[3]/u[0]};
  const FieldType ue_p = (u[e]+p);
  
  f[0][0]=u[1];          f[0][1]=u[2];           f[0][2]=u[3];
  f[1][0]=v[0]*u[1]+p;   f[1][1]=v[1]*u[1];      f[1][2]=v[2]*u[1];
  f[2][0]=v[0]*u[2];     f[2][1]=v[1]*u[2]+p;    f[2][2]=v[2]*u[2];
  f[3][0]=v[0]*u[3];     f[3][1]=v[1]*u[3];      f[3][2]=v[2]*u[3]+p;
  f[e][0]=v[0]*ue_p;     f[e][1]=v[1]*ue_p;      f[e][2]=v[2]*ue_p;
}

template <>
template <class RangeType, class FluxRangeType>
inline
void EulerAnalyticalFlux<2>::jacobian(const FieldType gamma,
                                      const RangeType& u, 
                                      const FluxRangeType& du,
                                      RangeType& A) const 
{
  assert(u[0]>1e-10);
  FieldType v[2] = {u[1]/u[0],u[2]/u[0]};
  
  A[0] = du[1][0] + du[1][1];
  A[1] = du[0][0]*((gamma-3.0)/2.0*v[0]*v[0] +(gamma-1.0)/2.0*v[1]*v[1])-
    du[0][1]*v[0]*v[1];
  A[1] += du[1][0]*(3.0-gamma)*v[0] + du[1][1]*v[1];
  A[1] += du[2][0]*(1.0-gamma)*v[1] + du[2][1]*v[0];
  A[1] += du[3][0]*(gamma-1.0);
  A[2] = du[0][1]*((gamma-3.0)/2.0*v[1]*v[1] +(gamma-1.0)/2.0*v[0]*v[0])-
    du[0][0]*v[0]*v[1];
  A[2] += du[1][1]*(1.0-gamma)*v[0] + du[1][0]*v[1];
  A[2] += du[2][1]*(3.0-gamma)*v[1] + du[2][0]*v[0];
  A[2] += du[3][1]*(gamma-1.0);
  A[3] = du[0][0]*(-gamma*u[3]*v[0]/u[0]+(gamma-1.0)*v[0]*(v[0]*v[0]+v[1]*v[1]))+
    du[0][1]*(-gamma*u[3]*v[1]/u[0]+(gamma-1.0)*v[1]*(v[0]*v[0]+v[1]*v[1]));
  A[3] += du[1][0]*(gamma*u[3]/u[0]-(gamma-1.0)/2.0*(3.0*v[0]*v[0]+v[1]*v[1]))-
    du[1][1]*(gamma-1.0)*v[0]*v[1];
  A[3] += -du[2][0]*(gamma-1.0)*v[0]*v[1]+
    du[2][1]*(gamma*u[3]/u[0]-(gamma-1.0)/2.0*(v[0]*v[0]+3.0*v[1]*v[1]));
  A[3] += du[3][0]*gamma*v[0] + du[3][1]*gamma*v[1];
}

template <>
template <class RangeType>
inline
EulerAnalyticalFlux<1>::FieldType 
EulerAnalyticalFlux<1>::maxSpeed(const FieldType gamma,
                                 const DomainType& n,
                                 const RangeType& u) const 
{
  assert(u[0]>1e-10);
  FieldType u_normal = u[1]*n[0] / u[0];
  FieldType p = pressure(gamma,u);
  const FieldType c2 = gamma*p/u[0]*n.two_norm2();
  assert(c2>1e-10);
  return std::abs(u_normal) + std::sqrt(c2);
}

template <>
template <class RangeType>
inline
EulerAnalyticalFlux<2>::FieldType 
EulerAnalyticalFlux<2>::maxSpeed(const FieldType gamma,
                                 const DomainType& n,
                                 const RangeType& u) const 
{
  assert( u[0] > 1e-10 );
  FieldType u_normal = (u[1]*n[0]+u[2]*n[1]) / u[0];
  FieldType p = pressure(gamma,u);
  FieldType c2 = gamma * p/ u[0] * n.two_norm2();
  assert( c2 > 1e-10 );
  return std::abs(u_normal) + std::sqrt(c2);
}

template <>
template <class RangeType>
inline
EulerAnalyticalFlux<3>::FieldType 
EulerAnalyticalFlux<3>::maxSpeed(const FieldType gamma,
                                 const DomainType& n,
                                 const RangeType& u) const 
{
  assert(u[0]>1e-10);
  const FieldType u_normal = (u[1]*n[0]+u[2]*n[1]+u[3]*n[2]) / u[0];
  const FieldType p = pressure(gamma,u);
  FieldType c2 = gamma*p/u[0]*n.two_norm2();
  assert(c2>1e-10);
  return std::abs(u_normal) + std::sqrt(c2);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// local includes
namespace EULERNUMFLUX
{
  ////////////////////////////////////////////////////////
  //
  // Implementation of Euler fluxes from Dennis Diehl.
  //
  ////////////////////////////////////////////////////////

  typedef enum {LLF, HLL, HLLC} EulerFluxType;
  const EulerFluxType std_flux_type = HLL;

  template< class Model, EulerFluxType flux_type=std_flux_type> 
  class EulerFlux
  {
  public:
    typedef Model  ModelType ;
    enum { dim = Model :: dimDomain };

    EulerFlux( const Model& model );

    void flux(const double U[dim+2], double *f[dim]) const;

    double num_flux(const double Uj[dim+2], const double Un[dim+2], 
                    const double normal[dim], double gj[dim+2]) const;

    const Model& model_;
    const double _gamma;

  private:
    double num_flux_LLF(const double Uj[dim+2], const double Un[dim+2], 
                        const double normal[dim], double gj[dim+2]) const;

    double num_flux_HLL(const double Uj[dim+2], const double Un[dim+2], 
                        const double normal[dim], double gj[dim+2]) const;
    
    double num_flux_HLLC(const double Uj[dim+2], const double Un[dim+2], 
                         const double normal[dim], double gj[dim+2]) const;


    static void rotate(const double normal[dim], 
                       const double u[dim], double u_rot[dim]);

    static void rotate_inv(const double normal[dim], 
                           const double u_rot[dim], double u[dim]);
  };



  // ======= class Euler inline implementation =================
  template<class Model, EulerFluxType flux_type>
  inline
  EulerFlux<Model,flux_type>::EulerFlux( const ModelType& model ) 
    : model_( model ), 
      _gamma( model.gamma() )
  {}

  template<class Model, EulerFluxType flux_type>
  inline
  void EulerFlux<Model,flux_type>::rotate(const double n[dim], 
                                        const double u[dim], double u_rot[dim])
  {
    if (dim == 1)
    {
      u_rot[0] = n[0] * u[0];
    }
    else if (dim == 2)
    {
      u_rot[0] = n[0]*u[0] + n[1]*u[1];
      u_rot[1] = -n[1]*u[0] + n[0]*u[1];
    }
    else if (dim == 3)
    {
      double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

      if (d > 1.0e-8) {
        double d_1 = 1.0/d;
        u_rot[0] = n[0]*u[0]           + n[1]*u[1]          + n[2]*u[2];
        u_rot[1] = -n[1]*d_1*u[0]      + n[0]*d_1* u[1];
        u_rot[2] = -n[0]*n[2]*d_1*u[0] - n[1]*n[2]*d_1*u[1] + d*u[2];
      } 
      else {
        u_rot[0] = n[2]*u[2];
        u_rot[1] = u[1];
        u_rot[2] = -n[2]*u[0];
      }
      //assert(0); // test it, not tested up to now
    }
    else
      DUNE_THROW( Dune::NotImplemented, "EulerFlux::rotate does not support dim>3." );
  }


  template<class Model, EulerFluxType flux_type>
  inline
  void EulerFlux<Model,flux_type>::rotate_inv(const double n[dim], 
                                           const double u_rot[dim], 
                                           double u[dim])
  {
    if (dim == 1){
      u[0] = n[0] * u_rot[0];
    }

    if (dim == 2){
      u[0] = n[0]*u_rot[0] - n[1]*u_rot[1];
      u[1] = n[1]*u_rot[0] + n[0]*u_rot[1];
    }

    if (dim == 3){
      double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

      if (d > 1.0e-8) {
        double d_1 = 1.0/d;
        u[0] = n[0]*u_rot[0] - n[1]*d_1*u_rot[1] - n[0]*n[2]*d_1*u_rot[2];
        u[1] = n[1]*u_rot[0] + n[0]*d_1*u_rot[1] - n[1]*n[2]*d_1*u_rot[2];
        u[2] = n[2]*u_rot[0]                     + d*u_rot[2];
      } 
      else {
        u[0] = -n[2]*u_rot[2];
        u[1] = u_rot[1];
        u[2] = n[2]*u_rot[0];
      }
      
      //assert(0); // test it, not tested up to now
    }

    if (dim > 3) assert(0);
  }


  // U[0] = rho, (U[1],...,U[dim]) = rho_\vect u, U[dim+1] = E
  template<class Model, EulerFluxType flux_type>
  inline
  void EulerFlux<Model,flux_type>::flux(const double U[dim+2], 
                                        double *f[dim]) const
  {
    const double rho = U[0];
    const double *rho_u = &U[1];
    const double E = U[dim+1];

    double u[dim], Ekin2 = 0.0;
    for(int i=0; i<dim; i++){
      u[i] = (1.0/rho) * rho_u[i];
      Ekin2 += rho_u[i] * u[i];
    }
   
    const double p = (_gamma-1.0)*(E - 0.5*Ekin2);

    for(int i=0; i<dim; i++){
      f[i][0] = rho_u[i];

      for(int j=0; j<dim; j++) f[i][1+j] = rho_u[i] * u[j];
      f[i][1+i] += p;

      f[i][dim+1] = (E+p) * u[i];
    }
  }



  // returns fastest wave speed
  // U[0] = rho, (U[1],...,U[dim]) = rho_\vect u, U[dim+1] = E
  template<class Model, EulerFluxType flux_type>
  inline
  double EulerFlux<Model,flux_type>::num_flux(const double Uj[dim+2], 
                                              const double Un[dim+2], 
                                              const double normal[dim], 
                                              double gj[dim+2]) const
  {
    if (flux_type == LLF) return num_flux_LLF(Uj, Un, normal, gj);
   
    if (flux_type == HLL) return num_flux_HLL(Uj, Un, normal, gj);

    if (flux_type == HLLC) return num_flux_HLLC(Uj, Un, normal, gj);

    DUNE_THROW( Dune::NotImplemented, "Numerical flux not implemented" );
  }


  // returns fastest wave speed
  // U[0] = rho, (U[1],...,U[dim]) = rho_\vect u, U[dim+1] = E
  template<class Model, EulerFluxType flux_type>
  inline
  double EulerFlux<Model,flux_type>::num_flux_LLF(const double Uj[dim+2], 
                                                  const double Un[dim+2], 
                                                  const double normal[dim], 
                                                  double gj[dim+2]) const
  {
    const double rhoj = Uj[0];
    const double *rho_uj = &Uj[1];
    const double Ej = Uj[dim+1];
    const double rhon = Un[0];
    const double *rho_un = &Un[1];
    const double En = Un[dim+1];

    double uj[dim], Ekin2j=0.0, un[dim], Ekin2n=0.0;
    double u_normal_j=0.0, u_normal_n=0.0;
    for(int i=0; i<dim; i++){
      uj[i] = (1.0/rhoj) * rho_uj[i];
      un[i] = (1.0/rhon) * rho_un[i];
      Ekin2j += rho_uj[i] * uj[i];
      Ekin2n += rho_un[i] * un[i];
      u_normal_j += uj[i] * normal[i];
      u_normal_n += un[i] * normal[i];
    }

    const double pj = (_gamma-1.0)*(Ej - 0.5*Ekin2j);
    const double cj = sqrt(_gamma*pj/rhoj);
    const double pn = (_gamma-1.0)*(En - 0.5*Ekin2n);
    const double cn = sqrt(_gamma*pn/rhon);

    assert(rhoj>0.0 && pj>0.0 && rhoj>0.0 && pj>0.0);

    const double alphaj = fabs(u_normal_j) + cj;
    const double alphan = fabs(u_normal_n) + cn;
    const double alpha = (alphaj > alphan)? alphaj : alphan;

    gj[0] = gj[dim+1] = 0.0;
    for(int i=0; i<dim; i++) gj[1 + i] = 0.0;

    for(int j=0; j<dim; j++){
      gj[0] += ( rho_uj[j] + rho_un[j] ) * normal[j];

      for(int i=0; i<dim; i++){
        gj[1 + i] += (rho_uj[i]*uj[j] + rho_un[i]*un[j]) * normal[j];
      }

      gj[dim+1] += ( (Ej+pj)*uj[j] + (En+pn)*un[j] ) * normal[j];
    }

    gj[0] = 0.5 * (gj[0] - alpha*(rhon - rhoj));
    for(int i=0; i<dim; i++){
      gj[1+i] = 0.5*(gj[1+i] + (pj+pn)*normal[i] - alpha*(rho_un[i]-rho_uj[i]));
    }
    gj[dim+1] = 0.5 * (gj[dim+1] - alpha*(En - Ej));

    return alpha;
  }



  // returns fastest wave speed
  // U[0] = rho, (U[1],...,U[dim]) = rho_\vect u, U[dim+1] = E
  template<class Model, EulerFluxType flux_type>
  inline
  double EulerFlux<Model,flux_type>::num_flux_HLL(const double Uj[dim+2], 
                                                  const double Un[dim+2], 
                                                  const double normal[dim], 
                                                  double gj[dim+2]) const

  {
    const double rhoj = Uj[0];
    double Ej = Uj[dim+1];
    const double rhon = Un[0];
    double En = Un[dim+1];
#ifdef POTTEMP
    double pressj, tempj;
    double pressn, tempn;
    model_.pressAndTemp( Uj, pressj, tempj );
    model_.pressAndTemp( Un, pressn, tempn );

    double Ekinj = 0;
    double Ekinn = 0;
    for(int i=1; i<dim+1; i++){
      Ekinj += (0.5/rhoj) * Uj[i] * Uj[i];
      Ekinn += (0.5/rhon) * Un[i] * Un[i];
    }
    Ej = pressj/(_gamma-1.) + Ekinj;
    En = pressn/(_gamma-1.) + Ekinn;
#endif

    double rho_uj[dim], rho_un[dim], uj[dim], un[dim];
    double Ekin2j=0.0, Ekin2n=0.0;
    rotate(normal, Uj+1, rho_uj);
    rotate(normal, Un+1, rho_un);
    for(int i=0; i<dim; i++){
      uj[i] = (1.0/rhoj) * rho_uj[i];
      un[i] = (1.0/rhon) * rho_un[i];
      Ekin2j += rho_uj[i] * uj[i];
      Ekin2n += rho_un[i] * un[i];
    }

    const double pj = (_gamma-1.0)*(Ej - 0.5*Ekin2j);
    const double pn = (_gamma-1.0)*(En - 0.5*Ekin2n);

    const double cj = sqrt(_gamma*pj/rhoj);
    const double cn = sqrt(_gamma*pn/rhon);

    assert(rhoj>0.0 && pj>0.0 && rhoj>0.0 && pj>0.0);

    const double rho_bar = 0.5 * (rhoj + rhon);
    const double c_bar = 0.5 * (cj + cn);
    const double p_star = 0.5 * ( (pj+pn) - (un[0]-uj[0])*rho_bar*c_bar );
    const double u_star = 0.5 * ( (uj[0]+un[0]) - (pn-pj)/(rho_bar*c_bar) );
    const double tmp = 0.5*(_gamma+1.0)/_gamma;
    const double qj = (p_star > pj)? sqrt( 1.0 + tmp*(p_star/pj - 1.0) ): 1.0;
    const double qn = (p_star > pn)? sqrt( 1.0 + tmp*(p_star/pn - 1.0) ): 1.0;

    const double sj = uj[0] - cj*qj;
    const double sn = un[0] + cn*qn;

    double guj[dim];

    if (u_star > 0.0){
      if (sj >= 0.0){
        gj[0] = rho_uj[0];

        for(int i=0; i<dim; i++) guj[i] = rho_uj[i]*uj[0];
        guj[0] += pj;

#ifdef POTTEMP
        gj[dim+1] = Uj[dim+1]*uj[0];
#else 
        gj[dim+1] = (Ej+pj)*uj[0];
#endif
      }
      else{
        const double tmp1 = sj * sn;
        const double tmp2 = 1.0/(sn - sj);      
        gj[0] = tmp2 * ( sn*rho_uj[0] - sj*rho_un[0] + tmp1*(rhon - rhoj) );

        for(int i=0; i<dim; i++){
          guj[i] = tmp2*((sn*uj[0]-tmp1)*rho_uj[i] - (sj*un[0]-tmp1)*rho_un[i]);
        }
        guj[0] += tmp2 * (sn*pj - sj*pn);

#ifdef POTTEMP
        const double Etmpj = Uj[dim+1]*uj[0];
        const double Etmpn = Un[dim+1]*un[0];
        gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(Un[dim+1] - Uj[dim+1]));
#else
        const double Etmpj = (Ej+pj)*uj[0];
        const double Etmpn = (En+pn)*un[0];
        gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(En - Ej));
#endif
      }
    }
    else{
      if (sn <= 0.0){
        gj[0] = rho_un[0];

        for(int i=0; i<dim; i++) guj[i] = rho_un[i]*un[0];
        guj[0] += pn;

#ifdef POTTEMP
        gj[dim+1] = Uj[dim+1]*un[0];
#else 
        gj[dim+1] = (En+pn)*un[0];
#endif
      }
      else{
        const double tmp1 = sj * sn;
        const double tmp2 = 1.0/(sn - sj);
        gj[0] = tmp2 * ( sn*rho_uj[0] - sj*rho_un[0] + tmp1*(rhon - rhoj) );

        for(int i=0; i<dim; i++){
          guj[i] = tmp2*((sn*uj[0]-tmp1)*rho_uj[i] - (sj*un[0]-tmp1)*rho_un[i]);
        }
        guj[0] += tmp2 * (sn*pj - sj*pn);

#ifdef POTTEMP
        const double Etmpj = Uj[dim+1]*uj[0];
        const double Etmpn = Un[dim+1]*un[0];
        gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(Un[dim+1] - Uj[dim+1]));
#else
        const double Etmpj = (Ej+pj)*uj[0];
        const double Etmpn = (En+pn)*un[0];
        gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(En - Ej));
#endif

      }
    }

    rotate_inv(normal, guj, gj+1);
    return (fabs(sj) > fabs(sn))? fabs(sj): fabs(sn);
  }




  template<class Model, EulerFluxType flux_type>
  inline
  double EulerFlux<Model,flux_type>::num_flux_HLLC(const double Um[dim+2], 
                                                   const double Up[dim+2], 
                                                   const double normal[dim], 
                                                   double g[dim+2]) const
  {
    const double rhom = Um[0];
    const double rhop = Up[0];
    const double Em = Um[dim+1];
    const double Ep = Up[dim+1];

    double rho_um[dim], rho_up[dim];
    rotate( normal, Um+1, rho_um );
    rotate( normal, Up+1, rho_up );

    double Ekinm = 0.;
    double Ekinp = 0.;
    double um[dim], up[dim];
    for( int i=0; i<dim; ++i )
    {
      um[i] = rho_um[i] / rhom;
      up[i] = rho_up[i] / rhop;
      Ekinm += rho_um[i] * um[i];
      Ekinp += rho_up[i] * up[i];
    }

    const double pm = (_gamma-1.0)*(Em - 0.5*Ekinm);
    const double pp = (_gamma-1.0)*(Ep - 0.5*Ekinp);

    assert( rhom>0.0 && pm>0.0 && rhop>0.0 && pp>0.0 );

    const double cm = sqrt(_gamma*pm/rhom);
    const double cp = sqrt(_gamma*pp/rhop);

    const double rho_bar = 0.5 * (rhom + rhop);
    const double c_bar = 0.5 * (cm + cp);
    const double p_star = 0.5 * ( (pm+pp) - (up[0]-um[0])*rho_bar*c_bar );
    const double u_star = 0.5 * ( (um[0]+up[0]) - (pp-pm)/(rho_bar*c_bar) );
    const double tmp = 0.5*(_gamma+1.0)/_gamma;
    const double qm = (p_star > pm) ? sqrt( 1.0 + tmp*(p_star/pm - 1.0) ) : 1.0;
    const double qp = (p_star > pp) ? sqrt( 1.0 + tmp*(p_star/pp - 1.0) ) : 1.0;

    const double sm = um[0] - cm*qm;
    const double sp = up[0] + cp*qp;

    double guj[dim];

    if (sm >= 0.0)
    {
      g[0] = rho_um[0];

      for(int i=0; i<dim; i++) 
        guj[i] = rho_um[i]*um[0];
      guj[0] += pm;

      g[dim+1] = (Em+pm)*um[0];
    }
    else if (sp <= 0.0)
    {
      g[0] = rho_up[0];

      for(int i=0; i<dim; i++) 
        guj[i] = rho_up[i]*up[0];
      guj[0] += pp;

      g[dim+1] = (Ep+pp)*up[0];
    }
    else
    {
      const double tmpm = sm*(sm-um[0])/(sm-u_star);
      const double tmpp = sp*(sp-up[0])/(sp-u_star);

      if (u_star >= 0.0)
      {
        g[0] = rho_um[0] + rhom*(tmpm-sm);

        for(int i=0; i<dim; i++) 
          guj[i] = rho_um[i]*um[0] + rhom*um[i]*tmpm - sm*rho_um[i];
        guj[0] += pm + rhom*(u_star-um[0])*tmpm;

        g[dim+1] = (Em+pm)*um[0] + Em*(tmpm-sm) 
          + tmpm*(u_star-um[0])*( rhom*u_star + pm/(sm-um[0]) );
      }
      else
      {
        g[0] = rho_up[0] + rhop*(tmpp-sp);

        for(int i=0; i<dim; i++) 
          guj[i] = rho_up[i]*up[0] + rhop*up[i]*tmpp - sp*rho_up[i];
        guj[0] += pp + rhop*(u_star-up[0])*tmpp;

        g[dim+1] = (Ep+pp)*up[0] + Ep*(tmpp-sp) 
          + tmpp*(u_star-up[0])*( rhop*u_star + pp/(sp-up[0]) );
      }
    }

    rotate_inv( normal, guj, g+1 );

    return (fabs(sm) > fabs(sp)) ? fabs(sm) : fabs(sp);
  }
}


////////////////////////////////////////////////////////
//
// Dune interface for Euler fluxes 
//
////////////////////////////////////////////////////////
// LLFFlux
//--------

template< class Model >
class LLFFlux
{
public:
  typedef Model                                       ModelType;
  enum { dimDomain = Model::dimDomain };
  enum { dimRange = Model::dimRange };
  typedef typename Model::Traits                      Traits;
  typedef typename Traits::GridType                   GridType;
  typedef typename GridType::ctype                    ctype;
  typedef typename Traits::EntityType                 EntityType;

  typedef typename Traits::DomainType                 DomainType;
  typedef typename Traits::FaceDomainType             FaceDomainType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::FluxRangeType              FluxRangeType;

  LLFFlux( const Model& mod )
    : model_(mod)
  {}

  static std::string name () { return "LLF"; }

  const Model& model() const { return model_; }

  // Return value: maximum wavespeed*length of integrationOuterNormal
  // gLeft,gRight are fluxed * length of integrationOuterNormal
  template< class Intersection, class QuadratureImp >
  inline double 
  numericalFlux( const Intersection& intersection,
                 const EntityType& inside,
                 const EntityType& outside,
                 const double time, 
                 const QuadratureImp& faceQuadInner,
                 const QuadratureImp& faceQuadOuter,
                 const int quadPoint,
                 const RangeType& uLeft, 
                 const RangeType& uRight,
                 RangeType& gLeft,
                 RangeType& gRight) const 
  {
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
    DomainType normal = intersection.integrationOuterNormal(x);  
    const double len = normal.two_norm();
    normal *= 1./len;
    
    RangeType visc;
    FluxRangeType anaflux;

    model_.advection( inside, time, faceQuadInner.point( quadPoint ),
                      uLeft, anaflux );

    // set gLeft 
    anaflux.mv( normal, gLeft );

    model_.advection( outside, time, faceQuadOuter.point( quadPoint ),
                      uRight, anaflux );
    anaflux.umv( normal, gLeft );

    double maxspeedl, maxspeedr, maxspeed;
    double viscparal, viscparar, viscpara;
    
    const DomainType xGlobal = intersection.geometry().global(x);

    model_.maxSpeed( normal, time, xGlobal, 
                     uLeft, viscparal, maxspeedl );
    model_.maxSpeed( normal, time, xGlobal,
                     uRight, viscparar, maxspeedr );

    maxspeed = (maxspeedl > maxspeedr) ? maxspeedl : maxspeedr;
    viscpara = (viscparal > viscparar) ? viscparal : viscparar;
    visc = uRight;
    visc -= uLeft;
    visc *= viscpara;
    gLeft -= visc;
    
    gLeft *= 0.5*len;
    gRight = gLeft;

#if WELLBALANCE
    const double g = model_.problem().g();

    // calculate geopotential in the grid elements sharing 'it'
    const double z = xGlobal[dimDomain-1];

    // calculate num. flux for pressure
    double pInside, TInside;
    double pOutside, TOutside;
    model_.pressAndTemp( uLeft, pInside, TInside );
    model_.pressAndTemp( uRight, pOutside, TOutside );
    const double pNumFlux = 0.5*(pInside + pOutside);
    
    const double rhoAvg = 0.5*(uLeft[0] + uRight[0]);
    const double rhoJump = uRight[0] - uLeft[0];

    // add well-balancing terms to vertical momentum flux, scaled with normal
    gLeft[dimDomain]  += 0.5*normal[dimDomain-1]*( -2.*(pNumFlux-pInside))*len;
    gRight[dimDomain] += 0.5*normal[dimDomain-1]*( -2.*(pNumFlux-pOutside))*len;

    gLeft[dimDomain]  -= 0.25*normal[dimDomain-1]*( g*z*rhoJump )*len;
    gRight[dimDomain] += 0.25*normal[dimDomain-1]*( g*z*rhoJump )*len;
#endif

    return maxspeed * len;
  }

 protected:
  const Model& model_;
};

/////////////////////////////////////////////////////////////
//
//  Dennis Flux implementations 
//
////////////////////////////////////////////////////////////
template <class Model, 
          EULERNUMFLUX::EulerFluxType fluxtype>
class NumFluxBase {
 public:
  typedef Model ModelType;
  typedef typename Model::Traits::GridType GridType;
  typedef typename GridType::ctype                    ctype;
  enum { dimDomain = GridType::dimensionworld };
  enum { dimRange = Model::dimRange };
  typedef typename Model::Traits Traits;
  typedef typename Traits::EntityType                 EntityType;

  typedef typename Traits::DomainType                 DomainType;
  typedef typename Traits::FaceDomainType             FaceDomainType;
  typedef typename Traits::RangeType                  RangeType;
  typedef typename Traits::FluxRangeType              FluxRangeType;

  // constructor 
  NumFluxBase(const Model& mod) : 
    model_(mod),
    numFlux_( mod )
  {}

  // Return value: maximum wavespeed*length of integrationOuterNormal
  // gLeft,gRight are fluxed * length of integrationOuterNormal
  template< class Intersection, class QuadratureImp >
  inline double 
  numericalFlux( const Intersection& intersection,
                 const EntityType& inside,
                 const EntityType& outside,
                 const double time, 
                 const QuadratureImp& faceQuadInner,
                 const QuadratureImp& faceQuadOuter,
                 const int quadPoint,
                 const RangeType& uLeft, 
                 const RangeType& uRight,
                 RangeType& gLeft,
                 RangeType& gRight) const;

  //! return reference to model 
  const Model& model() const { return model_; }

protected:
  const Model& model_;
  EULERNUMFLUX::EulerFlux<Model,fluxtype> numFlux_;
};

template<class Model, EULERNUMFLUX::EulerFluxType fluxtype>
template<class Intersection, class QuadratureImp>
inline double
NumFluxBase<Model, fluxtype> :: 
numericalFlux( const Intersection& intersection,
               const EntityType& inside,
               const EntityType& outside,
               const double time, 
               const QuadratureImp& faceQuadInner,
               const QuadratureImp& faceQuadOuter,
               const int quadPoint,
               const RangeType& uLeft, 
               const RangeType& uRight,
               RangeType& gLeft,
               RangeType& gRight) const
{
  DomainType normal = intersection.integrationOuterNormal( faceQuadInner.localPoint( quadPoint ) );
  const double len = normal.two_norm();
  normal *= 1./len;

  // for the sake of additional components we put...
  gLeft = 0.;

  double ldt = numFlux_.num_flux((&(uLeft [0])),
                                 (&(uRight[0])),
                                 (&(normal[0])),
                                 (&(gLeft [0])));

  // scaling and conservation
  gLeft *= len;
  gRight = gLeft;

  // return timestep restriction
  return ldt*len;
}


// LLFNumFlux
//-----------

template <class Model>
class LLFNumFlux : public NumFluxBase< Model, EULERNUMFLUX::LLF >
{
  typedef NumFluxBase< Model, EULERNUMFLUX::LLF > BaseType ;
public:
  LLFNumFlux(const Model& mod) : BaseType( mod ) {}
  static std::string name () { return "LLF (Dennis)"; }
};

// HLLNumFlux
//-----------
template <class Model> 
class HLLNumFlux : public NumFluxBase< Model, EULERNUMFLUX::HLL >
{
  typedef NumFluxBase< Model, EULERNUMFLUX::HLL >  BaseType;
public:
  HLLNumFlux(const Model& mod) : BaseType ( mod ) {}
  static std::string name () { return "HLL (Dennis)"; }
};

// HLLCNumFlux
//-----------
template <class Model> 
class HLLCNumFlux : public NumFluxBase< Model, EULERNUMFLUX::HLLC >
{
  typedef NumFluxBase< Model, EULERNUMFLUX::HLLC >  BaseType;
public:
  HLLCNumFlux(const Model& mod) : BaseType ( mod ) {}
  static std::string name () { return "HLLC (Dennis)"; }
};

#endif // file declaration
