#ifndef DUNE_FEM_DG_ANALYTICALEULERFLUXES_HH
#define DUNE_FEM_DG_ANALYTICALEULERFLUXES_HH

// system includes
#include <string>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#ifdef COUNT_FLOPS
#include <dune/fem/misc/double.hh>
#endif


/**
 *  \brief Analytical flux for the euler problem.
 *
 *  \ingroups AnalyticalFluxes
 */
template < int dimDomain >
class EulerAnalyticalFlux
{
 public:
  enum { e = dimDomain+1 };
#ifdef COUNT_FLOPS
  typedef Dune::Fem::Double FieldType;
#else
  typedef double FieldType;
#endif
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
    return (gamma-1.0)*rhoe;
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
  return (gamma-1.0)*rhoe;
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
  return (gamma-1.0)*re;
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
  assert( u[0] > FieldType( 1e-10 ) );
  FieldType u_normal = FieldType(u[1]*n[0]+u[2]*n[1]) / u[0];
  FieldType p = pressure(gamma,u);
  FieldType c2 = gamma * p/ u[0] * n.two_norm2();
  assert( c2 > FieldType( 1e-10 ) );
  FieldType maxspd = std::abs(double(u_normal)) + std::sqrt(double(c2));
  return maxspd;
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

#endif // file declaration
