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

namespace Dune
{
namespace Fem
{
  /**
   * \brief Analytical flux for the Euler problem.
   *
   * \ingroups AnalyticalFluxes
   *
   * The Euler problem we want to solve, is given in conservative form by
   * \f{eqnarray*}{
   *    \partial_t u + \nabla\cdot F(u) & = & 0   & \text{in } \Omega\times[t_\text{start},t_\text{end}], \\
   *                        u( 0,\cdot) & = & u_0 & \text{in } \Omega,
   * \f}
   * where the vector of conservative variables has the form
   *
   * \f[ u = \begin{pmatrix} \rho \\ \rho\mathbf{v} \\ \varepsilon \end{pmatrix}, \quad
   *     \rho\mathbf{v} = (\rho v_1,\ldots,\rho v_d)^\top,\quad
   *     \varepsilon = \rho \mathcal{E},\f]
   *
   * plus suitable boundary conditions.
   *
   * The analytical flux function is given (here: 3D case) by
   *
   * \f[ \begin{pmatrix}
   *     \rho v_1             & \rho v_2             & \rho v_3               \\
   *     \rho v_1^2 + p       & \rho v_2 v_1         & \rho v_3 v_1           \\
   *     \rho v_1 v_2         & \rho v_2^2 + p       & \rho v_3 v_2           \\
   *     \rho v_1 v_3         & \rho v_2 v_3         & \rho v_3^2 + p         \\
   *     (\varepsilon + p)v_1 & (\varepsilon + p)v_2 & (\varepsilon + p)v_3
   * \end{pmatrix} \f]
   *
   * Here, \f$ \rho \f$ denotes the density of the fluid,
   * \f$ \mathbf{v} \f$ the velocity of the fluid,
   * \f$ \varepsilon \f$ the internal energy and
   * \f$ \mathcal{E} \f$ the total energy.
   */

  template < int dimDomain, class Field = double >
  class EulerAnalyticalFlux
  {
  public:
    enum { e = dimDomain+1 };
    typedef Field  FieldType;
    typedef Dune :: FieldVector< FieldType, dimDomain> DomainType;

    /**
     * \brief Returns the analytical flux \f$ F \f$ of the Euler equations.
     *
     * \param[in] gamma the adiabatic constant \f$ \gamma \f$ from the Euler equations
     * \param[in] u evaluation of the solution vector \f$ u \f$
     * \param[out] f the analytical flux \f$ F(u) \f$.
     */
    template <class RangeType, class FluxRangeType>
    inline void analyticalFlux(const FieldType gamma,
                               const RangeType& u,
                               FluxRangeType& f) const
    {
      if( dimDomain == 1 )
      {
        assert(u[0]>1e-10);
        const FieldType p = pressure(gamma, u);
        f[0][0] = u[1];
        f[1][0] = u[1]/u[0]*u[1]+p;
        f[e][0] = u[1]/u[0]*(u[e]+p);
        return ;
      }
      else if ( dimDomain == 2 )
      {
        assert( u[0] >= 1e-10 );
        const FieldType v[2] = {u[1]/u[0],u[2]/u[0]};
        const FieldType p = pressure( gamma, u );
        const FieldType ue_p = (u[e]+p);
        f[0][0] = u[1];            f[0][1] = u[2];
        f[1][0] = v[0]*u[1]+p;     f[1][1] = v[1]*u[1];
        f[2][0] = v[0]*u[2];       f[2][1] = v[1]*u[2]+p;
        f[e][0] = v[0]*ue_p;       f[e][1] = v[1]*ue_p;
        return ;
      }
      else if ( dimDomain == 3 )
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
        return;
      }
      else
        std::abort();
    }

    /**
     * \brief Returns the analytical flux \f$ F \f$ of the Euler equations.
     *
     * \param[in] gamma the adiabatic constant \f$ \gamma \f$ from the Euler equations
     * \param[in] u evaluation of the solution vector \f$ u \f$
     * \param[in] du evaluation of jacobian of the solution vector \f$ \nabla u \f$
     * \param[out] f the analytical flux \f$ F(u) \f$.
     */
    template <class RangeType, class FluxRangeType>
    inline  void jacobian(const FieldType gamma,
                          const RangeType& u,
                          const FluxRangeType& du,
                          RangeType& A) const
    {
      if( dimDomain == 2 )
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
      else
      {
        std::abort();
      }
    };

    /**
     * \brief Return the total energy of the Euler equations.
     *
     * The pressure (for an ideal gas) is given by the equation
     *
     * \f[ \mathcal{E} = \varepsilon - \frac{\rho}{2} |\mathbf{v}|^2. \f]
     *
     * \param[in] u evaluation of the solution vector \f$ u \f$
     */
    template <class RangeType>
    inline FieldType rhoeps(const RangeType& u) const
    {
      assert( u[0] >= 1e-10 );
      double v = u[1]*u[1] ;
      for( int i=2; i<dimDomain+1; ++i )
        v += u[i]*u[i];

      return u[e]-0.5*(v)/u[0];
    }

    /**
     * \brief Returns the pressure \f$ p \f$ of the Euler equations.
     *
     * The pressure (for an ideal gas) is given by the equation
     *
     * \f[ p(u) = (\gamma - 1)\left( \varepsilon - \frac{\rho}{2}|\mathbf{v}|^2\right) = (\gamma - 1)\mathcal{E}. \f]
     *
     * \param[in] gamma the adiabatic constant \f$ \gamma \f$ from the Euler equations
     * \param[in] u evaluation of the solution vector \f$ u \f$
     */
    template <class RangeType>
    inline FieldType pressure(const FieldType gamma,const RangeType& u) const
    {
      assert(u[0]>1e-10);
      const FieldType rhoe = rhoeps(u);
      assert( rhoe>1e-10 );
      return (gamma-1.0)*rhoe;
    }

    /**
     * \brief Returns the maximal speed of the Euler equation, i.e. the maximal eigenvalue
     * of the flux jacobian (multiplied with the length of the normal of the intersection).
     *
     * For the Euler equations the maximal eigenvalue is given by the speed of sound, i.e.
     *
     * \f[ c_s (\rho,p) = \sqrt{\gamma\frac{p}{\rho}}. \f]
     */
    template <class RangeType>
    inline FieldType maxWaveSpeed(const FieldType gamma, const DomainType& n, const RangeType& u) const
    {
      assert( u[0] > FieldType( 1e-10 ) );
      // make sure that n is unit normal
      assert( std::abs( n.two_norm() - 1.0 ) < 1e12 );
      FieldType u_normal = u[1] * n[0];
      for (int i=1; i<dimDomain; ++i )
        u_normal += u[i+1] * n[i];
      u_normal /= u[0];

      FieldType p = pressure(gamma,u);
      FieldType c2 = gamma * p/ u[0];
      assert( c2 > FieldType( 1e-10 ) );
      FieldType maxspd = std::abs(double(u_normal)) + std::sqrt(double(c2));
      return maxspd;

    }
  };

} // end namespace Fem
} // end namespace Dune

#endif // file declaration
