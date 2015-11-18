#ifndef DUNE_FEM_DG_EULERFLUXES_HH
#define DUNE_FEM_DG_EULERFLUXES_HH

// system includes
#include <string>
#include <cmath>

#include "parameters.hh"
#include "../advection/fluxbase.hh"

// dune-grid includes
#if WELLBALANCE
#include <dune/grid/common/genericreferenceelements.hh>
#endif

#include <dune/fem-dg/operator/fluxes/rotator.hh>

// local includes
namespace EulerNumFlux
{
  ////////////////////////////////////////////////////////
  //
  // Implementation of Euler fluxes from Dennis Diehl.
  //
  ////////////////////////////////////////////////////////

  typedef enum {LLF, HLL, HLLC} EulerFluxType;
  const EulerFluxType std_flux_type = HLL;

  /**
  *  \brief Advection flux for the euler problem.
  *
  *  \ingroups AdvectionFluxes
  */
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

namespace Dune
{
namespace Fem
{
namespace Euler
{

  ////////////////////////////////////////////////////////
  //
  // Implementation of Euler fluxes from Dennis Diehl.
  //
  ////////////////////////////////////////////////////////

  /**
    *  \brief Advection flux base for the euler problem.
    *
    *  \ingroups AdvectionFluxes
    */
  template <class ModelImp,
            class FluxImp>
  class EulerFluxImpl
    : public DGAdvectionFluxBase< ModelImp, EulerFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, EulerFluxParameters > BaseType;

    typedef typename ModelImp::Traits             Traits;
    typedef typename ModelImp::Traits::GridType   GridType;
    enum { dimRange = ModelImp::dimRange };
    enum { dimDomain = GridType::dimensionworld };
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;
    typedef typename ModelImp::EntityType         EntityType;
    typedef typename ModelImp::IntersectionType   IntersectionType;

    typedef FluxImp                                     FluxType;

  public:
    typedef typename BaseType::MethodType         MethodType;
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    // constructor
    EulerFluxImpl (const ModelImp& mod,
                   const typename MethodType::id mthd,
                   const ParameterType& parameters )
      : BaseType( mod, mthd, parameters ),
        numFlux_( mod )
    {}

    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
    template< class LocalEvaluation >
    double numericalFlux( const LocalEvaluation& left,
                          const LocalEvaluation& right,
                          const RangeType& uLeft,
                          const RangeType& uRight,
                          const JacobianRangeType& jacLeft,
                          const JacobianRangeType& jacRight,
                          RangeType& gLeft,
                          RangeType& gRight) const
    {
      DomainType normal = left.intersection().integrationOuterNormal( left.localPoint() );
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

  protected:
    FluxType numFlux_;
  };



}
}
}

#endif // file declaration
