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

  typedef enum {LLF, HLL, HLL_PT, HLLC} EulerFluxType;
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
    //typedef double FieldType ;
    typedef typename Model :: RangeFieldType  FieldType;
    static const int dim = Model :: dimDomain;

    EulerFlux( const Model& model )
      : model_( model ),
        _gamma( model.gamma() )
    {}

    void flux(const FieldType U[dim+2], FieldType *f[dim]) const
    {
      const FieldType rho = U[0];
      const FieldType *rho_u = &U[1];
      const FieldType E = U[dim+1];

      FieldType u[dim], Ekin2 = 0.0;
      for(int i=0; i<dim; i++){
        u[i] = (1.0/rho) * rho_u[i];
        Ekin2 += rho_u[i] * u[i];
      }

      const FieldType p = (_gamma-1.0)*(E - 0.5*Ekin2);

      for(int i=0; i<dim; i++){
        f[i][0] = rho_u[i];

        for(int j=0; j<dim; j++) f[i][1+j] = rho_u[i] * u[j];
        f[i][1+i] += p;

        f[i][dim+1] = (E+p) * u[i];
      }

    }

    FieldType num_flux(const FieldType Uj[dim+2], const FieldType Un[dim+2],
                       const FieldType normal[dim], FieldType gj[dim+2]) const
    {
      if constexpr (flux_type == LLF) return num_flux_LLF(Uj, Un, normal, gj);

      if constexpr (flux_type == HLL) return num_flux_HLL(Uj, Un, normal, gj, std::false_type());

      // potential temperature
      if constexpr (flux_type == HLL_PT) return num_flux_HLL(Uj, Un, normal, gj, std::true_type());

      if constexpr (flux_type == HLLC) return num_flux_HLLC(Uj, Un, normal, gj);

      DUNE_THROW( Dune::NotImplemented, "Numerical flux not implemented" );
    }

    const Model& model_;
    const FieldType _gamma;

  private:
    FieldType num_flux_LLF(const FieldType Uj[dim+2], const FieldType Un[dim+2],
                           const FieldType normal[dim], FieldType gj[dim+2]) const
    {
      const FieldType rhoj = Uj[0];
      const FieldType *rho_uj = &Uj[1];
      const FieldType Ej = Uj[dim+1];
      const FieldType rhon = Un[0];
      const FieldType *rho_un = &Un[1];
      const FieldType En = Un[dim+1];

      FieldType uj[dim], Ekin2j=0.0, un[dim], Ekin2n=0.0;
      FieldType u_normal_j=0.0, u_normal_n=0.0;
      for(int i=0; i<dim; i++){
        uj[i] = (1.0/rhoj) * rho_uj[i];
        un[i] = (1.0/rhon) * rho_un[i];
        Ekin2j += rho_uj[i] * uj[i];
        Ekin2n += rho_un[i] * un[i];
        u_normal_j += uj[i] * normal[i];
        u_normal_n += un[i] * normal[i];
      }

      const FieldType pj = (_gamma-1.0)*(Ej - 0.5*Ekin2j);
      const FieldType cj = std::sqrt(_gamma*pj/rhoj);
      const FieldType pn = (_gamma-1.0)*(En - 0.5*Ekin2n);
      const FieldType cn = std::sqrt(_gamma*pn/rhon);

      assert(rhoj>0.0 && pj>0.0 && rhoj>0.0 && pj>0.0);

      const FieldType alphaj = std::abs(u_normal_j) + cj;
      const FieldType alphan = std::abs(u_normal_n) + cn;
      const FieldType alpha = (alphaj > alphan)? alphaj : alphan;

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


    template <bool pottemp>
    FieldType num_flux_HLL(const FieldType Uj[dim+2], const FieldType Un[dim+2],
                           const FieldType normal[dim], FieldType gj[dim+2],
                           const std::integral_constant<bool, pottemp> ) const
    {
      static constexpr bool potentialTemperature = pottemp;

      const FieldType rhoj = Uj[0];
      FieldType Ej = Uj[dim+1];
      const FieldType rhon = Un[0];
      FieldType En = Un[dim+1];

      if constexpr (potentialTemperature)
      {
        FieldType pressj, tempj;
        FieldType pressn, tempn;
        model_.pressureTemperature( Uj, pressj, tempj );
        model_.pressureTemperature( Un, pressn, tempn );

        FieldType Ekinj = 0;
        FieldType Ekinn = 0;
        for(int i=1; i<dim+1; i++){
          Ekinj += (0.5/rhoj) * Uj[i] * Uj[i];
          Ekinn += (0.5/rhon) * Un[i] * Un[i];
        }
        Ej = pressj/(_gamma-1.) + Ekinj;
        En = pressn/(_gamma-1.) + Ekinn;
      }

      FieldType rho_uj[dim], rho_un[dim], uj[dim], un[dim];
      FieldType Ekin2j=0.0, Ekin2n=0.0;
      rotate(normal, Uj+1, rho_uj);
      rotate(normal, Un+1, rho_un);
      for(int i=0; i<dim; i++){
        uj[i] = (1.0/rhoj) * rho_uj[i];
        un[i] = (1.0/rhon) * rho_un[i];
        Ekin2j += rho_uj[i] * uj[i];
        Ekin2n += rho_un[i] * un[i];
      }

      const FieldType pj = (_gamma-1.0)*(Ej - 0.5*Ekin2j);
      const FieldType pn = (_gamma-1.0)*(En - 0.5*Ekin2n);

      const FieldType cj = std::sqrt(_gamma*pj/rhoj);
      const FieldType cn = std::sqrt(_gamma*pn/rhon);

      assert(rhoj>0.0 && pj>0.0 && rhoj>0.0 && pj>0.0);

      const FieldType rho_bar = 0.5 * (rhoj + rhon);
      const FieldType c_bar = 0.5 * (cj + cn);
      const FieldType p_star = 0.5 * ( (pj+pn) - (un[0]-uj[0])*rho_bar*c_bar );
      const FieldType u_star = 0.5 * ( (uj[0]+un[0]) - (pn-pj)/(rho_bar*c_bar) );
      const FieldType tmp = 0.5*(_gamma+1.0)/_gamma;
      const FieldType qj = (p_star > pj)? std::sqrt( 1.0 + tmp*(p_star/pj - 1.0) ): 1.0;
      const FieldType qn = (p_star > pn)? std::sqrt( 1.0 + tmp*(p_star/pn - 1.0) ): 1.0;

      const FieldType sj = uj[0] - cj*qj;
      const FieldType sn = un[0] + cn*qn;

      FieldType guj[dim];

      if (u_star > 0.0){
        if (sj >= 0.0){
          gj[0] = rho_uj[0];

          for(int i=0; i<dim; i++) guj[i] = rho_uj[i]*uj[0];
          guj[0] += pj;

          if constexpr (potentialTemperature)
            gj[dim+1] = Uj[dim+1]*uj[0];
          else
            gj[dim+1] = (Ej+pj)*uj[0];
        }
        else{
          const FieldType tmp1 = sj * sn;
          const FieldType tmp2 = 1.0/(sn - sj);
          gj[0] = tmp2 * ( sn*rho_uj[0] - sj*rho_un[0] + tmp1*(rhon - rhoj) );

          for(int i=0; i<dim; i++){
            guj[i] = tmp2*((sn*uj[0]-tmp1)*rho_uj[i] - (sj*un[0]-tmp1)*rho_un[i]);
          }
          guj[0] += tmp2 * (sn*pj - sj*pn);

          if constexpr(potentialTemperature)
          {
            const FieldType Etmpj = Uj[dim+1]*uj[0];
            const FieldType Etmpn = Un[dim+1]*un[0];
            gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(Un[dim+1] - Uj[dim+1]));
          }
          else
          {
            const FieldType Etmpj = (Ej+pj)*uj[0];
            const FieldType Etmpn = (En+pn)*un[0];
            gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(En - Ej));
          }
        }
      }
      else{
        if (sn <= 0.0){
          gj[0] = rho_un[0];

          for(int i=0; i<dim; i++) guj[i] = rho_un[i]*un[0];
          guj[0] += pn;

          if constexpr(potentialTemperature)
            gj[dim+1] = Uj[dim+1]*un[0];
          else
            gj[dim+1] = (En+pn)*un[0];
        }
        else{
          const FieldType tmp1 = sj * sn;
          const FieldType tmp2 = 1.0/(sn - sj);
          gj[0] = tmp2 * ( sn*rho_uj[0] - sj*rho_un[0] + tmp1*(rhon - rhoj) );

          for(int i=0; i<dim; i++){
            guj[i] = tmp2*((sn*uj[0]-tmp1)*rho_uj[i] - (sj*un[0]-tmp1)*rho_un[i]);
          }
          guj[0] += tmp2 * (sn*pj - sj*pn);

          if constexpr(potentialTemperature)
          {
            const FieldType Etmpj = Uj[dim+1]*uj[0];
            const FieldType Etmpn = Un[dim+1]*un[0];
            gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(Un[dim+1] - Uj[dim+1]));
          }
          else
          {
            const FieldType Etmpj = (Ej+pj)*uj[0];
            const FieldType Etmpn = (En+pn)*un[0];
            gj[dim+1] = tmp2 * (sn*Etmpj-sj*Etmpn + tmp1*(En - Ej));
          }
        }
      }

      rotate_inv(normal, guj, gj+1);
      return (std::abs(sj) > std::abs(sn))? std::abs(sj): std::abs(sn);
    }


    FieldType num_flux_HLLC(const FieldType Um[dim+2], const FieldType Up[dim+2],
                            const FieldType normal[dim], FieldType g[dim+2]) const
    {
      const FieldType rhom = Um[0];
      const FieldType rhop = Up[0];
      const FieldType Em = Um[dim+1];
      const FieldType Ep = Up[dim+1];

      FieldType rho_um[dim], rho_up[dim];
      rotate( normal, Um+1, rho_um );
      rotate( normal, Up+1, rho_up );

      FieldType Ekinm = 0.;
      FieldType Ekinp = 0.;
      FieldType um[dim], up[dim];
      for( int i=0; i<dim; ++i )
      {
        um[i] = rho_um[i] / rhom;
        up[i] = rho_up[i] / rhop;
        Ekinm += rho_um[i] * um[i];
        Ekinp += rho_up[i] * up[i];
      }

      const FieldType pm = (_gamma-1.0)*(Em - 0.5*Ekinm);
      const FieldType pp = (_gamma-1.0)*(Ep - 0.5*Ekinp);

      assert( rhom>0.0 && pm>0.0 && rhop>0.0 && pp>0.0 );

      const FieldType cm = std::sqrt(_gamma*pm/rhom);
      const FieldType cp = std::sqrt(_gamma*pp/rhop);

      const FieldType rho_bar = 0.5 * (rhom + rhop);
      const FieldType c_bar = 0.5 * (cm + cp);
      const FieldType p_star = 0.5 * ( (pm+pp) - (up[0]-um[0])*rho_bar*c_bar );
      const FieldType u_star = 0.5 * ( (um[0]+up[0]) - (pp-pm)/(rho_bar*c_bar) );
      const FieldType tmp = 0.5*(_gamma+1.0)/_gamma;
      const FieldType qm = (p_star > pm) ? std::sqrt( 1.0 + tmp*(p_star/pm - 1.0) ) : 1.0;
      const FieldType qp = (p_star > pp) ? std::sqrt( 1.0 + tmp*(p_star/pp - 1.0) ) : 1.0;

      const FieldType sm = um[0] - cm*qm;
      const FieldType sp = up[0] + cp*qp;

      FieldType guj[dim];

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
        const FieldType tmpm = sm*(sm-um[0])/(sm-u_star);
        const FieldType tmpp = sp*(sp-up[0])/(sp-u_star);

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

      return (std::abs(sm) > std::abs(sp)) ? std::abs(sm) : std::abs(sp);
    }



    static void rotate(const FieldType n[dim],
                       const FieldType u[dim], FieldType u_rot[dim])
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
        FieldType d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

        if (d > 1.0e-8) {
          FieldType d_1 = 1.0/d;
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
    }

    static void rotate_inv(const FieldType n[dim],
                           const FieldType u_rot[dim], FieldType u[dim])
    {
      if (dim == 1){
        u[0] = n[0] * u_rot[0];
      }

      if (dim == 2){
        u[0] = n[0]*u_rot[0] - n[1]*u_rot[1];
        u[1] = n[1]*u_rot[0] + n[0]*u_rot[1];
      }

      if (dim == 3){
        FieldType d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

        if (d > 1.0e-8) {
          FieldType d_1 = 1.0/d;
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
  };

}

namespace Dune
{
namespace Fem
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
    : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters > BaseType;

    typedef typename ModelImp::Traits             Traits;
    typedef typename ModelImp::Traits::GridType   GridType;
    static const int dimRange = ModelImp::dimRange;
    static const int dimDomain = GridType::dimensionworld;
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::RangeFieldType     RangeFieldType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;

    typedef FluxImp                               FluxType;

  public:
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    /**
     * \copydoc DGAdvectionFluxBase::DGAdvectionFluxBase()
     */
    EulerFluxImpl (const ModelImp& mod, const ParameterType& parameter = ParameterType() )
      : BaseType( mod ),
        numFlux_( mod )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::numericalFlux()
     */
    template< class LocalEvaluation >
    RangeFieldType numericalFlux( const LocalEvaluation& left,
                                  const LocalEvaluation& right,
                                  const RangeType& uLeft,
                                  const RangeType& uRight,
                                  const JacobianRangeType& jacLeft,
                                  const JacobianRangeType& jacRight,
                                  RangeType& gLeft,
                                  RangeType& gRight) const
    {
      // use RangeFieldType here to avoid conflicts in flux calculation
      Dune::FieldVector< RangeFieldType, dimDomain > normal = left.intersection().integrationOuterNormal( left.localPosition() );
      const RangeFieldType len = normal.two_norm();
      normal *= 1./len;

      // for the sake of additional components we put...
      gLeft = 0.;

      RangeFieldType ldt = numFlux_.num_flux((&(uLeft [0])),
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


} // end namespace Fem
} // end namespace Dune

#endif // file declaration
