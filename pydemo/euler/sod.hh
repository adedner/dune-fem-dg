#ifndef DUNE_FEMDG_SOD_HH
#define DUNE_FEMDG_SOD_HH

#include <dune/common/fvector.hh>
#include <dune/fem-dg/examples/euler/chorjo.hh>

void chorin(const double (&Ulr)[6], double gamma,
            double t,double x,
            double& q_erg,double& u_erg,double& p_erg)
{
  EULERCHORIN::
    lsg(x,t,&q_erg,&u_erg,&p_erg,
        Ulr[0],Ulr[3],Ulr[1],Ulr[4],Ulr[2],Ulr[5],
        gamma);
}
template <class DomainType>
Dune::FieldVector<double, DomainType::dimension+2> sod(
     const std::vector<double> &UL,
     const std::vector<double> &UR,
     double gamma,
     const double t, const double x0, const DomainType& x )
{
  Dune::FieldVector<double, DomainType::dimension+2> res( 0 );
  double Ulr[6] = { UL[0],UL[1],UL[DomainType::dimension+1],
                    UR[0],UR[1],UR[DomainType::dimension+1] };
  chorin(Ulr,gamma,t,x[0]-x0,res[0],res[1],res[DomainType::dimension+1]);
  return res;
}

#endif
