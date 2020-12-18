#include <dune/common/fvector.hh>
#include <dune/fem-dg/examples/euler/chorjo.hh>
#include <dune/python/grid/function.hh>

void chorin(const double (&Ulr)[6], double gamma,
            double t,double x,
            double& q_erg,double& u_erg,double& p_erg)
{
  EULERCHORIN::
    lsg(x,t,&q_erg,&u_erg,&p_erg,
        Ulr[0],Ulr[3],Ulr[1],Ulr[4],Ulr[2],Ulr[5],
        gamma);
}
template <class GridView>
auto chorin(const std::vector<double> &UL, const std::vector<double> &UR, double gamma,
            double x0, double t)
{
  return [UL,UR,gamma,t,x0](const auto& en,const auto& x) -> auto
  {
    Dune::FieldVector<double,GridView::dimension+2> res( 0 );
    auto y = en.geometry().global(x);
    double Ulr[6] = { UL[0],UL[1],UL[GridView::dimension+1],
                      UR[0],UR[1],UR[GridView::dimension+1] };
    chorin(Ulr,gamma,t,y[0]-x0,res[0],res[1],res[GridView::dimension+1]);
    return res; // TODO convert to conservative using Model
  };
}
