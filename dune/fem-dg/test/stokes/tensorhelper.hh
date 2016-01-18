#ifndef TENSORHELPER_HH
#define TENSORHELPER_HH

namespace Dune
{
namespace Fem
{

  static int dimRange=2;
  class TensorHelper
  {
  public:
    template<class JacobianInverseType,class JacobianRangeType>
    static inline void makeStressTensor(const JacobianInverseType inv, const JacobianRangeType& sigma,JacobianRangeType& stress)
    {
      JacobianRangeType tau(0.0);

      for(int i = 0; i< dimRange; ++i)
      {
        inv.umv(sigma[i],tau[i]);
      }

     //  stress[0][0]=2*tau[0][0];
  //     stress[0][1]=tau[0][1]+tau[1][0];
  //     stress[1][0]=stress[0][1];
  //     stress[1][1]=2*tau[1][1];
      stress=tau;
    }

    template<class JacobianRangeType>
    static inline double dyadProduct(const JacobianRangeType a,const JacobianRangeType b)
    {
      double sum=0.;
      for(int i=0;i<dimRange;i++)
      {
        for(int j=0;j<dimRange;j++)
        {
          sum+=a[i][j]*b[i][j];
        }
      }
      return sum;
    }

    template<class RangeType,class DomainType,class JacobianRangeType>
    static inline void tensorProduct2d(const RangeType& v,const DomainType& n,JacobianRangeType& ret)
    {
      ret[0][0]= v[0]*n[0];
      ret[0][1]= v[0]*n[1];
      if( JacobianRangeType :: rows > 1 )
      {
        ret[1][0]= v[1]*n[0];
        ret[1][1]= v[1]*n[1];
      }

    }

    /*
      static inline void tensorProduct3d(const RangeType& v,const DomainType& n, GradientRangeType& ret)const
      {
      ret[0]= v[0]*n[0];
      ret[1]= v[0]*n[1];
      ret[2]= v[0]*n[2];
      ret[3]= v[1]*n[0];
      ret[4]= v[1]*n[1];
      ret[5]= v[1]*n[2];
      ret[6]= v[2]*n[0];
      ret[7]= v[2]*n[1];
      ret[8]= v[2]*n[2];
      }
    */
  };

}
}
#endif
