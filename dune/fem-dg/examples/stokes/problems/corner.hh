#ifndef CORNER_HH
#define CORNER_HH

#include <cassert>
#include <cmath>

#include <dune/fem-dg/models/stokesprobleminterfaces.hh>

namespace Dune
{
namespace Fem
{
namespace Stokes
{

  /**
   * \brief Interface class for Stokes problem.
   *
   * \ingroup StokesProblems
   */
  template< class GridImp>
  class ProblemCorner
  {
    typedef ProblemInterfaceTraits< GridImp >            Traits;

    typedef typename Traits::FunctionSpaceType         FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType        PoissonProblemBaseType;
    typedef typename Traits::StokesProblemType         StokesProblemBaseType;

  public:

    class PoissonProblem
      : public PoissonProblemBaseType
    {
    public:
      typedef PoissonProblemBaseType BaseType;

      static const int dimRange  = PoissonProblemBaseType::dimRange;
      static const int dimDomain = PoissonProblemBaseType::dimDomain;

      typedef typename PoissonProblemBaseType::DomainType          DomainType;
      typedef typename PoissonProblemBaseType::RangeType           RangeType;
      typedef typename PoissonProblemBaseType::JacobianRangeType   JacobianRangeType;
      typedef typename PoissonProblemBaseType::DomainFieldType     DomainFieldType;
      typedef typename PoissonProblemBaseType::RangeFieldType      RangeFieldType;

      typedef typename PoissonProblemBaseType::DiffusionMatrixType DiffusionMatrixType;

      //! the right hand side (i.e., the Laplace of u)
      void f ( const DomainType &x, RangeType &ret ) const
      {
        ret=0;
      }

      //! the exact solution
      void u ( const DomainType &x, RangeType &ret ) const
      {
        const double lambda=0.54448373678246;
        double r=sqrt(x[0]*x[0]+x[1]*x[1]);

        double phi=arg(std::complex<double>(x[0],x[1]));
        if (x[1]<0) phi+=2.*M_PI;

		  	ret[0]= pow(r,1.0*lambda)
		  		*(cos(phi)*
		  			(cos((1.0+lambda)*phi)
		  			 *cos(3.0/2.0*lambda*0.3141592653589793E1)
		  			 +sin((1.0+lambda)*phi)*(1.0+lambda)
		  			 -cos((1.0-lambda)*phi)*
		  			 cos(3.0/2.0*lambda*0.3141592653589793E1)
		  			 -sin((1.0-lambda)*phi)*(1.0-lambda))
		  			+(1.0+lambda)*sin(phi)
		  			*(sin((1.0+lambda)*phi)*cos(3.0/2.0*lambda*0.3141592653589793E1)/(1.0+lambda)
		  				-cos((1.0+lambda)*phi)-sin((1.0-lambda)*phi)
		  				*cos(3.0/2.0*lambda*0.3141592653589793E1)/(1.0-lambda)
		  				+cos((1.0-lambda)*phi)));

		  	ret[1]=pow(r,1.0*lambda)
		  		*(sin(phi)*
		  			(cos((1.0+lambda)*phi)
		  			 *cos(3.0/2.0*lambda*0.3141592653589793E1)
		  			 +sin((1.0+lambda)*phi)*(1.0+lambda)-cos((1.0-lambda)*phi)*
		  			 cos(3.0/2.0*lambda*0.3141592653589793E1)
		  			 -sin((1.0-lambda)*phi)*(1.0-lambda))
		  			-(1.0+lambda)*cos(phi)
		  			*(sin((1.0+lambda)*phi)*cos(3.0/2.0*lambda* 0.3141592653589793E1)/(1.0+lambda)
		  				-cos((1.0+lambda)*phi)-sin((1.0-lambda)*phi)
		  				*	cos(3.0/2.0*lambda*0.3141592653589793E1)/(1.0-lambda)+cos((1.0-lambda)*phi)));

      }

      //! the diffusion matrix
      void K ( const DomainType &x, DiffusionMatrixType &m ) const
      {
        m = 0;
        for( int i = 0; i < dimDomain; ++i )
          m[ i ][ i ] = 1.;
      }

      bool constantK () const
      {
        return true;
      }

      //! the gradient of the exact solution
      void gradient ( const DomainType &p, JacobianRangeType &grad ) const
      {
	  		double x=p[0];
	  		double y=p[1];
	  		double ur,uphi,vr,vphi,rx,ry,phix,phiy;
	  		double r=sqrt(x*x+y*y);
	  		double phi=arg(std::complex<double>(x,y));
        if (y<0)
          phi+=2.*M_PI;

	  		ur= dru(r,phi);
	  		uphi=dphiu(r,phi);
	  		vr= drv(r,phi);
	  		vphi= dphiv(r,phi);
	  		rx=dxr(x,y);
	  		ry=dyr(x,y);
	  		phix=dxphi(x,y);
	  		phiy=dyphi(x,y);

	  		grad[0][0]=ur*rx+uphi*phix;
	  		grad[0][1]=ur*ry+uphi*phiy;
	  		grad[1][0]=vr*rx+vphi*phix;
	  		grad[1][1]=vr*ry+vphi*phiy;
	  	}

    private:
      inline double dru(double r, double phi)const
  		{
  			return -0.1195311300e1 * pow(r, -0.4555162632e0) * (-0.4555162632e0 * cos(phi) * cos(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.7035374604e0 * cos(phi) * sin(0.1544483737e1 * phi) + 0.4555162632e0 * cos(phi) * cos(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) + 0.2074950656e0 * cos(phi) * sin(0.4555162632e0 * phi) - 0.4555162632e0 * sin(phi) * sin(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) + 0.7035374604e0 * sin(phi) * cos(0.1544483737e1 * phi) + 0.1544483737e1 * sin(phi) * sin(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.7035374604e0 * sin(phi) * cos(0.4555162632e0 * phi));


  		}
  		inline double dphiu(double r, double phi)const
  		{
  			return  -0.1195311300e1 * pow(r, 0.544483736782463929140876854601e0) * (-0.1118527593e1 * cos(phi) * cos(0.4555162632e0 * phi) + 0.207495066e0 * sin(phi) * sin(0.4555162632e0 * phi) - 0.7035374604e0 * cos(phi) * cos(0.1544483737e1 * phi) - 0.4555162632e0 * sin(phi) * cos(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) + 0.4555162632e0 * sin(phi) * cos(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) + 0.4555162632e0 * cos(phi) * sin(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) + 0.2455516263e1 * cos(phi) * sin(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.7035374604e0 * sin(phi) * sin(0.1544483737e1 * phi));


  		}
  		inline double drv(double r, double phi)const
  		{
  			return -0.1195311300e1 * pow(r, -0.4555162632e0) * (-0.4555162632e0 * sin(phi) * cos(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.7035374604e0 * sin(phi) * sin(0.1544483737e1 * phi) + 0.4555162632e0 * sin(phi) * cos(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) + 0.2074950656e0 * sin(phi) * sin(0.4555162632e0 * phi) + 0.4555162632e0 * cos(phi) * sin(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.7035374604e0 * cos(phi) * cos(0.1544483737e1 * phi) - 0.1544483737e1 * cos(phi) * sin(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) + 0.7035374604e0 * cos(phi) * cos(0.4555162632e0 * phi));
  		}

      inline double dphiv(double r, double phi) const
  		{
  			return	0.1195311300e1 * pow(r, 0.544483736782463929140876854601e0) * (0.1118527593e1 * sin(phi) * cos(0.4555162632e0 * phi) + 0.2074950656e0 * cos(phi) * sin(0.4555162632e0 * phi) + 0.4555162632e0 * cos(phi) * cos(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.2455516263e1 * sin(phi) * sin(0.4555162632e0 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.4555162632e0 * sin(phi) * sin(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.4555162632e0 * cos(phi) * cos(0.1544483737e1 * phi) * cos(0.8167256052e0 * 0.3141592654e1) - 0.7035374604e0 * cos(phi) * sin(0.1544483737e1 * phi) + 0.7035374604e0 * sin(phi) * cos(0.1544483737e1 * phi));
  		}

  		inline double dxr(double x, double y)const
  		{
  			return x * pow(x * x + y * y, -0.1e1 / 0.2e1);
  		}

  		inline double dxphi(double x, double y)const
  		{
  			double res= (pow(x * x + y * y, -0.1e1 / 0.2e1) - x * x * pow(x * x + y * y, -0.3e1 / 0.2e1))
  					* pow(0.1e1 - x * x / (x * x + y * y), -0.1e1 / 0.2e1);

  			if (y<0)
  				return res;
  			else
  				return -1.*res;
  		}

      inline double dyr(double x, double y)const
  		{
  			return pow(x * x + y * y, -0.1e1 / 0.2e1) * y;
  		}

  		inline double dyphi(double x, double y)const
  		{
  			double res=x * pow(x * x + y * y, -0.3e1 / 0.2e1) * y * pow(0.1e1 - x * x / (x * x + y * y), -0.1e1 / 0.2e1);

  			if (y<0)
  				return -1*res;
  			else
  				return res;
  		}
    };

    class StokesProblem
      : public StokesProblemBaseType
    {
    public:
      typedef StokesProblemBaseType BaseType;

      static const int dimRange  = StokesProblemBaseType::dimRange;
      static const int dimDomain = StokesProblemBaseType::dimDomain;

      typedef typename StokesProblemBaseType::DomainType        DomainType;
      typedef typename StokesProblemBaseType::RangeType         RangeType;
      typedef typename StokesProblemBaseType::JacobianRangeType JacobianRangeType;
      typedef typename StokesProblemBaseType::DomainFieldType   DomainFieldType;
      typedef typename StokesProblemBaseType::RangeFieldType    RangeFieldType;

      typedef typename StokesProblemBaseType::DiffusionMatrixType DiffusionMatrixType;

      //! the exact solution
      void u(const DomainType& x, RangeType& ret) const
      {
        const double lambda=0.54448373678246;
        double r=sqrt(x[0]*x[0]+x[1]*x[1]);
        double phi=arg(std::complex<double>(x[0],x[1]));
        if (x[1]<0) phi+=2.*M_PI;

  			ret[0]=-pow(r,1.0*lambda-1.0)*
  				(pow(1.0+lambda,2.0)
  				 *(cos((1.0+lambda)*phi)
  					 *cos(3.0/2.0*lambda*0.3141592653589793E1)
  					 +sin((1.0+lambda)*phi)*(1.0+lambda)
  					 -cos((1.0-lambda)*phi)*cos(3.0/2.0*lambda*0.3141592653589793E1)
  					 -sin((1.0-lambda)*phi)*(1.0-lambda))-cos((1.0+lambda)*phi)
  				 *pow(1.0+lambda,2.0)*cos(3.0/2.0*lambda*0.3141592653589793E1)
  				 -sin((1.0+lambda)*phi)*pow(1.0+lambda,3.0)
  				 +cos((1.0-lambda)*phi)*pow(1.0-lambda,2.0)
  				 *cos(3.0/2.0*lambda*0.3141592653589793E1)+sin((1.0-lambda)*phi)*pow(1.0-lambda,3.0))/(1.0-lambda);
      }
    };


    typedef PoissonProblem PoissonProblemType;
    typedef StokesProblem StokesProblemType;

  };


}
}
}
#endif // #ifndef PROBLEM_HH
