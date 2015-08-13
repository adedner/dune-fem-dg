#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cassert>
#include <cmath>

#include <dune/common/array.hh>
#include <dune/fem-dg/models/stokesprobleminterfaces.hh>
#include "corner.hh"

namespace Dune
{


  template< class GridImp>
  class StokesProblemDefault
    : public StokesProblemInterface<Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, GridImp :: dimension > ,
			                        	    Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, 1 >  >
  {
    typedef Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, GridImp :: dimension > FunctionSpaceType ;
    typedef Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, 1 > PressureFunctionSpaceType ;
    typedef StokesProblemInterface<FunctionSpaceType,PressureFunctionSpaceType> BaseType;

  public:

    static const int dimRange  = FunctionSpaceType::dimRange;
    static const int dimDomain = FunctionSpaceType::dimDomain;

    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType  RangeType;
    typedef typename BaseType::PressureRangeType  PressureRangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::DomainFieldType DomainFieldType;
    typedef typename BaseType::RangeFieldType  RangeFieldType;

    typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    typedef typename BaseType::DiffusionMatrixType DiffusionMatrixType;

   //  explicit StokesProblemDefault (  )


    //! the right hand side (i.e., the Laplace of u)
     void f ( const DomainType &x, RangeType &ret ) const
    {

      ret[0]=0;
      ret[1]=0;
#if 0
      ret[0] = sin(x[1]);
      ret[0] *=exp(x[0]);
      ret[0] *=-2.0;

      ret[1]=0;
      ret[1] = cos(x[1]);
      ret[1] *=exp(x[0]);
      ret[1] *=-2.0;
#endif

    }


    //! the exact solution
    void u ( const DomainType &p, RangeType &ret ) const
    {
      double x=p[0];
      double y=p[1];

      //u1
       ret[0]=cos(y);
       ret[0]*=y;
       ret[0]+=sin (y);
       ret[0]*=exp(x);
       ret[0]*=-1.0;

      //u2
       ret[1]=sin(y);
       ret[1]*=y;
       ret[1]*=exp(x);
    }
    //! the exact solution
    void p(const DomainType& x, PressureRangeType& ret) const
    {
      ret[0] = sin(x[1]);
      ret[0] *=exp(x[0]);
      ret[0] *=2.0;
     //ret=0;
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
     grad[0][0]=-exp(x)*(cos(y)*y+sin(y));
     grad[0][1]=-exp(x)*(2*cos(y)-sin(y)*y);
     grad[1][0]=sin(y)*y*exp(x);
     grad[1][1]=exp(x)*(cos(y)*y+sin(y));

    }

  private:

  };




  template< class GridImp>
  class StokesProblemPeriodic
    : public StokesProblemInterface<Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, GridImp :: dimension > ,
				    Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, 1 >  >
  {
    typedef Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, GridImp :: dimension > FunctionSpaceType ;
    typedef Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, 1 > PressureFunctionSpaceType ;
    typedef StokesProblemInterface<FunctionSpaceType,PressureFunctionSpaceType> BaseType;

  public:

    static const int dimRange = FunctionSpaceType::dimRange;
    static const int dimDomain = FunctionSpaceType::dimDomain;

    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType  RangeType;
    typedef typename BaseType::PressureRangeType  PressureRangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::DomainFieldType DomainFieldType;
    typedef typename BaseType::RangeFieldType  RangeFieldType;

    typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    typedef typename BaseType::DiffusionMatrixType DiffusionMatrixType;

   //  explicit StokesProblemDefault (  )


    //! the right hand side (i.e., the Laplace of u)
     void f ( const DomainType &p, RangeType &ret ) const
    {

      double x=p[0];
      double y=p[1];

      ret[0]=(12-24*y)*x*x*x*x;
      ret[0]+=(-24+48*y)*x*x*x;
      ret[0]+=(-48*y+72*y*y-48*y*y*y+12)*x*x;
      ret[0]+=(-2+24*y-72*y*y+48*y*y*y)*x+1-4*y+12*y*y-8*y*y*y;

      ret[1]=(8-48*y+48*y*y)*x*x*x;
      ret[1]+=(-12+72*y-72*y*y)*x*x;
      ret[1]+=(4-24*y+48*y*y-48*y*y*y+24*y*y*y*y)*x-12*y*y+24*y*y*y-12*y*y*y*y;

    }
    //! the exact solution
    void u ( const DomainType &p, RangeType &ret ) const
    {
      double x=p[0];
      double y=p[1];

      //u1
      ret[0]=x*x;
      ret[0]*=(1-x)*(1-x);
      ret[0]*=2*y-6*y*y+4*y*y*y;
			ret[0]+=1;
      //u2
      ret[1]=y*y;
      ret[1]*=-1.;
      ret[1]*=(1-y)*(1-y);
      ret[1]*=2*x-6*x*x+4*x*x*x;
			ret[1]+=1.;
    }

    //! the exact solution
    void p (const DomainType& x, PressureRangeType& ret) const
    {
      ret[0]=x[0]*(1-x[0]);
    }


    //! the diffusion matrix
    void K ( const DomainType &x, DiffusionMatrixType &m ) const
    {
      m = 0;
      for( int i = 0; i < dimDomain; ++i )
        m[ i ][ i ] = 1;
    }

    bool constantK () const
    {
      return true;
    }

    //! the gradient of the exact solution
     void gradient ( const DomainType &x, JacobianRangeType &grad ) const
    {
      grad=0.0;
    }

  private:

  };


  template< class GridImp>
  class GeneralizedStokesProblem
    : public StokesProblemInterface<Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, GridImp :: dimension > ,
				    Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, 1 >  >
  {
    typedef Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, GridImp :: dimension > FunctionSpaceType ;
    typedef Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, 1 > PressureFunctionSpaceType ;
    typedef StokesProblemInterface<FunctionSpaceType,PressureFunctionSpaceType> BaseType;

  public:

    static const int dimRange = FunctionSpaceType::dimRange;
    static const int dimDomain = FunctionSpaceType::dimDomain;

    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType  RangeType;
    typedef typename BaseType::PressureRangeType  PressureRangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::DomainFieldType DomainFieldType;
    typedef typename BaseType::RangeFieldType  RangeFieldType;

    typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    typedef typename BaseType::DiffusionMatrixType DiffusionMatrixType;

    GeneralizedStokesProblem()
      : mu_(Dune::Fem:: Parameter::getValue<double>( "mu", 1.0 ) ),
        alpha_(Dune::Fem:: Parameter::getValue<double>( "alpha", 1.0 ) )
    {}


    //! the right hand side (i.e., the Laplace of u)
     void f ( const DomainType &p, RangeType &ret ) const
    {
      double x=p[0];
      double y=p[1];
      ret[0] = cos(0.5*M_PI*(x+y)) * (-alpha_+0.5*mu_*M_PI*M_PI) + 0.5*M_PI*cos(0.5*M_PI*(x-y) );
      ret[1] = - ret[0];
    }
    //! the exact solution
    void u ( const DomainType &p, RangeType &ret ) const
    {
      double x=p[0];
      double y=p[1];
      //u1
      ret[0] = cos(0.5*M_PI*(x+y));
      ret[1] = -ret[0];
    }

    //! the exact solution
    void p (const DomainType& x, PressureRangeType& ret) const
    {
      ret[0] = sin(0.5*M_PI*(x[0]-x[1]));
    }

    //! the diffusion matrix
    void K ( const DomainType &x, DiffusionMatrixType &m ) const
    {
      m = 0;
      for( int i = 0; i < dimDomain; ++i )
        m[ i ][ i ] = mu_;
    }

    bool constantK () const
    {
      return true;
    }

    //! the gradient of the exact solution
     void gradient ( const DomainType &x, JacobianRangeType &grad ) const
    {
      grad[0][0] = -0.5*M_PI*sin(0.5*M_PI*(x[0]+x[1]));
//grad[0][1] = grad[0][0];
//grad[1][0] = 0.5*M_PI*sin(0.5*M_PI*(x[0]+x[1]));

      grad[1][0] = grad[0][0];
      grad[0][1] = 0.5*M_PI*sin(0.5*M_PI*(x[0]+x[1]));
      grad[1][1] = grad[1][0];
    }

    virtual double gamma() const { return alpha_; }

  private:
    double mu_;
    double alpha_;
  };


  template< class GridImp >
  static StokesProblemInterface<Dune::Fem::FunctionSpace< double, double, GridImp::dimensionworld,GridImp::dimensionworld  >,
    Dune::Fem::FunctionSpace< double, double, GridImp::dimensionworld, 1 > > *
  createProblem()
  {
    std::cout<<"CREATEPROBLEM\n";
    int problemFunction = 2; // default value

    switch( 0 )
      {
      case 0: return new StokesProblemDefault< GridImp >();
      case 1: return new StokesProblemPeriodic<GridImp> ( );
      case 2: return new GeneralizedStokesProblem<GridImp> ( );
      default: std::cerr << "Wrong problem value, bye, bye!" << std::endl;
	abort();
      }
    return 0;
  }

}

#endif // #ifndef PROBLEM_HH
