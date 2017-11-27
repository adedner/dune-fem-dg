#ifndef PROBLEM_HH
#define PROBLEM_HH

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
  class ProblemDefault
  {
    typedef ProblemInterfaceTraits< GridImp >            Traits;

    typedef typename Traits::FunctionSpaceType           FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType   PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType          PoissonProblemBaseType;
    typedef typename Traits::StokesProblemType           StokesProblemBaseType;

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

      //! the diffusion matrix
      void K( const DomainType &x, DiffusionMatrixType &m ) const
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
        ret[0] = sin(x[1]);
        ret[0] *=exp(x[0]);
        ret[0] *=2.0;
       //ret=0;
      }
    };

    typedef PoissonProblem PoissonProblemType;
    typedef StokesProblem StokesProblemType;

    ProblemDefault()
    {}

  };

  /**
   * \brief Stokes problem.
   *
   * \ingroup StokesProblems
   */
  template< class GridImp>
  class ProblemPeriodic
  {
    typedef ProblemInterfaceTraits< GridImp >            Traits;

    typedef typename Traits::FunctionSpaceType           FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType   PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType          PoissonProblemBaseType;
    typedef typename Traits::StokesProblemType           StokesProblemBaseType;

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

      PoissonProblem()
        : mu_(Dune::Fem:: Parameter::getValue<double>( "mu", 1.0 ) ),
          alpha_(Dune::Fem:: Parameter::getValue<double>( "alpha", 1.0 ) )
      {}

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

      //! the diffusion matrix
      void K( const DomainType &x, DiffusionMatrixType &m ) const
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
      void gradient ( const DomainType &p, JacobianRangeType &grad ) const
      {
        grad=0.0;
      }

      virtual double gamma() const { return alpha_; }

    private:
      double mu_;
      double alpha_;
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
        ret[0]=x[0]*(1-x[0]);
      }

    };

    typedef PoissonProblem PoissonProblemType;
    typedef StokesProblem StokesProblemType;

    ProblemPeriodic()
    {}

  };


  /**
   * \brief Stokes problem.
   *
   * \ingroup StokesProblems
   */
  template< class GridImp>
  class GeneralizedProblem
  {
    typedef ProblemInterfaceTraits< GridImp >            Traits;

    typedef typename Traits::FunctionSpaceType           FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType   PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType          PoissonProblemBaseType;
    typedef typename Traits::StokesProblemType           StokesProblemBaseType;
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

      PoissonProblem()
        : mu_(Dune::Fem:: Parameter::getValue<double>( "mu", 1.0 ) ),
          alpha_(Dune::Fem:: Parameter::getValue<double>( "alpha", 1.0 ) )
      {}

      //! the right hand side (i.e., the Laplace of u)
      void f ( const DomainType &p, RangeType &ret ) const
      {
        double x=p[0];
        double y=p[1];
        ret[0] = cos(2.0*M_PI*(x+y)) * (-alpha_ + 4.0 * 2.0*mu_*M_PI*M_PI) + 2.0*M_PI*cos(2.0*M_PI*(x-y) );
        ret[1] = - ret[0];
      }


      //! the exact solution
      void u ( const DomainType &p, RangeType &ret ) const
      {
        double x=p[0];
        double y=p[1];
        //u1
        ret[0] = cos(2.0*M_PI*(x+y));
        ret[1] = -ret[0];
      }

      //! the diffusion matrix
      void K( const DomainType &x, DiffusionMatrixType &m ) const
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
      void gradient ( const DomainType &p, JacobianRangeType &grad ) const
      {
        double x=p[0];
        double y=p[1];
        grad[0] = 1.0;
        grad[1] = -1.0;
        grad *= -2.0 * M_PI * sin( 2.0 * M_PI * ( x+y ) );
      }

      //! the Dirichlet boundary data function
      virtual void g(const DomainType& x, RangeType& ret) const
      {
        u( x, ret );
      }

      virtual std::string name() const
      {
        return "General Poisson";
      }
      virtual double gamma() const { return alpha_; }

    private:
      double mu_;
      double alpha_;
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
        ret[0] = sin(2.0*M_PI*(x[0]-x[1]));
      }

      //! the Dirichlet boundary data function
      virtual void g(const DomainType& x, RangeType& ret) const
      {
        u( x, ret );
      }

      virtual std::string name() const
      {
        return "General Stokes";
      }

    };

    typedef PoissonProblem PoissonProblemType;
    typedef StokesProblem StokesProblemType;

  };


  /**
   * \brief Stokes problem.
   *
   * \ingroup StokesProblems
   */
  template< class GridImp>
  class DrivenCavityProblem
  {
    typedef ProblemInterfaceTraits< GridImp >            Traits;

    typedef typename Traits::FunctionSpaceType           FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType   PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType          PoissonProblemBaseType;
    typedef typename Traits::StokesProblemType           StokesProblemBaseType;
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

      PoissonProblem()
        : mu_(Dune::Fem:: Parameter::getValue<double>( "mu", 1.0 ) ),
          alpha_(Dune::Fem:: Parameter::getValue<double>( "alpha", 1.0 ) )
      {}

      //! the right hand side (i.e., the Laplace of u)
      void f ( const DomainType &p, RangeType &ret ) const
      {
        ret = 0;
      }


      //! the exact solution
      void u ( const DomainType &p, RangeType &ret ) const
      {
        ret = 0;
      }

      //! the diffusion matrix
      void K( const DomainType &x, DiffusionMatrixType &m ) const
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
      void gradient ( const DomainType &p, JacobianRangeType &grad ) const
      {
        grad = 0;
      }

      //! the Dirichlet boundary data function
      virtual void g(const DomainType& x, RangeType& ret) const
      {
        ret = 0;
        ret[ 0 ] = ( x[ 1 ] == 1.0 ) * 2.0;
      }

      virtual std::string name() const
      {
        return "Driven Cavity Problem";
      }

      virtual double gamma() const { return alpha_; }

    private:
      double mu_;
      double alpha_;
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
        ret = 0;
      }

      //! the Dirichlet boundary data function
      virtual void g(const DomainType& x, RangeType& ret) const
      {
        ret = 0.;
      }

      virtual std::string name() const
      {
        return "General Stokes";
      }

    };

    typedef PoissonProblem PoissonProblemType;
    typedef StokesProblem StokesProblemType;

  };

 // template< class GridImp >
 // static StokesProblemInterface<Dune::Fem::FunctionSpace< double, double, GridImp::dimensionworld,GridImp::dimensionworld  >,
 //   Dune::Fem::FunctionSpace< double, double, GridImp::dimensionworld, 1 > > *
 // createProblem()
 // {
 //   std::cout<<"CREATEPROBLEM\n";
 //   int problemFunction = 2; // default value

 //   switch( 0 )
 //     {
 //     case 0: return new StokesProblemDefault< GridImp >();
 //     case 1: return new StokesProblemPeriodic<GridImp> ();
 //     case 2: return new GeneralizedStokesProblem<GridImp> ();
 //     default: std::cerr << "Wrong problem value, bye, bye!" << std::endl;
 // abort();
 //     }
 //   return 0;
 // }
}
}
}

#endif // #ifndef PROBLEM_HH
