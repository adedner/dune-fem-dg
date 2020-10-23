#ifndef FVCA6_BENCHMARKPROBLEMS_HH
#define FVCA6_BENCHMARKPROBLEMS_HH

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <set>
#include <cassert>
#include <complex>

namespace Dune
{
namespace Fem
{
namespace Poisson
{

  /**
   * \brief Data function interface for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template<int dim, class DomainField, class Field>
  class DataFunctionIF
  {
  protected:
    DataFunctionIF() {}
  public:
    // destructor
    virtual ~DataFunctionIF() {}

    // returns true if K is constant on one element
    virtual bool constantLocalK () const { return true; }

    virtual double advection () const { return 0.0; }

    // diffusion tensor
    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const = 0;
    // right hand side
    virtual Field rhs  (const DomainField x[dim]) const = 0;
    // right hand side
    virtual Field rhs  (const double time, const DomainField x[dim]) const
    {
      return rhs( x );
    }
    virtual void velocity(const DomainField x[dim], DomainField v[dim]) const
    {
      for(int i=0; i<dim; ++i)
        v[i] = 0.;
    }

    // exact solution
    virtual Field exact(const DomainField x[dim]) const = 0;

    // exact solution (time dependent)
    virtual Field exact(const double time, const DomainField x[dim]) const
    {
      return exact( x );
    }

    // exact gradient
    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const = 0;

    // boundary data
    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const = 0;

    // neumann boundary
    virtual void neumann(const DomainField x[dim], Field grad[dim]) const
    {
      Field tmp[dim];
      gradExact(x,tmp);
      Field k[dim][dim];
      K(x, k);
      for(int i=0; i<dim; ++i)
      {
        grad[i] = 0;
        for(int j=0; j<dim; ++j)
          grad[i] += tmp[j] * k[i][j];
      }
    }

    // boundary data
    virtual bool boundaryDataFunction(const double time,
                                      const DomainField x[dim],
                                      Field & val) const
    {
      return boundaryDataFunction(x, val);
    }

    // neumann boundary
    virtual void neumann(const double time,
                         const DomainField x[dim], Field grad[dim]) const
    {
      return neumann(x, grad);
    }

    virtual int getBoundaryId(const DomainField x[dim]) const
    {
      for(int i=0; i<dim; ++i)
      {
        const int id = 2 * i;
        if( x[i] <= 1e-10 )
        {
          return id;
        }
        else if ( x[i]  >= 0.99999999 )
        {
          return id + 1;
        }
      }

      //assert( false );
      //abort();
      return -1;
    }
    virtual int getBoundaryId(const DomainField x[dim],
                              const DomainField n[dim]) const
    {
      return getBoundaryId( x );
    }

    Field SQR( const Field& a ) const
    {
      return (a * a);
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class AlbertaProblem : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~AlbertaProblem() {}
    AlbertaProblem(Field globalShift, Field factor)
      : globalShift_(0.0)
      , factor_(1.0)
    {
      //assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j) k[i][j] = 0;
        k[i][i] = factor_;
      }
    }

    inline Field xSqr(const DomainField x[dim]) const
    {
      Field xsqr = 0.0;
      for(int i=0; i<dim; ++i) xsqr += x[i] * x[i];
      return xsqr;
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      const Field xsqr = xSqr( x );
      return -(400.0 * xsqr - 20.0 * dim) *  std :: exp( -10.0 * xsqr );
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return std :: exp( -10.0 * xSqr( x ) );
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      const Field factor = -20.0 * std :: exp( -10.0 * xSqr( x ) );
      for(int i=0; i<dim; ++i)
        grad[i] = x[i] * factor ;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      //if(x[0] <= 0.0) return false;
      //if(x[1] <= 0.0) return false;
      return true;
      //return false;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class SinSinSin : public DataFunctionIF<dim,DomainField,Field>
  {
    Field K_[dim][dim];
    const double advection_;
  public:
    virtual ~SinSinSin() {}
    SinSinSin(Field globalShift, Field factor, const double advection = 0.0 )
      : advection_( advection )
    {
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
        {
          if( i == j )
          {
            if( i == 0 )
              K_[i][j] = 10.;//factor;
            else
              K_[i][j] = 0.1;//factor;
          }
          /*
          else if( std::abs( i - j ) == 1 )
          {
            K_[i][j] = 0.5;
          }
          */
          else K_[i][j] = 0;
        }
      }

      /*
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
        {
          std::cout << K_[i][j] << " ";
        }
        std::cout << std::endl;
      }
      */
    }

    virtual double advection() const { return advection_; }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
        {
          k[i][j] = K_[i][j];
        }
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      Field sum = 0;
      int comp[ dim - 1 ] ;
      for(int i=0; i<dim; ++i)
      {
        comp[0] = i;
        for(int j=0; j<dim; ++j)
        {
          comp[ dim - 2 ] = j;
          sum -= K_[j][i] * laplace( x, comp );
        }
      }
      if( std::abs( advection_ ) > 0 )
      {
        for(int i=0; i<dim; ++i)
        {
          sum += advection_ * gradient( x, i );
        }
      }
      return sum;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      const Field pi = 2.0 * M_PI;
      Field val = 1.0;
      for(int i=0; i<dim; ++i)
      {
        val *= std::sin( pi * x[i] );
      }
      return val;
    }

    double laplace(const DomainField x[dim], const int comp[dim - 1] ) const
    {
      const Field pi = 2.0 * M_PI;
      Field val = pi * pi;

      std::set<int> comps ;
      for( int j=0; j<dim-1; ++j)
      {
        comps.insert( comp[j] );
      }

      if( comps.size() == 1 )
      {
        // add other components
        for(int i=0; i<dim; ++i)
        {
          val *= std::sin( pi * x[ i ] );
        }

        // minus because sin'' = -sin
        return -val;
      }
      else
      {
        for( int i=0; i<dim-1; ++i)
        {
          val *= std::cos( pi * x[ comp[i] ] );
        }

        for(int i=0; i<dim; ++i)
        {
          if( comps.find( i ) == comps.end() )
            val *= std::sin( pi * x[ i ] );
        }
        return val;
      }
    }

    double gradient(const DomainField x[dim], const int comp) const
    {
      const Field pi = 2.0 * M_PI;
      Field val = pi * std::cos( pi * x[ comp ] );
      // add other components
      for(int j=1; j<dim; ++j)
      {
        val *= std::sin( pi * x[ (comp + j) % dim ] );
      }
      return val;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        grad[i] = gradient( x, i );
      }
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      //if(x[0] <= 0.0) return false;
      return true;
      //return false;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class SinSin : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~SinSin() {}
    SinSin(Field globalShift, Field factor)
      : globalShift_(globalShift)
      , factor_(factor)
    {
      //assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      Field sin_x = std::sin(2.0*M_PI*x[0]);
      Field sin_y = std::sin(2.0*M_PI*x[1]);

      Field val = 8.0 * M_PI * M_PI * sin_x * sin_y ;
      val *= factor_;
      return val;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      Field val = std::sin(2.0*M_PI*x[0]) * std::sin(2.0*M_PI*x[1]);
      val += globalShift_;
      return val;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = 2.0*M_PI*std::cos(2.0*M_PI*x[0])*std::sin(2.0*M_PI*x[1]);
      grad[1] = 2.0*M_PI*std::sin(2.0*M_PI*x[0])*std::cos(2.0*M_PI*x[1]);

      // initial grad with zero for 3d
      for( int i=2; i<dim; ++i) grad[i] = 0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      //if(x[0] <= 0.0) return false;
      return true;
      //return false;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class CosCos : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~CosCos() {}
    CosCos(Field globalShift, Field factor)
      : globalShift_(globalShift)
      , factor_(factor)
    {
      //assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }


    virtual Field rhs  (const DomainField x[dim]) const
    {
      Field cos_x = std::cos(2.0*M_PI*x[0]);
      Field cos_y = std::cos(2.0*M_PI*x[1]);

      Field val = 8.0 * M_PI*M_PI* cos_x * cos_y ;
      val *= factor_;
      return val;
      //return 0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      Field val = std::cos(2.0*M_PI*x[0]) * std::cos(2.0*M_PI*x[1]);
      val += globalShift_;
      return val;
      //return x[1];
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = 2.0*M_PI*-std::sin(2.0*M_PI*x[0])*std::cos(2.0*M_PI*x[1]);
      grad[1] = 2.0*M_PI*std::cos(2.0*M_PI*x[0])*-std::sin(2.0*M_PI*x[1]);
      //grad[0] = 0.;
      //grad[1] = 1.;

      // initial grad with zero for 3d
      for( int i=2; i<dim; ++i) grad[i] = 0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      //if(x[0] <= 0.0) return false;
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * Problem from Castillo paper.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class CastilloProblem : public DataFunctionIF<dim,DomainField,Field>
  {
    using DataFunctionIF<dim,DomainField,Field> :: SQR;
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~CastilloProblem() {}
    CastilloProblem(Field globalShift, Field factor)
      : globalShift_(globalShift)
      , factor_(factor)
    {
      assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }


    virtual Field rhs  (const DomainField x[dim]) const
    {
      Field ret = 0.0;
      Field tmp =-23+7*SQR(x[1])-24*x[0]+24*x[0]*SQR(x[1])+7*SQR(x[0])+9*SQR(x[0])*SQR(x[1])-24
                  *x[1]+24*x[1]*SQR(x[0]);

      ret=-0.5*exp(0.75*(x[0]+x[1]));
      ret *= tmp;
      ret *= factor_;
      return ret;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      Field val = 4.0 *(1.-SQR(x[0]))*(1.-SQR(x[1]))*exp(0.75*(x[0]+x[1]));
      val += globalShift_;
      return val;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      // not implemented yet
      grad[0] = grad[1] = 0.0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

  /**
   * \brief Data functions "Reentrant corner" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <class DomainField, class Field>
  class ReentrantCorner2d : public DataFunctionIF<2,DomainField,Field>
  {
    enum { dim = 2 };
    const Field globalShift_;
    const Field factor_;
    const Field lambda_;
  public:
    virtual ~ReentrantCorner2d() {}
    ReentrantCorner2d(Field globalShift, Field factor)
      : globalShift_(0.0)
      , factor_(1.0)
      , lambda_(180./270.)
    {
      assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      return 0.0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      double r2 = radius(x[0],x[1]);
      double phi = argphi(x[0],x[1]);
      return pow(r2,lambda_*0.5)*std::sin(lambda_*phi);
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      double r2=radius(x[0],x[1]);
      double phi=argphi(x[0],x[1]);
      double r2dx=2.*x[0];
      double r2dy=2.*x[1];
      double phidx=-x[1]/r2;
      double phidy=x[0]/r2;
      double lambdaPow = (lambda_*0.5)*pow(r2,lambda_*0.5-1.);
      grad[0]= lambdaPow * r2dx * std::sin(lambda_ * phi)
               + lambda_ * std::cos( lambda_ * phi) * phidx * pow(r2,lambda_ * 0.5);
      grad[1]= lambdaPow * r2dy * std::sin(lambda_ * phi)
               + lambda_ * std::cos( lambda_ * phi) * phidy * pow(r2,lambda_*0.5);
      assert( grad[0] == grad[0] );
      assert( grad[1] == grad[1] );
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }

    private:
    /** \brief proper implementation of atan(x,y)
     */
    inline double argphi(double x,double y) const
    {
      double phi=arg(std::complex<double>(x,y));
      if (y<0) phi+=2.*M_PI;
      return phi;
    }

    /** \brief implementation for the radius squared
     * (^0.5 is part of function implementation)
     */
    inline double radius(double x, double y) const
    {
      double ret =0;
      ret = x*x +y*y;
      return ret;
    }
  };

  /**
   * \brief Data functions "Reentrant corner 3D" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class ReentrantCorner : public DataFunctionIF<dim,DomainField,Field>
  {
    typedef ReentrantCorner2d<DomainField, Field> Corner2dType ;
    Corner2dType corner2d_;
    const Field factor_;

  public:
    virtual ~ReentrantCorner() {}
    ReentrantCorner(Field globalShift, Field factor)
      : corner2d_(globalShift, factor),
        factor_( factor )
    {
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      return 0.0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      DomainField x2d[ 2 ];
      if ( dim == 2 )
      {
        for( int i=0; i<dim; ++i ) x2d[ i ] = x[ i ];
        return corner2d_.exact( x2d );
      }
      else
      {
        double sum = 0;
        // ( x, y)
        x2d[ 0 ] = -x[ 0 ]; x2d[ 1 ] =  x[ 1 ];
        sum += corner2d_.exact( x2d );
        // (-y,-z)
        //x2d[ 0 ] = -x[ 2 ]; x2d[ 1 ] = -x[ 1 ];
        x2d[ 0 ] = -x[ 1 ]; x2d[ 1 ] = -x[ 2 ];
        sum += corner2d_.exact( x2d );
        // ( z,-x)
        x2d[ 0 ] = -x[ 0 ]; x2d[ 1 ] = -x[ 2 ];
        sum += corner2d_.exact( x2d );

        //sum *= (x[ 0 ] + 1.0) / 2.0;
        return sum ;
      }
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      for( int i=0; i<dim; ++i ) grad[ i ] = 0;
      if ( dim == 2 )
      {
        DomainField x2d[ 2 ];
        for( int i=0; i<dim; ++i ) x2d[ i ] = x[ i ];
        corner2d_.gradExact( x2d, &grad[ 0 ] );
      }
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

  /**
   * \brief Data functions "Fichera corner" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class FicheraCorner : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;

  public:
    virtual ~FicheraCorner() {}
    FicheraCorner(Field globalShift, Field factor)
      : globalShift_( globalShift ),
        factor_( factor )
    {
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      const double u = exact( x );
      return -3.0/(4.0* u * u * u );
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      double xAbs = 0;
      for( int i=0; i<dim; ++i )
        xAbs += x[ i ] * x[ i ];
      return std::pow( xAbs, 0.25 );
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      const double u  = exact( x );
      const double u3 = 2.0 * u * u * u ;
      for( int i=0; i<dim; ++i ) grad[ i ] = x[ i ] / u3 ;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

  /**
   * \brief Data functions "Single Hump" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class Hump : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~Hump() {}
    Hump(Field globalShift, Field factor)
      : globalShift_(0.0)
      , factor_(1.0)
    {
      assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      double w=10.*x[0]*x[0]+10.*x[1];
      double v=(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
      double dwx = 20.*x[0];
      double dwy = 10.;
      double dwxx = 20.;
      double dwyy = 0.;
      double dvx = (1.-2.*x[0])*(x[1]-x[1]*x[1]);
      double dvy = (1.-2.*x[1])*(x[0]-x[0]*x[0]);
      double dvxx = -2.*(x[1]-x[1]*x[1]);
      double dvyy = -2.*(x[0]-x[0]*x[0]);
      Field grad[dim];
      grad[0] = exp(w)*(dwx*v*v + 2.*v*dvx);
      grad[1] = exp(w)*(dwy*v*v + 2.*v*dvy);
      double dxx = dwx*grad[0] + exp(w)*(dwxx*v*v+dwx*2.*v*dvx + 2.*dvx*dvx+2.*v*dvxx);
      double dyy = dwy*grad[1] + exp(w)*(dwyy*v*v+dwy*2.*v*dvy + 2.*dvy*dvy+2.*v*dvyy);
      return -dxx-dyy;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      double w=10.*x[0]*x[0]+10.*x[1];
      double v=(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
      return exp(w)*v*v;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      double w=10.*x[0]*x[0]+10.*x[1];
      double v=(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
      double dwx = 20.*x[0];
      double dwy = 10.;
      double dvx = (1.-2.*x[0])*(x[1]-x[1]*x[1]);
      double dvy = (1.-2.*x[1])*(x[0]-x[0]*x[0]);
      grad[0] = exp(w)*(dwx*v*v + 2.*v*dvx);
      grad[1] = exp(w)*(dwy*v*v + 2.*v*dvy);
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };


  /**
   * \brief Data functions "Rivere-Bastian" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class RiviereProblem : public DataFunctionIF<dim,DomainField,Field>
  {
    using DataFunctionIF<dim,DomainField,Field> :: SQR;
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~RiviereProblem() {}
    RiviereProblem(Field globalShift, Field factor)
      : globalShift_(0.0)
      , factor_(1.0)
    {
      assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      Field val = exact(x);
      Field x_part = val * SQR(-2.0 * (x[0] - 0.5)) - 2.0 * val;
      Field y_part = val * SQR(-2.0 * (x[1] - 0.5)) - 2.0 * val;
      return -(x_part + y_part);
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      Field power = -( SQR( x[0] - 0.5 ) + SQR( x[1] - 0.5 ) );
      return pow( M_E , power );
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      Field val = exact( x );
      grad[0] = val * ( -2.0*x[0] + 1.0 );
      grad[1] = val * ( -2.0*x[1] + 1.0 );
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };


  /**
   * \brief Data functions "Rivere-Bastian" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class HeatProblem : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~HeatProblem() {}
    HeatProblem(Field globalShift, Field factor)
      : globalShift_(0.0)
      , factor_(1.0)
    {
      assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    double scp(const DomainField x[dim]) const
    {
      double r2 = 0.0;
      for(int i=0; i<dim; ++i) r2 += x[i]*x[i];
      return r2;
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      return rhs(0.0, x);
    }

    virtual Field rhs  (const double time, const DomainField x[dim]) const
    {
      double r2 = scp( x );

      double  ux  = std::sin(M_PI*time)* std::exp(-10.0*r2);
      double  ut = M_PI*std::cos(M_PI*time)*std::exp(-10.0*r2);
      return(ut - (400.0*r2 - 20.0*dim)*ux);
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return exact(0.0, x);
    }

    virtual Field exact(const double time, const DomainField x[dim]) const
    {
      double r2 = scp( x );
      return(std::sin( M_PI * time) * std::exp( -10.0 * r2));
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = grad[1] = 0.0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      assert( false );
      val = exact( x );
      return true;
    }

    virtual bool boundaryDataFunction(const double time,
                                      const DomainField x[dim],
                                      Field & val) const
    {
      val = exact( time, x );
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BoundaryLayerProblem : public DataFunctionIF<dim,DomainField,Field>
  {
  public:
    virtual ~BoundaryLayerProblem() {}
    BoundaryLayerProblem(Field globalShift, Field factor)
      : eps_(factor)
    {
      assert( eps_ > 0 );
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j) k[i][j] = 0;
        k[i][i] = 1.; // eps_;
      }
    }

    void velocity(const DomainField x[dim], DomainField v[dim]) const
    {
      for(int i=0; i<dim; ++i)
        v[i] = 1./eps_;
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      return 0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      Field ret = 1;
      for (int i=0;i<dim;++i)
        ret *= u1( x[i] );
      return ret;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      for (int i=0;i<dim;++i)
      {
        grad[ i ] = 1.;
        for (int j=0;j<dim;++j)
          grad[ i ] *= (j==i)? du1( x[i] ):u1( x[i] );
      }
      return;
    }
    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  private:
    double u1(double x) const
    {
      return (exp(x/eps_)-1.) / (exp(1./eps_)-1.);
    }
    double du1(double x) const
    {
      return (exp(x/eps_)/eps_) / (exp(1./eps_)-1.);
    }
    const double eps_;
  };


  /**
   * \brief Data functions "CurvedRidges (deal.II step-14)" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class CurvedRidges : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;
  public:
    virtual ~CurvedRidges() {}
    CurvedRidges(Field globalShift, Field factor)
      : globalShift_(0.0)
      , factor_(1.0)
    {
      assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField p[dim]) const
    {
      double q = p[ 0 ];
      for (unsigned int i=1; i<dim; ++i)
      {
        q += std::sin(10*p[ i ]+5*p[ 0 ]*p[ 0 ]);
      }
      const double u = std::exp(q);
      double t1 = 1, t2 = 0, t3 = 0;
      for (unsigned int i=1; i<dim; ++i)
      {
        t1 += std::cos(10*p[ i ]+5*p[ 0 ]*p[ 0 ]) * 10 * p[ 0 ];
        t2 += 10*std::cos(10*p[ i ]+5*p[ 0 ]*p[ 0 ]) -
        100*std::sin(10*p[ i ]+5*p[ 0 ]*p[ 0 ]) * p[ 0 ]*p[ 0 ];
        t3 += 100*std::cos(10*p[ i ]+5*p[ 0 ]*p[ 0 ])*std::cos(10*p[ i ]+5*p[ 0 ]*p[ 0 ]) -
        100*std::sin(10*p[ i ]+5*p[ 0 ]*p[ 0 ]);
      }
      t1 = t1*t1;

      return -u*(t1+t2+t3);
    }

    virtual Field exact(const DomainField p[dim]) const
    {
      double q = p[ 0 ];
      for (unsigned int i=1; i<dim; ++i)
      {
        q += std::sin(10*p[ i ]+5*p[ 0 ]*p[ 0 ]);
      }
      const double exponential = std::exp(q);
      return exponential;
    }

    virtual void gradExact(const DomainField p[dim], Field grad[dim] ) const
    {
      double u = exact(p);
      grad[0] = 1.;
      for (unsigned int i=1; i<dim; ++i)
        grad[0] += std::cos(10*p[ i ]+5*p[ 0 ]*p[ 0 ]) * 10. * p[0];
      grad[0] *= u;
      for (int i=1;i<dim;++i)
      {
        grad[i] = std::cos(10*p[ i ]+5*p[ 0 ]*p[ 0 ]) * 10.;
        grad[i] *= u;
      }
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };


  /**
   * \brief Data functions "CurvedRidges (deal.II step-14)" for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class Excercise2_3 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field factor_;
    Field point_ [ dim ];
  public:
    virtual ~Excercise2_3() {}
    Excercise2_3(Field globalShift, Field factor)
      : globalShift_(0.0)
      , factor_(1.0)
    {
      for(int i=0; i<dim; ++i) point_[ i ] = 0.75;
      assert(dim == 2);
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_;
        for(int j=0; j<i; ++j)     k[i][j] = 0;
        for(int j=i+1; j<dim; ++j) k[i][j] = 0;
      }
    }

    virtual Field rhs  (const DomainField p[dim]) const
    {
      return 1.0;
    }

    virtual Field exact(const DomainField p[dim]) const
    {
      // exact solution not known
      return 0.0;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = 0;
      grad[1] = 0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };


  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_1 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    Field factor_[dim][dim];
  public:
    virtual ~BenchMark_1() {}
    BenchMark_1(Field globalShift, Field factor)
      : globalShift_( globalShift )
    {
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
        {
          if( i == j ) factor_[i][j] = 1.5;
          else if( std::abs( i - j ) == 1 ) factor_[i][j] = 0.5;
          else factor_[i][j] = 0;
        }
      }
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      // copy values
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
           k[i][j] = factor_[i][j];
      }
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      double x = arg[0];
      double y = arg[1];

      double uxx = -2.*y*(1-y)*16.;
      double uxy = (-2.*x+1)*(-2.*y+1)*16.;
      double uyy = -2.*x*(1-x)*16.;

      return -( factor_[0][0] * uxx + factor_[0][1] * uxy + factor_[1][0] * uxy + factor_[1][1] * uyy);
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      Field val = 16.;
      for(int i=0; i<dim; ++i)
        val *= (x[i] - (x[i]*x[i]));
      val += globalShift_ ;
      return val;
    }

    virtual void gradExact(const DomainField arg[dim], Field grad[dim] ) const
    {
      double x = arg[0];
      double y = arg[1];

      grad[0] = (-2.*x+1)*y*(1-y)*16.;
      grad[1] = (-2.*y+1)*x*(1-x)*16.;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_1_2 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    Field factor_[dim][dim];
  public:
    virtual ~BenchMark_1_2() {}
    BenchMark_1_2(Field globalShift, Field factor)
      : globalShift_(0.0)
    {
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
        {
          if( i == j ) factor_[i][j] = 1.5;
          else factor_[i][j] = 0.5;
        }
      }
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_[i][i];
        for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
        for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
      }
    }


    virtual Field rhs  (const DomainField arg[dim]) const
    {
      double x1 = 1.-arg[0];
      double y1 = 1.-arg[1];

      double uxx = -y1*y1*std::sin(x1*y1) + 6.*x1* y1*y1;
      double uxy = -x1*y1*std::sin(x1*y1) + std::cos(x1*y1) + 6.*y1* x1*x1;
      double uyy = -x1*x1*std::sin(x1*y1) + 2.* x1*x1*x1;
      return -(factor_[0][0] * uxx + factor_[0][1]*uxy + factor_[1][0]* uxy + factor_[1][1]*uyy);
    }

    virtual Field exact(const DomainField arg[dim]) const
    {
      double x1 = 1.-arg[0];
      double y1 = 1.-arg[1];
      return std::sin(x1*y1) + x1*x1*x1 * y1*y1;
    }

    virtual void gradExact(const DomainField arg[dim], Field grad[dim] ) const
    {
      double x1 = 1.-arg[0];
      double y1 = 1.-arg[1];

      grad[0] = -y1*std::cos(x1*y1) - 3.* (x1*y1)* (x1*y1);
      grad[1] = -x1*std::cos(x1*y1) - 2.*y1* x1*x1*x1;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_2 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    Field factor_[dim][dim];
    const Field delta_;
    const Field sqrtDelta_;
    const Field x1_;
    const Field x2_;
  public:
    virtual ~BenchMark_2() {}
    BenchMark_2(Field globalShift, Field factor)
      : globalShift_(0.0)
      , delta_(factor)
      , sqrtDelta_( sqrt(delta_) )
      , x1_ ( 8. * atan(1.) )
      , x2_ ( x1_ / sqrtDelta_ )
    {
      factor_[0][0] = 1;
      factor_[0][1] = factor_[1][0] = 0.0;
      factor_[1][1] = delta_;
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_[i][i];
        for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
        for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
      }
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      return 0.0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return std::sin(x1_* x[0])*exp(-x2_* x[1]);
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = x1_ * std::cos(x1_* x[0])*exp(-x2_ * x[1]);
      grad[1] = -x2_ * exact(x);
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact(x);
      // we have only neumann boundary here (see discretemodel.hh)
      return true;
      //return false;
    }
  };


  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_3 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    Field factor_[dim][dim];
    const Field delta_;
    const Field pi_;
    const Field cost_;
    const Field sint_;
  public:
    virtual ~BenchMark_3() {}
    BenchMark_3(Field globalShift, Field factor)
      : globalShift_(0.0)
      , delta_(1e-3)
      , pi_ ( 4. * atan(1.) )
      , cost_ ( std::cos( 40. * pi_ / 180. ) )
      , sint_ ( sqrt(1. - cost_*cost_) )
    {
      factor_[0][0] = cost_*cost_+delta_*sint_*sint_;
      factor_[1][0] = factor_[0][1] = cost_*sint_*(1.-delta_);
      factor_[1][1] = sint_*sint_+delta_*cost_*cost_;
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        k[i][i] = factor_[i][i];
        for(int j=0; j<i; ++j)     k[i][j] = factor_[i][j];
        for(int j=i+1; j<dim; ++j) k[i][j] = factor_[i][j];
      }
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      return 0.0;
    }

    double bndFunc (const double x,
                    const double lower,
                    const double upper,
                    const double lowval,
                    const double upval)  const
    {
      if( x <= lower ) return lowval;
      if( x >= upper ) return upval;

      const double scale = (x - lower)/(upper - lower);
      assert( scale >= 0 && scale <= 1 );

      return (1. - scale) * lowval + scale * upval;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      const int bndId = this->getBoundaryId( x , x );
      if( bndId == 0 )
      {
        return bndFunc(x[1],0.2,0.3,1,0.5);
      }
      else if( bndId == 1 )
      {
        return bndFunc(x[1],0.7,0.8,0.5,0);
      }
      else if( bndId == 2 )
      {
        return bndFunc(x[0],0.2,0.3,1,0.5);
      }
      else if( bndId == 3 )
      {
        return bndFunc(x[0],0.7,0.8,0.5,0);
      }
      return 0.5;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = 0.0;
      grad[1] = 0.0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      // only dirichlet for this problem
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_4 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
  public:
    virtual ~BenchMark_4() {}
    BenchMark_4(Field globalShift, Field factor)
      : globalShift_(0.0)
    {
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j )
          k[i][j] = 0;

      if( omega1(x) )
      {
        k[0][0] = 1e2;
        k[1][1] = 1e1;
        if( dim > 2 )
        {
          k[2][2] = 1e3;
          k[0][1] = 0.5;
          k[1][0] = 0.5;
          k[1][2] = 5;
          k[2][1] = 5;
        }
      }
      else
      {
        k[0][0] = 1e-2;
        k[1][1] = 1e-3;
        if( dim > 2 )
        {
          //k[0][1] = 1e-1;
          //k[1][0] = 1e-1;
          //k[1][2] = 1e-2;
          //k[2][1] = 1e-2;
          k[2][2] = 1e-2;
        }
      }
    }

    bool omega1(const DomainField x[dim]) const
    {
      if( dim == 2 )
      {
        if (x[0] <= 0.5)
        {
          int inty = int(10.0 * (x[1] + 0.15));
          // if even then omega1 else omega2
          return ((inty%2) == 0);
        }
        else
        {
          int inty = int(10.0 * x[1]);
          // if even then omega1 else omega2
          return ((inty%2) == 0);
        }
      }
      else if ( dim == 3 )
      {
        if (x[0] <= 0.5)
        {
          int inty = int(16.0 * (x[1] + 0.0625));
          // if even then omega1 else omega2
          return ((inty%2) == 0);
        }
        else
        {
          int inty = int(16.0 * x[1]);
          // if even then omega1 else omega2
          return ((inty%2) == 0);
        }
      }
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      return 0.0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      double val = (1.0 - x[0]);
      if( dim > 2 )
        val *= (1 - x[2]);
      return val;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = 0.0;
      grad[1] = 0.0;
      grad[dim-1] = 0.0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_5 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field delta_;
    const Field pi;
  public:
    virtual ~BenchMark_5() {}
    BenchMark_5(Field globalShift, Field factor)
      : globalShift_(0.0)
      , delta_(1e-3)
      , pi ( 4. * atan(1.) )
    {
    }

    virtual bool constantLocalK () const { return false; }

    virtual void K(const DomainField arg[dim], Field k[dim][dim] ) const
    {
      double x = arg[0];
      double y = arg[1];
      double rt = x*x+y*y;
      k[0][0] = (y*y+delta_*x*x)/rt;
      k[1][1] = (x*x+delta_*y*y)/rt;
      k[1][0] = k[0][1] = -(1-delta_)*x*y/rt;
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      Field k[dim][dim];
      K(arg,k);
      double x = arg[0];
      double y = arg[1];
      double rt = x*x+y*y;

      double ux = pi * std::cos(pi*x)*std::sin(pi*y);
      double uy = pi * std::cos(pi*y)*std::sin(pi*x);

      double f0 = std::sin(pi*x)*std::sin(pi*y)*pi*pi*(1+delta_)*(x*x+y*y)
                + std::cos(pi*x)*std::sin(pi*y)*pi*(1.-3.*delta_)*x
                + std::cos(pi*y)*std::sin(pi*x)*pi*(1.-3.*delta_)*y
                + std::cos(pi*y)*std::cos(pi*x)*2.*pi*pi*(1.-delta_)*x*y;
      double kxx = k[0][0];
      double kyy = k[1][1];
      double kxy = k[0][1];
      return (f0+2.*(x*(kxx*ux+kxy*uy)+y*(kxy*ux+kyy*uy)))/rt;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return std::sin(pi*x[0])*std::sin(pi*x[1]);
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = pi * std::cos(pi*x[0])*std::sin(pi*x[1]);
      grad[1] = pi * std::cos(pi*x[1])*std::sin(pi*x[0]);
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_6 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field delta_;
    const Field cost_;
    const Field sint_;
  public:
    virtual ~BenchMark_6() {}
    BenchMark_6(Field globalShift, Field factor)
      : globalShift_(0.0)
      , delta_(0.2)
      , cost_ ( 1./sqrt(1.+delta_*delta_) )
      , sint_ ( delta_*cost_ )
    {
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      double phi1 = x[1] - delta_ * (x[0] - .5) - .475;
      double phi2 = phi1 - .05;

      double alpha = 0.0;
      double beta  = 0.0;
      if (phi1<0 || phi2>0)
      {
         alpha = 1.0;
         beta  = 0.1;
      }
      else
      {
         alpha = 100.0;
         beta  = 10.0;
      }

      k[0][0] = alpha*cost_*cost_+beta*sint_*sint_;
      k[0][1] = k[1][0] = cost_*sint_*(alpha-beta);
      k[1][1] = alpha*sint_*sint_+beta*cost_*cost_;
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      return 0.0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return - x[0] - x[1] * delta_;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = -1.0;
      grad[1] = -delta_;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      // we have neumann boundary here
      return true;
    }
  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_7 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field delta_;
  public:
    virtual ~BenchMark_7() {}
    BenchMark_7(Field globalShift, Field factor)
      : globalShift_(0.0)
      , delta_(0.2)
    {
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      double phi1 = phi(x);
      double phi2 = phi1 - .05;

      int dom = domain( phi1, phi2 );
      if( dom == 1 || dom == 3 )
      {
        k[0][0] = k[1][1] = 1;
        k[1][0] = k[0][1] = 0;
      }
      else
      {
        k[0][0] = k[1][1] = 0.01;
        k[1][0] = k[0][1] = 0;
      }
    }

    double phi(const DomainField x[dim]) const
    {
      return x[1] - delta_ * (x[0] - .5) - .475;
    }

    int domain(const double phi1, const double phi2) const
    {
      if (phi1<0)
        return 1;
      else if (phi2<0)
        return 2;
      else
        return 3;
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      return 0.0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      double phi1 = phi(x);
      double phi2 = phi1 - .05;

      int dom = domain( phi1, phi2 );
      if( dom == 1 )
      {
        return -phi1;
      }
      else if( dom == 2 )
      {
        return -phi1/.01;
      }
      else
      {
        return -phi2 - 5.;
      }
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      double phi1 = phi(x);
      double phi2 = phi1 - .05;
      int dom = domain( phi1, phi2 );
      if( dom == 1 || dom == 3 )
      {
        grad[0] = delta_;
        grad[1] = -1.0;
      }
      else // if (dom == 2)
      {
        grad[0] = delta_/.01;
        grad[1] = - 1./.01;
      }
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

#if PROBLEM_8
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/fvspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

  typedef LeafGridPart<GridType> ParamGridPartType;
  typedef FunctionSpace<double,double,GridType::dimension,1> ParamFunctionSpaceType;
  typedef LagrangeDiscreteFunctionSpace<ParamFunctionSpaceType,ParamGridPartType,2>
           ParamDiscreteFunctionSpaceType;
  // typedef FiniteVolumeSpace<ParamFunctionSpaceType,ParamGridPartType,0>
  //         ParamDiscreteFunctionSpaceType;
  typedef AdaptiveDiscreteFunction<ParamDiscreteFunctionSpaceType>
          ParamDiscreteFunctionType;
  ParamDiscreteFunctionType* paramDiscreteFunction;

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_8 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field delta_;
  public:
    virtual ~BenchMark_8() {}
    BenchMark_8(Field globalShift, Field factor)
      : globalShift_(0.0)
      , delta_(0.)
    {
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      k[0][0] = k[1][1] = 1;
      k[1][0] = k[0][1] = 0;
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      FieldVector<double,dim> p(0);
      FieldVector<double,1> ret(0);
      p[0] = arg[0];
      p[1] = arg[1];

      const ParamGridPartType& gridPart = paramDiscreteFunction->space().gridPart();
      const typename ParamGridPartType::IndexSetType& index = gridPart.indexSet();
      HierarchicSearch<GridType,ParamGridPartType::IndexSetType>
        search(gridPart.grid(),index);
      typename ParamDiscreteFunctionType::LocalFunctionType lf =
        paramDiscreteFunction->localFunction((*(search.findEntity(p))));

      lf.evaluate(search.findEntity(p)->geometry().local(p),ret);
      /*
      std::cout << search.findEntity(p)->geometry().local(p)
                << " " << p
                << " " << ret[0] << std::endl;
                */
      return ret[0];
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return 0.0;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = 0.0;
      grad[1] = 0.0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = 0.0;
      return true;
    }
    virtual int getBoundaryId(const DomainField x[dim],
                              const DomainField n[dim]) const
    {
      if (std::abs(n[0]) > 1e-8) {
        if (n[1] < 0.) return 0;
        else return 1;
      }
      if (std::abs(n[0]) < 1e-8) {
        if (n[1] < 0.) return 2;
        else return 3;
      }
      abort();
      return 0;
    }
  };
#endif
  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark_9 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field globalShift_;
    const Field delta_;
  public:
    virtual ~BenchMark_9() {}
    BenchMark_9(Field globalShift, Field factor)
      : globalShift_(0.0)
      , delta_(0.)
    {
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      k[0][0] = k[1][1] = 1;
      k[1][0] = k[0][1] = 0;
    }


    virtual Field rhs  (const DomainField arg[dim]) const
    {
      return 0.;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return 0.0;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = 0.0;
      grad[1] = 0.0;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = (x[0]<0.6)?1.:0.;
      if (x[0]<=0 || x[0]>=1)
        return false;
      if (x[1]<=0 || x[1]>=1)
        return false;
      return true;
    }
  };


  /////////////////////////////////////////////////////////////////////
  //
  //  3D Benchmark Problems
  //
  /////////////////////////////////////////////////////////////////////
  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark3d_1 : public DataFunctionIF<dim,DomainField,Field>
  {

    Field K_[dim][dim];
  public:
    virtual ~BenchMark3d_1() {}
    BenchMark3d_1(Field globalShift, Field factor)
    {
      if( dim != 3 )
      {
        std::cerr << "Problem only implemented for dim=3" << std::endl;
        std::abort();
      }

      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
        {
          if( i == j ) K_[i][j] = 1.0;
          else if( std::abs( i - j ) == 1 )
          {
            K_[i][j] = 0.5;
          }
          else K_[i][j] = 0;
        }
      }
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      // just copy tensor
      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j)
          k[i][j] = K_[i][j];
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      return M_PI*M_PI*(3.0*std::sin(M_PI*x[0])*std::sin(M_PI*(x[1]+0.5))*std::sin(M_PI*(x[2]+(1.0/3.0)))
            -std::cos(M_PI*x[0])*std::cos(M_PI*(x[1]+0.5))*std::sin(M_PI*(x[2]+(1.0/3.0)))
            -std::sin(M_PI*x[0])*std::cos(M_PI*(x[1]+0.5))*std::cos(M_PI*(x[2]+(1.0/3.0))));
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return 1.0+std::sin(M_PI*x[0])*std::sin(M_PI*(x[1]+0.5))*std::sin(M_PI*(x[2]+(1.0/3.0)));
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      grad[0] = M_PI*std::cos(M_PI*x[0])*std::sin(M_PI*(x[1]+0.5))*std::sin(M_PI*(x[2]+(1.0/3.0)));
      grad[1] = M_PI*std::sin(M_PI*x[0])*std::cos(M_PI*(x[1]+0.5))*std::sin(M_PI*(x[2]+(1.0/3.0)));
      grad[2] = M_PI*std::sin(M_PI*x[0])*std::sin(M_PI*(x[1]+0.5))*std::cos(M_PI*(x[2]+(1.0/3.0)));
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }

  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark3d_3 : public DataFunctionIF<dim,DomainField,Field>
  {

    Field K_[dim][dim];
  public:
    virtual ~BenchMark3d_3() {}
    BenchMark3d_3(Field globalShift, Field factor)
    {
      if( dim != 3 )
      {
        std::cerr << "Problem only implemented for dim=3" << std::endl;
        std::abort();
      }

      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j)
          K_[ i ][ j ] = 0 ;

      // set diagonal
      K_[0][0] = 1;
      K_[1][1] = 1;
      K_[2][2] = 1e3 ;
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
      {
        for(int j=0; j<dim; ++j)
        {
          k[i][j] = K_[i][j];
        }
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      return 1002.0*4.0*M_PI*M_PI*std::sin(2.0*M_PI*x[0])*std::sin(2.0*M_PI*x[1])*std::sin(2.0*M_PI*x[2]);
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      return std::sin(2.0*M_PI*x[0])*std::sin(2.0*M_PI*x[1])*std::sin(2.0*M_PI*x[2]);
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      const Field pi2 = 2.0*M_PI;
      grad[0] = pi2 * std::cos( pi2 * x[0]) * std::sin( pi2 * x[1]) * std::sin( pi2 * x[2]);
      grad[1] = pi2 * std::sin( pi2 * x[0]) * std::cos( pi2 * x[1]) * std::sin( pi2 * x[2]);
      grad[2] = pi2 * std::sin( pi2 * x[0]) * std::sin( pi2 * x[1]) * std::cos( pi2 * x[2]);
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }

  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark3d_4 : public DataFunctionIF<dim,DomainField,Field>
  {
    const Field tau_ ;
  public:
    virtual ~BenchMark3d_4() {}
    BenchMark3d_4(Field globalShift, Field factor)
      : tau_( 0.2 )
    {
      if( dim != 3 )
      {
        std::cerr << "Problem only implemented for dim=3" << std::endl;
        std::abort();
      }
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j )
          k[i][j] = 0;

      for(int i=0; i<2; ++ i ) k[i][i] = 1;
      if( dim == 3 )
        k[dim-1][dim-1] = tau_ ;
    }

    int getDomain(const DomainField x[dim]) const
    {
      return 0;
    }

    virtual Field rhs  (const DomainField arg[dim]) const
    {
      return 0.0;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      Field u,u_x,u_y,u_z;
      calculsolwell(x[0],x[1],x[2],u,u_x,u_y,u_z,1.0,0.0,0.0,1.0,0.0,tau_);
      return u;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      Field u,u_x,u_y,u_z;
      calculsolwell(x[0],x[1],x[2],u,u_x,u_y,u_z,1.0,0.0,0.0,1.0,0.0,tau_);
      grad[0] = u_x;
      grad[1] = u_y;
      grad[dim-1] = u_z;
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }

    void calculsolwell (Field x, Field y, Field z,
                        Field& p, Field& px, Field& py, Field& pz,
                        Field lxx, Field lxy, Field lxz, Field lyy, Field lyz, Field lzz) const
    {
      Field alpha,mu0,
        A11,A12,A13,A21,A22,A23,A31,A32,A33,
        //XX,
        YY,ZZ,s2,a2,a1,a0,sh,frho,mu,ch,frho_p,
        GP_X,GP_Y,GP_Z;

      alpha = 1.0;
      mu0  = 1.85334252045523;

      lxx = 1.;
      lyy = 1.;
      lzz = 0.2;
      lxy = 0.;
      lxz = 0.;
      lyz = 0.;

      A11 = 0.122859140982213;
      A12 = 0.0;
      A13 = -1.68776357811111;
      A21 = 0.0;
      A22 = 0.76472449133173;
      A23 = 0.0;
      A31 = 0.754790818120945;
      A32 = 0.0;
      A33 = 0.274721390893459;

      a2 =  0.000603774740915493;


      //-------------------------------------------------
      // Eval solution
      //-------------------------------------------------
      //XX = A11*x+A12*y+A13*z;
      YY = A21*x+A22*y+A23*z;
      ZZ = A31*x+A32*y+A33*z;

      // computation of s2 = sh^2, solution to

      //  a2 S^2 + a1 S + a0 = 0

      a1 = a2 - ZZ*ZZ - YY*YY;

      a0 = - YY*YY;

      s2 = (-a1+sqrt(a1*a1 - 4.0*a0*a2))/(2.0*a2);

      sh = sqrt(s2);

      // argsh(w) = log(w+sqrt(w**2+1))

      frho = sh+sqrt(s2+1.0);

      mu = log(frho);

      p = alpha *(mu - mu0);

      //-------------------------------------------------
      // Eval gradient
      //-------------------------------------------------
      ch = (frho+1.0/frho)*0.5;

      frho_p = alpha /(  (ZZ*ZZ*sh)/(ch*ch*ch) + (YY*YY*ch)/(sh*sh*sh) );

      GP_X = 0.0;
      GP_Y = frho_p *  YY  / ( sh * sh );
      GP_Z = frho_p *  ZZ  / ( ch * ch );

      px = A11*GP_X + A21*GP_Y + A31*GP_Z;
      py = A12*GP_X + A22*GP_Y + A32*GP_Z;
      pz = A13*GP_X + A23*GP_Y + A33*GP_Z;
    }

  };

  /**
   * \brief Data functions for the Poisson problem.
   *
   * \ingroup PoissonDataFunctions
   */
  template <int dim, class DomainField, class Field>
  class BenchMark3d_5 : public DataFunctionIF<dim,DomainField,Field>
  {
    enum { numDomain = 4 };
    Field tensor_[ numDomain ][ dim ];
    Field alpha_[ numDomain ];
    Field trace_[ numDomain ];
    const Field pi2_ ;
  public:
    virtual ~BenchMark3d_5() {}
    BenchMark3d_5(Field globalShift, Field factor)
      : pi2_( 2.0 * M_PI )
    {
      if( dim != 3 )
      {
        std::cerr << "Problem only implemented for dim=3" << std::endl;
        std::abort();
      }

      {
        Field (&tensor)[dim] = tensor_[ 0 ];
        tensor[0] = 1.0;
        tensor[1] = 10.0;
        tensor[2] = 0.01;
      }

      {
        Field (&tensor)[dim] = tensor_[ 1 ];
        tensor[0] = 1.0;
        tensor[1] = 0.1;
        tensor[2] = 100.0;
      }

      {
        Field (&tensor)[dim] = tensor_[ 2 ];
        tensor[0] = 1.0;
        tensor[1] = 0.01;
        tensor[2] = 10.0;
      }

      {
        Field (&tensor)[dim] = tensor_[ 3 ];
        tensor[0] = 1.0;
        tensor[1] = 100.0;
        tensor[2] = 0.1;
      }

      alpha_[ 0 ] = 0.1;
      alpha_[ 1 ] = 10.0;
      alpha_[ 2 ] = 100.0;
      alpha_[ 3 ] = 0.01;

      for(int i=0; i<numDomain; ++i )
      {
        trace_[ i ] = 0;
        for( int j=0; j<dim; ++ j)
          trace_[ i ]  += tensor_[ i ][ j ];
      }
    }

    virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
    {
      const int domain = getDomain( x );
      assert( domain >= 0 && domain < numDomain );
      for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j )
          k[i][j] = 0;
      // set diagonal
      for(int j=0; j<dim; ++j )
        k[j][j] = tensor_[ domain ][ j ];
    }

    int getDomain(const DomainField x[dim]) const
    {
      if (x[1]<=0.5)
      {
        if (x[2]<=0.5) return 0; else return 3;
      }
      else
      {
        if (x[2]<=0.5) return 1; else return 2;
      }
    }

    virtual Field rhs  (const DomainField x[dim]) const
    {
      const int domain = getDomain(x);
      assert( domain >= 0 && domain < numDomain );
      Field result =  4.0*M_PI*M_PI*std::sin(pi2_*x[0])*std::sin(pi2_*x[1])*std::sin(pi2_*x[2]);
      result *= alpha_[ domain ] * trace_[ domain ];
      return result;
    }

    virtual Field exact(const DomainField x[dim]) const
    {
      const int domain = getDomain(x);
      assert( domain >= 0 && domain < numDomain );
      Field val = std::sin(pi2_*x[0])*std::sin(pi2_*x[1])*std::sin(pi2_*x[2]);
      val *= alpha_[ domain ];
      return val ;
    }

    virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const
    {
      const int domain = getDomain(x);
      assert( domain >= 0 && domain < numDomain );

      const Field x_pi = pi2_*x[0] ;
      const Field y_pi = pi2_*x[1] ;
      const Field z_pi = pi2_*x[2] ;

      grad[0] = pi2_ * std::cos( x_pi ) * std::sin( y_pi ) * std::sin( z_pi );
      grad[1] = pi2_ * std::sin( x_pi ) * std::cos( y_pi ) * std::sin( z_pi );
      grad[2] = pi2_ * std::sin( x_pi ) * std::sin( y_pi ) * std::cos( z_pi );

      for(int i=0; i<dim; ++i )
        grad[ i ] *= alpha_[ domain ];
    }

    virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
    {
      val = exact( x );
      return true;
    }
  };

}
}
}
#endif
