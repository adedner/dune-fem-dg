#ifndef ELLIPTPROBLEM_CC
#define ELLIPTPROBLEM_CC

#include <cmath>
#include <set>
#include <cassert>

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


template <int dim, class DomainField, class Field> 
class SinSinSin : public DataFunctionIF<dim,DomainField,Field>
{
  Field K_[dim][dim];
public:  
  virtual ~SinSinSin() {}
  SinSinSin(Field globalShift, Field factor)
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
        sum += K_[j][i] * laplace( x, comp );
      }
    }
    return -sum;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    const Field pi = 2.0 * M_PI;
    Field val = 1.0;
    for(int i=0; i<dim; ++i) 
    {
      val *= sin( pi * x[i] );
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
        val *= sin( pi * x[ i ] );
      }

      // minus because sin'' = -sin 
      return -val;
    }
    else 
    {
      for( int i=0; i<dim-1; ++i) 
      {
        val *= cos( pi * x[ comp[i] ] );
      }

      for(int i=0; i<dim; ++i) 
      {
        if( comps.find( i ) == comps.end() )
          val *= sin( pi * x[ i ] );
      }
      return val;
    }
  }

  double gradient(const DomainField x[dim], const int comp) const 
  {
    const Field pi = 2.0 * M_PI; 
    Field val = pi * cos( pi * x[ comp ] );
    // add other components 
    for(int j=1; j<dim; ++j) 
    {
      val *= sin( pi * x[ (comp + j) % dim ] );
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
    Field sin_x = sin(2.0*M_PI*x[0]);
    Field sin_y = sin(2.0*M_PI*x[1]);

    Field val = 8.0 * M_PI * M_PI * sin_x * sin_y ;
    val *= factor_;
    return val;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field val = sin(2.0*M_PI*x[0]) * sin(2.0*M_PI*x[1]);
    val += globalShift_;
    return val;
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 2.0*M_PI*cos(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1]);
    grad[1] = 2.0*M_PI*sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);

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
    Field cos_x = cos(2.0*M_PI*x[0]);
    Field cos_y = cos(2.0*M_PI*x[1]);
       
    Field val = 8.0 * M_PI*M_PI* cos_x * cos_y ;
    val *= factor_;
    return val;
    //return 0;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    Field val = cos(2.0*M_PI*x[0]) * cos(2.0*M_PI*x[1]);
    val += globalShift_;
    return val;
    //return x[1];
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = 2.0*M_PI*-sin(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1]);
    grad[1] = 2.0*M_PI*cos(2.0*M_PI*x[0])*-sin(2.0*M_PI*x[1]);
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

//! problem from Castillo paper 
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

//! problem Einstringende Ecke 
template <class DomainField, class Field> 
class InSpringingCorner2d : public DataFunctionIF<2,DomainField,Field>
{
  enum { dim = 2 };
  const Field globalShift_;
  const Field factor_;
  const Field lambda_;
public:  
  virtual ~InSpringingCorner2d() {}
  InSpringingCorner2d(Field globalShift, Field factor)
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
    return pow(r2,lambda_*0.5)*sin(lambda_*phi);
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
    grad[0]= lambdaPow * r2dx * sin(lambda_ * phi)
             + lambda_ * cos( lambda_ * phi) * phidx * pow(r2,lambda_ * 0.5);
    grad[1]= lambdaPow * r2dy * sin(lambda_ * phi)
             + lambda_ * cos( lambda_ * phi) * phidy * pow(r2,lambda_*0.5);
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

//! problem Einstringende Ecke 3d
template <int dim, class DomainField, class Field> 
class InSpringingCorner : public DataFunctionIF<dim,DomainField,Field>
{
  typedef InSpringingCorner2d<DomainField, Field> Corner2dType ;
  Corner2dType corner2d_;
  const Field factor_;

public:  
  virtual ~InSpringingCorner() {}
  InSpringingCorner(Field globalShift, Field factor)
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

//! problem Einstringende Ecke 3d
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

//! problem single Hump
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


//! problem Riviere-Bastian 
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


//! problem Riviere-Bastian 
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


//! problem CurvedRidges (deal.II step-14) 
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


//! problem CurvedRidges (deal.II step-14) 
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




#include "benchmark.cc"
#endif
