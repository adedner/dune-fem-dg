#ifndef FVCA5_BENCHMARK_CC
#define FVCA5_BENCHMARK_CC

#include <cmath>
#include <dune/common/exceptions.hh>
#include "benchmark.hh"

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

    double uxx = -y1*y1*sin(x1*y1) + 6.*x1* y1*y1;
    double uxy = -x1*y1*sin(x1*y1) + cos(x1*y1) + 6.*y1* x1*x1;
    double uyy = -x1*x1*sin(x1*y1) + 2.* x1*x1*x1;
    return -(factor_[0][0] * uxx + factor_[0][1]*uxy + factor_[1][0]* uxy + factor_[1][1]*uyy);
  }

  virtual Field exact(const DomainField arg[dim]) const
  {
    double x1 = 1.-arg[0];
    double y1 = 1.-arg[1];
    return sin(x1*y1) + x1*x1*x1 * y1*y1;
  }
  
  virtual void gradExact(const DomainField arg[dim], Field grad[dim] ) const 
  {
    double x1 = 1.-arg[0];
    double y1 = 1.-arg[1];

    grad[0] = -y1*cos(x1*y1) - 3.* (x1*y1)* (x1*y1);
    grad[1] = -x1*cos(x1*y1) - 2.*y1* x1*x1*x1;
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};


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
    return sin(x1_* x[0])*exp(-x2_* x[1]);
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = x1_ * cos(x1_* x[0])*exp(-x2_ * x[1]);
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
    , cost_ ( cos( 40. * pi_ / 180. ) )
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

    double ux = pi * cos(pi*x)*sin(pi*y);
    double uy = pi * cos(pi*y)*sin(pi*x);

    double f0 = sin(pi*x)*sin(pi*y)*pi*pi*(1+delta_)*(x*x+y*y) 
              + cos(pi*x)*sin(pi*y)*pi*(1.-3.*delta_)*x
              + cos(pi*y)*sin(pi*x)*pi*(1.-3.*delta_)*y
              + cos(pi*y)*cos(pi*x)*2.*pi*pi*(1.-delta_)*x*y;
    double kxx = k[0][0];
    double kyy = k[1][1];
    double kxy = k[0][1];
    return (f0+2.*(x*(kxx*ux+kxy*uy)+y*(kxy*ux+kyy*uy)))/rt;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    return sin(pi*x[0])*sin(pi*x[1]);
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = pi * cos(pi*x[0])*sin(pi*x[1]);
    grad[1] = pi * cos(pi*x[1])*sin(pi*x[0]);
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }
};


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
      DUNE_THROW(Dune::NotImplemented,"Problem only implemented for dim=3");
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
    return M_PI*M_PI*(3.0*sin(M_PI*x[0])*sin(M_PI*(x[1]+0.5))*sin(M_PI*(x[2]+(1.0/3.0)))
          -cos(M_PI*x[0])*cos(M_PI*(x[1]+0.5))*sin(M_PI*(x[2]+(1.0/3.0)))
          -sin(M_PI*x[0])*cos(M_PI*(x[1]+0.5))*cos(M_PI*(x[2]+(1.0/3.0))));
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    return 1.0+sin(M_PI*x[0])*sin(M_PI*(x[1]+0.5))*sin(M_PI*(x[2]+(1.0/3.0)));
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    grad[0] = M_PI*cos(M_PI*x[0])*sin(M_PI*(x[1]+0.5))*sin(M_PI*(x[2]+(1.0/3.0)));
    grad[1] = M_PI*sin(M_PI*x[0])*cos(M_PI*(x[1]+0.5))*sin(M_PI*(x[2]+(1.0/3.0)));
    grad[2] = M_PI*sin(M_PI*x[0])*sin(M_PI*(x[1]+0.5))*cos(M_PI*(x[2]+(1.0/3.0)));
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }

};


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
      DUNE_THROW(Dune::NotImplemented,"Problem only implemented for dim=3");
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
    return 1002.0*4.0*M_PI*M_PI*sin(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1])*sin(2.0*M_PI*x[2]);
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    return sin(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1])*sin(2.0*M_PI*x[2]);
  }
  
  virtual void gradExact(const DomainField x[dim], Field grad[dim] ) const 
  {
    const Field pi2 = 2.0*M_PI;
    grad[0] = pi2 * cos( pi2 * x[0]) * sin( pi2 * x[1]) * sin( pi2 * x[2]);
    grad[1] = pi2 * sin( pi2 * x[0]) * cos( pi2 * x[1]) * sin( pi2 * x[2]);
    grad[2] = pi2 * sin( pi2 * x[0]) * sin( pi2 * x[1]) * cos( pi2 * x[2]);
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x ); 
    return true; 
  }

};


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
      DUNE_THROW(Dune::NotImplemented,"Problem only implemented for dim=3");
    }
  }

  virtual void K(const DomainField x[dim], Field k[dim][dim] ) const
  {
    for(int i=0; i<dim; ++i) 
      for(int j=0; j<dim; ++j )
        k[i][j] = 0;

    for(int i=0; i<2; ++ i ) k[i][i] = 1;
    k[2][2] = tau_ ;
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
      DUNE_THROW(Dune::NotImplemented,"Problem only implemented for dim=3");
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
    Field result =  4.0*M_PI*M_PI*sin(pi2_*x[0])*sin(pi2_*x[1])*sin(pi2_*x[2]);
    result *= alpha_[ domain ] * trace_[ domain ];
    return result;
  }

  virtual Field exact(const DomainField x[dim]) const
  {
    const int domain = getDomain(x);
    assert( domain >= 0 && domain < numDomain );
    Field val = sin(pi2_*x[0])*sin(pi2_*x[1])*sin(pi2_*x[2]);
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

    grad[0] = pi2_ * cos( x_pi ) * sin( y_pi ) * sin( z_pi );
    grad[1] = pi2_ * sin( x_pi ) * cos( y_pi ) * sin( z_pi );
    grad[2] = pi2_ * sin( x_pi ) * sin( y_pi ) * cos( z_pi );

    for(int i=0; i<dim; ++i ) 
      grad[ i ] *= alpha_[ domain ];
  }
  
  virtual bool boundaryDataFunction(const DomainField x[dim], Field & val) const
  {
    val = exact( x );
    return true; 
  }
};


#endif
