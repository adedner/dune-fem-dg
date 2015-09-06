#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <tgmath.h>

#include "fluxprojection.hh"
#include "polynomial.hh"

template <class Range, class NumFlux>
struct Hermite
{
  struct Op;

  typedef Range ResultType;
  typedef typename Range::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Hermite<Range,NumFlux> This;
  typedef std::vector<Range*> Vector;
  static typename Op::Approximate getApprox()
  {
    const std::string approxNames [] = { "fd3", "fd4", "fd6" };
    const typename Op::Approximate approx [] = { Op::fd3, Op::fd4, Op::fd6 };
    const int approxNumber = Dune::Fem::Parameter::getEnum("hermite.approx", approxNames );
    return approx[approxNumber];
  }
    
  Hermite( const DiscreteFunctionSpaceType &space,
           int q, int q1, typename Op::Approximate approx)
  : space_(space),
    numflux_(0),
    q_(q), q1_(q1), 
    tau_( Dune::Fem::Parameter::getValue<double>("fixedTimeStep" ) ),
    x_(q+1), stuetz_(0),
    operator_(tau_, q, approx)
  {
  }
  Hermite( const DiscreteFunctionSpaceType &space, const NumFlux &numflux)
  : Hermite( space,numflux,
             Fem::Parameter::getValue<int>("hermite.q"),
             Fem::Parameter::getValue<int>("hermite.q1"),
             getApprox() )
  {
  }
  bool add( const Range &u, const Range &fu)
  {
    return operator_.add(u,fu);
  }
  template <class DF>
  bool add( const DF &v, const DF &fv)
  {
    Range u("u-tmp",space_),
          fu("fu-tmp",space_);
    Dune::Fem::WeightDefault<GridPartType> weight;
#if HERMITEDG==2
    Dune::Fem::FullProjectionImpl::project(v,u, weight, *numflux_);
#else
    Dune::Fem::FluxProjectionImpl::project(v,u, weight, *numflux_);
#endif
    Dune::Fem::VtxProjectionImpl::project(fv,fu,weight);
    return add(u,fu);
  }
  const typename Range::DiscreteFunctionSpaceType &space() const
  {
    return space_;
  }
  void setup(const NumFlux &numflux)
  {
    numflux_ = &numflux;
    std::cout << "#" << " q=" << q_ << " q1=" << q1_;
    int i=0;
    for (;i<q1_;++i)
    {
      x_[i] = 0;
      std::cout << " " << x_[i];
    }
    for (;i<=q_;++i)
    {
      x_[i] = 1.; // tau;
      std::cout << " " << x_[i];
    }
    std::cout << std::endl;
  }
  void init()
  {
    if (stuetz_.size()==0)
      for (int i=0;i<q_+1;++i)
        stuetz_.push_back( new Range("tmp", operator_.func().space() ) );
    for (int i=0;i<q1_;++i)
    {
      stuetz_[i]->assign( *( operator_(i,0) ) );
      *(stuetz_[i]) *= pow(tau_,i);
    }
    for (int i=q1_;i<=q_;++i)
    {
      stuetz_[i]->assing( *(operator_(i-q1_,1) ) );
      *(stuetz_[i]) *= pow(tau_,i-q1_);
    }
    diviDiff( );
  }
  void evaluate(double point, Range &val) const
  {
    val.clear();
    for (int i=0;i<stuetz_.size();++i)
    {
      Polynomial p(i,x_);
      val.axpy(p.evaluate(point), *(stuetz_[i]));
    }
  }
  template <class DF>
  void evaluate(double point, DF &val) const
  {
    Range v("tmp",space_);
    evaluate(point,v);
#if HERMITEDG==2
    Dune::Fem::FullProjectionImpl::project(v,val, *numflux_);
#else
    Dune::Fem::FluxProjectionImpl::project(v,val, *numflux_);
#endif
  }
  void derivative(double point, Range &val) const
  {
    val.clear();
    for (int i=1;i<stuetz_.size();++i)
    {
      Polynomial p(i,x_);
      val.axpy(p.derivative(1,point), *(stuetz_[i]));
    }
    val /= tau_;
  }
  /*
  Range derivative(int k,double point) const
  {
    Range val(0);
    for (int i=0;i<stuetz_.size();++i)
    {
      Polynomial p(i,x_);
      val.axpy(p.derivative(k,point), *(stuetz_[i]));
    }
    val /= pow(tau_,k);
    return val;
  }
  Range derivative_end(int k) const
  {
    Range val(0);
    for (int i=0;i<stuetz_.size();++i)
    {
      Polynomial p(i,x_);
      val.axpy(p.derivative_end(k), stuetz_[i]);
    }
    val /= pow(tau_,k);
    return val;
  }
  */

  protected:

  void diviDiff() 
  {
    int n = x_.size();
    int i, j, k, l, m;
    // std::vector<Range> koeffSchema( n*n, *(stuetz_[0]) );
    std::vector<Range> koeffSchema( n*n, {"tmp",stuetz_[0]->space()} );

    for(i=n-1; i>=0; i--) 
    {
      for(j=i; j<n; j++) 
      {
        if(i == j) 
        {
          // Ersten Index finden, an dem x[i] vorkommt
          k = 0; 
          while(x_[i] != x_[k]) 
            k++;
          koeffSchema[i*n + j].assign( *(stuetz_[k]) );
        }
        else if(x_[i] == x_[j]) 
        {
          // Fakultät von j-i bestimmen
          l = 1; m = j - i;
          while(m > 0) 
            l *= m--;
          // Wiederum den ersten Index finden
          k = 0; 
          while(x_[i] != x_[k]) 
            k++;
          koeffSchema[i*n + j].assign( *(stuetz_[j - i + k]) );
          koeffSchema[i*n + j] /= l;
        }
        else 
        {
          // Das ist wie im klassischen Fall
          koeffSchema[i*n + j].assign( koeffSchema[i*n + n + j] );  
          koeffSchema[i*n + j] -= koeffSchema[i*n + j - 1];
          koeffSchema[i*n + j] /= (x_[j] - x_[i]);
        }
      }
    }

    for(i=0; i<n; i++) 
      stuetz_[i]->assign( koeffSchema[i] );
  }
  const typename Range::DiscreteFunctionSpaceType &space_;
  const NumFlux *numflux_;
  int q_, q1_;
  double tau_;
  std::vector<double> x_;
  Vector stuetz_;
  Op operator_;
};


template <class Range,class NumFlux>
struct Hermite2 : public Hermite<Range,NumFlux>
{
  typedef Hermite<Range,NumFlux> Base;
  typedef Hermite2<Range,NumFlux> This;
  typedef std::vector<Range> Vector;
  typedef typename Base::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  Hermite2( const DiscreteFunctionSpaceType &space, int q)
  : Base(space,q,2,Base::Op::Approximate::fd3)
  {
  }
  Hermite2( const DiscreteFunctionSpaceType &space)
  : Hermite2(space, Fem::Parameter::getValue<int>("hermite.q"))
  {
  }
  void setup(const NumFlux &numflux)
  {
    Base::numflux_ = &numflux;
    std::cout << "#" << " q=" << q_ << " q1=" << q1_;
    double point = 1.; // tau;
    for (int i=0;i<=q_;++i)
    {
      x_[i] = point;
      if (i%2==1) point -= 1.; // tau;
      std::cout << " " << x_[i];
    }
    std::cout << std::endl;
  }
  void init()
  {
    if (stuetz_.size()==0)
      for (int i=0;i<=q_;++i)
        stuetz_.push_back( new Range("tmp", operator_.func().space() ) );
    int val = 0;
    for (int i=0;i<=q_;++i)
    {
      // assert( u.size() > val );
      stuetz_[i]->assign( *(operator_(i%2,1-val) ) );
      *(stuetz_[i]) *= std::pow(tau_,i%2);
      if (i%2==1) ++val;
    }
    Base::diviDiff( );
  }
  protected:
  using Base::q_;
  using Base::q1_;
  using Base::tau_;
  using Base::x_;
  using Base::stuetz_;
  using Base::operator_;
};


template <class Range,class NumFlux>
struct Hermite<Range,NumFlux>::Op
{
  std::vector<Range*> dtf, dttf;
  // fd3:   finite difference of order 3 at t^n using t^{n+1} and one sided of order 4 at t^{n+1}
  // fd4:   one sided finite difference of order 4 at both t^n and t^{n+1}
  enum class Approximate
  { fd3, fd4, fd6};
  Op(double ptau,int pq, Approximate papprox)
  : dtf(0)
  , dttf(0)
  , tau(ptau)
  , q(pq)
  , approx(papprox)
  , values(0)
  , fvalues(0)
  , steps(0)
  , first(0)
  {}
  bool add(const Range &u, const Range &fu)
  {
    if (dtf.size()==0)
    {
      dtf.push_back(new Range(u));
      dtf.push_back(new Range(u));
      dtf.push_back(new Range(u));
      dttf.push_back(new Range(u));
      dttf.push_back(new Range(u));
      dttf.push_back(new Range(u));
    }
    const int values_needed = std::max(8, (int)(std::floor(q/2.)+1));

    if (fvalues.size()>values_needed)
    {
      Range* newu = values.back();
      newu->assign(u);
      values.push_front(newu);
      values.pop_back();

      Range* newfu = fvalues.back();
      newfu->assign(fu);
      fvalues.push_front(newfu);
      fvalues.pop_back();

      fd_approx_forward();
      return false;
    }
    else
    {
      Range* newu = new Range("tmp",u.space());
      newu->assign(u);
      values.push_front(newu);
      Range* newfu = new Range("ftmp",u.space());
      newfu->assign(fu);
      fvalues.push_front(newfu);
    }
    return true;
  }
  const Range *operator()(int i,int step) 
  {
    int index = 1-step;
    switch (i)
    {
    case 0: return values[index]; break;
    case 1: return fvalues[index]; break;
    case 2: {
            abort();
            Range *dfval_approx = dtf[0];
            switch (approx)
            {
              case Approximate::fd3:
              case Approximate::fd4:
              {
                assert( values.size() > 3);
                if (step == 0)
                {
                  if (approx==Approximate::fd3)
                  {
                    // (f^(n-2) - 6f^(n-1) + 3 f^n + 2 f^(n+1))/(6k)
                    dfval_approx->assign( *(fvalues[3]) );
                    dfval_approx->axpy( -6., *(fvalues[2]) );
                    dfval_approx->axpy(  3., *(fvalues[1]) );
                    dfval_approx->axpy(  2., *(fvalues[0]) );
                    *(dfval_approx) /= 6.*tau;
                  }
                  else
                  {
                    // −25/12 	4 	−3 	4/3 	−1/4
                    dfval_approx->assign( *(fvalues[1]) );
                    *(dfval_approx) *= 25./12.;
                    dfval_approx->axpy( -4., *(fvalues[2]) );
                    dfval_approx->axpy( 3., *(fvalues[3]) );
                    dfval_approx->axpy( -4./3., *(fvalues[4]) );
                    dfval_approx->axpy( 1./4., *(fvalues[5]) );
                    *(dfval_approx) /= tau;
                  }
                }
                else
                {
                  // −25/12 	4 	−3 	4/3 	−1/4
                  dfval_approx->assign( *(fvalues[0]) );
                  *(dfval_approx) *= 25./12.;
                  dfval_approx->axpy( -4., *(fvalues[1]) );
                  dfval_approx->axpy( 3., *(fvalues[2]) );
                  dfval_approx->axpy( -4./3., *(fvalues[3]) );
                  dfval_approx->axpy( 1./4., *(fvalues[4]) );
                  (*dfval_approx) /= tau;
                }
                return dfval_approx;
              }

              case Approximate::fd6:
              { 
                return dtf[index];
              }
            }
            break;
         }
     case 3: {
            return dttf[index];
         }
    default: std::cout << "No approximation of d^4/dt^4 y available" << std::endl; abort();
    }
    return values[index];
  }
  const Range& func() const
  {
    assert(dtf.size()>0);
    return *(dtf[0]);
  }
  private:
  
  void fd_approx_backward()
  {
    /*
    // f' order 6: −49/20 	6 	−15/2 	20/3 	−15/4 	6/5 	−1/6
    dtf[0] = fvalues[0] ;
    dtf[0] *= -1./6.;
    dtf[0].axpy( fvalues[1], 6./5. );
    dtf[0].axpy( fvalues[2], -15./4. );
    dtf[0].axpy( fvalues[3], 20./3. );
    dtf[0].axpy( fvalues[4], -15./2. );
    dtf[0].axpy( fvalues[5], 6. );
    dtf[0].axpy( fvalues[6], -49./20. ); 
    dtf[0] /= tau;
    dtf[1] = fvalues[1] ;
    dtf[1] *= -1./6.;
    dtf[1].axpy( fvalues[2], 6./5. );
    dtf[1].axpy( fvalues[3], -15./4. );
    dtf[1].axpy( fvalues[4], 20./3. );
    dtf[1].axpy( fvalues[5], -15./2. );
    dtf[1].axpy( fvalues[6], 6. );
    dtf[1].axpy( fvalues[7], -49./20. ); // boundary node
    dtf[1] /= tau;
    // f'' order 5: 203/45 	−87/5 	117/4 	−254/9 	33/2 	−27/5 	137/180
    dttf[0]= fvalues[0];
    dttf[0] *= 137./180.;
    dttf[0].axpy( fvalues[1], -27./5. );
    dttf[0].axpy( fvalues[2], 33./2. );
    dttf[0].axpy( fvalues[3], -254./9. );
    dttf[0].axpy( fvalues[4], 117./4. );
    dttf[0].axpy( fvalues[5], -87./5. );
    dttf[0].axpy( fvalues[6], 203./45. );
    dttf[0] /= tau*tau;
    dttf[1]= fvalues[1];
    dttf[1] *= 137./180.;
    dttf[1].axpy( fvalues[2], -27./5. );
    dttf[1].axpy( fvalues[3], 33./2. );
    dttf[1].axpy( fvalues[4], -254./9. );
    dttf[1].axpy( fvalues[5], 117./4. );
    dttf[1].axpy( fvalues[6], -87./5. );
    dttf[1].axpy( fvalues[7], 203./45. );
    dttf[1] /= tau*tau;
    */
  }
  void fd_approx_forward()
  {
    /*
    // f' order 6: −49/20 	6 	−15/2 	20/3 	−15/4 	6/5 	−1/6
    dtf[0] = fvalues[6] ;
    dtf[0] *= -1./6.;
    dtf[0].axpy( fvalues[5], 6./5. );
    dtf[0].axpy( fvalues[4], -15./4. );
    dtf[0].axpy( fvalues[3], 20./3. );
    dtf[0].axpy( fvalues[2], -15./2. );
    dtf[0].axpy( fvalues[1], 6. );
    dtf[0].axpy( fvalues[0], -49./20. ); 
    dtf[0] /= -tau;
    dtf[1] = fvalues[7] ;
    dtf[1] *= -1./6.;
    dtf[1].axpy( fvalues[6], 6./5. );
    dtf[1].axpy( fvalues[5], -15./4. );
    dtf[1].axpy( fvalues[4], 20./3. );
    dtf[1].axpy( fvalues[3], -15./2. );
    dtf[1].axpy( fvalues[2], 6. );
    dtf[1].axpy( fvalues[1], -49./20. ); // boundary node
    dtf[1] /= -tau;
    // f'' order 5: 203/45 	−87/5 	117/4 	−254/9 	33/2 	−27/5 	137/180
    dttf[0]= fvalues[6];
    dttf[0] *= 137./180.;
    dttf[0].axpy( fvalues[5], -27./5. );
    dttf[0].axpy( fvalues[4], 33./2. );
    dttf[0].axpy( fvalues[3], -254./9. );
    dttf[0].axpy( fvalues[2], 117./4. );
    dttf[0].axpy( fvalues[1], -87./5. );
    dttf[0].axpy( fvalues[0], 203./45. );
    dttf[0] /= tau*tau;
    dttf[1]= fvalues[7];
    dttf[1] *= 137./180.;
    dttf[1].axpy( fvalues[6], -27./5. );
    dttf[1].axpy( fvalues[5], 33./2. );
    dttf[1].axpy( fvalues[4], -254./9. );
    dttf[1].axpy( fvalues[3], 117./4. );
    dttf[1].axpy( fvalues[2], -87./5. );
    dttf[1].axpy( fvalues[1], 203./45. );
    dttf[1] /= tau*tau;
    */
  }

  double tau;
  int q;
  Approximate approx;
  std::deque<Range*>  values, fvalues;
  std::deque<double> steps;
  int first;
};

