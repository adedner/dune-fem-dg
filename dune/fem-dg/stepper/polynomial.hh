#include <set>
#include <vector>
#include <iostream>

struct Polynomial
{
  std::vector<double> x_;
  int n_;
  Polynomial(int n, const std::vector<double> &x) : x_(x), n_(n) {}
  double evaluate(double point, std::set<int> ignore={})
  {
    // if (ignore.size()==n_)
    //   return 1;
    double phi = 1;
    for (int j=0;j<n_;++j)
    {
      if ( ignore.find(j) == ignore.end() )
      {
        if (std::abs(point-x_[j])<1e-10)
          phi = 0;
        else
          phi *= (point-x_[j]);
      }
    }
    return phi;
  }
  double derivative(int k,double point, std::set<int> ignore={}, int newignore=-1)
  {
    if (newignore >= 0)
      if ( ! ignore.insert(newignore).second )
        return 0;
    if (k==0)
      return evaluate(point,ignore);
    else
    {
      double ret=0;
      for (int i=0;i<n_;++i)
        ret += derivative(k-1,point,ignore,i);
      return ret;
    }
  }
  double evaluate_end(std::set<int> ignore={})
  {
    // if (ignore.size()==n_)
    //   return 1;
    double phi = 1;
    for (int j=0;j<n_;++j)
    {
      if ( ignore.find(j) == ignore.end() )
        phi *= (1.-x_[j]);
    }
    return phi;
  }
  double derivative_end(int k,std::set<int> ignore={}, int newignore=-1)
  {
    if (newignore >= 0)
      if ( ! ignore.insert(newignore).second )
        return 0;
    if (k==0)
      return evaluate_end(ignore);
    else
    {
      double ret=0;
      for (int i=0;i<n_;++i)
        ret += derivative_end(k-1,ignore,i);
      return ret;
    }
  }
};

#if 0
int main()
{
  std::vector<double> x={-1,0,1};
  Polynomial p(3,x);
  std::cout << "val  = " << p.evaluate(0.5) << std::endl;
  std::cout << "d/dt = " << p.derivative(1,0.5) << std::endl;
  std::cout << "d^2/dt^2 = " << p.derivative(2,0.5) << std::endl;
  std::cout << "d^3/dt^3 = " << p.derivative(3,0.5) << std::endl;
  std::cout << "d^4/dt^4 = " << p.derivative(4,0.5) << std::endl;
}
#endif
