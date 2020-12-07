#include <cmath>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/space/shapefunctionset/orthonormal.hh>
#include <dune/fem-dg/operator/limiter/indicatorbase.hh>

template <class LF>
double smoothnessIndicator(const LF& uLocal)
{
  using ONB = Dune::Fem::OrthonormalShapeFunctions<LF::dimDomain>;
  const std::size_t R = LF::dimRange; // will be using density
  const std::size_t P = uLocal.order();
  assert(uLocal.size()/double(R) == ONB::size(P));
  // actually some inverse mass matrix needed here
  const double area   = uLocal.entity().geometry().volume();
  const double factor = 1/std::sqrt( area );
  // 1D monoments vector
  double q[P+1];
  double b2[P+1];
  double f = 0;
  std::size_t k = ONB::size(0);
  q[0] = uLocal[0]*uLocal[0]; // constant part not used
  double l2norm2 = 0; // q[0];
  for (std::size_t i=1;i<=P;++i)
  {
    q[i]  = 0;
    b2[i] = 0;
    double nofMoments = ONB::size(i)-k;
    assert(nofMoments == 1.);
    for (;k<ONB::size(i);++k)
    {
      q[i]    += uLocal[k*R]*uLocal[k*R] / nofMoments;
      l2norm2 += uLocal[k*R]*uLocal[k*R] / nofMoments;
      b2[i]   += pow(1/double(i),2*P) / nofMoments;
      f       += pow(1/double(i),2*P) / nofMoments;
    }
  }
#if 0
  for (std::size_t i=0;i<=P;++i)
    q[i] = std::sqrt( q[i] ) / factor;
#endif
  for (std::size_t i=1;i<=P;++i)
  {
    /*
    std::cout << "====: " << i << " " << q[i] << " "
              << l2norm2 << " " << b2[i]/f << "     " << l2norm2*b2[i]/f
              << " -> " << std::sqrt( q[i] + l2norm2*b2[i]/f ) / factor
              << std::endl; */
    q[i] = std::sqrt( q[i] + l2norm2*b2[i]/f ) / factor;
  }
  double maxQ = std::max( q[P], q[P-1] );
  // find first 'significant' mode
  std::size_t significant = 0;
  for (std::size_t i=P;i>=1;--i)
  {
    maxQ = std::max(maxQ, q[i]);
    if (maxQ>1e-14)
    {
      significant = i;
      break;
    }
  }
  if (significant==0) // constant, i.e., very smooth indeed
    return 100;
  if (significant==1) // not quite clear if this needs to be fixed -
    // solution is linear and we do not have enough info to fit
    return 10;

  Dune::DynamicMatrix<double> matrix(significant,2);
  Dune::DynamicVector<double> rhs(significant);
  for (std::size_t r=significant; r-->0; )
  {
    maxQ = std::max(maxQ, q[r+1]);
    rhs[r]       = std::log( maxQ );
    matrix[r][0] = 1;
    matrix[r][1] = -std::log(double(r+1));
  }
  Dune::FieldMatrix<double,2,2> A;
  Dune::FieldVector<double,2> b;
  for (std::size_t r=0;r<2;++r)
  {
    for (std::size_t c=0;c<2;++c)
    {
      A[r][c] = 0;
      for (std::size_t k=0;k<significant;++k)
        A[r][c] += matrix[k][r]*matrix[k][c];
    }
    b[r] = 0;
    for (std::size_t k=0;k<significant;++k)
      b[r] += matrix[k][r]*rhs[k];
  }
  Dune::FieldVector<double,2> x;
  A.solve(x,b);
  double s = x[1];
  assert(s==s);
  // std::cout << "       indicator= " << s << std::endl;
  return s;
}
template <class DiscreteFunction>
struct ModalSmoothnessIndicator
: public Dune::Fem::TroubledCellIndicatorBase<DiscreteFunction>
{
  typedef typename  Dune::Fem::TroubledCellIndicatorBase<DiscreteFunction>::LocalFunctionType
    LocalFunctionType;
  virtual ~ModalSmoothnessIndicator() {}

  double operator()( const DiscreteFunction& U, const LocalFunctionType& uEn) const override
  {
    double modalInd = smoothnessIndicator( uEn );
    if( std::abs( modalInd ) > 1e-14 )
      return 1.0 / modalInd ;
    else
      return 0.0;
  }
};
