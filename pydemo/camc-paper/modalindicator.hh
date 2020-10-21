#pragma once

#include <cmath>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/space/shapefunctionset/orthonormal.hh>

#include <dune/fem-dg/operator/limiter/indicatorbase.hh>

template <class Model, class DiscreteFunction>
struct ModalIndicator
: public Dune::Fem::TroubledCellIndicatorBase<DiscreteFunction>
{
  typedef Model ModelType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DFS;
  typedef typename DFS::template ToNewDimRange<1>::Type ScalarDFS;
  typedef typename DFS::GridPartType GridPartType;

  ModalIndicator (const ModelType& model, const double factor, const DFS &dfs,
                  const std::string &component)
  : model_(model), factor_(factor),
    scalarDF_(dfs.gridPart()), pEn_(scalarDF_),
    comp_(0)
  {
    if (component == "density") comp_ = 1;
    if (component == "entropy") comp_ = 2;
    if (component == "pressure") comp_ = 3;
  }

  template <class LF>
  double rho( const LF& uEn) const
  {
    return smoothnessIndicator( uEn );
  }
  template <class LF>
  double pressure( const LF& uEn) const
  {
    model_.init(uEn.entity());
    pEn_.init(uEn.entity());
    auto interpolate = scalarDF_.interpolation(uEn.entity());
    interpolate(TCIValue(model_,scalarDF_.gridPart(),uEn), pEn_);
    return smoothnessIndicator( pEn_ );
  }
  template <class LF>
  double entropy( const LF& uEn) const
  {
    model_.init(uEn.entity());
    pEn_.init(uEn.entity());
    auto interpolate = scalarDF_.interpolation(uEn.entity());
    interpolate(TCIValue(model_,scalarDF_.gridPart(),uEn,false), pEn_);
    return smoothnessIndicator( pEn_ );
  }
  double operator()( const DiscreteFunction& U, const LocalFunctionType& uEn) const override
  {
    double modalInd = 0;
    if (comp_==1)
      modalInd = rho(uEn);
    else if (comp_==2)
      modalInd = entropy(uEn);
    else if (comp_==3)
      modalInd = pressure(uEn);
    if( std::abs( modalInd ) > 1e-14 )
      return factor_ / modalInd ;
    else
      return 0.0;
  }
private:
  struct TCIValue : public Dune::Fem::BindableGridFunction< GridPartType, Dune::Dim<1> >
  {
    typedef Dune::Fem::BindableGridFunction<GridPartType, Dune::Dim<1> > Base;
    TCIValue(const ModelType &model, const GridPartType &gridPart,
             const LocalFunctionType &uEn, bool pressure=true)
    : Base(gridPart), model_(model), uEn_(uEn), pressure_(pressure)
    {Base::bind(uEn.entity());}
    template <class Quadrature>
    void evaluateQuadrature(const Quadrature &quadrature, std::vector<typename Base::RangeType> &values) const
    {
      typename LocalFunctionType::RangeType u;
      for(auto qp : quadrature )
      {
        uEn_.evaluate(qp, u);
        values[qp.index()][0] = 0.3999999999999999 * ( u[ 3 ] -
                   (( u[ 0 ] * (u[ 1 ] / u[ 0 ] * u[ 1 ] / u[ 0 ] +
                                u[ 2 ] / u[ 0 ] * u[ 2 ] / u[ 0 ] )) / 2));
        if (!pressure_)
          values[qp.index()][0] = u[0] * std::log(values[qp.index()][0]/pow(u[0],1.4));
      }
    }
    int order() const {return uEn_.order();}
    const Model &model_;
    const LocalFunctionType &uEn_;
    bool pressure_;
  };

  const Model &model_;
  double factor_;
  ScalarDFS scalarDF_;
  mutable Dune::Fem::TemporaryLocalFunction<ScalarDFS> pEn_;
  int comp_;

  template <class LF>
  static double smoothnessIndicator(const LF& uLocal)
  {
    using ONB = Dune::Fem::OrthonormalShapeFunctions<LF::dimDomain>;
    const std::size_t R = LF::dimRange; // will be using first component if vector
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
      for (;k<ONB::size(i);++k)
      {
        q[i]    += uLocal[k*R]*uLocal[k*R] / nofMoments;
        l2norm2 += uLocal[k*R]*uLocal[k*R] / nofMoments;
        b2[i]   += pow(1/double(i),2*P) / nofMoments;
        f       += pow(1/double(i),2*P) / nofMoments;
      }
    }
    for (std::size_t i=1;i<=P;++i)
      q[i] = std::sqrt( q[i] + l2norm2*b2[i]/f ) / factor;
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
      return 1000;
    if (significant==1) // not quite clear if this needs to be fixed -
      // solution is linear and we do not have enough info to fit
      return 100;

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
    return s;
  }
};

template <class GridView, class GF, class Model>
auto smoothnessGF(const Model &model, const GF &gf) {
  typedef ModalIndicator<Model,GF> Indicator;
  return [&model, &gf, lgf=gf.localFunction()]
  (const auto& en,const auto& xLocal) mutable -> auto {
    Indicator indicator(model,1,gf.space());
    lgf.bind(en);
    return Dune::FieldVector<double,3>{indicator.rho(lgf),
                                       indicator.pressure(lgf),
                                       indicator.entropy(lgf)};
  };
}
