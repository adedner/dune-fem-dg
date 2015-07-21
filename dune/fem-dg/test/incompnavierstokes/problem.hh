#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cassert>
#include <cmath>

#include <dune/common/array.hh>
#include <dune/fem-dg/models/stokesprobleminterfaces.hh>

namespace Dune
{


  template< class GridImp>
  class NavierStokesProblemDefault
    : public StokesProblemInterface<Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, GridImp :: dimension > ,
			                        	    Dune::Fem::FunctionSpace< double, double, GridImp :: dimension, 1 >  >,
      public EvolutionProblemInterface< Dune::Fem::FunctionSpace< double, double,  GridImp :: dimension, GridImp :: dimension >, false >
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

    NavierStokesProblemDefault()
      : time_( 0 ),
        mu_( 1 ),
        alpha_( 1 ),
        beta_( 1 ),
        theta_( 1.0 - 1.0/M_SQRT2 )
    {}

    double theta () const { return theta_; }

    void setTime( const double time ) const
    {
      time_ = time ;
    }

    //! the right hand side (i.e., the Laplace of u)
    void f ( const DomainType &x, RangeType &ret ) const
    {
      ret = 0;

      const double t = time_;
      const double t2 = t*t;
      const double t5 = t2*t2*t;

      ret[ 0 ] = 3.0*t2*x[1]*x[0] - mu_*2.0*t2*t + t + 2.0*t5*x[0];
      ret[ 1 ] = 2.0*t*x[0] + 1.0 + t5*x[1]*x[1];
    }


    //! the exact solution
    void u ( const DomainType &x, RangeType &ret ) const
    {
      const double t = time_;

      ret[ 0 ] = t*t*t * x[1]*x[1];
      ret[ 1 ] = t*t* x[0];
    }

    virtual void evaluate(const DomainType& arg,
                          const double t, RangeType& res) const
    {
      setTime( t );
      u( arg, res );
    }

    //! the exact solution
    void p(const DomainType& x, PressureRangeType& ret) const
    {
      const double t = time_;
      ret = t * x[ 0 ] + x[ 1 ] - (t + 1.0) * 0.5;
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
    void gradient ( const DomainType &p, JacobianRangeType &grad ) const
    {
    }

    double betaMu() const { return beta_ * mu_; }
    double alphaMu() const { return alpha_ * mu_; }

    using BaseType :: dataPrefix;
  protected:
    mutable double time_;
    const double mu_;
    const double alpha_;
    const double beta_;
    const double theta_;
  };

}

#endif // #ifndef PROBLEM_HH
