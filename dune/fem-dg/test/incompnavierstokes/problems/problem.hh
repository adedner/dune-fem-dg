#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cassert>
#include <cmath>

#include <dune/common/array.hh>
#include <dune/fem-dg/models/stokesprobleminterfaces.hh>

namespace Dune
{
namespace Fem
{

  // ProblemInterface
  //-----------------
  /**
   * \brief problem interface for a poisson problem
   *
   * \ingroup Problems
   */
  template <class FunctionSpaceImp>
  class ThetaProblemInterface
    : public ProblemInterface< FunctionSpaceImp >
  {
    typedef ProblemInterface< FunctionSpaceImp >                  BaseType;
  public:
    typedef FunctionSpaceImp                                      FunctionSpaceType;
    typedef ThetaProblemInterface< FunctionSpaceType >            ThisType;

    enum { dimDomain = FunctionSpaceType :: dimDomain };
    enum { dimRange  = FunctionSpaceType :: dimRange  };

    typedef typename FunctionSpaceType :: DomainType              DomainType;
    typedef typename FunctionSpaceType :: RangeType               RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType       JacobianRangeType;
    typedef typename FunctionSpaceType :: DomainFieldType         DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType          RangeFieldType;

    typedef FieldMatrix< RangeFieldType, dimDomain, dimDomain >   DiffusionMatrixType;

  public:
    static constexpr double theta_ = 1.0 - 1.0/M_SQRT2 ;
    static constexpr double alpha_ = (1.0 - 2.0 * theta_ )/(1.0 - theta_);
    static constexpr double beta_  = theta_/ ( 1.0 - theta_ );

    ThetaProblemInterface()
      : time_( 0 ),
        mu_( 1 )
    {}

    virtual void setTime( const double time ) const
    {
      time_ = time;
    };

    virtual void evaluate(const DomainType& arg,
                          const double t, RangeType& res) const
    {
      setTime( t );
      BaseType::u( arg, res );
    }

    virtual double theta () const { return theta_; }
    virtual double betaMu() const { return beta_ * mu_; }
    virtual double alphaMu() const { return alpha_ * mu_; }

  protected:
    mutable double time_;
    const double mu_;
  };


  template< class GridImp>
  class NavierStokesProblemDefault
    : public StokesProblemInterface< ThetaProblemInterface< Dune::Fem::FunctionSpace< double, double, GridImp::dimension, GridImp::dimension > > ,
			                        	     ThetaProblemInterface< Dune::Fem::FunctionSpace< double, double, GridImp::dimension, 1 > > >,
      public EvolutionProblemInterface< Dune::Fem::FunctionSpace< double, double, GridImp::dimension, GridImp::dimension >, false >
  {
  public:
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, GridImp::dimension > FunctionSpaceType;
  private:
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, 1 > PressureFunctionSpaceType;

    typedef ThetaProblemInterface< FunctionSpaceType >         PoissonProblemBaseType;
    typedef ThetaProblemInterface< PressureFunctionSpaceType > StokesProblemBaseType;

    typedef StokesProblemInterface< PoissonProblemBaseType, StokesProblemBaseType > BaseType;
    typedef EvolutionProblemInterface< FunctionSpaceType, false >                   EvolBaseType;


  public:

    class PoissonProblem
      : public PoissonProblemBaseType
    {
    public:
      static const int dimRange  = PoissonProblemBaseType::dimRange;
      static const int dimDomain = PoissonProblemBaseType::dimDomain;

      typedef typename PoissonProblemBaseType::DomainType          DomainType;
      typedef typename PoissonProblemBaseType::RangeType           RangeType;
      typedef typename PoissonProblemBaseType::JacobianRangeType   JacobianRangeType;
      typedef typename PoissonProblemBaseType::DomainFieldType     DomainFieldType;
      typedef typename PoissonProblemBaseType::RangeFieldType      RangeFieldType;

      typedef typename PoissonProblemBaseType::DiffusionMatrixType DiffusionMatrixType;

      using PoissonProblemBaseType::time_;
      using PoissonProblemBaseType::mu_;

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
      }
    };

    class StokesProblem
      : public StokesProblemBaseType
    {
    public:
      static const int dimRange  = StokesProblemBaseType::dimRange;
      static const int dimDomain = StokesProblemBaseType::dimDomain;

      typedef typename StokesProblemBaseType::DomainType        DomainType;
      typedef typename StokesProblemBaseType::RangeType         RangeType;
      typedef typename StokesProblemBaseType::JacobianRangeType JacobianRangeType;
      typedef typename StokesProblemBaseType::DomainFieldType   DomainFieldType;
      typedef typename StokesProblemBaseType::RangeFieldType    RangeFieldType;

      typedef typename StokesProblemBaseType::DiffusionMatrixType DiffusionMatrixType;

      using StokesProblemBaseType::time_;
      using StokesProblemBaseType::mu_;

      //! the exact solution
      void u(const DomainType& x, RangeType& ret) const
      {
        const double t = time_;
        ret = t * x[ 0 ] + x[ 1 ] - (t + 1.0) * 0.5;
      }
    };

    using EvolBaseType::evaluate;

    static const int dimRange  = FunctionSpaceType::dimRange;
    static const int dimDomain = FunctionSpaceType::dimDomain;

    typedef typename FunctionSpaceType::DomainType          DomainType;
    typedef typename FunctionSpaceType::RangeType           RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType   JacobianRangeType;
    typedef typename FunctionSpaceType::DomainFieldType     DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType      RangeFieldType;


    NavierStokesProblemDefault()
      : BaseType( std::make_tuple( PoissonProblem(), StokesProblem() ) ),
        mu_( 1 )
    {}

    virtual void evaluate(const DomainType& x,
                          const double t, RangeType& res) const
    {
      //todo: same as poisson problem -> improve...
      res[ 0 ] = t*t*t * x[1]*x[1];
      res[ 1 ] = t*t* x[0];
    }


    double theta () const { return PoissonProblemBaseType::theta_; }
    double betaMu() const { return PoissonProblemBaseType::beta_ * mu_; }
    double alphaMu() const { return PoissonProblemBaseType::alpha_ * mu_; }

    using EvolBaseType::dataPrefix;

  protected:
    const double mu_;

  };

}
}
#endif // #ifndef PROBLEM_HH
