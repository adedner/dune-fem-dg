#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cassert>
#include <cmath>

#include <dune/fem-dg/models/stokesprobleminterfaces.hh>

#include "operatorsplitting.hh"

namespace Dune
{
namespace Fem
{

  struct ProblemToEvolutionInterface
  {

    ProblemToEvolutionInterface()
      : time_( 0 ),
        deltaT_( 0 ),
        unique_( 0 )
    {
      static int counter = 0;
      unique_ = counter;
      std::cout << "00000CREATED: " << unique_ << std::endl;
      counter++;
    }

    //set time for stationary problems to fake time provider
    virtual void setTime( const double time ) const
    {
      time_ = time;
    };

    //set time for stationary problems to fake time provider
    virtual double time() const
    {
      return time_;
    };

    //set time step size for mass matrix scaling
    virtual void setDeltaT( const double deltaT ) const
    {
      std::cout << "0000setDeltaT: " << unique_ << std::endl;
      deltaT_ = deltaT;
    };

    double deltaT() const
    {
      std::cout << "0000deltaT: " << unique_ << std::endl;
      return deltaT_;
    }

  protected:
    mutable double time_;
    mutable double deltaT_;
    mutable double unique_;
  };

  // ProblemInterface
  //-----------------
  /**
   * \brief problem interface for a poisson problem
   *
   * \ingroup Problems
   */
  template <class FunctionSpaceImp>
  class ThetaProblemInterface
    : public ProblemInterface< FunctionSpaceImp >,
      public ProblemToEvolutionInterface
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

    static constexpr double theta_ = OperatorSplittingScheme::theta_;
    static constexpr double alpha_ = OperatorSplittingScheme::alpha_;
    static constexpr double beta_  = OperatorSplittingScheme::beta_;

  public:
    ThetaProblemInterface()
      :  mu_( 1 )
    {}

    //! mass factor gamma
    virtual double gamma() const
    {
      return 1.0/ProblemToEvolutionInterface::deltaT();
    }

    virtual void evaluate(const DomainType& arg,
                          const double t, RangeType& res) const
    {
      setTime( t );
      BaseType::u( arg, res );
    }

    double theta () const { return theta_; }
    double beta() const { return beta_; }
    double alpha() const { return alpha_; }
    double mu() const { return mu_; }

  protected:
    mutable double time_;
    const double mu_;
  };

  template< class FunctionSpaceImp, bool constVel = false >
  class ThetaNavierStokesProblemInterface
    : public EvolutionProblemInterface< FunctionSpaceImp, constVel >
  {
    typedef EvolutionProblemInterface< FunctionSpaceImp, constVel > BaseType;

    static constexpr double theta_ = OperatorSplittingScheme::theta_;
    static constexpr double alpha_ = OperatorSplittingScheme::alpha_;
    static constexpr double beta_  = OperatorSplittingScheme::beta_;
  public:

    static const int dimRange  = BaseType::dimRange;
    static const int dimDomain = BaseType::dimDomain;

    typedef typename BaseType::DomainType        DomainType;
    typedef typename BaseType::RangeType         RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::DomainFieldType   DomainFieldType;
    typedef typename BaseType::RangeFieldType    RangeFieldType;

    ThetaNavierStokesProblemInterface()
      : BaseType(),
        mu_( 1 )
    {}

    virtual void evaluate( const DomainType& x,
                           const double t, RangeType& res) const
    {}

    double theta () const { return theta_; }
    double beta() const { return beta_; }
    double alpha() const { return alpha_; }
    double mu() const { return mu_; }

  private:
    const double mu_;

  };

  template< class GridImp>
  class NavierStokesProblemInterfaceTraits
  {
  public:
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, GridImp::dimension > FunctionSpaceType;
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, 1 > PressureFunctionSpaceType;

    typedef ThetaProblemInterface< FunctionSpaceType >            PoissonProblemType;
    typedef ThetaProblemInterface< PressureFunctionSpaceType >    StokesProblemType;
    typedef ThetaNavierStokesProblemInterface< FunctionSpaceType, false > NavierStokesProblemType;

  };


  template< class GridImp >
  class NavierStokesProblemInterface
  {
    typedef NavierStokesProblemInterfaceTraits< GridImp > Traits;
  public:

    typedef typename Traits::FunctionSpaceType             FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType     PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType            PoissonProblemType;
    typedef typename Traits::StokesProblemType             StokesProblemType;
    typedef typename Traits::NavierStokesProblemType       NavierStokesProblemType;

    typedef std::tuple< PoissonProblemType*, StokesProblemType*, NavierStokesProblemType* >         ProblemTupleType;

    /**
     *  \brief constructor constructing a combined problem of the interface sub problems,
     *  i.e. the poisson and the stokes problem.
     *
     *  \note Use the StokesProblem class to create derived objects.
     */
    NavierStokesProblemInterface()
      : problems_( std::make_tuple( new PoissonProblemType(), new StokesProblemType(), new NavierStokesProblemType() ) )
    {}

    NavierStokesProblemInterface( PoissonProblemType* poisson, StokesProblemType* stokes, NavierStokesProblemType* navier )
      : problems_( std::make_tuple( poisson, stokes, navier ) )
    {}

    template< int i >
    const std::remove_pointer_t< std::tuple_element_t<i,ProblemTupleType> >& get() const
    {
      return *(std::get<i>( problems_) );
    }

    template< int i >
    std::remove_pointer_t< std::tuple_element_t<i,ProblemTupleType> >& get()
    {
      return *(std::get<i>( problems_) );
    }

    virtual ~NavierStokesProblemInterface()
    {
      delete std::get<0>(problems_);
      delete std::get<1>(problems_);
      delete std::get<2>(problems_);
    }


    typename NavierStokesProblemType::ExactSolutionType exactSolution( const double time=0.0 ) const
    {
      return std::get<2>( problems_ )->exactSolution(time);
    }
  private:
    mutable ProblemTupleType   problems_;
  };



  /**
   * \brief helper class which helps for the correct (virtual) construction
   * of the problem tuple.
   *
   * \tparam GridImp type of the unterlying grid
   * \tparam StokesProblemImp type of the stokes problem
   *
   * \ingroup Problems
   */
  template< class GridImp,
            template<class> class NavierStokesProblemImp >
  class NavierStokesProblem
    : public NavierStokesProblemInterface< GridImp >
  {
    typedef NavierStokesProblemInterface< GridImp > BaseType;

    typedef typename BaseType::PoissonProblemType                             PoissonProblemBaseType;
    typedef typename BaseType::StokesProblemType                              StokesProblemBaseType;
    typedef typename BaseType::NavierStokesProblemType                        NavierStokesProblemBaseType;

  public:
    typedef typename NavierStokesProblemImp<GridImp>::PoissonProblemType      PoissonProblemType;
    typedef typename NavierStokesProblemImp<GridImp>::StokesProblemType       StokesProblemType;
    typedef typename NavierStokesProblemImp<GridImp>::NavierStokesProblemType NavierStokesProblemType;

    NavierStokesProblem()
      : BaseType( new PoissonProblemType(), new StokesProblemType(), new NavierStokesProblemType() )
    {}

  };


  template< class GridImp>
  class NavierStokesProblemDefault
  {
    typedef NavierStokesProblemInterfaceTraits< GridImp >  Traits;

    typedef typename Traits::FunctionSpaceType             FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType     PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType            PoissonProblemBaseType;
    typedef typename Traits::StokesProblemType             StokesProblemBaseType;
    typedef typename Traits::NavierStokesProblemType       NavierStokesProblemBaseType;

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

        ret = { 3.0*t2*x[1]*x[0] - mu_*2.0*t2*t + t + 2.0*t5*x[0], 2.0*t*x[0] + 1.0 + t5*x[1]*x[1] };
      }


      //! the exact solution
      void u ( const DomainType &x, RangeType &ret ) const
      {
        const double t = time_;

        ret = {t*t*t * x[1]*x[1], t*t* x[0] };
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
        ret = t * x[0] + x[1] - (t + 1.0) * 0.5;
      }
    };

    class NavierStokesProblem
      : public NavierStokesProblemBaseType
    {
    public:

      static const int dimRange  = StokesProblemBaseType::dimRange;
      static const int dimDomain = StokesProblemBaseType::dimDomain;

      typedef typename NavierStokesProblemBaseType::DomainType        DomainType;
      typedef typename NavierStokesProblemBaseType::RangeType         RangeType;
      typedef typename NavierStokesProblemBaseType::JacobianRangeType JacobianRangeType;
      typedef typename NavierStokesProblemBaseType::DomainFieldType   DomainFieldType;
      typedef typename NavierStokesProblemBaseType::RangeFieldType    RangeFieldType;

      virtual void evaluate(const DomainType& x,
                            const double t, RangeType& res) const
      {
        //todo: same as poisson problem -> improve...
        res = { t*t*t * x[1]*x[1], t*t* x[0]};
      }
    };

 public:
    typedef PoissonProblem       PoissonProblemType;
    typedef StokesProblem        StokesProblemType;
    typedef NavierStokesProblem  NavierStokesProblemType;
  };

}
}
#endif // #ifndef PROBLEM_HH
