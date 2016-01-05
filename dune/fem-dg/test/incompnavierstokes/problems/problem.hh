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

  struct OperatorSplittingScheme
  {
    static constexpr double theta_ = 1.0 - 1.0/M_SQRT2 ;
    static constexpr double alpha_ = (1.0 - 2.0 * theta_ )/(1.0 - theta_);
    static constexpr double beta_  = theta_/ ( 1.0 - theta_ );

    static const bool hasAdvection = true;
    static const bool hasDiffusion = true;
    static const bool hasSource    = true;
    static const bool hasGrad      = true;
  };

  template< int Step, bool RightHandSide >
  class FractionalStepThetaScheme;

  //step 0
  template<>
  class FractionalStepThetaScheme<0,false> //new timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return 1.0; }
    static constexpr double diffusion(){ return Split::alpha_; }
    static constexpr double advection(){ return 0.0; }
    static constexpr double grad(){ return 1.0; }
    static constexpr double mass(){ return Split::theta_;}

    static const bool hasAdvection = false;

    template< class DF >
    static void velocity( const DF& un, DF& velocity  ){ /*no velocity needed*/ }
  };
  template<>
  class FractionalStepThetaScheme<0,true> //old timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return 0.0;}
    static constexpr double diffusion(){ return Split::beta_;}
    static constexpr double advection(){ return 1.0;}
    static constexpr double grad(){ return 0.0;}
    static constexpr double mass(){ return Split::theta_;}

    static const bool hasSource    = false;
    static const bool hasGrad      = false;

    template< class DF >
    static void velocity( const DF& un, DF& velocity  )
    {
      velocity.assign( un );
    }
  };


  //step 1
  template<>
  class FractionalStepThetaScheme<1,false> //new timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return 1.0;}
    static constexpr double diffusion(){ return Split::beta_;}
    static constexpr double advection(){ return 1.0;}
    static constexpr double grad(){ return 0.0;}
    static constexpr double mass(){ return 1.0 - 2.0 * Split::theta_;}

    static const bool hasGrad      = false;

    template< class DF >
    static void velocity( const DF& un, const DF& untheta, DF& velocity )
    {
      velocity.assign( un );
      velocity *= ( 2.0*theta_ - 1.0 ) / theta_ ;
      velocity.axpy( (1.0-theta_)/theta_, untheta );
    }
  };
  template<>
  class FractionalStepThetaScheme<1,true> //old timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return  0.0;}
    static constexpr double diffusion(){ return  Split::alpha_;}
    static constexpr double advection(){ return  0.0;}
    static constexpr double grad(){ return  1.0;}
    static constexpr double mass(){ return 1.0 - 2.0 * Split::theta_;}

    static const bool hasAdvection = false;
    static const bool hasSource    = false;

    template< class DF >
    static void velocity( const DF& un, const DF& untheta, DF& velocity ){ /*no velocity needed*/ }
  };

  //step 2
  template<>
  class FractionalStepThetaScheme<2,false> //new timestep
    : public FractionalStepThetaScheme<0,0>
  {};
  template<>
  class FractionalStepThetaScheme<2,true> //old timestep
    : public FractionalStepThetaScheme<0,1>
  {};



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
    static constexpr double theta_ = OperatorSplittingScheme::theta_;
    static constexpr double alpha_ = OperatorSplittingScheme::alpha_;
    static constexpr double beta_  = OperatorSplittingScheme::beta_;

    ThetaProblemInterface()
      : time_( 0 ),
        mu_( 1 ),
        deltaT_( 0 )
    {}

    //set time for stationary problems to fake time provider
    virtual void setTime( const double time ) const
    {
      time_ = time;
    };

    //set time step size for mass matrix scaling
    virtual void setDeltaT( const double deltaT ) const
    {
      deltaT_ = deltaT;
    };

    //! mass factor gamma
    virtual double gamma() const
    {
      return 1.0/deltaT_;
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
    double deltaT() const { return deltaT_; }

  protected:
    mutable double time_;
    mutable double deltaT_;
    const double mu_;
  };

  template< class FunctionSpaceImp, bool constVel = false >
  class ThetaNavierStokesProblemInterface
    : public EvolutionProblemInterface< FunctionSpaceImp, constVel >
  {
    typedef EvolutionProblemInterface< FunctionSpaceImp, constVel > BaseType;
  public:

    static constexpr double theta_ = OperatorSplittingScheme::theta_;
    static constexpr double alpha_ = OperatorSplittingScheme::alpha_;
    static constexpr double beta_  = OperatorSplittingScheme::beta_;


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
  class NavierStokesProblemInterfaceBase
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
    typedef NavierStokesProblemInterfaceBase< GridImp >        BaseType;
  public:

    typedef typename BaseType::FunctionSpaceType         FunctionSpaceType;
    typedef typename BaseType::PressureFunctionSpaceType PressureFunctionSpaceType;

    typedef typename BaseType::PoissonProblemType        PoissonProblemType;
    typedef typename BaseType::StokesProblemType         StokesProblemType;
    typedef typename BaseType::NavierStokesProblemType   NavierStokesProblemType;

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

    template< int i >
    const typename std::remove_pointer< typename std::tuple_element<i,ProblemTupleType>::type >::type& get() const
    {
      return *(std::get<i>( problems_) );
    }

    template< int i >
    typename std::remove_pointer< typename std::tuple_element<i,ProblemTupleType>::type >::type& get()
    {
      return *(std::get<i>( problems_) );
    }

    // return prefix for data loops
    virtual std::string dataPrefix() const
    {
      return get<0>().dataPrefix();
    }

  protected:
    template< class PoissonProblemPtrImp, class StokesProblemPtrImp, class NavierStokesProblemPtrImp >
    void create( const PoissonProblemPtrImp& poisson, const StokesProblemPtrImp& stokes,
                 const NavierStokesProblemPtrImp& navier )
    {
      std::get<0>( problems_ ) = poisson;
      std::get<1>( problems_ ) = stokes;
      std::get<2>( problems_ ) = navier;
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
      : BaseType(),
        poisson_( new PoissonProblemType() ),
        stokes_( new StokesProblemType() ),
        navier_( new NavierStokesProblemType() )
    {
      BaseType::create( poisson_, stokes_, navier_ );
    }

  private:
    mutable PoissonProblemBaseType* poisson_;
    mutable StokesProblemBaseType*  stokes_;
    mutable NavierStokesProblemBaseType* navier_;

  };


  template< class GridImp>
  class NavierStokesProblemDefault
    : public NavierStokesProblemInterfaceBase< GridImp >
  {
    typedef NavierStokesProblemInterfaceBase< GridImp >    BaseType;

    typedef typename BaseType::FunctionSpaceType           FunctionSpaceType;
    typedef typename BaseType::PressureFunctionSpaceType   PressureFunctionSpaceType;

    typedef typename BaseType::PoissonProblemType          PoissonProblemBaseType;
    typedef typename BaseType::StokesProblemType           StokesProblemBaseType;
    typedef typename BaseType::NavierStokesProblemType     NavierStokesProblemBaseType;

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
        ret = t * x[ 0 ] + x[ 1 ] - (t + 1.0) * 0.5;
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
        res[ 0 ] = t*t*t * x[1]*x[1];
        res[ 1 ] = t*t* x[0];
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
