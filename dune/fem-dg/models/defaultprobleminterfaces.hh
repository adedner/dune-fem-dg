#ifndef DUNE_DEFAULTPROBLEMINTERFACE_HH
#define DUNE_DEFAULTPROBLEMINTERFACE_HH

#include <dune/common/version.hh>
#include <dune/common/exceptions.hh>

#include <dune/fem/function/common/instationary.hh>
#include <dune/fem/misc/gridsolution.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Interface class for
   * initial and exact solution of the advection-diffusion model
   *
   * \ingroup Problems
   */
  template< class FunctionSpaceImp, bool constantVelocity = false >
  class EvolutionProblemInterface
  {
    typedef EvolutionProblemInterface< FunctionSpaceImp,
                                       constantVelocity >            ThisType;

  public:
    typedef FunctionSpaceImp                                         FunctionSpaceType;

    enum { ConstantVelocity = constantVelocity };
    enum { dimDomain = FunctionSpaceType::dimDomain };
    enum { dimRange  = FunctionSpaceType::dimRange  };

    typedef typename FunctionSpaceType::DomainType                   DomainType;
    typedef typename FunctionSpaceType::RangeType                    RangeType;
    typedef typename FunctionSpaceType::DomainFieldType              DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType               RangeFieldType;
    typedef typename FunctionSpaceType::JacobianRangeType            JacobianRangeType;

    typedef Fem::Parameter ParameterType;

  protected:
    /**
     * \brief define problem parameters
     */
    EvolutionProblemInterface()
    {}

  public:
    typedef Fem::InstationaryFunction< ThisType, Fem::__InstationaryFunction::HoldReference > InstationaryFunctionType;

    typedef InstationaryFunctionType ExactSolutionType;


    //! turn timedependent function into function by fixing time
    InstationaryFunctionType fixedTimeFunction( const double time ) const
    {
      return InstationaryFunctionType( *this, time );
    }

    ExactSolutionType exactSolution( const double time ) const
    {
      return InstationaryFunctionType( *this, time );
    }

    //! destructor
    virtual ~EvolutionProblemInterface() {}

    //! return true is source term is available
    virtual inline bool hasStiffSource() const { return true ; }
    virtual inline bool hasNonStiffSource() const { return false; }

    //! return true if mass term is not the identity
    virtual inline bool hasMass() const { return false; }

    //! stiff source term
    virtual inline double stiffSource (const DomainType& arg,
                                       const double time,
                                       const RangeType& u,
                                       RangeType& res) const
    {
      res = 0;
      return 0.0;
    }

    //! non stiff source term
    virtual inline double nonStiffSource (const DomainType& arg,
                                          const double time,
                                          const RangeType& u,
                                          RangeType& res) const
    {
      res = 0;
      return 0.0;
    }

    //! diagonal of the mass term
    virtual inline void mass (const DomainType& arg,
                              const double time,
                              const RangeType& u,
                              RangeType& diag ) const
    {
      assert( hasMass() );
      diag = 1;
    }

    /** \brief decide if the refinement/coarsening should be allowed in certain regions of the
     *    computational domain
     *  \param[in] x Global coordinates
     *  \return whether or not the refinement is allowed
     */
    bool allowsRefinement( const DomainType& x ) const
    {
      // by default refinement is allowed in the whole computational domain
      return true;
    }

    //! use both indicator for a grid adaptation
    inline bool twoIndicators() const
    {
      return false;
    }

    /** \brief return a first indicator for an adaptation of the grid
     *  \param[in] u Value of the numerical solution in a point
     *
     *  \note NaN does no adaptation
     *
     *  \return scalar indicator of importance for a grid adaptation
     *    like pot. temperature, density etc
     */
    inline double indicator1( const DomainType& x, const RangeType& u ) const
    {
      return std::numeric_limits<float>::quiet_NaN();
    }

    /** \brief return a second indicator for an adaptation of the grid
     *  \param[in] u Value of the numerical solution in a point
     *
     *  \note NaN does no adaptation
     *
     *  \return scalar indicator of importance for a grid adaptation
     *    like pot. temperature, density etc
     */
    inline double indicator2( const DomainType& x, const RangeType& u ) const
    {
      return std::numeric_limits<float>::quiet_NaN();
    }


    //! return diffusion coefficient (default returns epsilon)
    virtual inline double diffusion ( const RangeType& u, const JacobianRangeType& gradU ) const
    {
      return epsilon();
    }

    /** \brief start time of problem */
    virtual double startTime() const { return 0.0; }

    /** \brief diffusion coefficient of problem */
    virtual double epsilon() const { return 0.0; }

    /** \brief problem velocity (for advection-diffusion) */
    virtual void velocity(const DomainType& x, const double, DomainType& v) const
    {
      v = 0;
    }

    /**
     * \brief old version of the exact solution
     *
     * old version of evaluate(const DomainType& arg, double t, RangeType& res),
     * which is still needed by the DataWriter
     */
    virtual inline void evaluate(const double t,
                                 const DomainType& arg, RangeType& res) const
    {
      evaluate(arg, t, res);
    }

    /**
     * \brief evaluate exact solution, to be implemented in derived classes
     */
    virtual void evaluate(const DomainType& arg,
                          const double t, RangeType& res) const = 0 ;


    /**
     * \brief evaluate exact solution, to be implemented in derived classes
     */
    virtual double boundaryFlux(const DomainType& arg,
                                const double t,
                                const RangeType& u, RangeType& flux) const
    {
      flux = 0;
      return 0.0;
    }


    /**
     * \brief latex output for EocOutput, default is empty
     */
    virtual std::string description() const
    {
      return std::string("");
    }

  };



  // ProblemInterface
  //-----------------
  /**
   * \brief Interface class for a Poisson problem
   *
   * \ingroup Problems
   */
  template <class FunctionSpaceImp>
  class ProblemInterface
  {
  public:
    typedef FunctionSpaceImp                                         FunctionSpaceType;
    typedef ProblemInterface< FunctionSpaceType >                    ThisType;

    enum { dimDomain = FunctionSpaceType::dimDomain };
    enum { dimRange  = FunctionSpaceType::dimRange  };

    typedef typename FunctionSpaceType::DomainType                   DomainType;
    typedef typename FunctionSpaceType::RangeType                    RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType            JacobianRangeType;
    typedef typename FunctionSpaceType::DomainFieldType              DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType               RangeFieldType;

    typedef FieldMatrix< RangeFieldType, dimDomain, dimDomain >      DiffusionMatrixType;

  public:

    //! destructor
    virtual ~ProblemInterface() {}

    virtual std::string name () const { return "ProblemInterface"; }

    //! evaluates the right hand side (Function like bahavior)
    inline void evaluate(const DomainType& x, RangeType& ret) const
    {
      u(x , ret);
    }

    //! advection factor a if available
    virtual double constantAdvection() const { return 0.0; }

    //! the right hand side
    virtual void f(const DomainType& x, RangeType& ret) const { std::abort(); }// = 0;

    //! the exact solution
    virtual void u(const DomainType& x, RangeType& ret) const { std::abort(); }//= 0;

    //! the diffusion matrix
    virtual void K(const DomainType& x, DiffusionMatrixType& m) const { std::abort(); }//= 0;

    //! returns true if diffusion coefficient is constant
    virtual bool constantK() const { std::abort(); return false; }//= 0;

    //! the Dirichlet boundary data function
    virtual void g(const DomainType& x, RangeType& ret) const
    {
      u( x, ret );
    }

    //! mass factor gamma
    virtual double gamma() const { return 0.0; }


    //! the Neumann boundary data function
    virtual void psi(const DomainType& x,
                     JacobianRangeType& dn) const
    {
      abort();
      // gradient( x, dn );
    }

    /**
     * \brief getter for the velocity
     */
    virtual void velocity(const DomainType& x, DomainType& v) const
    {
      v = 0;
    }

    //! the gradient of the exact solution (default is empty)
    virtual void gradient(const DomainType& x,
                          JacobianRangeType& grad) const
    {
      assert( false );
      DUNE_THROW( NotImplemented, "ProblemInterface::gradient not overloaded but called!" );
    }

    //! return whether boundary is Dirichlet (true) or Neumann (false)
    virtual bool dirichletBoundary(const int bndId, const DomainType& x) const
    {
      return true;
    }

    /**
     * \brief latex output for EocOutput, default is empty
     */
    virtual std::string description() const
    {
      return std::string("");
    }

    ProblemInterface() : exactSolution_( *this ) {}

  protected:
    //! the exact solution to the problem for EOC calculation
    class ExactSolution
    : public Fem:: Function< FunctionSpaceType, ExactSolution >
    {
    private:
      typedef Fem:: Function< FunctionSpaceType, ExactSolution >      BaseType;

      typedef ProblemInterface< FunctionSpaceType>   DataType;
    protected:
      FunctionSpaceType  functionSpace_;
      const DataType &data_;

    public:
      inline ExactSolution ( const ThisType& data )
      : BaseType( ),
        functionSpace_(),
        data_( data )
      {
      }

      inline void evaluate ( const DomainType &x, RangeType &ret ) const
      {
        data_.u( x, ret );
      }

      inline void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
      {
        data_.gradient( x, ret );
      }

      inline void evaluate (const DomainType &x,
                            const double time, RangeType &phi ) const
      {
        evaluate( x, phi );
      }
    }; // end class ExactSolution

  public:
    //! type of function converter for exact solution and gradient
    typedef ExactSolution ExactSolutionType;

  protected:
    ExactSolutionType exactSolution_;

  public:
    const ExactSolutionType& exactSolution( const double time=0.0 ) const
    {
      return exactSolution_;
    }
  };

}
}
#endif  /*DUNE_PROBLEMINTERFACE_HH*/
