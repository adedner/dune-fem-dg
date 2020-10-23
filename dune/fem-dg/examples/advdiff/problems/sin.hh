#ifndef  DUNE_PROBLEM_SIN_HH
#define  DUNE_PROBLEM_SIN_HH

// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Describes a Sin problem of the advection-diffusion equation.
   *
   * \ingroup AdvDiffProblems
   */
  template <class GridType, int dimRange>
  class U0Sin : public EvolutionProblemInterface<
                  Fem::FunctionSpace< typename GridType :: ctype, typename GridType::ctype,
                                      GridType::dimension, dimRange> >
  {
  public:
    typedef EvolutionProblemInterface<
                   Fem::FunctionSpace< typename GridType :: ctype, typename GridType::ctype,
                                       GridType::dimension, dimRange > >  BaseType;
    using BaseType::evaluate;

    enum{ dimDomain = BaseType :: dimDomain };
    typedef typename BaseType :: DomainType            DomainType;
    typedef typename BaseType :: RangeType             RangeType;
    typedef typename BaseType :: JacobianRangeType     JacobianRangeType;

    typedef Fem::Parameter  ParameterType;

    /**
     * \brief define problem parameters
     */
    U0Sin () :
      BaseType () ,
      velocity_( 1 ),
      startTime_( ParameterType::getValue<double>("femdg.stepper.starttime",0.0) ),
      epsilon_( ParameterType::getValue<double>("epsilon",0.1) ),
      rhsFactor_( epsilon_ * 2.0 * std:: pow( 2.0, (double) dimDomain) * M_PI * M_PI ),
      massFactor_( 0.5 ),
      myName_("sin")
    {
      std::cout <<"Problem: "<<myName_<< ", epsilon " << epsilon_ << "\n";
    }

    //! this problem has no source term
    bool hasStiffSource() const override { return false; }
    bool hasNonStiffSource() const override { return true; }
    bool hasMass() const override { return true; }

    double stiffSource(const DomainType& arg,
                       const double t,
                       const RangeType& u,
                       RangeType& res) const override
    {
      return 0.0;
    }

    double nonStiffSource(const DomainType& arg,
                          const double t,
                          const RangeType& u,
                          RangeType& res) const override
    {
      // eval solution
      evaluate( arg, t, res );
      // apply factor
      res *= rhsFactor_ * massFactor_;
      return 0.0;
    }

    //! diagonal of the mass term
    virtual inline void mass (const DomainType& arg,
                              const double time,
                              const RangeType& u,
                              RangeType& diag ) const override
    {
      diag = massFactor_;
    }


    double diffusion( const RangeType& u, const JacobianRangeType& gradU ) const override
    {
      return epsilon() * massFactor_;
    }

    //! return start time
    double startTime() const override { return startTime_; }

    //! return start time
    double epsilon() const override { return epsilon_; }

    /**
     * \brief getter for the velocity
     */
    void velocity(const DomainType& x, const double time, DomainType& v) const override
    {
      v = velocity_;
    }

    /**
     * \brief evaluates \f$ u_0(x) \f$
     */
    void evaluate(const DomainType& arg, RangeType& res) const
    {
      evaluate(arg, startTime_, res);
    }

    /**
     * \brief evaluate exact solution
     */
    void evaluate(const DomainType& arg, const double t, RangeType& res) const override
    {
      res = 1.0;
      const double pi = 2.0 * M_PI ;
      for( int i=0; i< DomainType :: dimension; ++i)
      {
        res *= std :: sin( pi * (arg[i] - t) );
      }
    }

    /**
     * \brief latex output for EocOutput
     */
    std::string description() const override
    {
      std::ostringstream ofs;

      ofs << "Problem: " << myName_
        << ", Epsilon: " << epsilon_ ;

      ofs << ", End time: " << ParameterType::template getValue<double>("femdg.stepper.endtime");

      return ofs.str();
    }

  protected:
    const DomainType velocity_;
    const double  startTime_;
    const double  epsilon_;
    const double  rhsFactor_;
    const double  massFactor_;
    std::string myName_;
  };

}
}
#endif  /*DUNE_PROBLEM_HH__*/

