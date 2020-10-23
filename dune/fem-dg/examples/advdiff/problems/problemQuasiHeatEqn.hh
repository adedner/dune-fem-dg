#ifndef  DUNE_PROBLEM__QUASI_HH__
#define  DUNE_PROBLEM__QUASI_HH__

#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune
{
namespace Fem
{

  /**
   * \brief describes the initial and exact solution of the advection-diffusion model
   *
   * \ingroup AdvDiffProblems
   */
  template <class GridType, int dimRange >
  struct QuasiHeatEqnSolution : public EvolutionProblemInterface<
                    Dune::Fem::FunctionSpace< typename GridType::ctype, typename GridType::ctype,
                                              GridType::dimension, dimRange> >
  {
  public:
    typedef EvolutionProblemInterface<
                   Dune::Fem::FunctionSpace< typename GridType::ctype, typename GridType::ctype,
                                             GridType::dimension, dimRange > >   BaseType;
    using BaseType::evaluate;

    enum{ dimDomain = BaseType::dimDomain };
    typedef typename BaseType::DomainType            DomainType;
    typedef typename BaseType::RangeType             RangeType;
    typedef typename BaseType::JacobianRangeType     JacobianRangeType;

    typedef Fem :: Parameter  ParameterType ;

    /**
     * \brief define problem parameters
     */
    QuasiHeatEqnSolution () :
      BaseType () ,
      velocity_( 0 ),
      startTime_( ParameterType::getValue<double>("femdg.stepper.starttime",0.0) ),
      epsilon_( ParameterType::template getValue<double>("epsilon") )
    {
      if ( (dimRange != 1) || (dimDomain != 2) )
      {
        std::cout <<"QuasiHeatEqn only supports dimRange=1 and dimDomain=2\n";
        abort();
      }
      std::cout <<"Problem: QuasiHeatEqn" <<epsilon_ << " epsilon \n";

      myName = "quasi";
    }

    double diffusion( const RangeType& u, const JacobianRangeType& gradU ) const override
    {
      return epsilon_ * u.infinity_norm();
    }

    double startTime() const override { return startTime_; }

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
      res = std::sin(arg[0]) * std::sin(arg[1]) * std::exp( -epsilon_*t );
    }

    bool hasStiffSource() const override { return true; }
    bool hasNonStiffSource() const override { return false; }

    /**
     * \brief evaluate stiff source function
     */
    double nonStiffSource( const DomainType& arg,
                           const double t,
                           const RangeType& u,
                           RangeType& res ) const override
    {
      res = 0.;

      // time step restriction from the stiff source term
      return 0.;
    }

    /**
     * \brief evaluate non stiff source function
     */
    double stiffSource( const DomainType& arg,
                        const double t,
                        const RangeType& u,
                        RangeType& res ) const override
    {
      const double x = arg[0];
      const double y = arg[1];
      const double cosx2 = cos(x)*cos(x);
      const double sinx2 = sin(x)*sin(x);
      const double cosy2 = cos(y)*cos(y);
      const double siny2 = sin(y)*sin(y);

      const double expVal = -std::exp(-2.*epsilon_*t ) * (cosx2*siny2 + cosy2*sinx2);
      for( int r=0; r<dimRange; ++r )
      {
        res[ r ] = expVal + 2.0 * u[ r ] * u[ r ] - u[ r ];
      }
      res *= epsilon_;

      // time step restriction from the non stiff source term
      return 0.0;
    }

    /**
     * \brief latex output for EocOutput
     */
    std::string description() const override
    {
      std::ostringstream ofs;

      ofs << "Problem: " << myName
        << ", Epsilon: " << epsilon_
        << ", End time: " << ParameterType::template getValue<double>("femdg.stepper.endtime");

      return ofs.str();
    }



  private:
    DomainType velocity_;
    double     startTime_;

  public:
    double      epsilon_;
    std::string myName;
  };

}
}
#endif  /*DUNE_PROBLEM_HH__*/

