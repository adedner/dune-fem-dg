#ifndef  DUNE_PROBLEM_PULSE_HH
#define  DUNE_PROBLEM_PULSE_HH

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
   * \brief Describes the initial and exact solution of the advection-diffusion model
   * described in:
   *
   * P. Bastian. Higher Order Discontinuous Galerkin Methods for Flow and Transport in Porous Media
   * Challenges in Scientific Computing - CISC 2002, Volume 35 of the series
   * Lecture Notes in Computational Science and Engineering pp 1-22
   *
   */
  template <class GridType, int dimRange>
  class Pulse : public EvolutionProblemInterface<
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
    Pulse () :
      BaseType () ,
      startTime_( ParameterType::getValue<double>("femdg.stepper.starttime",0.0) ),
      epsilon_( ParameterType::getValue<double>("epsilon",0.1) ),
      spotmid_( 0 ),
      center_( 0.5 ), // assume unit square
      myName_("pulse")
    {
      spotmid_[0] = -0.25;
      std::cout <<"Problem: "<<myName_<< ", epsilon " << epsilon_ << "\n";
    }

    //! this problem has no source term
    bool hasStiffSource() const override { return false; }
    bool hasNonStiffSource() const override { return false; }

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
      return 0.0;
    }

    double diffusion( const RangeType& u, const JacobianRangeType& gradU ) const override
    {
      return epsilon();
    }

    //! return start time
    double startTime() const override { return startTime_; }

    //! return start time
    double epsilon() const override { return epsilon_; }

    /**
     * \brief problem velocity
     */
    void velocity(const DomainType& x, const double time, DomainType& v) const override
    {
      // rotation in 2d
      v[0] = -4.0*(x[1] - center_[ 1 ]);
      v[1] =  4.0*(x[0] - center_[ 0 ]);
      for(int i=2; i<DomainType :: dimension; ++i) v[i] = 0;
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
      const double x = arg[0] - center_[ 0 ];
      const double y = arg[1] - center_[ 1 ];

      const double sig2 = 0.004; /* Siehe Paper P.Bastian Gl. 30 */
      const double sig2PlusDt4 = sig2+(4.0*epsilon_*t);
      const double xq = ( x*cos(4.0*t) + y*sin(4.0*t)) - spotmid_[0];
      const double yq = (-x*sin(4.0*t) + y*cos(4.0*t)) - spotmid_[1];

      res = (sig2/ (sig2PlusDt4) ) * exp (-( xq*xq + yq*yq ) / sig2PlusDt4 );
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
    const double  startTime_;
    const double  epsilon_;
    DomainType spotmid_;
    DomainType center_;
    std::string myName_;
  };

}
}
#endif  /*DUNE_PROBLEM_HH__*/

