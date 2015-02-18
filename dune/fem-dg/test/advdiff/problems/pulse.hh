#ifndef  DUNE_PROBLEM_PULSE_HH
#define  DUNE_PROBLEM_PULSE_HH

// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

namespace Dune {

/**
 * @brief describes the initial and exact solution of the advection-diffusion model
 * for given constant velocity vector v=(v1,v2)
 *
 * \f[u(x,y,z,t):=\displaystyle{\sum_{i=0}^{1}} T_i(t) \cdot X_i(x) \cdot
 * Y_i(y) \cdot Z_i(z)\f]
 *
 * with
 *
 * \f{eqnarray*}{
 * T_0(t) &:=&  e^{-\varepsilon t \pi^2 (2^2 + 1^2 + 1.3^2 )} \\
 * X_0(x) &:=&  0.6\cdot \cos(2\pi (x-v_1t)) + 0.8\cdot \sin(2\pi (x-v_1t)) \\
 * Y_0(y) &:=&  1.2\cdot \cos(1\pi (y-v_2t)) + 0.4\cdot \sin(1\pi (y-v_2t)) \\
 * Z_0(z) &:=&  0.1\cdot \cos(1.3\pi (z-v_3t)) - 0.4\cdot \sin(1.3\pi (z-v_3t)) \\
 * T_1(t) &:=&  e^{-\varepsilon t \pi^2 (0.7^2 + 0.5^2 + 0.1^2 )} \\
 * X_1(x) &:=&  0.9\cdot \cos(0.7\pi (x-v_1t)) + 0.2\cdot \sin(0.7\pi (x-v_1t)) \\
 * Y_1(y) &:=&  0.3\cdot \cos(0.5\pi (y-v_2t)) + 0.1\cdot \sin(0.5\pi (y-v_2t))
 * Z_1(z) &:=&  -0.3\cdot \cos(0.1\pi (z-v_3t)) + 0.2\cdot \sin(0.1\pi (z-v_3t)) \\
 * \f}
 *
 * This is a solution of the AdvectionDiffusionModel for \f$g_D = u|_{\partial
 * \Omega}\f$.
 *
 */
template <class GridType, int dimRange>
class Pulse : public EvolutionProblemInterface<
                  Fem::FunctionSpace< double, double, GridType::dimension, dimRange>,
                  false >
{
public:
  typedef EvolutionProblemInterface<
                 Fem::FunctionSpace< double, double,
                                     GridType::dimension, dimRange >,
                 false >                                              BaseType;

  enum{ dimDomain = BaseType :: dimDomain };
  typedef typename BaseType :: DomainType            DomainType;
  typedef typename BaseType :: RangeType             RangeType;
  typedef typename BaseType :: JacobianRangeType     JacobianRangeType;

  typedef Fem::Parameter  ParameterType;

  /**
   * @brief define problem parameters
   */
  Pulse () :
    BaseType () ,
    startTime_( ParameterType::getValue<double>("femhowto.startTime",0.0) ),
    epsilon_( ParameterType::getValue<double>("femhowto.epsilon",0.1) ),
    spotmid_( 0 ),
    myName_("AdvDiff")
  {
    spotmid_[0] = -0.25;
    std::cout <<"Problem: "<<myName_<< ", epsilon " << epsilon_ << "\n";
    //std::cout <<"Problem: HeatEqnWithAdvection, epsilon_" <<  epsilon_ << "\n";
  }

  //! this problem has no source term
  bool hasStiffSource() const { return false; }
  bool hasNonStiffSource() const { return false; }

  double stiffSource(const DomainType& arg,
                const double t,
                const RangeType& u,
                RangeType& res) const
  {
    return 0.0;
  }

  double nonStiffSource(const DomainType& arg,
                const double t,
                const RangeType& u,
                RangeType& res) const
  {
    return 0.0;
  }

  double diffusion( const RangeType& u, const JacobianRangeType& gradU ) const
  {
    return epsilon();
  }

  //! return start time
  double startTime() const { return startTime_; }

  //! return start time
  double epsilon() const { return epsilon_; }

  /**
   * @brief getter for the velocity
   */
  void velocity(const DomainType& x, DomainType& v) const
  {
    // rotation in 2d
    v[0] = -4.0*x[1];
    v[1] =  4.0*x[0];
    for(int i=2; i<DomainType :: dimension; ++i) v[i] = 0;
  }

  /**
   * @brief evaluates \f$ u_0(x) \f$
   */
  void evaluate(const DomainType& arg, RangeType& res) const
  {
    evaluate(arg, startTime_, res);
  }

  /**
   * @brief evaluate exact solution
   */
  void evaluate(const DomainType& arg, const double t, RangeType& res) const
  {
    const double x = arg[0];
    const double y = arg[1];

    const double sig2 = 0.004; /* Siehe Paper P.Bastian Gl. 30 */
    const double sig2PlusDt4 = sig2+(4.0*epsilon_*t);
    const double xq = ( x*cos(4.0*t) + y*sin(4.0*t)) - spotmid_[0];
    const double yq = (-x*sin(4.0*t) + y*cos(4.0*t)) - spotmid_[1];

    res = (sig2/ (sig2PlusDt4) ) * exp (-( xq*xq + yq*yq ) / sig2PlusDt4 );
  }

  /**
   * @brief latex output for EocOutput
   */
  std::string description() const
  {
    std::ostringstream ofs;

    ofs << "Problem: " << myName_
      << ", Epsilon: " << epsilon_ ;

    ofs << ", End time: " << ParameterType::template getValue<double>("femhowto.endTime");

    return ofs.str();
  }

  /*  \brief finalize the simulation using the calculated numerical
   *  solution u for this problem
   *
   *  \param[in] variablesToOutput Numerical solution in the suitably chosen variables
   *  \param[in] eocloop Specific EOC loop
   */
  template< class DiscreteFunctionType >
  void finalizeSimulation( DiscreteFunctionType& variablesToOutput,
                           const int eocloop) const
  {}

protected:
  const double  startTime_;
  const double  epsilon_;
  DomainType spotmid_;
  std::string myName_;
};

}
#endif  /*DUNE_PROBLEM_HH__*/

