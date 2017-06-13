#ifndef DUNE_DEFORMATIONALFLOW_HH
#define DUNE_DEFORMATIONALFLOW_HH

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
   * \brief describes the initial and exact solution of the advection-diffusion model
   * for given constant velocity vector v=(v1,v2)
   *
   * \ingroup AdvDiffProblems
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
   * Y_1(y) &:=&  0.3\cdot \cos(0.5\pi (y-v_2t)) + 0.1\cdot \sin(0.5\pi (y-v_2t)) \\
   * Z_1(z) &:=&  -0.3\cdot \cos(0.1\pi (z-v_3t)) + 0.2\cdot \sin(0.1\pi (z-v_3t))
   * \f}
   *
   * This is a solution of the AdvectionDiffusionModel for \f$g_D = u|_{\partial
   * \Omega}\f$.
   *
   */
  template <class GridType, int dimRange>
  struct DeformationalFlow : public EvolutionProblemInterface<
                    Fem::FunctionSpace< typename GridType :: ctype, typename GridType::ctype,
                                        GridType::dimension, dimRange> >
  {
  public:
    typedef EvolutionProblemInterface<
                   Fem::FunctionSpace< typename GridType :: ctype, typename GridType::ctype,
                                       GridType::dimension, dimRange > >  BaseType;

    enum{ dimDomain = BaseType :: dimDomain };
    typedef typename BaseType :: DomainType            DomainType;
    typedef typename BaseType :: RangeType             RangeType;
    typedef typename BaseType :: JacobianRangeType     JacobianRangeType;

    typedef Fem::Parameter  ParameterType;

    /**
     * \brief define problem parameters
     */
    DeformationalFlow () :
      BaseType () ,
      center_( 0.5 ),
      startTime_( ParameterType::getValue<double>("femdg.stepper.starttime",0.0) ),
      endTime_( ParameterType::template getValue<double>("femdg.stepper.endtime") ),
      epsilon_( ParameterType::getValue<double>("epsilon", 0.0 ) )
    {
      myName = "deform";
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
      return epsilon_;
    }

    //! return start time
    double startTime() const { return startTime_; }

    //! return start time
    double epsilon() const { return epsilon_; }

    /**
     * \brief getter for the velocity
     */
    void velocity(const DomainType& x, const double time, DomainType& v) const
    {
      DomainType p( x );
      p -= center_;

      const double r = p.two_norm();
      const double f1 = std::pow( (4.0*r), 6.0 );
      const double factor = (1.0 - f1) / (1.0 + f1);

      const double theta  = std::atan( p[ 1 ]/p[ 0 ] );
      const double utheta = 4.0 * M_PI * r / endTime_ * ( 1.0 - ( std::cos( 2.0 * M_PI * time / endTime_ ) * factor )) ;

      v[ 0 ] =  utheta * std::sin( theta );
      v[ 1 ] = -utheta * std::cos( theta );
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
    void evaluate(const DomainType& x, const double t, RangeType& res) const
    {
      DomainType c2( 0.5 );
      c2[ 0 ] = 0.3;

      DomainType p( x );
      p -= c2;

      const double rtilde = p.two_norm() * 5.0;
      if( rtilde <= 1.0 )
      {
        double rc = 0.5 * ( 1.0 + std::cos( M_PI * rtilde) );
        res = rc * rc ;
      }
      else
      {
        res = 0.0;
      }
    }

    /**
     * \brief latex output for EocOutput
     */
    std::string description() const
    {
      std::ostringstream ofs;

      ofs << "Problem: " << myName
        << ", Epsilon: " << epsilon_ ;
      return ofs.str();
    }

  private:
    DomainType center_;
    const double  startTime_;
    const double  endTime_;
  public:
    const double  epsilon_;
    std::string myName;
  };

} // end namespace
} // end namespace Dune
#endif

