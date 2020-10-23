#ifndef  DUNE_PROBLEM_HH__
#define  DUNE_PROBLEM_HH__

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
  struct U0 : public EvolutionProblemInterface<
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
    U0 () :
      BaseType () ,
      velocity_( 0 ),
      startTime_( ParameterType::getValue<double>("femdg.stepper.starttime",0.0) ),
      epsilon_( ParameterType::getValue<double>("epsilon",0.1) )
    {
      std::cout <<"Problem: HeatEqnWithAdvection, epsilon " << epsilon_ << "\n";
      //std::cout <<"Problem: HeatEqnWithAdvection, epsilon_" <<  epsilon_ << "\n";

      // an advection vector, default to 0
      velocity_[0] = ParameterType::getValue<double>( "xvelocity", 0. );
      if (dimDomain > 1)
        velocity_[1] = ParameterType::getValue<double>( "yvelocity", 0. );
      if (dimDomain > 2)
        velocity_[2] = ParameterType::getValue<double>( "zvelocity", 0. );

      max_n_of_coefs_ = 2;

      //x coordinate
      common_coef_x_[0] = 2.0;
      sin_coef_x_[0]    = 0.8;
      cos_coef_x_[0]    = 0.6;

      //x coordinate
      common_coef_x_[1] = 0.7;
      sin_coef_x_[1]    = 0.2;
      cos_coef_x_[1]    = 0.9;

      //y coordinate
      common_coef_y_[0] = 1.0;
      sin_coef_y_[0]    = 0.4;
      cos_coef_y_[0]    = 1.2;


      //y coordinate
      common_coef_y_[1] = 0.5;
      sin_coef_y_[1]    = 0.1;
      cos_coef_y_[1]    = 0.3;


      //z coordinate
      common_coef_z_[0] = 1.3;
      sin_coef_z_[0]    = -0.4;
      cos_coef_z_[0]    = 0.1;

      //z coordinate
      common_coef_z_[1] = 0.1;
      sin_coef_z_[1]    = 0.2;
      cos_coef_z_[1]    = -0.3;

      myName = "heat";
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

    template <class T>
    double SQR( const T& a ) const
    {
      return ( a * a );
    }


    /**
     * \brief evaluate exact solution
     */
    void evaluate(const DomainType& arg, const double t, RangeType& res) const override
    {

      res = 0;
      DomainType x(arg);
      x -= DomainType(0.5);
      x.axpy(-t,velocity_);

      for(int i=0;i<max_n_of_coefs_;++i)
        {
          if(dimDomain == 1)
            res += exp(-epsilon_*t*(SQR(common_coef_x_[i]*M_PI)))\
                   *((cos_coef_x_[i]*cos(common_coef_x_[i]*M_PI*x[0])\
                      +  sin_coef_x_[i]*sin(common_coef_x_[i]*M_PI*x[0])));
          else if(dimDomain == 2)
            res += exp(-epsilon_*t*(SQR(common_coef_x_[i]*M_PI)+\
                                   SQR(common_coef_y_[i]*M_PI)))\
                   *((cos_coef_x_[i]*cos(common_coef_x_[i]*M_PI*x[0])\
                      +  sin_coef_x_[i]*sin(common_coef_x_[i]*M_PI*x[0]))\
                     *(cos_coef_y_[i]*cos(common_coef_y_[i]*M_PI*x[1])\
                       + sin_coef_y_[i]*sin(common_coef_y_[i]*M_PI*x[1])));
          else if(dimDomain == 3)
            res += exp(-epsilon_*t*(SQR(common_coef_x_[i]*M_PI)+\
                                   SQR(common_coef_y_[i]*M_PI)+\
                                   SQR(common_coef_z_[i]*M_PI)))\
                   *((cos_coef_x_[i]*cos(common_coef_x_[i]*M_PI*x[0])\
                      +  sin_coef_x_[i]*sin(common_coef_x_[i]*M_PI*x[0]))\
                     *(cos_coef_y_[i]*cos(common_coef_y_[i]*M_PI*x[1])\
                       + sin_coef_y_[i]*sin(common_coef_y_[i]*M_PI*x[1]))\
                     *(cos_coef_z_[i]*cos(common_coef_z_[i]*M_PI*x[2])\
                       + sin_coef_z_[i]*sin(common_coef_z_[i]*M_PI*x[2])));
        }
    }

    /**
     * \brief latex output for EocOutput
     */
    std::string description() const override
    {
      std::ostringstream ofs;

      ofs << "Problem: " << myName
        << ", Epsilon: " << epsilon_
        << ", Advection vector: (" <<velocity_[0];

      if (dimDomain > 1)
        ofs <<"," <<velocity_[1];
      if (dimDomain > 2)
        ofs <<"," <<velocity_[2];
      ofs <<")";

      ofs << ", End time: " << ParameterType::template getValue<double>("femdg.stepper.endtime");

      return ofs.str();
    }

  private:
    DomainType velocity_;
    int        max_n_of_coefs_;
    double     common_coef_x_[2];
    double     sin_coef_x_[2];
    double     cos_coef_x_[2];
    double     common_coef_y_[2];
    double     sin_coef_y_[2];
    double     cos_coef_y_[2];
    double     common_coef_z_[2];
    double     sin_coef_z_[2];
    double     cos_coef_z_[2];
    const double  startTime_;
  public:
    const double  epsilon_;
    std::string myName;
  };

}
}
#endif  /*DUNE_PROBLEM_HH__*/

