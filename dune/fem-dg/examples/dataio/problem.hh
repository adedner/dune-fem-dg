#ifndef  DUNE_PROBLEM_IMPLEMENTATION_HH
#define  DUNE_PROBLEM_IMPLEMENTATION_HH

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief describes the initial and exact solution of the advection-diffusion model
   *
   * \f[u(x,y,z,t):=\displaystyle{\sum_{i=0}^{1}} v_i(t) \cdot \mu_i(x) \cdot
   * \nu_i(y) \cdot \omega_i(z)\f]
   *
   * with
   *
   * \f{eqnarray*}{
   * v_0(t)&:=&e^{-\varepsilon t \pi^2 (2^2 + 1^2 + 1.3^2 )} \\
   * \mu_0(x)&:=&0.6\cdot \cos(2\pi (x-at)) + 0.8\cdot \sin(2\pi (x-at)) \\
   * \nu_0(y)&:=&1.2\cdot \cos(1\pi (y-at)) + 0.4\cdot \sin(1\pi (y-at)) \\
   * v_1(t)&:=&e^{-\varepsilon t \pi^2 (0.7^2 + 0.5^2 + 0.1^2 )} \\
   * \mu_1(x)&:=&0.9\cdot \cos(0.7\pi (x-at)) + 0.2\cdot \sin(0.7\pi (x-at)) \\
   * \nu_1(y)&:=&0.3\cdot \cos(0.5\pi (y-at)) + 0.1\cdot \sin(0.5\pi (y-at))
   * \f}
   *
   * This is a solution of the AdvectionDiffusionModel for \f$g_D = u|_{\partial
   * \Omega}\f$.
   *
   */
  template <class GridType>
  struct U0 : public EvolutionProblemInterface<
                        Dune::Fem::FunctionSpace< double, double,
                        GridType::dimensionworld, GridType::dimensionworld+2 >,
                        true >
  {
    double SQR( const double a ) const
    {
      return (a * a);
    }

  public:
    typedef EvolutionProblemInterface<
                   Dune::Fem::FunctionSpace< double, double,
                   GridType::dimensionworld, GridType::dimensionworld+2 >,
                   true >                                              BaseType;
    using BaseType::evaluate;
    using BaseType::velocity;

    enum{ dimDomain = BaseType :: dimDomain };
    typedef typename BaseType :: DomainType                            DomainType;
    typedef typename BaseType :: RangeType                             RangeType;

    typedef Dune::Fem::Parameter  ParameterType ;

    /**
     * \brief define problem parameters
     */
    U0 () :
      BaseType () ,
      velocity_( 0 ),
      startTime_( ParameterType::getValue<double>("femdg.stepper.starttime",0.0) ),
      epsilon_  ( ParameterType::getValue<double>("epsilon"  ,0.1) )
    {
        velocity_[0]=0.8;
        velocity_[1]=0.6;

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

        myName = "AdvectDiff";
      }

    //! return start time
    double startTime() const override { return startTime_; }

    //! return start time
    double epsilon() const override { return epsilon_; }

    /**
     * \brief getter for the velocity
     */
    void velocity(const DomainType& x, DomainType& v) const
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
      ofs << "Problem: " << myName << "\n\n"
        << "Epsilon = " << epsilon_ << "\n\n"
        << "Exact solution: $u(x,y,z,t):=\\displaystyle{\\sum_{i=0}^{"
        << max_n_of_coefs_-1
        << "}} v_i(t) \\cdot \\mu_i(x) \\cdot \\nu_i(y) \\cdot \\omega_i(z)$\n\n";

      for(int i=0;i<max_n_of_coefs_;++i)
        {
          std::ostringstream temp;
          if(dimDomain > 1) {
            temp << common_coef_y_[i] << "^2 ";
          }
          if(dimDomain > 2) {
            temp << "+ " << common_coef_z_[i] << "^2 ";
          }
          ofs << "$v_" << i << "(t):=e^{-\\varepsilon t \\pi^2 ("
            << common_coef_x_[i] << "^2 + "
            << temp.str()
            << ")} $\n\n"
            << "$\\mu_" << i << "(x):=" << cos_coef_x_[i]
            << "\\cdot \\cos(" << common_coef_x_[i] << "\\pi (x-at)) + " << sin_coef_x_[i]
            << "\\cdot \\sin(" << common_coef_x_[i] << "\\pi (x-at)) $";
          if(dimDomain > 1)
          {
            ofs << "\n\n"
              << "$\\nu_" << i
              << "(y):=" << cos_coef_y_[i] << "\\cdot \\cos("
              << common_coef_y_[i] << "\\pi (y-at)) + " << sin_coef_y_[i]
              << "\\cdot \\sin(" << common_coef_y_[i] << "\\pi (y-at)) $";
          }
          if(dimDomain >2)
          {
            ofs << "\n\n"
              << "$\\omega_" << i << "(z):=" << cos_coef_z_[i]
              << "\\cdot \\cos(" << common_coef_z_[i] << "\\pi (z-at)) + " << sin_coef_z_[i]
              << "\\cdot \\sin(" << common_coef_z_[i] << "\\pi (z-at)) $";
          }
          ofs << "\n\n";
        }
      ofs << "\n\n";
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
    double     startTime_;
  public:
    double      epsilon_;
    std::string myName;
  };

} // end namespace
} // end namespace Dune
#endif

