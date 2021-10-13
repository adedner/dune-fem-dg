#ifndef EULERPROBLEMS_HH
#define EULERPROBLEMS_HH

#include <cmath>

// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/l1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/analyticaleulerflux.hh>
#include <dune/fem/io/parameter.hh>

// local includes
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>
#include <dune/fem-dg/examples/euler/chorjo.hh>


namespace Dune
{
namespace Fem
{

  /**
   * \brief Interface class for Euler problems
   *
   * \ingroup EulerProblems
   */

  template <class GridType>
  class ProblemBase :
    public EvolutionProblemInterface<
                Dune::Fem::FunctionSpace< typename GridType::ctype,
                                          typename GridType::ctype,
                                          GridType::dimensionworld, GridType::dimensionworld+2>,
                false >
  {
    typedef EvolutionProblemInterface<
                Dune::Fem::FunctionSpace< typename GridType::ctype,
                                          typename GridType::ctype,
                                          GridType::dimensionworld, GridType::dimensionworld+2>,
                false > BaseType ;

  public:
    using BaseType :: evaluate ;
    using BaseType::fixedTimeFunction;

    typedef Dune::Fem::Parameter ParameterType;

    enum { Inflow = 1, Outflow = 2, Reflection = 3 , Slip = 4 };
    enum { MaxBnd = Slip };

    typedef typename BaseType :: FunctionSpaceType  FunctionSpaceType;

    enum { dimDomain = GridType::dimensionworld };
    typedef typename FunctionSpaceType :: RangeType  RangeType ;
    typedef typename FunctionSpaceType :: DomainType DomainType ;

    typedef typename FunctionSpaceType :: RangeFieldType  RangeFieldType ;
    typedef RangeFieldType FieldType ;

    const FieldType  gamma_;
    const FieldType  V_;

    ProblemBase() : gamma_( 1.4 ), V_( 1.0 )  {}
    ProblemBase( const FieldType gamma, const FieldType  V )
      : gamma_( gamma ), V_( V ) {}

    FieldType gamma() const { return gamma_ ; }
    FieldType V () const { return V_; }

    virtual double endTime () const { return 0.1 ; }

    void bg ( const DomainType&, RangeType& ) const {}

    //! methods for gradient based indicator
    bool twoIndicators() const { return false ; }

    //! methods for gradient based indicator
    double indicator1( const DomainType& xgl, const RangeType& u ) const
    {
      // use density as indicator
      return u[ 0 ];
    }

    virtual int boundaryId ( const int id ) const
    {
      return 2;
    }
  };


  /**
   * \brief Euler problem.
   *
   * \ingroup EulerProblems
   */
  template <class GridType>
  class U0Smooth1D : public ProblemBase < GridType >
  {
  public:
    typedef ProblemBase < GridType > BaseType ;
    typedef typename BaseType :: RangeType   RangeType ;
    typedef typename BaseType :: DomainType  DomainType ;

    using BaseType :: gamma ;
    using BaseType :: evaluate ;

    DomainType u_;
    std::string myName;

  public:
    U0Smooth1D() : BaseType( 1.4, 1.0 ),
      u_(0), myName("Advection")
    {
      init();
    }

    void init()
    {
      enum { dim = DomainType :: dimension  };
      u_[ 0 ] = std::cos(M_PI/5.0);
      for(int i=1; i<dim; ++i)
      {
        u_[ i ] = std::sin(M_PI/5.0);
      }
    }

    double endtime() {
      return 0.2;
    }

    void evaluate(const DomainType& x, const double time, RangeType& res) const
    {
      enum { dimR = RangeType :: dimension  };
      enum { dim = DomainType :: dimension  };

      res = 0;
      DomainType y(0.5);
      const double p = 0.3;

      double tmp_c[dim], tmp=0.0, u_sq=0.0;
      for(int i=0; i<dim; ++i)
      {
        tmp_c[i] = x[i] - y[i] - time * u_[i];
        tmp += tmp_c[i]*tmp_c[i];
        u_sq += u_[i] * u_[i];
      }
      tmp *= 16.0;

      const double cos_tmp = cos(tmp * M_PI) + 1.0;
      const double rho = (tmp > 1.0)? 0.5: 0.25 * cos_tmp*cos_tmp + 0.5;

      res[0] = rho;
      for(int i=0; i<dim; i++) res[1+i] = u_[i] * rho;
      res[dim+1] = p/(gamma() - 1.0) + 0.5*rho*u_sq;
    }

  };


  /**
   * \brief Euler problem.
   *
   * \ingroup EulerProblems
   */
  template <class GridType>
  class U0IsentropicVortex : public ProblemBase < GridType >
  {
  public:
    typedef ProblemBase < GridType > BaseType ;
    typedef typename BaseType :: RangeType   RangeType ;
    typedef typename BaseType :: DomainType  DomainType ;

    using BaseType :: gamma ;
    using BaseType :: evaluate ;

    std::string myName;

  public:
    U0IsentropicVortex() :
      BaseType( 1.4, 1.0 ),
      myName("Advection")
    {
      init();
    }

    void init()
    {
    }

    double endtime() {
      return 100.;
    }

    void evaluate(const DomainType& point, const double time, RangeType& res) const
    {
      enum { dimR = RangeType :: dimension  };
      enum { dim = DomainType :: dimension  };

      const double S = 5.0;
      const double M = 0.5;
      const double gamma_1 = gamma() - 1.0;
      // const double R = 1.0; // radius, not used yet
      const double x = point[ 0 ] - time;
      const double y = point[ 1 ] - time;

      const double f = (1.0 - x*x - y*y);

      const double rho = std::pow((1. - S*S*gamma_1*M*M*std::exp(f)/(8.*M_PI*M_PI)), (1./gamma_1) );
      DomainType u( 0 );
      u[ 0 ] = 1. - S*y*std::exp(0.5*f)/(2.*M_PI);
      u[ 1 ] = S*x*std::exp(0.5*f)/(2.*M_PI);
      const double p = std::pow( rho, gamma())/(gamma()*M*M);

      res[0] = rho;
      for(int i=0; i<dim; ++i) res[1+i] = u[i] * rho;
      res[dim+1] = p/(gamma_1) + (u*u)*0.5/rho;

      //std::cout << "Initial " << res << std::endl;
    }

  };


  /**
   * \brief Euler problem.
   *
   * \ingroup EulerProblems
   */
  template <class GridType>
  class U0FFS : public ProblemBase< GridType > {
  public:
    typedef ProblemBase< GridType > BaseType ;
    typedef typename BaseType :: RangeType   RangeType ;
    typedef typename BaseType :: DomainType  DomainType ;

    using BaseType :: gamma ;
    using BaseType :: V ;

    U0FFS() : BaseType( 1.4, 3.0 ),
              myName("Forward Facing Step") {}

    const std::string myName;

    double endtime()
    {
      return 5.0;
    }

    void evaluate(const DomainType& arg, double t, RangeType& res) const
    {
      evaluate(t,arg,res);
    }

    template <class DomainType, class RangeType>
    void evaluate(double t,const DomainType& arg, RangeType& res) const
    {
      enum { dimR = RangeType :: dimension };
      // reset all values
      res = 0;

      // initial values
      res[0] = gamma(); // density
      res[1] = V() * res[0]; // velocity-x
      res[dimR-1] = 8.8;

      //std::cout << "Initial value : "<< res << "\n";
      //RangeType prim;
      //cons2prim( res, prim );
      //std::cout << "Initial value : "<< prim << "\n";
    }

    void cons2prim(const RangeType& u,RangeType& v) const {
      v[0] = u[0];
      enum { dimDomain = BaseType :: dimDomain };
      for (int i=1;i<=dimDomain;i++)
      {
        v[i] = u[i]/u[0];
      }
      v[dimDomain+1] = EulerAnalyticalFlux<dimDomain>().pressure(gamma,u);
    }

    int boundaryId( const int id ) const
    {
      return (id > BaseType::MaxBnd) ? BaseType::MaxBnd : id;
    }

    int determineBndId(const DomainType& x) const
    {
      if( x[0] <= 1e-8     ) return 1; // Inflow
      if( x[0] >= 3.0-1e-8 ) return 2; // Outflow
      // reflection elsewhere
      return 3;
    }

  };

  /**
   * \brief Euler problem.
   *
   * \ingroup EulerProblems
   */
  template <class GridType>
  class U0Sod : public ProblemBase< GridType > {
    double T,startTime;
    Dune::FieldVector<double,6> Ulr;
    enum { dim = GridType :: dimension };
  public:
    typedef ProblemBase< GridType > BaseType ;
    typedef typename BaseType :: RangeType   RangeType ;
    typedef typename BaseType :: DomainType  DomainType ;

    typedef typename BaseType :: ParameterType ParameterType;

    using BaseType :: gamma ;
    using BaseType :: V ;

    enum { dimDomain = GridType::dimensionworld,
           dimRange = dimDomain+2,
           energ = dimRange-1};

    U0Sod() : T(0.4), startTime(0), flag(0) {
      myName = "RP-Sod";

      // default is sod's rp
      Ulr[0] = 1.0;
      Ulr[3] = 0.125;
      Ulr[1] = 0.;
      Ulr[4] = 0.;
      Ulr[2] = 1.0;
      Ulr[5] = 0.1;

      //ParameterType::get("euler.riemanndata", Ulr, Ulr );
      T = ParameterType::template getValue<double>("femdg.stepper.endtime"/*, T*/);
      ParameterType::get("femdg.stepper.starttime", startTime, startTime );
      flag = ParameterType::getValue("problemflag", flag);
    }

    int boundaryId( const int id ) const
    {
      return ( id > 2 ) ? BaseType::Reflection : id;
      //return BaseType::Inflow;
      //return ( id == 1 ) ? BaseType::Inflow : BaseType::Outflow ;
      //return BaseType::OutFlow;
      //return (id > BaseType::Inflow) ? BaseType::OutFlow : id;
    }

    double endtime() {
      return T;
    }
    double saveinterval() {
      return 0.01;
    }
    bool useReflectBndLR() const {
      return false;
    }
    bool useReflectBndTB() const {
      return flag==0;
    }

    void evaluate(const DomainType& arg, const double t, RangeType& res) const
    {
      if (flag==0) {
        res[2] = 0.;
        if (t>1e-8) {
          chorin(t,arg[0]-0.5,res[0],res[1],res[3]);
        }
        else {
          if (arg[0]<0.5) {
            res[0]=Ulr[0];
            res[1]=Ulr[1];
            res[3]=Ulr[2];
          } else {
            res[0]=Ulr[3];
            res[1]=Ulr[4];
            res[3]=Ulr[5];
          }
        }
      }
      else {
        DomainType c(0.075);
        DomainType velo(0.5);
        c[0] = 0.1;
        velo[0] = 4.5;

        DomainType x(arg);
        DomainType v(velo);
        x -= c;
        v *= t;
        x -= v;
        double r2=0.004;
        double z=x*x/r2;

        res[0] = 0.1;
        if (flag<=3) {
          if (z < 1.) {
            if (flag==1) {
              double f = 0.5*(cos(z*M_PI)+1.);
              double pf = pow(f,4.);
              res[0] += pf;
            }
            else if (flag==2)
              res[0] += (1.-z);
            else
              res[0] += 1.0;
          }
        } else if (flag>=4) {
          double s = 1.;
          for (int i=0;i<dimDomain;++i)
            s *= sin(M_PI*x[i]);
          if (flag==4)
            res[0] = 0.6+0.5*s;
          else if (flag==5 || flag==6 || flag==7 || flag==8) {
            while (x[0]<0)    x[0]+=1.;
            while (x[0]>1)    x[0]-=1.;
            while (x[1]<0)    x[1]+=0.25;
            while (x[1]>0.25) x[1]-=0.25;
            res[0] = 0.6+0.5*std::abs(sin(M_PI*x[0]))*std::abs(sin(2.*M_PI*x[1]));
            if (flag==6)
              if (x[0]-M_PI*x[1]>0)
                res[0]+=1.;
            if (flag==7)
              if (x[0]<0.5)
                res[0]+=1.0;
            if (flag==8)
              if (x[1]<0.1)
                res[0]+=1.0;
          }
        }
        for (int i=0;i<dimDomain;++i)
          res[i+1] = velo[i];
        res[energ] = 0.4;
      }
      double kinEnerg = 0;
      for (int i=0;i<dimDomain;++i) {
        kinEnerg += res[i+1]*res[i+1];
        res[i+1] *= res[0];
      }
      kinEnerg *= 0.5*res[0];
      res[energ]  = res[energ]/(gamma()-1.0)+kinEnerg;
    }

    void chorin(double t,double x,
                double& q_erg,double& u_erg,double& p_erg) const
    {
      EULERCHORIN::
        lsg(x,t,&q_erg,&u_erg,&p_erg,
            Ulr[0],Ulr[3],Ulr[1],Ulr[4],Ulr[2],Ulr[5],
            double(gamma()));
    }

    template <class Field>
    void chorin(double t, double x,
                Field& pq_erg, Field& pu_erg, Field& pp_erg) const
    {
      double q_erg = double(pq_erg);
      double u_erg = double(pu_erg);
      double p_erg = double(pp_erg);

      chorin( t, x, q_erg, u_erg, p_erg );

      pq_erg = q_erg;
      pu_erg = u_erg;
      pp_erg = p_erg;
    }

    void printmyInfo(std::string filename)
    {
      std::ostringstream filestream;
      filestream << filename;

      std::ofstream ofs(filestream.str().c_str(), std::ios::app);

      ofs << "Problem: " << myName << "\n\n"
          << "gamma = " << gamma() << "\n\n";
      ofs << "\n\n";

      ofs.close();

    }
    std::string myName;
    int flag;
  };

  /**
   * \brief Euler problem.
   *
   * \ingroup EulerProblems
   */
  template< class Grid >
  class RiemannProblem
  : public ProblemBase< Grid >
  {
    typedef ProblemBase< Grid > BaseType;

    static const int dimension = Grid::dimension;

    double T, startTime;

    typedef Dune::FieldVector< double, 6 > RiemannDataType;
    RiemannDataType Ulr;

  public:
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType :: ParameterType ParameterType;

    using BaseType::gamma;
    using BaseType::V;

    static const int dimDomain = Grid::dimensionworld;
    static const int dimRange = dimDomain + 2;
    static const int energ = dimRange - 1;

    RiemannProblem ()
      : T( 0.15 ), startTime( 0 )
    {
      myName = "RP";

      FieldVector<double,6> data;
      data[ 0 ] =  1.0;
      data[ 3 ] =  1.0;
      data[ 1 ] = -2.0;
      data[ 4 ] =  2.0;
      data[ 2 ] =  0.4;
      data[ 5 ] =  0.4;

      T = ParameterType::template getValue<double>( "femdg.stepper.endtime"/*, T*/ );
      ParameterType::get( "femdg.stepper.starttime", startTime, startTime );
      data = ParameterType::getValue("riemanndata", data );
      for(int i=0; i<6; ++i ) Ulr[ i ] = data[ i ];
    }

    int boundaryId ( const int id ) const
    {
      return ( id > 2 ) ? BaseType::Reflection : id;
      //return BaseType::Inflow;
      //return ( id == 1 ) ? BaseType::Inflow : BaseType::Outflow ;
      //return BaseType::OutFlow;
      //return (id > BaseType::Inflow) ? BaseType::OutFlow : id;
    }

    double endtime () { return T; }
    double saveinterval () { return 0.01; }
    bool useReflectBndLR() const { return false; }
    bool useReflectBndTB() const { return true; }

    void evaluate ( const DomainType &arg, double t, RangeType &res ) const
    {
      res[2] = 0.;
      if( t > 1e-8 )
        chorin( t, arg[ 0 ] - 0.5, res[ 0 ], res[ 1 ], res[ 3 ] );
      else
      {
        if( arg[ 0 ] < 0.5 )
        {
          res[0]= Ulr[ 0 ];
          res[1]= Ulr[ 1 ];
          res[3]= Ulr[ 2 ];
        }
        else
        {
          res[ 0 ] = Ulr[ 3 ];
          res[ 1 ] = Ulr[ 4 ];
          res[ 3 ] = Ulr[ 5 ];
        }
      }
      double kinEnerg = 0;
      for( int i = 0; i < dimDomain; ++i )
      {
        kinEnerg += res[ i+1 ]*res[ i+1 ];
        res[ i+1 ] *= res[ 0 ];
      }
      kinEnerg *= 0.5*res[ 0 ];
      res[energ]  = res[ energ ] / (gamma()-1.0) + kinEnerg;
    }

    void chorin ( double t, double x, double &q_erg, double &u_erg, double &p_erg ) const
    {
      EULERCHORIN::lsg( x, t, &q_erg, &u_erg, &p_erg, Ulr[ 0 ], Ulr[ 3 ], Ulr[ 1 ], Ulr[ 4 ], Ulr[ 2 ], Ulr[ 5 ], gamma() );
    }

    template <class Field>
    void chorin(double t, double x,
                Field& pq_erg, Field& pu_erg, Field& pp_erg) const
    {
      double q_erg = double(pq_erg);
      double u_erg = double(pu_erg);
      double p_erg = double(pp_erg);

      chorin( t, x, q_erg, u_erg, p_erg );

      pq_erg = q_erg;
      pu_erg = u_erg;
      pp_erg = p_erg;
    }

    void printmyInfo( const std::string &filename )
    {
      std::ostringstream filestream;
      filestream << filename;

      std::ofstream ofs(filestream.str().c_str(), std::ios::app);

      ofs << "Problem: " << myName << "\n\n"
          << "gamma = " << gamma() << "\n\n";
      ofs << "\n\n";

      ofs.close();
    }

    std::string myName;
  };


  #if 0
  /*****************************************************************/
  // Diffraction
  template <class GridType>
  class U0Diffraction : public ProblemBase< GridType > {
  public:
    enum { dimDomain = GridType::dimensionworld };
    typedef FunctionSpace<double,double,dimDomain,dimDomain+2> FunctionSpaceType;
    U0Diffraction() :
      gamma(1.4),myName("Shock Diffraction Problem")
    , pAB_( pAlphaBeta( V() ) )
    , rhoAB_( rhoAlphaBeta( pAB_ ) )
    // c = std::sqrt( 1.4 * 1 / 1.4 ) = 1
    , us_ ( u_s(1., pAB_ ))
    {}
    U0Diffraction(std::string,double,bool diff_timestep=true) :
      gamma(1.4),myName("Shock Diffraction Problem")
    , pAB_( pAlphaBeta( V() ) )
    , rhoAB_( rhoAlphaBeta( pAB_ ) )
    , us_ ( u_s(1., pAB_ ))
    {}

    // public member, Andreas .....
    const double gamma;
    const std::string myName;
    const double pAB_;
    const double rhoAB_;
    const double us_;

    void printmyInfo(std::string filename) {}
    double endtime()
    {
      return 2.3;
    }

    double V() const { return 5.09; }

    double pAlphaBeta(const double machNumber ) const
    {
      return ( 1. + ( (2*gamma)/(gamma+1) * ( (machNumber*machNumber) - 1.) ) );
    }

    double rhoAlphaBeta(const double pAB ) const
    {
      const double gammaPlusMinus = ((gamma+1.)/(gamma-1.));
      return (1. + gammaPlusMinus * pAB ) /( gammaPlusMinus + pAB );
    }

    double u_s (const double c, double pAB ) const
    {
      return c/gamma * (pAB - 1.) * std::sqrt(((2*gamma)/(gamma+1))/(pAB + ((gamma-1.)/(gamma+1.)) ));
    }

    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg, RangeType& res) const {
      evaluate(0,arg,res);
    }
    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg,double t, RangeType& res) const
    {
      evaluate(t,arg,res);
    }

    template <class DomainType, class RangeType>
    void evaluate(double t,const DomainType& arg, RangeType& res) const
    {
      enum { dimR = RangeType :: dimension };
      enum { e = dimR - 1 };

      // reset all values
      res = 0;

      // density
      res[0] = gamma;

      // pressure
      res[e] = 1;

      // check area
      double x = arg[0];

      if( x <= 0.5 )
      {
        res[0] *= rhoAB_;
        res[1] = us_;
        res[e] *= pAB_;
      }

      const double sqrVelo = res[1] * res[1];

      // get impuls
      res[1] *= res[0];

      // calculate energy
      res[e] = res[e]/(gamma-1.) + 0.5 * sqrVelo * res[0];

      //std::cout << "Initial value : "<< res << "\n";
      //RangeType prim;
      //cons2prim( res, prim );
      //std::cout << "Initial value : "<< prim << "\n";
    }

    template <class RangeType>
    void cons2prim(const RangeType& u,RangeType& v) const {
      v[0] = u[0];
      for (int i=1;i<=dimDomain;i++)
      {
        v[i] = u[i]/u[0];
      }
      v[dimDomain+1] = EulerFlux<dimDomain>().pressure(gamma,u);
    }

    template <class DomainType>
    int determineBndId(const DomainType& x) const
    {
      if( x[0] <= 0.) return 1;
      if( x[1] >= 11.0 ) return 4;
      if( x[0] >= 13.0 ) return 2;

      return 3;
    }
  };
  /*****************************************************************/
  /*****************************************************************/
  // Diffraction
  template <class GridType>
  class U0DoubleMachReflect : public ProblemBase {
  public:
    enum { dimDomain = GridType::dimensionworld };
    typedef FunctionSpace<double,double,dimDomain,dimDomain+2> FunctionSpaceType;
    U0DoubleMachReflect() :
      gamma(1.4),myName("Double Mach Reflection")
    , pAB_( pAlphaBeta( V() ) )
    , rhoAB_( rhoAlphaBeta( pAB_ ) )
    // c = sqrt( 1.4 * 1 / 1.4 ) = 1
    , us_ ( u_s(1., pAB_ ))
    {}
    U0DoubleMachReflect(std::string,double,bool diff_timestep=true) :
      gamma(1.4),myName("Double Mach Reflection")
    , pAB_( pAlphaBeta( V() ) )
    , rhoAB_( rhoAlphaBeta( pAB_ ) )
    , us_ ( u_s(1., pAB_ ))
    {}

    // public member, Andreas .....
    const double gamma;
    const std::string myName;
    const double pAB_;
    const double rhoAB_;
    const double us_;

    void printmyInfo(std::string filename) {}
    double endtime()
    {
      return 1.0;
    }

    double V() const { return 10.; }

    double pAlphaBeta(const double machNumber ) const
    {
      return ( 1. + ( (2*gamma)/(gamma+1) * ( (machNumber*machNumber) - 1.) ) );
    }

    double rhoAlphaBeta(const double pAB ) const
    {
      const double gammaPlusMinus = ((gamma+1.)/(gamma-1.));
      return (1. + gammaPlusMinus * pAB ) /( gammaPlusMinus + pAB );
    }

    double u_s (const double c, double pAB ) const
    {
      return c/gamma * (pAB - 1.) * std::sqrt(((2*gamma)/(gamma+1))/(pAB + ((gamma+1.)/(gamma-1.))));
    }

    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg, RangeType& res) const {
      evaluate(time(),arg,res);
    }
    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg,double t, RangeType& res) const
    {
      evaluate(t,arg,res);
    }

    template <class DomainType, class RangeType>
    void evaluate(double t,const DomainType& arg, RangeType& res) const
    {
      enum { dimR = RangeType :: dimension };
      enum { e = dimR - 1 };

      // reset all values
      res = 0;

      // density
      res[0] = gamma;

      // pressure
      res[e] = 1;

      // check area
      double x = arg[0];

      if( x <= 0. )
      {
        res[0] *= rhoAB_;
        res[1] = us_;
        res[e] *= pAB_;
      }

      const double sqrVelo = res[1] * res[1];

      // get impuls
      res[1] *= res[0];

      // calculate energy
      res[e] = res[e]/(gamma-1.) + 0.5 * sqrVelo * res[0];

      //std::cout << "Conservative value : "<< res << "\n";
      //RangeType prim;
      //cons2prim( res, prim );
      //std::cout << "Primitve value : "<< prim << "\n";
    }

    template <class RangeType>
    void cons2prim(const RangeType& u,RangeType& v) const {
      v[0] = u[0];
      for (int i=1;i<=dimDomain;i++)
      {
        v[i] = u[i]/u[0];
      }
      v[dimDomain+1] = EulerFlux<dimDomain>().pressure(gamma,u);
    }
  };

  /*****************************************************************/
  // ShockBubble
  template <class GridType>
  class U0ShockBubble : public ProblemBase {
  public:
    enum { dimDomain = GridType::dimensionworld };
    typedef FunctionSpace<double,double,dimDomain,dimDomain+2> FunctionSpaceType;

    typedef typename FunctionSpaceType :: DomainType DomainType;

    // public member, Andreas .....
    const double gamma;
    const std::string myName;
    DomainType center_;
    const double radius2_;

    U0ShockBubble() :
      gamma(1.4),myName("Shock Bubble Problem")
     , center_(0.5) , radius2_( 0.2 * 0.2 )
    {
      center_[dimDomain-1] = 0;
    }

    U0ShockBubble(std::string,double,bool diff_timestep=true) :
      gamma(1.4),myName("Shock Bubble Problem")
     , center_(0) , radius2_( 0.2 * 0.2 )
    {
      center_[0] = 0.5;
    }

    void printmyInfo(std::string filename) {}
    double endtime()
    {
      return 0.3;
    }

    double V() const { return 2.95; }

    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg, RangeType& res) const {
      evaluate(time(),arg,res);
    }
    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg,double t, RangeType& res) const
    {
      evaluate(t,arg,res);
    }

    template <class DomainType, class RangeType>
    void evaluate(double t, const DomainType& arg, RangeType& res) const
    {
      enum { dimR = RangeType :: dimension };

      // reset all values
      res = 0;

      // behind shock
      if ( arg[0] <= 0.2 )
      {
        /*
        const double gamma1 = gamma-1.;
        // pressure left of shock
        const double pinf = 10.0;
        const double rinf = ( gamma1 + (gamma+1)*pinf )/( (gamma+1) + gamma1*pinf );
        const double vinf = V();
        //(1.0/std::sqrt(gamma)) * (pinf - 1.)/
        //        std::sqrt( 0.5*((gamma+1)/gamma) * pinf + 0.5*gamma1/gamma);

        res[0] = rinf;
        res[dimR-1] = 0.5*rinf*vinf*vinf + pinf/gamma1;

        //res[1] = vinf * rinf;
        res[1] = V() * rinf;
        */

        const double gamma1 = gamma-1.;
        // pressure left of shock
        const double pinf = 5;
        const double rinf = ( gamma1 + (gamma+1)*pinf )/( (gamma+1) + gamma1*pinf );
        const double vinf = (1.0/std::sqrt(gamma)) * (pinf - 1.)/
                std::sqrt( 0.5*((gamma+1)/gamma) * pinf + 0.5*gamma1/gamma);

        res[0] = rinf;
        res[dimR-1] = 0.5*rinf*vinf*vinf + pinf/gamma1;
        res[1] = vinf * rinf;

        //RangeType prim;
        //cons2prim( res, prim );
        //std::cout << "Primitve behind: "<< prim << "\n";
      }
      else if( (arg - center_).two_norm2() <= radius2_ )
      {
        res[0] = 0.1;
        // pressure in bubble
        res[dimR-1] = 2.5;

        //RangeType prim;
        //cons2prim( res, prim );
        //std::cout << "Primitve inside : "<< prim << "\n";
      }
      // elsewhere
      else
      {
        res[0] = 1;
        res[dimR-1] = 2.5;

        //RangeType prim;
        //cons2prim( res, prim );
        //std::cout << "Primitve outside : "<< prim << "\n";
      }
    }

    template <class RangeType>
    void cons2prim(const RangeType& u,RangeType& v) const {
      v[0] = u[0];
      for (int i=1;i<=dimDomain;i++)
      {
        v[i] = u[i]/u[0];
      }
      v[dimDomain+1] = EulerFlux<dimDomain>().pressure(gamma,u);
    }
  };
  template <class GridType>
  class U0Sin : public ProblemBase {
  public:
    enum { dimDomain = GridType::dimensionworld };
    typedef FunctionSpace<double,double,dimDomain,dimDomain+2> FunctionSpaceType;
    U0Sin() :
      gamma(), startTime(0) {
      ParameterType::get("femdg.stepper.starttime",startTime,startTime);
    }
    U0Sin(std::string,double,bool diff_timestep=true) :
      gamma(1.4), startTime(0) {
      ParameterType::get("femdg.stepper.starttime",startTime,startTime);
    }
    double endtime() {
      return 2.;
    }
    double saveinterval() {
      return 0.1;
    }
    bool useReflectBndLR() const {
      return false;
    }
    bool useReflectBndTB() const {
      return true;
    }

    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg, RangeType& res) const {
      evaluate(time(),arg,res);
    }
    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg,double t, RangeType& res) const
    {
      evaluate(t,arg,res);
    }
    template <class DomainType, class RangeType>
    void evaluate(double t,const DomainType& arg, RangeType& res) const {
      int e = DomainType::size + 1;
      double x = (arg[0]-10.);
      res=0.;
      if (x<0) {
        res[0]=3.857143;
        res[1]=-0.920279;      // 2.629369;
        res[e]=10.33333;
      } else if (x<10) {
        res[0]= 1. + 0.2 * sin(5.*x);
        res[1]=-3.549648;     // 0.0
        res[e]=1.0;
      } else {
        res[0]= 1.;
        res[1]=-3.549648;     // 0.0
        res[e]=1.0;
      }
      res[1] *= res[0];
      res[e] = res[e]/(gamma-1.0)+
        0.5*(res[1]*res[1])/res[0];
    }
    void printmyInfo(std::string filename)
    {
      std::ostringstream filestream;
      filestream << filename;

      std::ofstream ofs(filestream.str().c_str(), std::ios::app);

      ofs << "Problem: " << myName << "\n\n"
          << "gamma = " << gamma << "\n\n";
      ofs << "\n\n";

      ofs.close();
    }
    double gamma,startTime;
    std::string myName;
  };
  template <class GridType>
  class U0RotatingCone : public ProblemBase {
  public:
    enum { dimDomain = GridType::dimensionworld };
    typedef FunctionSpace<double,double,dimDomain,dimDomain+2> FunctionSpaceType;
    U0RotatingCone() :
      gamma(), startTime(0) {
      ParameterType::get("femdg.stepper.starttime",startTime,startTime);
    }
    U0RotatingCone(std::string,double,bool diff_timestep=true) :
      gamma(1.4), startTime(0) {
      ParameterType::get("femdg.stepper.starttime",startTime,startTime);
    }
    double endtime() {
      return 0.5;
    }
    double saveinterval() {
      return 0.01;
    }
    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg, RangeType& res) const {
      evaluate(startTime,arg,res);
    }
    template <class DomainType, class RangeType>
    void evaluate(const DomainType& arg,double t, RangeType& res) const
    {
      evaluate(t,arg,res);
    }
    template <class DomainType, class RangeType>
    void evaluate(double t,const DomainType& arg, RangeType& res) const {
      res=0.;
      DomainType c(0.5);
      DomainType x=arg;
      x-=c;
      double r2=0.04;
      if (x*x < r2) {
        res[0] =cos(x*x/r2*M_PI)+2;
      }
      else {
        res[0] = 1.0;
      }
      x=arg;
      x-=DomainType(1.0);
      if (DomainType::size>1) {
        res[1] = x[1]*res[0];
        res[2] = -x[0]*res[0];
      } else {
        res[1] = -1.*res[0];
      }
      res[DomainType::size+1] = 2.;
      if (arg.size>1) {
        res[DomainType::size+1] +=
    0.5*(res[1]*res[1]+res[2]*res[2])/res[0];
      } else {
        res[DomainType::size+1] +=
    0.5*(res[1]*res[1])/res[0];
      }
    }
    void printmyInfo(std::string filename)
    {
      std::ostringstream filestream;
      filestream << filename;

      std::ofstream ofs(filestream.str().c_str(), std::ios::app);

      ofs << "Problem: " << myName << "\n\n"
          << "gamma = " << gamma << "\n\n";
      ofs << "\n\n";

      ofs.close();

    }
    double gamma,startTime;
    std::string myName;
  };

#endif

} // end namespace Fem
} // end namespace Dune

#endif

/* CLAWPACK example
 *
 * Uin value in bubble
 * Uout values outside bubble rigth of shock
 * Uinf values left of shock
 *
 * 0.0              t0          = initial time
 * 0.0              xlower      = left edge of computational domain
 * 1.2              xupper      = right edge of computational domain
 * 0.0              ylower      = bottom edge of computational domain
 * 0.5              yupper      = top edge of computational domain
 * 0.0              zlower      = front edge of computational domain
 * 0.5              zupper      = back edge of computational domain
 *
 * 1.4                       gamma
 * 0.5    0.5    0.2         x0, y0, r0:  center and radius of bubble
 * 0.1                       rhoin: density in bubble
 * 5.0                       pinf:  pressure behind shock
 *
 * # density outside bubble and pressure ahead of shock are fixed:
 * rhoout = 1.d0
 * pout   = 1.d0
 * pin    = 1.d0
 *
 * # Compute density and velocity behind shock from Hugoniot relations:
 * # gamma1 = gamma-1
 * rinf = ( gamma1 + (gamma+1)*pinf )/( (gamma+1) + gamma1*pinf )
 * vinf = (1.0d0/sqrt(gamma)) * (pinf - 1.d0)/
 *        sqrt( 0.5*((gamma+1)/gamma) * pinf + 0.5*gamma1/gamma )
 * einf = 0.5*rinf*vinf*vinf + pinf/gamma1
 */

