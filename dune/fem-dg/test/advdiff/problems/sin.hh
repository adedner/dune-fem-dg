#ifndef  DUNE_PROBLEM_SIN_HH
#define  DUNE_PROBLEM_SIN_HH
 
// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

namespace Dune {

template <class GridType>
class U0Sin : public EvolutionProblemInterface<
                Fem::FunctionSpace< double, double, GridType::dimension, DIMRANGE>,
                false >
{
public:
  typedef EvolutionProblemInterface<
                 Fem::FunctionSpace< double, double,
                                     GridType::dimension, DIMRANGE >,
                 false >                                BaseType;

  enum{ dimDomain = BaseType :: dimDomain };
  typedef typename BaseType :: DomainType            DomainType;
  typedef typename BaseType :: RangeType             RangeType;
  typedef typename BaseType :: JacobianRangeType     JacobianRangeType;

  typedef Fem::Parameter  ParameterType;

  /**
   * @brief define problem parameters
   */
  U0Sin () :
    BaseType () ,
    velocity_( 1 ), 
    startTime_( ParameterType::getValue<double>("femhowto.startTime",0.0) ),
    epsilon_( ParameterType::getValue<double>("femhowto.epsilon",0.1) ),
    rhsFactor_( epsilon_ * 2.0 * std:: pow( 2.0, (double) dimDomain) * M_PI * M_PI ),
    myName_("U0Sin")
  {
    std::cout <<"Problem: "<<myName_<< ", epsilon " << epsilon_ << "\n";
  }

  //! this problem has no source term
  bool hasStiffSource() const { return false; }
  bool hasNonStiffSource() const { return true; }

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
    // eval solution 
    evaluate( arg, t, res );
    // apply factor
    res *= rhsFactor_;
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
    v = velocity_;
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
    res = 1.0;
    const double pi = 2.0 * M_PI ;
    for( int i=0; i< DomainType :: dimension; ++i)
    {
      res *= std :: sin( pi * (arg[i] - t) );
    }
  }

  /**
   * @brief latex output for EocOutput
   */
  std::string description() const
  {
    std::ostringstream ofs;

    ofs << "Problem: " << myName_
      << ", Epsilon: " << epsilon_ ;

    ofs << ", End time: " << ParameterType::getValue<double>("femhowto.endTime");

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
  const DomainType velocity_;
  const double  startTime_;
  const double  epsilon_;
  const double  rhsFactor_;
  std::string myName_;
};

}
#endif  /*DUNE_PROBLEM_HH__*/

