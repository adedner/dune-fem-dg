#ifndef DUNE_FEMDG_OPERATORSPLITTING_HH
#define DUNE_FEMDG_OPERATORSPLITTING_HH

#include <cmath>

namespace Dune
{
namespace Fem
{

  struct OperatorSplittingScheme
  {
    static constexpr double theta_ = 1.0 - 1.0/M_SQRT2 ;
    static constexpr double alpha_ = (1.0 - 2.0 * theta_ )/(1.0 - theta_);
    static constexpr double beta_  = theta_/ ( 1.0 - theta_ );

    static const bool hasAdvection = true;
    static const bool hasDiffusion = true;
    static const bool hasSource    = true;
    static const bool hasGrad      = true;
  };

  template< int Step, bool RightHandSide >
  class FractionalStepThetaScheme;

  //step 0
  template<>
  class FractionalStepThetaScheme<0,false> //new timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return 1.0; }
    static constexpr double diffusion(){ return Split::alpha_; }
    static constexpr double advection(){ return 0.0; }
    static constexpr double grad(){ return 1.0; }
    static constexpr double mass(){ return Split::theta_;}

    static const bool hasAdvection = false;

    template< class DF >
    static void velocity( const DF& un, DF& velocity  ){ /*no velocity needed*/ }
  };
  template<>
  class FractionalStepThetaScheme<0,true> //old timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return 0.0;}
    static constexpr double diffusion(){ return Split::beta_;}
    static constexpr double advection(){ return 1.0;}
    static constexpr double grad(){ return 0.0;}
    static constexpr double mass(){ return Split::theta_;}

    static const bool hasSource    = false;
    static const bool hasGrad      = false;

    template< class DF >
    static void velocity( const DF& un, DF& velocity  )
    {
      velocity.assign( un );
    }
  };


  //step 1
  template<>
  class FractionalStepThetaScheme<1,false> //new timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return 1.0;}
    static constexpr double diffusion(){ return Split::beta_;}
    static constexpr double advection(){ return 1.0;}
    static constexpr double grad(){ return 0.0;}
    static constexpr double mass(){ return 1.0 - 2.0 * Split::theta_;}

    static const bool hasGrad      = false;

    template< class DF >
    static void velocity( const DF& un, const DF& untheta, DF& velocity )
    {
      velocity.assign( un );
      velocity *= ( 2.0*theta_ - 1.0 ) / theta_ ;
      velocity.axpy( (1.0-theta_)/theta_, untheta );
    }
  };
  template<>
  class FractionalStepThetaScheme<1,true> //old timestep
    : public OperatorSplittingScheme
  {
    using Split = OperatorSplittingScheme;
  public:
    static constexpr double source(){ return  0.0;}
    static constexpr double diffusion(){ return  Split::alpha_;}
    static constexpr double advection(){ return  0.0;}
    static constexpr double grad(){ return  1.0;}
    static constexpr double mass(){ return 1.0 - 2.0 * Split::theta_;}

    static const bool hasAdvection = false;
    static const bool hasSource    = false;

    template< class DF >
    static void velocity( const DF& un, const DF& untheta, DF& velocity ){ /*no velocity needed*/ }
  };

  //step 2
  template<>
  class FractionalStepThetaScheme<2,false> //new timestep
    : public FractionalStepThetaScheme<0,0>
  {};
  template<>
  class FractionalStepThetaScheme<2,true> //old timestep
    : public FractionalStepThetaScheme<0,1>
  {};

}
}
#endif // #ifndef OPERATORSPLITTING_HH
