#ifndef DUNE_FEM_LDG_FLUXOPERATOR_HH
#define DUNE_FEM_LDG_FLUXOPERATOR_HH

#include <string>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#include <dune/fem-dg/operator/limiter/limitpass.hh>

// note to me: it doesn't make sense to include primaldiscretemodel.hh
//             but it was design in that way. to be removed later!!!!!!
// local includes
#include <dune/fem-dg/operator/dg/primaldiscretemodel.hh>
#include <dune/fem-dg/operator/dg/fluxdiscretemodel.hh>
#include <dune/fem-dg/operator/dg/operatorbase.hh>
#include <dune/fem-dg/pass/dgpass.hh>


namespace Dune
{
namespace Fem
{

  // LDGAdvectionDiffusionOperator
  //-----------------------------

  template< class Traits, bool advection = true , bool diffusion = true >
  class LDGAdvectionDiffusionOperator :
    public Fem::SpaceOperatorInterface< typename Traits::DestinationType >
  {
    // choose ids for the three passes (including start pass) different to the tuple entries
    enum PassIdType { u = Traits::ModelType::modelParameterSize, gradPass = u+1, advectPass = u + 2 };

    struct GradientTraits
      : public Traits
    {
      // overload discrete function space
      typedef typename Traits::DiscreteFunctionSpaceType::template
        ToNewDimRange< Traits::ModelType::Traits::dimGradRange >::Type DiscreteFunctionSpaceType;

      template < class DF >
      struct ToNewSpace;

      template < class Space, template <class> class DF >
      struct ToNewSpace< DF< Space > >
      {
        typedef DF< DiscreteFunctionSpaceType > Type;
      };

      template < class Space, typename... Args, template < class, typename... > class DF >
      struct ToNewSpace< DF< Space, Args... > >
      {
        typedef DF< DiscreteFunctionSpaceType, Args... > Type;
      };

      typedef typename ToNewSpace< typename Traits::DestinationType >::Type DestinationType;
    };

    typedef Fem::SpaceOperatorInterface< typename Traits::DestinationType > BaseType;

  public:
    typedef typename Traits::AdvectionFluxType                        AdvectionFluxType;
    typedef typename Traits::ModelType                                ModelType;
    typedef typename ModelType::ProblemType                           ProblemType;

    enum { dimRange  = ModelType::dimRange };
    enum { dimDomain = ModelType::Traits::dimDomain };
    enum { dimGradRange = ModelType::Traits::dimGradRange };

    typedef typename Traits::LimiterIndicatorType                     LimiterIndicatorType;
    typedef typename LimiterIndicatorType::DiscreteFunctionSpaceType  LimiterIndicatorSpaceType;

    // Pass 2 Model (advection)
    typedef AdvectionDiffusionLDGModel< Traits, u, gradPass, advection, diffusion >
                                                                      DiscreteModel2Type;

    typedef typename DiscreteModel2Type::DiffusionFluxType            DiffusionFluxType;

    // Pass 1 Model (gradient)
    typedef GradientModel< GradientTraits, u >                        DiscreteModel1Type;

    typedef typename DiscreteModel1Type::Traits                       Traits1;
    typedef typename DiscreteModel2Type::Traits                       Traits2;

    typedef typename Traits::GridType                                 GridType;

    //typedef typename Traits2::Traits::DomainType                              DomainType;
    typedef typename Traits2::DestinationType                         DiscreteFunction2Type;

    typedef typename Traits1::DiscreteFunctionSpaceType               Space1Type;
    typedef typename Traits2::DiscreteFunctionSpaceType               Space2Type;
    typedef typename Traits1::DestinationType                         Destination1Type;
    typedef typename Traits2::DestinationType                         Destination2Type;
    typedef Destination2Type                                          DestinationType;
    typedef Space2Type                                                SpaceType;

    typedef typename Traits1::GridPartType                            GridPartType;

    typedef Fem::StartPass< DiscreteFunction2Type, u >                Pass0Type;
    typedef LocalCDGPass< DiscreteModel1Type, Pass0Type, gradPass >   Pass1Type;
    typedef LocalCDGPass< DiscreteModel2Type, Pass1Type, advectPass > Pass2Type;

    typedef typename Traits::ExtraParameterTupleType                  ExtraParameterTupleType;

  public:
    template< class ExtraParameterTupleImp >
    LDGAdvectionDiffusionOperator( GridPartType& gridPart,
                                   const ModelType& model,
                                   ExtraParameterTupleImp tuple,
                                   const std::string keyPrefix = "" ) :
      model_( model ),
      numflux_( model_ ),
      gridPart_( gridPart ),
      space1_( gridPart_ ),
      space2_( gridPart_ ),
      diffFlux_( gridPart_, model_ ),
      discModel1_(model_, diffFlux_ ),
      discModel2_(model_, numflux_, diffFlux_),
      pass0_ (),
      pass1_(discModel1_, pass0_, space1_),
      pass2_(discModel2_, pass1_, space2_)
    {}

    void setTime(const double time) {
	    pass2_.setTime( time );
    }

    double timeStepEstimate() const {
	    return pass2_.timeStepEstimate();
    }

    void operator()( const DestinationType& arg, DestinationType& dest ) const {
	    pass2_( arg, dest );
    }

    inline const SpaceType& space() const {
	    return space2_;
    }

    inline void switchupwind()
    {
      diffFlux_.switchUpwind();
    }

    inline double maxAdvectionTimeStep() const
    {
      return discModel2_.maxAdvectionTimeStep();
    }
    inline double maxDiffusionTimeStep() const
    {
      return discModel2_.maxDiffusionTimeStep();
    }

    inline double computeTime() const
    {
      return pass2_.computeTime();
    }

    inline size_t numberOfElements () const
    {
      return pass2_.numberOfElements();
    }

    void printmyInfo(std::string filename) const {
	    std::ostringstream filestream;
            filestream << filename;
            std::ofstream ofs(filestream.str().c_str(), std::ios::app);
            ofs << "LDG Op., polynomial order: " << Traits::polynomialOrder << "\\\\\n\n";
            ofs.close();
    }

    std::string description() const
    {
      std::stringstream stream;
      stream <<" {\\bf LDG Diff. Op.}, flux formulation, order: " << Traits::polynomialOrder+1
             <<", $\\eta = ";
      diffFlux_.diffusionFluxPenalty( stream );
      stream <<"$, {\\bf Adv. Flux:} ";
      //TODO has to be implemented
      stream <<",\\\\\n";
      return stream.str();
    }

    using BaseType::discreteModel;

  private:
    ModelType           model_;
    AdvectionFluxType   numflux_;
    GridPartType&       gridPart_;
    Space1Type          space1_;
    Space2Type          space2_;

  protected:
    DiffusionFluxType   diffFlux_;

  private:
    DiscreteModel1Type  discModel1_;
    DiscreteModel2Type  discModel2_;
    Pass0Type           pass0_;
    Pass1Type           pass1_;
    Pass2Type           pass2_;
  };


  // LDGAdvectionTraits
  //-------------------

  template <class Traits,
            bool advection>
  struct LDGAdvectionTraits : public Traits
  {
    // choose ids for the two passes (including start pass) different to the tuple entries
    enum PassIdType { u = Traits::ModelType::modelParameterSize, cdgpass = u + 1 };

    typedef AdvectionDiffusionDGPrimalModel
      // put a method_none here to avoid diffusion
      < Traits, u, advection, false> DiscreteModelType;
  };


  // DGAdvectionOperator
  //--------------------

  template< class OpTraits >
  class LDGAdvectionOperator : public
    DGAdvectionDiffusionOperatorBase< LDGAdvectionTraits< OpTraits, true > >
  {
    typedef LDGAdvectionTraits< OpTraits, true> Traits;
    typedef DGAdvectionDiffusionOperatorBase< Traits > BaseType;
  public:
    typedef typename BaseType::GridPartType            GridPartType;
    typedef typename BaseType::ProblemType             ProblemType;
    typedef typename BaseType::ExtraParameterTupleType ExtraParameterTupleType;

    // constructor: do not touch/delegate everything
    template< class ... Args>
    LDGAdvectionOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

    std::string description() const
    {
      std::stringstream stream;
      stream <<"{\\bf Adv. Op.}, flux formulation, order: " << Traits::polynomialOrder+1
             <<", {\\bf Adv. Flux:} ";
      /*if (FLUX==1)
        stream <<"LLF";
      else if (FLUX==2)
        stream <<"HLL";*/
      stream <<",\\\\\n";
      return stream.str();
    }
  };


  // LDGDiffusionOperator
  //--------------------

  template< class Traits >
  class LDGDiffusionOperator : public
    LDGAdvectionDiffusionOperator< Traits, false >
  {
    typedef LDGAdvectionDiffusionOperator< Traits, false >  BaseType;

  public:
    typedef typename BaseType::GridPartType            GridPartType;
    typedef typename BaseType::ProblemType             ProblemType;
    typedef typename BaseType::ExtraParameterTupleType ExtraParameterTupleType;

    // constructor: do not touch/delegate everything
    template< class ... Args>
    LDGDiffusionOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

    std::string description() const
    {
      std::stringstream stream;
      stream <<"{\\bf LDG Diff. Op.}, flux formulation, order: " << Traits::polynomialOrder+1
             <<", $\\eta = ";
      diffFlux_.diffusionFluxPenalty( stream );
      stream <<"$, {\\bf Adv. Flux:} ";
      stream <<"None";
      stream <<",\\\\\n";
      return stream.str();
    }

  private:
    using BaseType::diffFlux_;
  };


  // LDGLimitedAdvectionDiffusionOperator
  //------------------------------------

  /** \class DGLimitedAdvectionDiffusionOperator
   *  \brief Dual operator for NS equtions with a limiting
   *         of the numerical solution
   *
   *  \tparam ModelType Analytical model
   *  \tparam NumFlux Numerical flux
   *  \tparam polOrd Polynomial degree
   *  \tparam advection Advection
   */
  template< class Traits,
            bool advection = true >
  class LDGLimitedAdvectionDiffusionOperator
   : public Fem::SpaceOperatorInterface< typename Traits::DestinationType >
  {
    // choose ids for the three passes (including start pass) different to the tuple entries
    enum PassIdType { u = Traits::ModelType::modelParameterSize,
                      limitPassId  = u + 1,
                      gradPassId   = u + 2,
                      advectPassId = u + 3 };

    struct GradientTraits : public Traits
    {
      static const int dimRange = Traits::ModelType::Traits::dimGradRange ;
      // overload discrete function space
      typedef typename Traits :: DiscreteFunctionSpaceType ::template ToNewDimRange< dimRange > :: Type  DiscreteFunctionSpaceType;
      // set new discrete function type
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >  DestinationType;
    };

    typedef Fem::SpaceOperatorInterface< typename Traits::DestinationType > BaseType;

  public:
    typedef typename Traits::ModelType                                  ModelType;
    typedef typename Traits::AvectionFluxType                           AdvectionFluxType;
    typedef typename ModelType::ProblemType                             ProblemType;

    enum { dimRange  = ModelType::dimRange };
    enum { dimDomain = ModelType::Traits::dimDomain };
    enum { dimGradRange = ModelType::Traits::dimGradRange };

    // Pass 2 Model (advectPassId)
    typedef AdvectionDiffusionLDGModel < Traits, limitPassId, gradPassId, advection, true >
                                                                        DiscreteModel3Type;

    // Pass 1 Model (gradPassId)
    typedef typename DiscreteModel3Type::DiffusionFluxType              DiffusionFluxType;
    typedef GradientModel< GradientTraits, limitPassId >                DiscreteModel2Type;
    // The model of the limiter pass (limitPassId)
    typedef Fem::StandardLimiterDiscreteModel< Traits, ModelType, u >   LimiterDiscreteModelType;

    // Pass 0 Model (limitPassId)
    typedef LimiterDiscreteModelType                                    DiscreteModel1Type;


    typedef typename DiscreteModel1Type::Traits                         Traits1;
    typedef typename DiscreteModel2Type::Traits                         Traits2;
    typedef typename DiscreteModel3Type::Traits                         Traits3;

    typedef typename Traits::GridType                                   GridType;

    //typedef typename Traits3::DomainType                                DomainType;
    typedef typename Traits3::DiscreteFunctionType                      DiscreteFunction3Type;

    typedef typename Traits1::DiscreteFunctionSpaceType                 Space1Type;
    typedef typename Traits2::DiscreteFunctionSpaceType                 Space2Type;
    typedef typename Traits3::DiscreteFunctionSpaceType                 Space3Type;
    typedef typename Traits2::DestinationType                           Destination2Type;
    typedef typename Traits3::DestinationType                           Destination3Type;
    typedef Destination3Type                                            DestinationType;
    typedef Space3Type                                                  SpaceType;

    typedef typename Traits2::GridPartType                              GridPartType;

    typedef Fem::StartPass< DiscreteFunction3Type, u >                  Pass0Type;
    typedef LimitDGPass< DiscreteModel1Type, Pass0Type, limitPassId >   Pass1Type;
    typedef LocalCDGPass< DiscreteModel2Type, Pass1Type, gradPassId >   Pass2Type;
    typedef LocalCDGPass< DiscreteModel3Type, Pass2Type, advectPassId > Pass3Type;

    typedef typename Traits::LimiterIndicatorType                       LimiterIndicatorType;
    typedef typename LimiterIndicatorType::DiscreteFunctionSpaceType    LimiterIndicatorSpaceType;
    typedef typename Traits::ExtraParameterTupleType                    ExtraParameterTupleType;

    template <class Limiter, int pOrd>
    struct LimiterCall
    {
      template <class ArgumentType, class DestinationType>
      static inline void limit(const Limiter& limiter,
                               ArgumentType* arg,
                               DestinationType& dest)
      {
        limiter.enableFirstCall();
        assert( arg );
        arg->assign(dest);
        limiter(*arg,dest);
        limiter.disableFirstCall();
      }
    };

    template <class Limiter>
    struct LimiterCall<Limiter,0>
    {
      template <class ArgumentType, class DestinationType>
      static inline void limit(const Limiter& limiter,
                               const ArgumentType* arg,
                               DestinationType& dest)
      {
      }
    };

  public:
    template< class ExtraParameterTupleImp >
    LDGLimitedAdvectionDiffusionOperator( GridPartType& gridPart, const ModelType& model,
                                          ExtraParameterTupleImp tuple,
                                          const std::string keyPrefix = "" )
      : model_( model )
      , numflux_( model_ )
      , gridPart_( gridPart )
      , space1_( gridPart_ )
      , space2_( gridPart_ )
      , space3_( gridPart_ )
      , uTmp_( (Traits::polynomialOrder > 0) ? (new DestinationType("limitTmp", space3_)) : 0 )
      , fvSpc_( gridPart_ )
      , indicator_( "Indicator", fvSpc_ )
      , diffFlux_( gridPart_, model_ )
      , discreteModel1_( model_ )
      , discreteModel2_( model_, diffFlux_ )
      , discreteModel3_( model_, numflux_, diffFlux_ )
      , pass0_()
      , pass1_(discreteModel1_, pass0_, space1_)
      , pass2_(discreteModel2_, pass1_, space2_)
      , pass3_(discreteModel3_, pass2_, space3_)
    {
      discreteModel1_.setIndicator( &indicator_ );
    }

    ~LDGLimitedAdvectionDiffusionOperator() { delete uTmp_; }

    void setTime(const double time) {
	    pass3_.setTime( time );
    }

    double timeStepEstimate() const {
	    return pass3_.timeStepEstimate();
    }

    void operator()( const DestinationType& arg, DestinationType& dest ) const {
	    pass3_( arg, dest );
      pass1_.enableFirstCall();
    }

    inline const SpaceType& space() const {
	    return space3_;
    }

    inline void switchupwind()
    {
      diffFlux_.switchUpwind();
    }
    double limitTime() const
    {
      return pass1_.computeTime();
    }
    std::vector<double> limitSteps() const
    {
      return pass1_.computeTimeSteps();
    }
    const Pass1Type& limitPass() const
    {
      return pass1_;
    }

    inline void limit( const DestinationType& arg, DestinationType& dest) const
    {
      pass1_.enableFirstCall();
      LimiterCall< Pass1Type, Traits::polynomialOrder >::limit( pass1_, uTmp_, dest );
    }

    void printmyInfo(std::string filename) const {
	    std::ostringstream filestream;
            filestream << filename;
            std::ofstream ofs(filestream.str().c_str(), std::ios::app);
            ofs << "Limited LDG Op., polynomial order: " << Traits::poolynomialOrder << "\\\\\n\n";
            ofs.close();
    }

    std::string description() const
    {
      std::stringstream stream;
      stream <<" {\\bf LDG Diff. Op.}, dual form, order: " << Traits::poolynomialOrder+1
             <<", penalty: ";
      diffFlux_.diffusionFluxPenalty( stream );
      stream <<", {\\bf Adv. Flux:} ";
      /*if (FLUX==1)
        stream <<"LLF";
      else if (FLUX==2)
        stream <<"HLL";*/
      stream <<",\\\\\n";
      return stream.str();
    }

    using BaseType::discreteModel;

  private:
    ModelType                   model_;
    AdvectionFluxType           numflux_;
    GridPartType&               gridPart_;
    Space1Type                  space1_;
    Space2Type                  space2_;
    Space3Type                  space3_;
    mutable DestinationType*    uTmp_;
    LimiterIndicatorSpaceType   fvSpc_;
    LimiterIndicatorType        indicator_;

  protected:
    DiffusionFluxType   diffFlux_;

  private:
    DiscreteModel1Type  discreteModel1_;
    DiscreteModel2Type  discreteModel2_;
    DiscreteModel3Type  discreteModel3_;
    Pass0Type           pass0_;
    Pass1Type           pass1_;
    Pass2Type           pass2_;
    Pass3Type           pass3_;
  };

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  template <class Traits,
            bool advection, bool diffusion >
  struct AdaptationIndicatorTraits : public Traits
  {
    // choose ids for the two passes (including start pass) different to the tuple entries
    enum PassIdType { u = Traits::ModelType::modelParameterSize, cdgpass = u + 1 };

    typedef AdaptiveAdvectionDiffusionDGPrimalModel
      < Traits, u, advection, diffusion> DiscreteModelType;
  };

  // LDGAdaptationIndicatorOperator
  //------------------------------

  template< class OpTraits,
            bool advection, bool diffusion = false >
  struct LDGAdaptationIndicatorOperator : public
    DGAdvectionDiffusionOperatorBase<
       AdaptationIndicatorTraits< OpTraits, advection, diffusion > >
  {
    typedef AdaptationIndicatorTraits< OpTraits, advection, diffusion > Traits;
    typedef DGAdvectionDiffusionOperatorBase< Traits >                  BaseType;
    typedef typename BaseType::GridPartType                             GridPartType;
    typedef typename BaseType::ProblemType                              ProblemType;
    typedef typename BaseType::ExtraParameterTupleType                  ExtraParameterTupleType;

    // constructor: do not touch/delegate everything
    template< class ... Args>
    LDGAdaptationIndicatorOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {
      if ( Fem::Parameter::verbose() )
      {
        std::cerr <<"\nWARNING: The adaptation indicator based on Ohlberger's a-posteriori\n";
        std::cerr <<"         error estimator is not ment to be used with flux formulation.\n\n";
      }
    }

    std::string description() const
    {
      std::stringstream stream;
      discreteModel_.diffusionFlux().diffusionFluxName( stream );
      stream <<" {\\bf Adv. Op.} in primal formulation, order: " << Traits::polynomialOrder+1
             <<", $\\eta = ";
      discreteModel_.diffusionFlux().diffusionFluxPenalty( stream );
      stream <<"$, $\\chi = ";
      discreteModel_.diffusionFlux().diffusionFluxLiftFactor( stream );
      stream << "$, {\\bf Adv. Flux:} " << numflux_.name() << ",\\\\" << std::endl;
      return stream.str();
    }

  protected:
    using BaseType::discreteModel_;
    using BaseType::numflux_;
  };

}
}
#endif
