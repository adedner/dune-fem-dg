#ifndef DUNE_FEM_DG_PRIMALOPERATOR_HH
#define DUNE_FEM_DG_PRIMALOPERATOR_HH

// system includes
#include <string>

// dune-fem includes
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/limiter/limiterdiscretemodel.hh>
#include <dune/fem-dg/operator/limiter/limitpass.hh>
#include <dune/fem-dg/operator/dg/primaldiscretemodel.hh>
#include <dune/fem-dg/operator/dg/operatorbase.hh>
#include <dune/fem-dg/pass/dgpass.hh>
#include <dune/fem-dg/misc/parameterkey.hh>

namespace Dune
{
namespace Fem
{

  // CDGAdvectionDiffusionTraits
  //----------------------------

  template <class Traits,
            bool advection, bool diffusion >
  struct CDGAdvectionDiffusionTraits : public Traits
  {
    // choose ids different to the tuple entries
    enum { u = Traits::ModelType::modelParameterSize , cdgpass = u + 1 };

    typedef AdvectionDiffusionDGPrimalModel< Traits, advection, diffusion, u> DiscreteModelType;
  };

  /**
   * \brief advection diffusion operator for CDG
   *
   * \note This operator is based on the Pass-Concept
   *
   * \ingroup PassBased
   * \ingroup PassOperator
   */
  template< class OpTraits >
  struct DGAdvectionDiffusionOperator
  : public DGAdvectionDiffusionOperatorBase< CDGAdvectionDiffusionTraits< OpTraits, true, true > >
  {
    typedef CDGAdvectionDiffusionTraits< OpTraits, true, true > Traits;
    typedef DGAdvectionDiffusionOperatorBase< Traits >          BaseType;
    typedef typename BaseType::GridPartType                     GridPartType;
    typedef typename BaseType::ProblemType                      ProblemType;
    typedef typename BaseType::ExtraParameterTupleType          ExtraParameterTupleType;

    // constructor: do not touch/delegate everything
    template< class ... Args>
    DGAdvectionDiffusionOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

    std::string description() const
    {
      std::stringstream stream;
      discreteModel_.diffusionFlux().diffusionFluxName( stream );
      stream <<" {\\bf Diff. Op.} in primal formulation, order: " << Traits::polynomialOrder+1
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


  // DGAdvectionOperator
  //--------------------

  /**
   * \brief advection operator for CDG
   *
   * \note This operator is based on the Pass-Concept
   *
   * \ingroup PassBased
   * \ingroup PassOperator
   */
  template< class OpTraits >
  struct DGAdvectionOperator : public
    DGAdvectionDiffusionOperatorBase< CDGAdvectionDiffusionTraits< OpTraits, true, false > >
  {
    typedef CDGAdvectionDiffusionTraits< OpTraits, true, false > Traits;
    typedef DGAdvectionDiffusionOperatorBase< Traits >           BaseType;
    typedef typename BaseType::GridPartType                      GridPartType;
    typedef typename BaseType::ProblemType                       ProblemType;
    typedef typename BaseType::ExtraParameterTupleType           ExtraParameterTupleType;

    // constructor: do not touch/delegate everything
    template< class ... Args>
    DGAdvectionOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

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


  // DGDiffusionOperator
  //--------------------

  /**
   * \brief diffusion operator for CDG
   *
   * \note This operator is based on the Pass-Concept
   *
   * \ingroup PassBased
   * \ingroup PassOperator
   */
  template< class OpTraits >
  class DGDiffusionOperator : public
    DGAdvectionDiffusionOperatorBase< CDGAdvectionDiffusionTraits< OpTraits, false, true > >
  {
  public:
    typedef CDGAdvectionDiffusionTraits< OpTraits, false, true > Traits;
    typedef DGAdvectionDiffusionOperatorBase< Traits >           BaseType;
    typedef typename BaseType::GridPartType                      GridPartType;
    typedef typename BaseType::ProblemType                       ProblemType;
    typedef typename BaseType::ExtraParameterTupleType           ExtraParameterTupleType;

  private:
    using BaseType::discreteModel_;

  public:
    // constructor: do not touch/delegate everything
    template< class ... Args>
    DGDiffusionOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

    std::string description() const
    {
      std::stringstream stream;
      discreteModel_.diffusionFlux().diffusionFluxName( stream );
      stream <<" {\\bf Diff. Op.} in primal formulation, order: " << Traits::polynomialOrder+1
             <<", $\\eta = ";
      discreteModel_.diffusionFlux().diffusionFluxPenalty( stream );
      stream <<"$, $\\chi = ";
      discreteModel_.diffusionFlux().diffusionFluxLiftFactor( stream );
      stream <<"$, {\\bf Adv. Flux:} ";
      stream <<"None";
      stream <<",\\\\\n";
      return stream.str();
    }
  };

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  template <class Traits, bool advection, bool diffusion >
  struct CDGAdaptationIndicatorTraits : public Traits
  {
    // choose ids different to the tuple entries
    enum { u = Traits::ModelType::modelParameterSize, cdgpass = u + 1 };

    typedef AdaptiveAdvectionDiffusionDGPrimalModel< Traits, advection, diffusion, u> DiscreteModelType;
  };

  // DGAdaptationIndicatorOperator
  //------------------------------

  /**
   * \brief diffusion operator for CDG
   *
   * This is a helper operator class for indicating error estimates
   *
   * \note This operator is based on the Pass-Concept
   *
   * \ingroup PassBased
   * \ingroup PassOperator
   */
  template< class OpTraits,
            bool advection = OpTraits::ModelType::hasAdvection,
            bool diffusion = OpTraits::ModelType::hasDiffusion >
  struct DGAdaptationIndicatorOperator : public
    DGAdvectionDiffusionOperatorBase< CDGAdaptationIndicatorTraits< OpTraits, advection, diffusion > >
  {
    typedef CDGAdaptationIndicatorTraits< OpTraits, advection, diffusion > Traits;
    typedef DGAdvectionDiffusionOperatorBase< Traits >  BaseType;
    typedef typename BaseType::GridPartType             GridPartType;
    typedef typename BaseType::ProblemType              ProblemType ;
    typedef typename BaseType::ExtraParameterTupleType  ExtraParameterTupleType;

  public:
    // constructor: do not touch/delegate everything
    template< class ... Args>
    DGAdaptationIndicatorOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

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


  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////


  // DGLimitedAdvectionOperator
  //------------------------------------

  /**
   * \brief Advection operator for CDG with a limiting
   *         of the numerical solution
   *
   * \ingroup PassBased
   * \ingroup PassOperator
   *
   * \tparam Traits
   * \tparam advection Advection
   * \tparam diffusion Diffusion
   */
  template< class Traits,
            bool advection = true, bool diffusion = false>
  class DGLimitedAdvectionOperator :
    public Fem::SpaceOperatorInterface< typename Traits::DestinationType >
  {
    enum PassIdType { u, limitPassId, advectPassId };
    enum { polOrd = Traits::polynomialOrder };
    enum { limiterPolOrd = Traits::limiterPolynomialOrder };

    typedef Fem::SpaceOperatorInterface< typename Traits::DestinationType >       BaseType;

  public:
    typedef typename Traits::ExtraParameterTupleType                              ExtraParameterTupleType;

    typedef typename Traits::ModelType                                            ModelType;
    typedef typename Traits::AdvectionFluxType                                    AdvectionFluxType;
    enum { dimRange  = ModelType::dimRange };
    enum { dimDomain = ModelType::Traits::dimDomain };


    // The model of the advection pass (advectPassId)
    typedef AdvectionDiffusionDGPrimalModel< Traits, advection, diffusion, limitPassId > DiscreteModel1Type;

    typedef typename DiscreteModel1Type::DiffusionFluxType                        DiffusionFluxType;
    typedef typename DiscreteModel1Type::AdaptationType                           AdaptationType;


    typedef typename ModelType::ProblemType                                       ProblemType;
    typedef typename ModelType::Traits::GridType                                  GridType;

    typedef typename ModelType::DomainType                                        DomainType;

    typedef typename Traits::DiscreteFunctionSpaceType                            SpaceType;
    typedef typename Traits::DestinationType                                      DestinationType;
    typedef typename DestinationType::GridPartType                                GridPartType;


    //typedef Traits                                                              PassTraitsType;
    typedef PassTraits< Traits, limiterPolOrd, dimRange, GridPartType  >          LimiterTraitsType;
    // The model of the limiter pass (limitPassId)
    typedef Fem::StandardLimiterDiscreteModel< LimiterTraitsType, ModelType, u >  LimiterDiscreteModelType;

    typedef typename LimiterTraitsType::DestinationType                           LimiterDestinationType ;
    typedef typename LimiterDestinationType::DiscreteFunctionSpaceType            LimiterSpaceType;

#ifdef USE_SMP_PARALLEL
    typedef Fem::StartPass < DestinationType, u, NonBlockingCommHandle< DestinationType > > Pass0Type;
    typedef LimitDGPass    < LimiterDiscreteModelType, Pass0Type, limitPassId >             InnerPass1Type;
    typedef ThreadPass     < InnerPass1Type, Fem::ThreadIterator< GridPartType >, true >    Pass1Type;
    typedef LocalCDGPass   < DiscreteModel1Type, Pass1Type, advectPassId >                  InnerPass2Type;
    typedef ThreadPass     < InnerPass2Type,
              //Fem::DomainDecomposedIteratorStorage<GridPartType >,
              Fem::ThreadIterator< GridPartType >,
              true > Pass2Type;
#else
    typedef Fem::StartPass < DestinationType, u >                                 Pass0Type;
    typedef LimitDGPass    < LimiterDiscreteModelType, Pass0Type, limitPassId >   Pass1Type;
    typedef LocalCDGPass   < DiscreteModel1Type, Pass1Type, advectPassId >        Pass2Type;
#endif

    typedef Limiter< DestinationType > LimiterOperator;

    typedef typename LimiterDiscreteModelType::IndicatorType                      LimiterIndicatorType;
    typedef typename LimiterIndicatorType::DiscreteFunctionSpaceType              LimiterIndicatorSpaceType;

    template< class Limiter, int pO >
    struct LimiterCall
    {
      template <class ArgumentType, class DestinationType>
      static inline void limit(const Limiter& limiter,
                               ArgumentType& arg,
                               DestinationType& dest)
      {
        limiter(arg, dest);
      }
    };

    template <class Limiter>
    struct LimiterCall<Limiter,0>
    {
      template <class ArgumentType, class DestinationType>
      static inline void limit(const Limiter& limiter,
                               const ArgumentType& arg,
                               DestinationType& dest)
      {
      }
    };

    void createIndicator()
    {
      // if indicator output is enabled create objects
      if( Fem::Parameter::getValue<bool> ("femdg.limiter.indicatoroutput", false ) )
      {
        fvSpc_     = new LimiterIndicatorSpaceType( gridPart_ );
        indicator_ = new LimiterIndicatorType( "SE", *fvSpc_ );
        limitDiscreteModel_.setIndicator( indicator_ );
      }
    }

  public:
    template< class ExtraParameterTupleImp >
    DGLimitedAdvectionOperator( GridPartType& gridPart,
                                const ModelType& model,
                                ExtraParameterTupleImp tuple,
                                const std::string name = "" )
      : model_( model )
      , advflux_( model_ )
      , gridPart_( gridPart )
      , space_( gridPart_ )
      , limiterSpace_( gridPart_ )
      , uTmp_( 0 )
      , fvSpc_( 0 )
      , indicator_( 0 )
      , diffFlux_( gridPart_, model_, DGPrimalDiffusionFluxParameters( ParameterKey::generate( name, "dgdiffusionflux." ) ) )
      , discreteModel1_( model_, advflux_, diffFlux_ )
      , limitDiscreteModel_( model_ , space_.order() )
      , pass0_()
      , pass1_( limitDiscreteModel_, pass0_, limiterSpace_ )
      , pass2_( discreteModel1_, pass1_, space_ )
      , limOp_( space_ )
    {
      // create indicator if enabled
      createIndicator();
    }

    ~DGLimitedAdvectionOperator()
    {
      delete uTmp_; uTmp_ = 0;
      if( indicator_ )
      {
        delete indicator_; indicator_ = 0;
        delete fvSpc_;     fvSpc_     = 0;
      }
    }

    void activateLinear() const {
      limitPass().disable();
    }
    void deactivateLinear() const {
      limitPass().enable();
    }

    void setAdaptationHandler( AdaptationType& adHandle, double weight = 1 )
    {
#ifdef USE_SMP_PARALLEL
      pass2_.setAdaptation( adHandle, weight );
#else
      discreteModel1_.setAdaptation( adHandle, weight );
#endif
    }

    void setTime(const double time)
    {
	    pass2_.setTime( time );
    }

    double timeStepEstimate() const
    {
	    return pass2_.timeStepEstimate();
    }

    void operator()( const DestinationType& arg, DestinationType& dest ) const
    {
	    pass2_( arg, dest );
    }

    inline const SpaceType& space() const {
	    return space_;
    }

    inline void switchupwind()
    {}

    double limitTime() const
    {
      return limitPass().computeTime();
    }

    std::vector<double> limitSteps() const
    {
      return limitPass().computeTimeSteps();
    }

    inline double computeTime() const
    {
      return pass2_.computeTime();
    }

    inline size_t numberOfElements () const
    {
      return pass2_.numberOfElements();
    }

    const Pass1Type& limitPass() const
    {
      return pass1_;
    }

    // return pointer to indicator function
    LimiterIndicatorType* indicator() { return indicator_ ; }

    inline void limit( DestinationType& U ) const
    {
      // copy U to uTmp_
      if( polOrd > 0 )
      {
        if( ! uTmp_ )
          uTmp_ = new LimiterDestinationType("limitTmp", limiterSpace_);

        assert( uTmp_ );
        uTmp_->assign( U );

        limit( *uTmp_, U );
      }
    }

    inline void limit( const DestinationType& arg, DestinationType& U ) const
    {
      LimiterCall< Pass1Type, polOrd >::limit( limitPass(), arg, U );
    }

    void printmyInfo(std::string filename) const
    {
	    std::ostringstream filestream;
            filestream << filename;
            std::ofstream ofs(filestream.str().c_str(), std::ios::app);
            ofs << "Limited Adv. Op., polynomial order: " << polOrd << "\\\\\n\n";
            ofs.close();
    }

    std::string description() const
    {
      std::cerr <<"DGLimitedAdvectionOperator::description() not implemented" <<std::endl;
      abort();

      /*
      std::stringstream stream;
      stream <<", {\\bf Adv. Flux:} ";
      if (FLUX==1)
        stream <<"LLF";
      else if (FLUX==2)
        stream <<"HLL";
      stream <<",\\\\\n";
      return stream.str();
      */
    }

    const DiscreteModel1Type& discreteModel() const
    {
      return discreteModel1_;
    }

  private:
    ModelType                  model_;
    AdvectionFluxType          advflux_;
    GridPartType&              gridPart_;
    SpaceType                  space_;
    LimiterSpaceType           limiterSpace_;
    mutable LimiterDestinationType* uTmp_;

    LimiterIndicatorSpaceType* fvSpc_;
    LimiterIndicatorType*      indicator_;


  protected:
    DiffusionFluxType   diffFlux_;

  private:
    DiscreteModel1Type  discreteModel1_;
    LimiterDiscreteModelType  limitDiscreteModel_;
    Pass0Type           pass0_;
    Pass1Type           pass1_;
    Pass2Type           pass2_;

    LimiterOperator  limOp_;
  };


  /**
   * \brief Advection diffusion operator for CDG with a limiting
   *         of the numerical solution
   *
   * \ingroup PassBased
   * \ingroup PassOperator
   *
   * \tparam Traits
   * \tparam advection Advection
   * \tparam diffusion Diffusion
   */
  template< class Traits, bool advection = true, bool diffusion = true >
  class DGLimitedAdvectionDiffusionOperator
  : public DGLimitedAdvectionOperator< Traits, advection, diffusion >
  {
    typedef DGLimitedAdvectionOperator< Traits, advection, diffusion > BaseType;

    typedef typename BaseType::GridPartType                            GridPartType;
    typedef typename BaseType::ProblemType                             ProblemType;

  public:
    // constructor: do not touch/delegate everything
    template< class ... Args>
    DGLimitedAdvectionDiffusionOperator( Args&&... args )
    : BaseType( std::forward<Args>(args)... )
    {}

    void printmyInfo(std::string filename) const
    {
	    std::ostringstream filestream;
            filestream << filename;
            std::ofstream ofs(filestream.str().c_str(), std::ios::app);
            ofs << "Limited Adv. Diff. Op., polynomial order: " << Traits::polynomialOrder << "\\\\\n\n";
            ofs.close();
    }

    std::string description() const
    {
      std::cerr <<"DGLimitedAdvectionDiffusionOperator::description() not implemented" <<std::endl;
      abort();

      /*
      std::stringstream stream;
      stream <<", {\\bf Adv. Flux:} ";
      if (FLUX==1)
        stream <<"LLF";
      else if (FLUX==2)
        stream <<"HLL";
      stream <<",\\\\\n";
      return stream.str();
      */
    }
  };

}
}
#endif
