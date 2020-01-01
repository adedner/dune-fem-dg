#ifndef DUNE_FEM_DG_PY_OPERATOR_HH
#define DUNE_FEM_DG_PY_OPERATOR_HH

// system includes
#include <string>

// dune-fem includes
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#include <dune/fem-dg/algorithm/evolution.hh>

// dune-fem-dg includes
#include <dune/fem-dg/algorithm/evolution.hh>
#include <dune/fem-dg/operator/fluxes/advection/fluxes.hh>
#include <dune/fem-dg/operator/fluxes/euler/fluxes.hh>
#include <dune/fem-dg/operator/dg/operatortraits.hh>
#include <dune/fem-dg/operator/dg/primaloperator.hh>
#include <dune/fem-dg/models/modelwrapper.hh>
#include <dune/fem-dg/misc/algorithmcreatorselector.hh>

#ifdef EULER_WRAPPER_TEST
#include <dune/fem-dg/models/additional.hh>
#endif

#if HAVE_DUNE_FEMPY
#include <dune/fempy/quadrature/fempyquadratures.hh>
#endif

namespace Dune
{
namespace Fem
{

  // DG solver operator
  //---------------------

  template < class DestinationImp,
             class AdvectionModel,
             class DiffusionModel,
             class Additional>
#ifdef EULER_WRAPPER_TEST
#error
  class DGOperator : public DuneODE :: OdeSolverInterface< DestinationImp >
#else
  class DGOperator : public Fem::SpaceOperatorInterface< DestinationImp >
#endif
  {
  public:
    static const Solver::Enum solverId             = Additional::solverId;
    static const Formulation::Enum formId          = Additional::formId;
    static const AdvectionLimiter::Enum limiterId  = Additional::limiterId;
    static const AdvectionLimiterFunction::Enum limiterFunctionId = Additional::limiterFunctionId;
    // for valid advection fluxes see dune/fem-dg/operator/fluxes/advection/parameters.hh
    static const AdvectionFlux::Enum advFluxId     = Additional::advFluxId;
    // for valid diffusion fluxes see dune/fem-dg/operator/fluxes/diffusion/parameters.hh
    static const DiffusionFlux::Enum diffFluxId    = Additional::diffFluxId;
    typedef DestinationImp   DestinationType;
    typedef typename DestinationType :: DiscreteFunctionSpaceType    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType  FunctionSpaceType;
    static const int polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder ;

    typedef typename DiscreteFunctionSpaceType :: GridPartType    GridPartType;
    typedef typename GridPartType::GridType                       GridType;

    typedef typename GridType :: CollectiveCommunication          CollectiveCommunicationType;
    typedef TimeProvider< CollectiveCommunicationType >           TimeProviderType;

#ifdef EULER_WRAPPER_TEST
    typedef U0Sod< GridType > Problem;
#endif

    typedef typename AdvectionLimiterFunctionSelector<
      typename FunctionSpaceType::DomainFieldType, limiterFunctionId > :: type
      LimiterFunctionType;

    typedef ModelWrapper< GridType, AdvectionModel, DiffusionModel, Additional,
                          LimiterFunctionType
#ifdef EULER_WRAPPER_TEST
        , Problem
#endif
      >  ModelType;

    static constexpr bool symmetric  =  false ;
    static constexpr bool matrixfree =  true  ;
    static constexpr bool threading  = Additional::threading;

    typedef DGAdvectionFlux< ModelType, advFluxId >       AdvectionFluxType;
    typedef typename DiffusionFluxSelector< ModelType, DiscreteFunctionSpaceType, diffFluxId, formId >::type  DiffusionFluxType;

    typedef DefaultOperatorTraits< ModelType, DestinationType, AdvectionFluxType, DiffusionFluxType,
                std::tuple<>, typename DiscreteFunctionSpaceType::FunctionSpaceType,
#if HAVE_DUNE_FEMPY
                Dune::FemPy::FempyQuadratureTraits, // use quadratures from dune-fempy
#else
                Dune::Fem::DefaultQuadratureTraits,
#endif
                threading >  OpTraits;

    typedef AdvectionDiffusionOperatorSelector< OpTraits, formId, limiterId > OperatorSelectorType ;

    typedef typename OperatorSelectorType :: FullOperatorType      FullOperatorType;
    typedef typename OperatorSelectorType :: ExplicitOperatorType  ExplicitOperatorType;
    typedef typename OperatorSelectorType :: ImplicitOperatorType  ImplicitOperatorType;

    // solver selection, available fem, istl, petsc, ...
    typedef typename MatrixFreeSolverSelector< solverId, symmetric > :: template LinearInverseOperatorType< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType >  LinearSolverType ;

    static std::string name() { return std::string(""); }

    DGOperator( const DiscreteFunctionSpaceType& space,
                const AdvectionModel &advectionModel,
                const DiffusionModel &diffusionModel,
                const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : space_( space ),
        extra_(),
        model_( advectionModel, diffusionModel ),
        fullOperator_( space.gridPart(), model_, extra_, name(), parameter ),
        explOperator_( space.gridPart(), model_, extra_, name(), parameter ),
        implOperator_( space.gridPart(), model_, extra_, name(), parameter )
    {
    }

    virtual void description( std::ostream&) const {}

    const DiscreteFunctionSpaceType& space ()       const { return space_; }
    const DiscreteFunctionSpaceType& domainSpace () const { return space(); }
    const DiscreteFunctionSpaceType& rangeSpace ()  const { return space(); }

    FullOperatorType&     fullOperator()     const { return fullOperator_; }
    ExplicitOperatorType& explicitOperator() const { return explOperator_; }
    ImplicitOperatorType& implicitOperator() const { return implOperator_; }
    std::tuple<int,int,int> counter() const
    {
      return {fullOperator_.counter(),explOperator_.counter(),implOperator_.counter()};
    }

    //! evaluate the operator, which always referrers to the fullOperator here
    void operator()( const DestinationType& arg, DestinationType& dest ) const
    {
      fullOperator_( arg, dest );
    }

    /// Methods from SpaceOperatorInterface ////

    bool hasLimiter() const { return explOperator_.hasLimiter(); }

    /** \copydoc SpaceOperatorInterface::limit */
    void limit (const DestinationType& arg, DestinationType& dest) const
    {
      if( hasLimiter() )
      {
        explOperator_.limit( arg, dest );
      }
    }

    /** \copydoc SpaceOperatorInterface::setTime */
    void setTime( const double time )
    {
      fullOperator_.setTime( time );
    }

    double timeStepEstimate() const { return fullOperator_.timeStepEstimate(); }

    //// End Methods from SpaceOperatorInterface /////

  protected:
    const DiscreteFunctionSpaceType&            space_;

    std::tuple<>                          extra_;
    ModelType                             model_;
    mutable FullOperatorType              fullOperator_;
    mutable ExplicitOperatorType          explOperator_;
    mutable ImplicitOperatorType          implOperator_;
    mutable double                        fixedTimeStep_ ;
    mutable bool                          initialized_;
  };

} // end namespace Fem
} // end namespace Dune
#endif
