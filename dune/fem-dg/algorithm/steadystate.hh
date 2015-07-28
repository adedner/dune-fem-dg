#ifndef DUNE_FEMDG_ALGORITHM_STEADYSTATEALGORITHM_HH
#define DUNE_FEMDG_ALGORITHM_STEADYSTATEALGORITHM_HH

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>

#include <dune/fem-dg/algorithm/base.hh>

namespace Dune
{

namespace Fem
{

  template< class Grid,
            class ProblemTraits,
            int polOrd,
            ExtraParameterTuple = std::tuple<> >
  struct SteadyStateTraits
  {
    enum { polynomialOrder = polOrd };

    // type of Grid
    typedef Grid                                      GridType;

    // Choose a suitable GridView
    typedef AdaptiveLeafGridPart< GridType >          HostGridPartType;

    typedef typename ProblemTraits::GridPartType      GridPartType;

    typedef typename ProblemTraits::AnalyticalTraits  ModelTraits;

    // obtain the problem dependent types, analytical context
    typedef typename ModelTraits::ModelType           ModelType;
    typedef typename ModelTraits::ProblemType         ProblemType;

    struct OperatorTraits :
      public Dune::PassTraits< ModelTraits, polynomialOrder, ModelType::dimRange >
    {
      static const int limiterPolynomialOrder = polynomialOrder;
      typedef ExtraParameterTuple ExtraParameterTupleType;
    };


    // type of discrete function space and discrete function
    typedef typename OperatorTraits::DestinationType                  DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;


    typedef typename OperatorTraits::LinearOperatorType               LinearOperatorType;
    typedef typename OperatorTraits::LinearInverseOperatorType        LinearInverseOperatorType;

    typedef Dune::DGAdvectionDiffusionOperator< OperatorTraits >      FullOperatorType;
    typedef DGPrimalMatrixAssembly< FullOperatorType >                AssembledOperatorType;

    // tpye of jacobian operator used in the nested newton loops
    typedef typename FullOperatorType::JacobianOperatorType           JacobianOperatorType;

    // type of IOTuple
    typedef Dune::tuple< DiscreteFunctionType * >                     IOTupleType;

    // type of restriction/prolongation projection for adaptive simulations
    typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
  };


};



    template< class Grid, class ProblemTraits, int polynomialOrder, class ExtraParameterTuple = std::tuple< > >
    class SteadyStateAlgorithm
      : public AlgorithmBase< SteadyStateTraits< Grid, ProblemTraits, polynomialOrder, ExtraParameterTuple > >
    {
      // my traits class
      typedef SteadyStateTraits< Grid, ProblemTraits, polynomialOrder, ExtraParameterTuple > Traits;

      typedef AlgorithmBase< Traits > BaseType;
      typedef SteadyStateAlgorithm< Grid, ProblemTraits, polynomialOrder > This;

    public:
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IOTupleType IOTupleType;
      typedef typename BaseType::SolverMonitorType SolverMonitorType;

      typedef typename Traits::HostGridPartType HostGridPartType;

      // initial data type
      typedef typename Traits::ProblemType ProblemType;

      // An analytical version of our model
      typedef typename Traits::ModelType ModelType;

      typedef typename Traits::GridPartType GridPartType;
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;

      // The DG space operator
      typedef typename Traits::FullOperatorType FullOperatorType;

      // type of steady state solver
      typedef typename Traits::InverseOperatorType InverseOperatorType;

      // type of analytical traits
      typedef typename Traits::AnalyticalTraits AnalyticalTraits;

      // type of discrete traits
      typedef typename Traits::DiscreteTraits DiscreteTraits;

      // Error Handling
      typedef typename AnalyticalTraits::EOCErrorIDs EOCErrorIDs;

      SteadyStateAlgorithm ( GridType &grid )
        : BaseType( grid ),
          problem_( AnalyticalTraits::problem() ),
          model_( problem() ),
          gridPart_( grid ),
          space_( gridPart_ ),
          solution_( "solution", space_ ),
          dgOperator_( space_, model(), DiscreteTraits::quadOrder ),
          eocIds_( AnalyticalTraits::initEoc() ),
          description_( "Steady state scheme,\nOperator: " + dgOperator_.description() +"\nProblem: " + problem().description() )
      {}

      //! return reference to discrete space
      DiscreteFunctionSpaceType &space () { return space_; }

      //! returns data prefix for EOC loops ( default is loop )
      virtual std::string dataPrefix () const
      {
        return problem().dataPrefix();
      }

      IOTupleType dataTuple ()
      {
        return IOTupleType( &solution_ );
      }

      SolverMonitorType solve ( int step )
      {
        SolverMonitorType monitor;
        DiscreteFunctionType rhs( "rhs", space_ );

        rhs.clear();
        solution_.clear();

        InverseOperatorType invOp( dgOperator_, NestedNewtonParameter() );

        auto callBack = [ this ] ( const double f ) { dgOperator_.setRelaxFactor( f ); };
        invOp( rhs, solution_, callBack );

        monitor.gridwidth_ = GridWidth::calcGridWidth( gridPart_ );
        monitor.totalNewtonIterations_ = invOp.iterations();
        monitor.totalLinearSolverIterations_ = invOp.linearIterations();

        monitor.finalize();

        return monitor;
      }


      //! finalize computation by calculating errors and EOCs
      void finalize ( const int eocloop )
      {
        // submit error to the FEM EOC calculator
        AnalyticalTraits::addEOCErrors( eocIds_, 0., solution_, model(), problem() );
      }

      std::string description () const { return description_; }

      const ProblemType &problem () const
      {
        assert( problem_ );
        return *problem_;
      }

      ProblemType &problem ()
      {
        assert( problem_ );
        return *problem_;
      }

      ModelType &model () { return model_; }
      const ModelType &model () const { return model_; }

    protected:
      ProblemType *problem_;
      ModelType model_;

      GridPartType gridPart_;      // reference to grid part, i.e. the leaf grid
      DiscreteFunctionSpaceType space_;    // the discrete function space
      DiscreteFunctionType solution_;

      FullOperatorType dgOperator_;
      EOCErrorIDs eocIds_;

      std::string description_;
    };

  }  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALGORITHM_STEADYSTATEALGORITHM_HH
