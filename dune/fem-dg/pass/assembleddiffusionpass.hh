#ifndef DUNE_ASSEMBLEDDIFFUSIONPASS_HH
#define DUNE_ASSEMBLEDDIFFUSIONPASS_HH

#ifndef HEADERCHECK
#error "Deprecated header, check usage"

#include <dune/fem-dg/pass/dgpass.hh>
#include <dune/fem-dg/operator/dg/assembled/cdgrowwise.hh>

namespace Dune {

  //! Concrete implementation of Pass for first hyperbolic systems using
  //! LDG
  template <class AdvectionDiscreteModelImp,
            class DiffusionDiscreteModelImp,
            class PreviousPassImp,
            bool advection,
            bool assembled,
            int passId >
  class AssembledAdvectionDiffusionDGPass :
    public LocalCDGPass<AdvectionDiscreteModelImp, PreviousPassImp, passId>
  {
  public:
    typedef AdvectionDiscreteModelImp AdvectionDiscreteModelType ;
    typedef DiffusionDiscreteModelImp DiffusionDiscreteModelType ;

    //- Typedefs and enums
    //! Base class for advection
    typedef LocalCDGPass<AdvectionDiscreteModelImp, PreviousPassImp, passId
      > AdvectionBaseType;

    //! Base class for diffusion
    typedef CDGRowwiseOperator< DiffusionDiscreteModelImp, PreviousPassImp, assembled, passId >
                DiffusionOperatorType;

    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    typedef AssembledAdvectionDiffusionDGPass<AdvectionDiscreteModelType, DiffusionDiscreteModelType,
                                     PreviousPassType, advection, assembled, passId> ThisType;

    // Types from the base class
    typedef typename AdvectionBaseType::Entity EntityType;
    typedef const EntityType ConstEntityType;

    typedef typename AdvectionBaseType::ArgumentType ArgumentType;

    typedef typename AdvectionBaseType :: SelectorType SelectorType ;

    // Types from the traits
    typedef typename AdvectionDiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename AdvectionDiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  protected:
    using AdvectionBaseType :: dest_ ;
    using AdvectionBaseType :: spc_ ;

    mutable DiffusionOperatorType diffusionOperator_;
    mutable const DestinationType* U_;
    mutable DestinationType* rightHandSide_;
    const bool constantDiffusion_;
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order
    AssembledAdvectionDiffusionDGPass(AdvectionDiscreteModelType& advProblem,
                             DiffusionDiscreteModelType& difProblem,
                             PreviousPassType& pass,
                             const DiscreteFunctionSpaceType& spc,
                             const int volumeQuadOrd =-1,
                             const int faceQuadOrd=-1)
    : AdvectionBaseType(advProblem, pass, spc, volumeQuadOrd, faceQuadOrd),
      diffusionOperator_(difProblem, pass, spc, "" ),
      U_( 0 ),
      rightHandSide_( ( assembled ) ? 0 :
          new DestinationType("AdvDiffPass::rhs", spc_ )),
      constantDiffusion_( difProblem.constantCoefficient() )
    {
    }

    //! Destructor
    virtual ~AssembledAdvectionDiffusionDGPass()
    {
      delete rightHandSide_;
    }

    //! print tex info
    void printTexInfo(std::ostream& out) const {
      AdvectionBaseType::printTexInfo(out);
      diffusionOperator_.printTexInfo(out);
      out << "AssembledAdvectionDiffusionDGPass: "
          << "\\\\ \n";
    }

    void switchUpwind()
    {
      diffusionOperator_.switchUpwind();
    }

    //! set time
    void setTime(const double time)
    {
      AdvectionBaseType::setTime( time );
      diffusionOperator_.setTime( time );
    }

    //! Estimate for the timestep size
    double timeStepEstimateImpl() const
    {
      const double advdt = (advection) ? AdvectionBaseType::timeStepEstimateImpl () :
                                         std::numeric_limits<double>::max() ;
      const double difdt = diffusionOperator_.timeStepEstimate ();

      // return minimal dt
      return std :: min( advdt , difdt );
    }

  protected:
    //! In the preparations, store pointers to the actual arguments and
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // prepare operator
      diffusionOperator_.prepare( arg, dest );

      // compute matrix if grid has been changed
      diffusionOperator_.computeMatrix( arg, ! constantDiffusion_ );

      // also clears dest
      AdvectionBaseType::prepare( arg, dest );

      // get pointer to U
      U_ = Fem :: Element<0> :: get(arg);

      // prepare again (compute is only done the first time)
      if( constantDiffusion_ )
      {
        diffusionOperator_.prepare( arg, dest );
      }

      // apply matrix multiplication
      if( assembled )
      {
        // apply diffusion part
        // !!! overwrites dest !!!
        assert( U_ );
        diffusionOperator_.applyGlobal( *U_ , dest );

        // we need sign -1
        dest *= -1.0 ;
      }
      else
      {
        assert( rightHandSide_ );
        assert( U_ );
        // store U because U changes during operator application
        rightHandSide_->assign( *U_ );
        // we need sign -1
        (*rightHandSide_) *= -1.0;
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      AdvectionBaseType::finalize( arg, dest );
      diffusionOperator_.finalize( arg, dest );

      // std::cout << diffusionOperator_.computeTime() << "Diffusion Time  matrix " << std::endl ;

      // reset pointer
      U_ = 0 ;
    }

    void applyLocal(ConstEntityType& entity) const
    {
      if( assembled )
      {
        // apply diffusion boundary conditions
        diffusionOperator_.applyBoundary( entity );
      }
      else
      {
        assert( dest_ );
        assert( rightHandSide_ );

        // apply diffusion part
        diffusionOperator_.applyLocal( entity, *rightHandSide_ , *dest_ );
      }

      // if advection is enabled
      if( advection )
      {
        // apply advection part
        AdvectionBaseType :: applyLocal( entity );
      }
      else
      {
        // apply only inverse mass part
        AdvectionBaseType :: applyLocalMass( entity );
      }
    }

  private:
    AssembledAdvectionDiffusionDGPass();
    AssembledAdvectionDiffusionDGPass(const AssembledAdvectionDiffusionDGPass&);
    AssembledAdvectionDiffusionDGPass& operator=(const AssembledAdvectionDiffusionDGPass&);
  };
//! @}
} // end namespace Dune

#endif

#endif
