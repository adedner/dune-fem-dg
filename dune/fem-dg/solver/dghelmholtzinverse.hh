#ifndef DUNE_FEM_DG_SOLVER_DGHELMHOLTZ_INV_HH
#define DUNE_FEM_DG_SOLVER_DGHELMHOLTZ_INV_HH

#include <iostream>
#include <utility>

#include <dune/fem/solver/rungekutta/timestepcontrol.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/operator/dghelmholtz.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

#include <dune/fem/io/parameter.hh>

namespace Dune
{
  namespace Fem
  {
    namespace detail
    {
      template <class DiscreteFunction>
      class UpdatePreconditionerWrapper
      {
      public:
        typedef DiscreteFunction DestinationType ;
        typedef std::reference_wrapper< const DestinationType >  ConstReferenceType;
        typedef std::function< void( ConstReferenceType& ) > UpdatePreconditionerFunctionType;

        UpdatePreconditionerWrapper( const UpdatePreconditionerFunctionType& upPre )
          : updatePrecond_( upPre )
        {}

        void update( const DestinationType& ubar ) const
        {
          ConstReferenceType u( ubar );
          // call Python side update function
          updatePrecond_( u );
        }

      protected:
        const UpdatePreconditionerFunctionType& updatePrecond_;

      };
    }

    template< class SpaceOperator >
    class DGHelmholtzOperatorWithUpdate
    : public DGHelmholtzOperator< SpaceOperator >
    {
      typedef DGHelmholtzOperatorWithUpdate< SpaceOperator > ThisType;
      typedef DGHelmholtzOperator< SpaceOperator >           BaseType;

    public:
      typedef typename BaseType :: DomainFunctionType    DomainFunctionType;
      typedef typename BaseType :: JacobianOperatorType  JacobianOperatorType;
      typedef typename BaseType :: SpaceOperatorType     SpaceOperatorType ;

      typedef detail::UpdatePreconditionerWrapper< DomainFunctionType > UpdatePreconditionerWrapperType;

      using BaseType :: spaceOperator;
      using BaseType :: lambda;

      explicit DGHelmholtzOperatorWithUpdate ( SpaceOperatorType &spaceOp )
        : BaseType( spaceOp )
      {}

      // overload jacobian to call update function
      void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
      {
        // update preconditioner to new linearization point
        if( updatePreconditioner_ )
          updatePreconditioner_->update( u );

        spaceOperator().jacobian( u, jOp );
        jOp.setLambda( lambda() );
      }

      void bind( const UpdatePreconditionerWrapperType& upPre )
      {
        updatePreconditioner_ = &upPre;
      }

      void unbind() { updatePreconditioner_ = nullptr; }

    protected:
      const UpdatePreconditionerWrapperType* updatePreconditioner_ = nullptr;
    };


    /** \class DGHelmholtzInverseOperator
     *  \brief Operator implementing the solution of
     *
     *  @f[
     *  L[\bar{u} + \lambda u ] = rhs
     *  @f]
     *
     *  \tparam SpaceOperator operator implementing L.
     *  \tparam DiscreteFunction type of discrete function representing u.
     *
     */
    template <class SpaceOperator, class DiscreteFunction>
    class DGHelmholtzInverseOperator
      : public virtual Operator< DiscreteFunction, DiscreteFunction >
    {
    public:
      typedef SpaceOperator SpaceOperatorType;
      typedef DiscreteFunction DestinationType;
      typedef typename DestinationType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef std::reference_wrapper< const DestinationType >  ConstReferenceType;
      typedef std::reference_wrapper< DestinationType >        ReferenceType;
      typedef std::function< void( ConstReferenceType& , ReferenceType& ) > PreconditionerType;

      typedef detail::UpdatePreconditionerWrapper< DestinationType > UpdatePreconditionerWrapperType;
      typedef typename UpdatePreconditionerWrapperType :: UpdatePreconditionerFunctionType UpdatePreconditionerType;

    protected:
      typedef Dune::Fem::GmresInverseOperator< DestinationType >                LinearInverseOperatorType;

      typedef Dune::Fem::DGHelmholtzOperatorWithUpdate< SpaceOperatorType >     HelmholtzOperatorType;
      typedef Dune::Fem::NewtonInverseOperator< typename HelmholtzOperatorType::JacobianOperatorType,
                                                LinearInverseOperatorType >     NonlinearInverseOperatorType;

      struct SolverInfo
      {
        SolverInfo ( bool converged, int linearIterations, int nonlinearIterations, double residualNorm )
          : converged( converged ), linearIterations( linearIterations ), nonlinearIterations( nonlinearIterations ), residualNorm( residualNorm )
        {}

        bool converged;
        int linearIterations, nonlinearIterations;
        double residualNorm;
      };

      // class wrapping a Python function to look like a Fem::Operator
      class PreconditionerWrapper : public virtual Operator< DestinationType, DestinationType >
      {
        const PreconditionerType& pre_; // function given from Python side
      public:
        PreconditionerWrapper( const PreconditionerType& pre )
          : pre_( pre ) {}

        virtual void operator() ( const DestinationType &u, DestinationType &v ) const final override
        {
          // convert to reference_wrapper to avoid copying
          ConstReferenceType uR( u );
          ReferenceType vR( v );

          // callback to python applying preconditioner
          pre_( uR, vR );
        }
      };
      typedef PreconditionerWrapper PreconditionerWrapperType;

    public:
      DGHelmholtzInverseOperator( SpaceOperatorType& op,
                                  const Dune::Fem::ParameterReader& parameter = Dune::Fem::Parameter::container() )
        : op_( op ),
          helmholtzOp_( op ),
          invOp_( parameter )
      {}

      /** \brief set lambda */
      void setLambda( const double lambda )
      {
        helmholtzOp_.setLambda( lambda );
      }

      /** \brief set Newton tolerance */
      void setTolerance( const double tol )
      {
        typedef DestinationType RangeFunctionType;
        auto finished = [ tol ] ( const RangeFunctionType &w, const RangeFunctionType &dw, double res ) { return res < tol; };
        invOp_.setErrorMeasure( finished );
        invOp_.eisenstatWalker().setTolerance( tol );
      }

      /** solve
       *
       * @f[
       * L[\bar{u} + \lambda u ] = rhs
       * @f]
       *
       */
      virtual void operator() ( const DestinationType &rhs, DestinationType &u ) const final override
      {
        // lambda needs to be set beforehand
        // bind operator (this will not overwrite an already set preconditioner)
        invOp_.bind( helmholtzOp_ );
        invOp_( rhs, u );
        invOp_.unbind();
      }

      /** \brief Solve
       *
       * @f[
       * L[\bar{u} + \lambda u ] = rhs
       * @f]
       *
       *  \param[in]     rhs     right hand side of F(w) = rhs
       *  \param[inout]  u       initial guess and returned solution
       *  \param[in]     lambda  lambda for Helmholtz operator
       */
      SolverInfo solve( const DestinationType& rhs, DestinationType &u, const double lambda ) const
      {
        helmholtzOp_.setLambda( lambda );

        // apply inv op
        (*this)( rhs, u );

        helmholtzOp_.unbind();
        return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations(), invOp_.residual() );
      }

      /** \brief Preconditioned solve
       *
       * @f[
       * L[\bar{u} + \lambda u ] = rhs
       * @f]
       *
       *  \param[in]     pre     preconditioner passed from Python side
       *  \param[in]     upPre   update preconditioner function passed from Python side
       *  \param[in]     rhs     right hand side of F(w) = rhs
       *  \param[inout]  u       initial guess and returned solution
       *  \param[in]     lambda  lambda for Helmholtz operator
       */
      SolverInfo preconditionedSolve( const PreconditionerType& pre,
                                      const UpdatePreconditionerType& upPre,
                                      const DestinationType& rhs,
                                      DestinationType &u,
                                      const double lambda ) const
      {
        PreconditionerWrapperType p( pre );
        UpdatePreconditionerWrapperType up( upPre );

        helmholtzOp_.bind( up );

        // bind operator and preconditioner
        invOp_.bind( helmholtzOp_, p );
        return solve( rhs, u, lambda );
      }

      std::pair< int, int > iterations() const
      {
        return std::make_pair( invOp_.iterations(), invOp_.linearIterations() );
      }

    protected:
      SpaceOperatorType& op_;
      mutable HelmholtzOperatorType helmholtzOp_;
      mutable NonlinearInverseOperatorType invOp_;
    };

  } // end namespace Fem
} // end namespace Dune
#endif
