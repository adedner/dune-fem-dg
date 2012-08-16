#ifndef DUNE_CDGROWWISEOPERATOR_HH
#define DUNE_CDGROWWISEOPERATOR_HH

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem-dg/operator/dg/assembled/rowwisematrix.hh>
#include <dune/fem-dg/operator/dg/assembled/cdgprimaloperator.hh>

namespace Dune {

/*! @ingroup EllipticOperator
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   div(A(x)\nabla u) &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/

  //! non-assembled version 
  template <class DiscreteFunctionSpaceType, bool assembled> 
  struct CDGRowwiseTraits 
  {
    typedef RowwiseMatrixTraits< 
                  DiscreteFunctionSpaceType, 
                  DiscreteFunctionSpaceType > MatrixTraits; 

  };

  //! assembled version 
  template <class DiscreteFunctionSpaceType> 
  struct CDGRowwiseTraits<DiscreteFunctionSpaceType, true> 
  {
    typedef Fem::SparseRowMatrixTraits< 
                  DiscreteFunctionSpaceType, 
                  DiscreteFunctionSpaceType > MatrixTraits; 
  };

  ////////////////////////////////////////////////////////////
  //
  //  --CDGPrimalOperator 
  //
  ////////////////////////////////////////////////////////////
  /** \brief Operator assembling matrix for DG methods for elliptic problems. 
      Currently implemented are:
        - Interior Penalty
        - Baumann-Oden
        - NIPG 
        - Compact LDG (CDG) 

        References:
          The first 4 methods can be found for example in:
            D.N. Arnold, F. Brezzi, B. Cockburn, L.D. Marini: Unified
            Analysis of Discontinuous Galerkin Methods for Elliptic
            Problems SIAM J. Num. Anal, 39  (2002), 1749-1779.
            http://www.imati.cnr.it/~marini/reports/dgm_siam.ps.gz

          The Compact LDG method is described in detail in:
            J. Peraire and P.-O. Persson,
            The Compact Discontinuous Galerkin (CDG) Method for Elliptic Problems.
            SIAM J. Sci. Comput., to appear.
            http://www.mit.edu/~persson/pub/peraire07cdg.pdf
  */
  template <class DiscreteModelImp,
            class PreviousPassImp, 
            bool assembled , 
            int passId >
  class CDGRowwiseOperator 
    : public CDGPrimalOperatorImpl< DiscreteModelImp, 
                                    PreviousPassImp, 
                                    // matrix traits 
                                    typename CDGRowwiseTraits< 
                                      typename DiscreteModelImp::Traits::DiscreteFunctionSpaceType, 
                                      assembled 
                                    >::MatrixTraits, 
                                    passId > 
  {
  public:
    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;

    //! I need to switch PreviousPassType
    typedef PreviousPassImp PreviousPassType;

    //! type of discrete function space 
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // matrix traits 
    typedef typename CDGRowwiseTraits< 
        DiscreteFunctionSpaceType, assembled >::MatrixTraits  MatrixObjectTraits; 
    
    //- Typedefs and enums
    //! Base class
    typedef CDGPrimalOperatorImpl< DiscreteModelType, 
                                   PreviousPassType, 
                                   MatrixObjectTraits, 
                                   passId >  BaseType ;

    //! type of this class 
    typedef CDGRowwiseOperator<DiscreteModelType,
                PreviousPassType, assembled, passId > ThisType;

  private:  
    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef const EntityType ConstEntityType;

    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    // needed for InverseOperator
    typedef DestinationType DomainType ;
    typedef DestinationType RangeType ;

    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef Fem::CombinedSelector< ThisType , SelectorType >  CombinedSelectorType;
    typedef EllipticDiscreteModelCaller< DiscreteModelType, ArgumentType,
              CombinedSelectorType> DiscreteModelCallerType;

    //my Typedefs
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    
    using BaseType :: spc_;
    using BaseType :: matrixObj_;

    typedef typename BaseType :: MatrixObjectType MatrixObjectType; 

    mutable double computeTime_;
  public:
    using BaseType :: applyBoundary ;

    //- Public methods
    /** \copydoc CDGPrimalOperator */ 
    CDGRowwiseOperator(DiscreteModelType& problem, 
                       PreviousPassType& pass, 
                       const DiscreteFunctionSpaceType& spc,
                       const std::string paramFile = "")
      : BaseType(problem, pass, spc, paramFile )
    {
    }

    //! Destructor
    virtual ~CDGRowwiseOperator() 
    {
    }

    //////////////////////////////////////////////////////
    //
    //  OEM Solver interface 
    //
    //////////////////////////////////////////////////////
    const ThisType& systemMatrix() const 
    {
      return *this;
    }

    //! mult method for OEM solver 
    virtual void multOEM(const double* arg, double* dest) const 
    {
      //std :: cout << "CDG-Rowwise :: multOEM \n";
      typedef Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
        TemporaryDiscreteFunctionType;

      TemporaryDiscreteFunctionType Arg ("CDG-Row::Arg" , spc_ , arg  );
      TemporaryDiscreteFunctionType Dest("CDG-Row::Dest", spc_ , dest );

      // call application operator 
      multiply( Arg, Dest );
    }

    //! ddot method for OEM solver 
    virtual double ddotOEM(const double* x, const double* y) const 
    {
      //std::cout << "CDG-Rowwise :: ddotOEM \n";
      typedef Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
        TemporaryDiscreteFunctionType;

      double dot = 0.0 ;
#if HAVE_MPI 
      {
        TemporaryDiscreteFunctionType X ("CDG-Row::X", spc_ , x );
        TemporaryDiscreteFunctionType Y ("CDG-Row::Y", spc_ , y );
        dot = X.scalarProductDofs( Y );
      }
#else 
      const int spcSize = spc_.size();
      for( int i = 0; i<spcSize; ++i) 
      {
        dot += x[i] * y[i]; 
      }
#endif

      //std::cout << "Finished CDG-Rowwise :: ddotOEM \n";
      return dot ;
    }

    //////////////////////////////////////////////////////////
    //  
    //  InverseOperator interface 
    //
    //////////////////////////////////////////////////////////

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg,
                       const bool rebuild )
    {
      if( assembled ) 
      {
        BaseType :: computeMatrix( arg, rebuild );
      }

      computeTime_ = 0.0 ;
    }

  protected:  
    //! application operator does multiplication with matrix 
    //! row wise with grid traversal 
    void multiply(const DestinationType& arg, DestinationType& dest) const 
    {
      // clear destination 
      dest.clear();

      // apply mult
      multiplyAdd( arg, dest );
    }

    //! application operator does multiplication with matrix 
    //! row wise with grid traversal 
    void multiplyAdd(const DestinationType& arg, DestinationType& dest) const 
    {
      if( assembled ) 
      {
        // apply matrix multiplication 
        matrixObj_.apply( arg, dest );
      }
      else 
      {
        //std::cout << "Call CDGRowwise operator () \n"; 
        typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
        const IteratorType end = spc_.end();
        for(IteratorType it = spc_.begin(); it != end; ++it)
        {
          applyLocal( *it , arg, dest );
        }
      }
    }

  public:  
    double computeTime() const { return computeTime_; }

    void applyGlobal(const DestinationType& arg,
                     DestinationType& dest ) const
    {
      if( assembled ) 
      {
        multiplyAdd( arg, dest ); 
      }
    }

    /////////////////////////////////
    //! apply operator on entity 
    void applyLocal(ConstEntityType& entity, 
                    const DestinationType& arg,
                    DestinationType& dest ) const
    {
      if( ! assembled )
      {
        Timer timer; 

        // call applyLocal of CDGPrimalOperator 
        // this assembles a row of the system matrix 
        BaseType :: applyLocal( entity );

        computeTime_ += timer.elapsed();

        // multiply row   
        MObject<MatrixObjectType, assembled> ::
            applyLocal(matrixObj_, entity, arg, dest );

        // clear matrices 
        matrixObj_.clear();
      }
    } // end apply local 

  protected:
    template <class MatrixObject, bool> 
    struct MObject
    {
      static inline void applyLocal(MatrixObject& mobj,
                                    ConstEntityType& entity,
                                    const DestinationType& arg,
                                    DestinationType& dest) 
      {
        assert( false );
        abort();
      }
    };

    template <class MatrixObject> 
    struct MObject<MatrixObject, false >
    {
      static inline void applyLocal(MatrixObject& mObj,
                                    ConstEntityType& entity,
                                    const DestinationType& arg,
                                    DestinationType& dest) 
      {
        // multiply row   
        mObj.applyLocal( entity, arg, dest );
      }
    };
  };
} // end namespace Dune
#endif
