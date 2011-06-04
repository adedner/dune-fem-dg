#ifndef DUNE_CDGPRIMALOPERATOR_HH
#define DUNE_CDGPRIMALOPERATOR_HH

//- Dune includes 
#include <dune/common/typetraits.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

//- local includes 
#include <dune/fem/pass/pass.hh>
// use model caller from fem-dg 
#include <dune/fem-dg/pass/ellipticmodelcaller.hh>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>

#include <dune/fem/space/combinedspace.hh>

#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/arrays.hh>

#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/operator/2order/dgmatrixsetup.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

// eigen value calculation 
#include <dune/common/fmatrixev.hh>

#define USE_COMPACT_LDG 1 

namespace Dune {

// double feature only works in serial runs 
//#if HAVE_MPI == 0
//#define DG_DOUBLE_FEATURE 
//#endif

/*! @ingroup EllipticOperator
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   div(A(x)\nabla u) &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/
  ////////////////////////////////////////////////////////////
  //
  //  --CDGPrimalOperatorImpl 
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
            class MatrixObjectTraits, 
            int passId >
  class CDGPrimalOperatorImpl 
    : public LocalPass<DiscreteModelImp, PreviousPassImp, passId
      > 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp, passId
      > BaseType;

    typedef CDGPrimalOperatorImpl<DiscreteModelImp,
            PreviousPassImp,MatrixObjectTraits, passId > ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;

    //! I need to switch PreviousPassType
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef const EntityType ConstEntityType;

    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef FieldVector<double,DomainType::dimension-1> FaceDomainType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    enum { dim = GridType :: dimension };
    
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType::polynomialOrder };



    typedef FunctionSpace < DomainFieldType, 
                            RangeFieldType, 
                            dimDomain, 
                            dimGradRange > ContainedFunctionSpaceType ;

    typedef typename DiscreteFunctionSpaceType :: template ToNewFunctionSpace<
            ContainedFunctionSpaceType > :: Type DiscreteGradientSpaceType ;

    enum { elementMassSize = DiscreteFunctionSpaceType :: localBlockSize };
    enum { massSize = DiscreteGradientSpaceType :: localBlockSize };

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;

    typedef typename DiscreteModelType::SelectorType SelectorType;

    typedef CombinedSelector< ThisType , SelectorType >  CombinedSelectorType;
    typedef EllipticDiscreteModelCaller< DiscreteModelType, ArgumentType,
              CombinedSelectorType> DiscreteModelCallerType;

    typedef typename GridType :: ctype ctype;
    typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;
    
    //my Typedefs
    typedef DiscreteFunctionSpaceType SingleDiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType:: IteratorType IteratorType ;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    typedef typename DiscreteGradientSpaceType::RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType::JacobianRangeType GradJacobianRangeType;
    typedef GradJacobianRangeType GradientJacobianRangeType;
    enum { GradDimRange = GradientRangeType :: dimension };
    
    typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;
     // type of temporary local function belonging to lower space 
    typedef TemporaryLocalFunction< DiscreteGradientSpaceType > TemporaryLocalFunctionType;
    typedef MutableArray< std::auto_ptr< TemporaryLocalFunctionType > > TemporaryLocalFunctionArrayType;

    typedef DGMatrixTraits< MatrixObjectTraits > MyOperatorTraits;
    //! type of underlying matrix implementation 
    typedef typename MatrixObjectTraits :: template MatrixObject<
      MyOperatorTraits > :: MatrixObjectType MatrixObjectType;

    typedef typename MatrixObjectType::LocalMatrixType LocalMatrixType;
    typedef typename MatrixObjectType::MatrixType MatrixType;
    typedef typename MatrixObjectType::PreconditionMatrixType PreconditionMatrixType;
    
    typedef typename DiscreteModelType :: BoundaryIdentifierType BoundaryIdentifierType;    

    typedef FieldMatrix< RangeFieldType, dimDomain, dimDomain> FluxRangeType;
    typedef FieldVector< FluxRangeType , dimRange > CoeffMatrixType;

    // typedef GradJacobianRangeType FluxRangeType; 
    typedef TemporaryLocalFunction< DiscreteFunctionSpaceType > SingleLFType;
    //typedef typename DestinationType :: LocalFunctionType SingleLFType; 

    //! singleton list , key type is const pointer to grid 
    typedef GridWidthProvider< GridType > GridWidthType;
    typedef typename GridWidthType :: ProviderType GridWidthProviderType;

    typedef typename DiscreteGradientSpaceType :: RangeType GradRangeType;

    typedef FieldMatrix<double, massSize , massSize > MassMatrixType; 
    typedef FieldVector<double, massSize > MassVectorType; 

    //! type of local mass matrix 
    typedef LocalDGMassMatrix< DiscreteFunctionSpaceType, VolumeQuadratureType > RhsMassMatrixType;
    typedef LocalDGMassMatrix< DiscreteGradientSpaceType, VolumeQuadratureType > LocalMassMatrixType;

    ////////////////////////////////////////////////////
    //
    //  Coefficient and RHS caller 
    //
    ////////////////////////////////////////////////////
    template <class CallerType> 
    class CoefficientCallerTrue
    {
    public:  
      template <class QuadratureType, class CoeffType> 
      void evaluateCoefficient(CallerType& caller, 
                               ConstEntityType& en, 
                               QuadratureType& quad, 
                               const int l,
                               CoeffType& coeff) const 
      {
        caller.evaluateCoefficient(en, quad, l, coeff );
      }         

      // version for GradRangeType 
      void applyCoefficient(const CoeffMatrixType& coeffEn, 
                            const GradRangeType& psi,
                            GradRangeType& coeffPsi) const 
      {
        // reset value
        coeffPsi = 0;
        // apply coefficient matrix 
        for(int i=0; i<dimRange; ++i)
        {
          const FluxRangeType& matrix = coeffEn[i];
          const int idx = i * dimDomain;

          for(int j=0; j<dimDomain; ++j) 
          {
            for(int k=0; k<dimDomain; ++k) 
            {
              coeffPsi[ idx + j ] += matrix[ j ][ k ] * psi[ idx + k ];
            }
          }
        }
      }

      // version for JacobianRangeType 
      void applyCoefficient(const CoeffMatrixType& coeffEn, 
                            const JacobianRangeType& psi,
                            JacobianRangeType& coeffPsi) const 
      {
        for(int i=0; i<dimRange; ++i)
        {
          coeffEn[i].mv(psi[i], coeffPsi[i]);
        }
      }
    };

    template <class CallerType> 
    struct CoefficientCallerFalse
    {
      template <class QuadratureType, class CoeffType> 
      void evaluateCoefficient(const CallerType& caller, 
                               const EntityType& en, 
                               const QuadratureType& quad, 
                               const int l,
                               CoeffType& coeff) const 
      {
      }         

      // just copy in default case 
      template <class CoeffType, class PsiType> 
      void applyCoefficient(const CoeffType& coeffEn, 
                            const PsiType& psi,
                            PsiType& coeffPsi) const 
      {
        coeffPsi = psi;
      }
    };


    //! if right hand side available 
    template <class CallerType> 
    class CoefficientCallerRHS 
    {
      mutable SingleLFType& singleRhs_;
      const BaseFunctionSetType& bsetEn_;
      mutable RangeType rhsval_;

    public:
      CoefficientCallerRHS(SingleLFType& singleRhs)
        : singleRhs_ ( singleRhs )
        , bsetEn_( singleRhs_.baseFunctionSet()) 
        , rhsval_ (0.0)  
      {}

      template <class QuadratureType> 
      void rightHandSide(CallerType& caller, 
                         ConstEntityType& en, 
                         QuadratureType& volQuad, 
                         const int l, 
                         const double intel) const 
      {
        // eval rightHandSide function 
        // if empty, rhs stays 0.0
        caller.rightHandSide(en, volQuad, l, rhsval_ );

        // scale with intel 
        rhsval_ *= intel;

        // add to right hand side 
        singleRhs_.axpy(volQuad[l] , rhsval_);
      }
    };

    //! no rhs 
    template <class CallerType> 
    class CoefficientCallerNoRHS 
    {
    public:  
      template <class QuadratureType> 
      void rightHandSide(const CallerType& caller, 
                         const EntityType& en, 
                         const QuadratureType& volQuad, 
                         const int l,
                         const double intel) const 
      {
      }
    };

    template <class CallerType, bool hasCoeff, bool hasRHS> 
    class CoefficientCaller : public CoefficientCallerTrue<CallerType> ,
                              public CoefficientCallerRHS<CallerType> 
    {
      typedef CoefficientCallerRHS<CallerType> BaseType;
    public:
      CoefficientCaller(SingleLFType& singleRhs) : BaseType( singleRhs ) {}
    };

    template <class CallerType> 
    class CoefficientCaller<CallerType,false,true>
      : public CoefficientCallerFalse<CallerType> ,
        public CoefficientCallerRHS<CallerType> 
    {
      typedef CoefficientCallerRHS<CallerType> BaseType;
    public:
      CoefficientCaller(SingleLFType& singleRhs) : BaseType( singleRhs ) {}
    };

    template <class CallerType> 
    class CoefficientCaller<CallerType,true,false> 
      : public CoefficientCallerTrue<CallerType>
      , public CoefficientCallerNoRHS<CallerType>
    {
      public:
    };

    template <class CallerType> 
    class CoefficientCaller<CallerType,false,false> 
      : public CoefficientCallerFalse<CallerType>
      , public CoefficientCallerNoRHS<CallerType>
    {
      public:
    };

    bool hasLifting () const { return method_ < method_ip ; } 

    double getLiftFactor() const 
    {
      double liftFac = Parameter::getValue<double>("dgprimal.liftfactor");
      if( theoryParams_ ) 
      {
        //if( spc_.multipleGeometryTypes() ) 
        //  DUNE_THROW(NotImplemented,"theory parameter for hybrid grids not implemented");

        if( theoryParams_ == 1 ) 
        {
          const IteratorType end = spc_.end();
          double maxFaces = 0;
          for(IteratorType it = spc_.begin(); it != end; ++ it ) 
          {
            double faces =  0;
            const EntityType& en = *it ;
            const IntersectionIteratorType endnit = gridPart_.iend(en); 
            for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
            { 
              faces += 1;
            }

            maxFaces = std::max( faces, maxFaces );
          }
          liftFac = 0.5 * maxFaces ;
        }
        else
        {
          IteratorType begin = spc_.begin();
          if( begin != spc_.end() )
          {
            double numFaces = (*begin).template count< 1 >();
            liftFac = 0.5 * numFaces ;
          }
        }
      }
      return liftFac;
    }

  public:
    //- Public methods
    /**  \brief Constructor
     \param problem Actual problem definition (see problem.hh)
     \param[in]  gradPass  types and discrete model for the gradient pass
     \param pass Previous pass
     \param spc Space belonging to the discrete function local to this pass
     \param paramFile parameter file to read necessary parameters, if empty 
             default parameters will be applied 
    */         
    CDGPrimalOperatorImpl(DiscreteModelType& problem, 
                          PreviousPassType& pass, 
                          const DiscreteFunctionSpaceType& spc,
                          const std::string paramFile = "")
    : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      arg_(0),
      rhs_(0),
      uh_(0),
      spc_(spc),
      gridPart_(const_cast<GridPartType &> (spc_.gridPart())),
      gradientSpace_(gridPart_),
      gridWidth_ ( GridWidthProviderType :: getObject( &spc_.grid())),
      upwind_ ( M_PI ),
      localGradientMass_( gradientSpace_, (2 * spc_.order()) ),
      applyInverseMass_( ! localGradientMass_.affine () ),
      rhsFunc_( spc_ ),
      //volumeQuadOrd_( ( applyInverseMass_ ) ? (2 * spc_.order() + 2 ) : (2 * spc_.order())  ),
      //faceQuadOrd_( ( applyInverseMass_ ) ? (2 * spc_.order() + 3) : (2 * spc_.order() + 1) ),
      volumeQuadOrd_( (2 * spc_.order())  ),
      faceQuadOrd_( (2 * spc_.order() + 1) ),
      matrixObj_(spc_,spc_, paramFile ),
      coeffEn_(),
      coeffNb_(),
      inverseMassEn_(0),
      inverseMassNb_(0),
      matrixAssembled_(false),
      power_( 2*spc_.order() ),
      dtMin_ (std::numeric_limits<double>::max()),
      minLimit_(2.0*std::numeric_limits<double>::min()),
      sequence_ ( -1 ),
      numberOfElements_( 0 ),
      method_( DGPrimalMethodNames::getMethod() ),
      penalty_( Parameter::getValue<double>("dgprimal.penalty") ),
      penaltyNotZero_( std::abs( penalty_ ) > 0 ),
      bilinearPlus_( (method_ == method_nipg) || (method_ == method_bo) ? true : false ),
      areaSwitch_( method_ == method_cdg2 ),
      theoryParams_( Parameter::getValue< int >("dgprimal.theoryparams", 1 ) ),
      liftFactor_( getLiftFactor() ),
      affineGeoms_( true )
    {
      for(int i=0; i<dimRange; ++i) 
      {
        // set to unit matrix 
        setToUnitMatrix(coeffEn_[i], coeffNb_[i]);
      }

      // set upwind to some strange numbers 
      if( dimDomain > 1 ) upwind_[1] = M_LN2 ;
      if( dimDomain > 2 ) upwind_[2] = M_E ;

      if( ! (spc_.order() > 0))
      {
        std::cerr << "ERROR: DG Primal operator only working for spaces with polynomial order > 0! \n";
        assert(false);
        abort();
      }

      const bool output = (gridPart_.grid().comm().rank() == 0);

      if( output ) 
      {
        std::cout << "CDGPrimalOperatorImpl: using k=" << spc_.order() << "  "<< DGPrimalMethodNames::methodNames( method_ );
        if( theoryParams_ ) 
        {
          if( penaltyNotZero_ ) std::cout  <<" theory: (eta = " << penalty_ << ") "; 
        }
        else 
          if( penaltyNotZero_ ) std::cout  <<"   eta = " << penalty_; 

        if( liftFactor_ > 0 ) std::cout << "  chi = " << liftFactor_ ;
        std::cout << std::endl;
      }

      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );

      if(problem_.hasSource())
      {
        std::cerr << "Source for DGElliptPass not supported yet. \n";
        abort();
      }
      switchUpwind() ;
    }

    // set matrices to unit matrices 
    void setToUnitMatrix(FluxRangeType& coeffEn, FluxRangeType& coeffNb) const
    {
      for( int i=0; i<FluxRangeType :: rows ; ++i) 
      {
        // set diagonal to 1 
        coeffEn[i][i] = coeffNb[i][i] = 1;
        // all others to zero 
        for(int j=i+1; j< FluxRangeType :: cols; ++j) 
        {
          // set default value to fMat 
          coeffEn[i][j] = coeffNb[i][j] = 0;
          coeffEn[j][i] = coeffNb[j][i] = 0;
        }
      }
    }

    //! Destructor
    virtual ~CDGPrimalOperatorImpl() 
    {
    }

    // switch sign of upwind vector 
    void switchUpwind() 
    {
      upwind_ *= -1; 
    }

    // return number of elements 
    size_t numberOfElements () const 
    {
      return numberOfElements_ ; 
    }

    // set tau and theta 
    void setTauAndTheta(const double tau, const double theta)
    {
      tau_1_ = 1.0/tau;
      theta_ = theta;
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, 
                       const DestinationType &uh, 
                       DestinationType & rhs)
    {
      // store uh_ ;
      uh_ = &uh;

      computeMatrix( arg, rhs );
    }

    void buildRhs(DestinationType & rhs)
    {
      // rhs.clear();
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      const IteratorType end = spc_.end();
      rhs_ = & rhs ;
      for(IteratorType it = spc_.begin(); it != end; ++it)
      {
        ConstEntityType& entity = *it ; 
        if( entity.hasBoundaryIntersections() ) 
        {
          applyBoundary( *it );
        }
      }
      rhs_ = 0;
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, 
                       const bool rebuild )
    {
      if( rebuild || sequence_ != spc_.sequence() )
      {
        Timer timer ;

        // time initialisation to max value 
        dtMin_ = std::numeric_limits<double>::max();

        // prepare operator  
        prepare( arg );

        // if grid has changed, then matrix structure is re-build
        matrixObj_.reserve();

        // clear matrix 
        matrixObj_.clear();

        // reset number of elements 
        numberOfElements_ = 0; 
        // build matrix 
        const IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it) 
        {
          // assemble local matrices  
          applyLocal( *it );

          // increase element counter 
          ++numberOfElements_ ; 
        }

        /*
        // also calculate matrix for overlap (needed for YaspGrid)
        if( spc_.grid().overlapSize(0) > 0 )
        {
          typedef typename GridPartNewPartitionType<GridPartType,Overlap_Partition>:: NewGridPartType NewGridPartType;
          typedef typename NewGridPartType :: template Codim<0> :: IteratorType IteratorType;

          NewGridPartType gridPart( const_cast<GridPartType&> (gridPart_).grid() );
          IteratorType endit = gridPart. template end<0>();
          for(IteratorType it = gridPart. template begin<0>();
              it != endit; ++it)
          {
            applyLocal( *it );
          }
        }
        */

        // build matrix and rhs 
        matrixAssembled_ = true;
        
        // create pre-condition matrix if activated 
        matrixObj_.createPreconditionMatrix();
        
        // finalize 
        finalize( );

        // store current sequence 
        sequence_ = spc_.sequence();

        // store compute time 
        this->computeTime_ = timer.elapsed();

        // in verbose mode make output 
        if( Parameter :: verbose () )
        {
          std::cout << "Matrix assembly took " << this->computeTime_ << " seconds! \n";
        }
      }
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, 
                       DestinationType & rhs, const bool rebuild = true )
    {
      if( rebuild || sequence_ != spc_.sequence() )
      {
        //set right hand side 
        setRightHandSide( rhs );

        // compute matrix 
        computeMatrix( arg, rebuild );
      }
    }

  public:  
    //! do matrix vector multiplication, used by InverseOp  
    void operator () (const DestinationType & arg, DestinationType& dest) const 
    {
      matrixObj_.multOEM( arg.leakPointer(), dest.leakPointer() ); 
    }
    
  public:
    //set right hand side 
    void setRightHandSide( DestinationType& rhs ) const 
    {
      rhs.clear ();
      rhs_ = & rhs ;
    }

    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    void prepare(const ArgumentType& arg) const
    {
      // set argument 
      arg_ = const_cast<ArgumentType*>(&arg);
      caller_.setArgument(*arg_);

      // set time to caller 
      caller_.setTime( this->time() );

      // only calculate in case of beta not zero 
      /*
      if( penaltyNotZero_ )
      {
        // calculate beta = O(1/h)
        const double h = gridWidth_.gridWidth();
        globalBeta_ = penalty_ / h; 
      }
      */

      affineGeoms_ = true ;
    }

    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // prepare operator 
      prepare( arg );

      //set right hand side 
      setRightHandSide( dest );
    }

    //! Some timestep size management.
    virtual void finalize() const
    {
      caller_.finalize();
      rhs_ = 0;
      uh_  = 0;
      //std::cout << "CDGPrimal: geometries affine " << affineGeoms_ << std::endl;
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      finalize ();
    }

    //! Estimate for the timestep size 
    double timeStepEstimateImpl() const
    {
      // factor for LDG  Discretization 
      const double p  = 2 * spc_.order() + 1;
      return dtMin_ / p;
    }

    void printTexInfo(std::ostream& out) const 
    {
      out << "CDGPrimalOperator: ";
      out << " method = " << DGPrimalMethodNames::methodNames( method_ ) ; 
      if( penaltyNotZero_ ) out << " ( eta = " << penalty_ << " )"; 
      if( liftFactor_>0 ) out << " ( chi  = " << liftFactor_ << " )"; 
      out  << "\\\\ \n";
      matrixObj_.printTexInfo( out );
    }

  public:
    //! return reference to function space 
    const DiscreteFunctionSpaceType & space() const 
    {
      return spc_;
    }

    //! set all data belonging to this entity to zero 
    void restrictLocal(const EntityType& father, const EntityType & son, bool firstCall ) const
    {
    }

    //! set all data belonging to this entity to zero 
    void prolongLocal(const EntityType& father, const EntityType & son, bool ) const
    {
    }

protected:
    void resizeCaches(const size_t numDofs) const
    {
      // resize caches 
      if(tau_.size() != numDofs) 
      { 
        tau_.resize(numDofs);
        tauNeigh_.resize(numDofs);
        phi_.resize(numDofs);
        phiNeigh_.resize(numDofs);
        psi_.resize(numDofs);
        coeffPsi_.resize(numDofs);
      }
    }

    //! assemble right hand side modifications 
    void applyBoundary(ConstEntityType& en) const
    {
      // this method should not be called for ghost entities 
      assert( en.partitionType() != GhostEntity );
      
      // make entities known in callers
      caller_.setEntity(en);

      // get geometry
      const Geometry & geo = en.geometry();

      if( ! geo.affine() ) 
      {
        affineGeoms_ = false ;
      }

      // create volume quadrature  
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      // get base function set of single space 
      const BaseFunctionSetType bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      assert( numDofs > 0 );
      // resize caches 
      resizeCaches(numDofs);

      // calculate inverse mass matrix for gradient space 
      calculateInverseMass(en, geo, volQuad);

      rhsFunc_.init( en );
      rhsFunc_.clear();

      /////////////////////////////////
      // Surface integral part
      /////////////////////////////////
      const IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      { 
        // neighbor volume  
        double wspeedS = 0.0;

        const IntersectionType& inter = *nit;
  
        // if intersection with boundary 
        // --boundary
        if( inter.boundary() ) 
        { 
          applyLocalBoundary(inter, 
              en, geo, volQuad, numDofs, 
              bsetEn, (LocalMatrixType* ) 0, 
              rhsFunc_, 
              wspeedS ); 
        } // end if boundary
      } // end intersection iterator 

      if( rhs_ ) 
      {
        // local function for right hand side 
        rhs_->localFunction(en) += rhsFunc_ ;
      }
    } // end applyBoundary  

    template<class QuadratureType, class CoeffCallerType> 
    void volumetricPart(ConstEntityType& en, 
                        const Geometry& geo,
                        QuadratureType& volQuad,
                        const CoeffCallerType& coeffCaller,
                        const BaseFunctionSetType bsetEn,
                        const int numDofs, 
                        LocalMatrixType& matrixEn) const
    {
      // number of different dofs 
      const int numDiffDofs = numDofs / dimRange ;

      const int quadNop = volQuad.nop();
      // loop over all quadrature points 
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calc integration element 
        const double intel = volQuad.weight(l)
            *geo.integrationElement( volQuad.point(l) );

        const JacobianInverseType& inv =
          geo.jacobianInverseTransposed( volQuad.point(l) );

        ////////////////////////////////////
        // create rightHandSide
        ////////////////////////////////////
        coeffCaller.rightHandSide(caller_, en, volQuad, l, intel);
        
        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        
        // call coefficient of discrete model 
        coeffCaller.evaluateCoefficient(caller_, en, volQuad, l, coeffEn_ ); 

        ///////////////////////////////////////////
        // apply jacobian inverse and coefficient 
        ///////////////////////////////////////////
        for(int k = 0; k < numDofs; ++k)
        {
          JacobianRangeType& psi = psi_[k]; 

          // eval grad psi on reference element
          bsetEn.jacobian( k, volQuad[l], psitmp_ );
  
          // apply inverse jacobian 
          for(int i=0; i<dimRange; ++i) 
          {
            inv.mv(psitmp_[i], psi[i]);
          }

          // apply coefficient 
          coeffCaller.applyCoefficient(coeffEn_, psi, coeffPsi_[ k ] ); 
        }

        ///////////////////////////////////////////
        // fill element matrix 
        ///////////////////////////////////////////
        for(int dk = 0; dk < numDiffDofs; ++dk)
        {
          // add diagonal entry
          {
            int k = dk * dimRange ;
            // loop over dimRange 
            for (int i = 0; i <dimRange; ++i, ++k ) 
            {
              RangeFieldType val = coeffPsi_[k][i] * psi_[k][i];
              val *= intel;
              matrixEn.add( k , k , val );
            }
          }

          // add secondary diagonal
          // assume matrix is symectric
          // entry (k,j) == entry (j,k)
          //for (int dj = 0; dj < numDiffDofs; ++dj) 
          for (int dj = dk+1; dj < numDiffDofs; ++dj)
          {
            int k = dk * dimRange ;
            int j = dj * dimRange ;

            // loop over dimRange 
            for (int i = 0; i <dimRange; ++i, ++k, ++j ) 
            {
              RangeFieldType val = intel * (coeffPsi_[k][i] * psi_[j][i]);
              
              // add k,j 
              matrixEn.add( k , j , val );

              // add j,k
              matrixEn.add( j , k , val );
            }
          }
        }
      } // end element integral 
    }
 
    // calculate inverse mass matrix for gradient space 
    void calculateInverseMass(ConstEntityType& en,
                              const Geometry & geo,
                              VolumeQuadratureType& volQuad) const 
    {
#if USE_COMPACT_LDG 
      // for compact LDG calculate inverse mass matrix 
      if( applyInverseMass_ && hasLifting() )
      {
        typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType GradBaseFunctionSetType;
        const GradBaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
        // get number of base functions for gradient space 
        const size_t numGradBase = enSet.numBaseFunctions();

        // for en 
        inverseMassEn_ = 0;
        if( rRets_.size() < numGradBase ) 
        {
          rRets_.resize( numGradBase );
        }
        getMassMatrix( geo, volQuad, enSet, numGradBase, rRets_, inverseMassEn_);

        // invert matrix 
        inverseMassEn_.invert();
      }
#endif
    }

    ///////////////////////////////////////////
    //! --apply operator on entity 
    ///////////////////////////////////////////
    void applyLocal(ConstEntityType& en) const
    {
      //std::cout << "Call CDGPrimal :: applyLocal for entity " <<
      //  spc_.indexSet().index( en ) << std::endl; 

      // this method should not be called for ghost entities 
      assert( en.partitionType() != GhostEntity );
      
      // get local element matrix 
      LocalMatrixType matrixEn = matrixObj_.localMatrix(en,en); 

      // make entities known in callers
      caller_.setEntity(en);

      // create volume quadrature  
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      // get geometry
      const Geometry & geo = en.geometry();

      // get base function set of single space 
      const BaseFunctionSetType bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      assert( numDofs > 0 );

      const bool rightHandSide = problem_.hasRHS() && rhs_; 

      rhsFunc_.init( en );
      rhsFunc_.clear();

      // local function for right hand side 
      SingleLFType& singleRhs = rhsFunc_; 
      
      // resize caches 
      resizeCaches(numDofs);

      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////
      if(problem_.hasCoefficient() && rightHandSide )
      {
        CoefficientCaller<DiscreteModelCallerType,true,true> coeffCaller( singleRhs ); 
        volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else if( problem_.hasCoefficient() )
      {
        CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
        volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else 
      if ( rightHandSide )
      {
        CoefficientCaller<DiscreteModelCallerType,false,true> coeffCaller( singleRhs ); 
        volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else 
      {
        CoefficientCaller<DiscreteModelCallerType,false,false> coeffCaller; 
        volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }

      // get elements volume 
      const double enVolume = geo.volume();

      // calculate inverse mass matrix for gradient space 
      calculateInverseMass(en, geo, volQuad);

#ifdef DG_DOUBLE_FEATURE
      //const int enIdx = spc_.indexSet().index( en );
#endif

      /////////////////////////////////
      // Surface integral part
      /////////////////////////////////
      IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      { 
        // neighbor volume  
        double nbVolume = enVolume;
        double wspeedS = 0.0;

        // get intersetion 
        const IntersectionType& inter = *nit ;
        // when space is conitnuous then only apply intersection fluxes 
        // on non-conforming intersections 
        const bool applyNeighborFlux = ( spc_.continuous() ) 
          ? ( ! inter.conforming()) : true ;
  
        // if neighbor exists 
        if( inter.neighbor() && applyNeighborFlux ) 
        {
          EntityPointerType neighEp = inter.outside();
          ConstEntityType& nb = *neighEp;

          // get partition type 
          const bool ghostEntity = 
            //( nb.partitionType() == GhostEntity );
            ( nb.partitionType() != InteriorEntity );

#ifdef DG_DOUBLE_FEATURE
          /*
          const int nbIdx = spc_.indexSet().index( nb );

          // only once per intersection or when outside is not interior 
          if( (enIdx < nbIdx) 
              || ghostEntity
            )
            */
#endif
          {
            // check conformity 
            if( inter.conforming() ) 
            {
              FaceQuadratureType faceQuadInner(gridPart_, inter, faceQuadOrd_,
                                               FaceQuadratureType::INSIDE);

        
              FaceQuadratureType faceQuadOuter(gridPart_, inter, faceQuadOrd_,
                                               FaceQuadratureType::OUTSIDE);

              // apply neighbor part 
              nbVolume = applyLocalNeighbor(inter,en,geo,nb,volQuad,
                              faceQuadInner,faceQuadOuter, 
                              bsetEn,matrixEn, singleRhs, 
                              ! ghostEntity ,
                              wspeedS  
                           );
            }
            else 
            {
              // we only should get here whne a non-conforming situation 
              // occurs in a non-conforming grid 
              assert( GridPartType :: conforming == false );
              
              typedef typename FaceQuadratureType :: NonConformingQuadratureType 
                NonConformingFaceQuadratureType;
              
              NonConformingFaceQuadratureType 
                nonConformingFaceQuadInner(gridPart_, inter, faceQuadOrd_,
                                           NonConformingFaceQuadratureType::INSIDE);
          
              NonConformingFaceQuadratureType 
                nonConformingFaceQuadOuter(gridPart_, inter, faceQuadOrd_,
                                           NonConformingFaceQuadratureType::OUTSIDE);

              // apply neighbor part 
              nbVolume = applyLocalNeighbor(inter, en,geo,nb,volQuad,
                            nonConformingFaceQuadInner,
                            nonConformingFaceQuadOuter, 
                            bsetEn,matrixEn, singleRhs, 
                            ! ghostEntity,
                            wspeedS 
                          );
            }
          }
        } // end if neighbor 

        // if intersection with boundary 
        // --boundary
        if( inter.boundary() ) 
        { 
          applyLocalBoundary(inter, 
              en, geo, volQuad, numDofs, 
              bsetEn, &matrixEn, singleRhs, 
              wspeedS ); 

        } // end if boundary

        // check timestep size 
        if ( wspeedS > minLimit_ )
        {
          const double minvolS = std::min(enVolume, nbVolume);
          const double p1 = 8.0 * (spc_.order() + 1);
          wspeedS *= p1 ;
          dtMin_ = std::min(dtMin_,minvolS/wspeedS);
        }

      } // end intersection iterator 

      if( rhs_ ) 
      {
        // local function for right hand side 
        rhs_->localFunction(en) += rhsFunc_ ;
      }

      // resort corresponding matrix rows for ascending numbering 
      matrixEn.resort(); 

    } // end apply local 

    void multLocal(LocalMatrixType& matrix, SingleLFType& rhsLF, 
                   const SingleLFType& uhLf) const 
    {
      const double factor = 1.0 - theta_ ;
      typedef FieldVector<double, elementMassSize > MassVectorType; 
      MassVectorType rhs (0);

      matrix.multiplyAdd( uhLf, rhs );

      const int numDofs = rhsLF.numDofs();
      // add to system matrix and adjust right hand side 
      for(int i=0; i<numDofs; ++i) 
      {
        rhsLF[i] += factor * rhs[ i ];
      }
    }

    void multiplyVec(const DomainType& unitNormal,
                     const GradRangeType& value,
                     RangeType& ret) const 
    {
      //if( dimRange == 1 ) 
      //  ret = unitNormal * value ;
      //else 
      {
        for(int i=0; i<dimRange; ++i) 
        {
          ret[i] = 0;
          // multiply with unitNormal
          int idx = i * dimDomain; 
          for(int j=0; j<dimDomain; ++j, ++idx) 
            ret[i] += unitNormal[j] * value[ idx ];
        }
      }
    }

    RangeFieldType multiplyVec(const GradRangeType& x,
                               const GradRangeType& y,
                               const int i) const 
    {
      //  returns x * y (only using components != 0)
      RangeFieldType ret = 0;
      int idx = i * dimDomain; 
      for(int j = 0 ; j<dimDomain; ++j, ++idx) 
        ret += x[ idx ] * y[ idx ];
      return ret;
    }

    //! apply boundary integrals to matrix and right hand side 
    template <class LFType>
    void applyLocalBoundary(const IntersectionType& inter, 
                            ConstEntityType& en,
                            const Geometry& geo, 
                            VolumeQuadratureType& volQuad,
                            const int numDofs,
                            const BaseFunctionSetType& bsetEn, 
                            LocalMatrixType* matrixEnPtr, 
                            LFType& singleRhs,
                            double& wspeedS) const 
    {
      {
        // create quadrature 
        FaceQuadratureType faceQuadInner(gridPart_, inter, faceQuadOrd_,
                                         FaceQuadratureType::INSIDE);

        applyLocalBoundary(inter, 
            en, geo, volQuad, faceQuadInner, numDofs, 
            bsetEn, matrixEnPtr, 
            singleRhs, 
            wspeedS ); 
      }
    }

    //! apply boundary integrals to matrix and right hand side 
    template <class FaceQuadType, class LFType>
    void applyLocalBoundary(const IntersectionType& nit, 
                            ConstEntityType& en,
                            const Geometry& geo, 
                            VolumeQuadratureType& volQuad,
                            FaceQuadType& faceQuadInner,
                            const int numDofs,
                            const BaseFunctionSetType& bsetEn, 
                            LocalMatrixType* matrixEnPtr, 
                            LFType& singleRhs,
                            double& wspeedS) const 
    {
      LocalMatrixType& matrixEn = *matrixEnPtr;

#if USE_COMPACT_LDG 
      typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType BaseFunctionSetType;
      const BaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );

      // get number of base functions for gradient space 
      const int numGradBase = enSet.numBaseFunctions();

      bool addCompactDG = false;

      // apply local inverse (multiply with vol/refVol)
      const double massInv = (applyInverseMass_) ? 
                1.0 : localGradientMass_.getAffineMassFactor( geo );
      
      if( hasLifting() ) 
      {
        // resize and reset temporary functions 
        resizeTemporaryFunctions( en, r_e_, numDofs, numGradBase );
        resizeTemporaryFunctions( en, r_e_neigh_ , numDofs, numGradBase );
      } // end compact LDG 
#endif

      const bool affineIntersection = nit.geometry().affine();

      // number of different dofs 
      const int numDiffDofs = numDofs / dimRange ;

      // loop over quadrature points 
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calculate normal 
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
        unitNormal *= 1.0/faceVol;

        assert( std::abs( unitNormal.two_norm() - 1.0 ) < 1e-8 ) ; 
        assert( faceVol > 0 );
        
        const double intelFactor = faceQuadInner.weight(l);
        
        // integration element factor 
        const double intel = intelFactor * faceVol;

        // intel switching between bilinear from B_+ and B_-  
        const double bilinIntel = (bilinearPlus_) ? intel : -intel;

        // get boundary value 
        RangeType boundaryValue(0.0);

        // call boundary value function 
        BoundaryIdentifierType bndType = 
          caller_.boundaryValue(nit, faceQuadInner, l, boundaryValue);

        // only Dirichlet and Neumann Boundary supported right now 
        assert( bndType.isDirichletType() || bndType.isNeumannType() );

        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        assert( psi_.size() > 0 );
        JacobianRangeType& norm = psi_[0];
        if( problem_.hasCoefficient() )
        {
          // evaluate coefficient on boundary
          caller_.evaluateCoefficientBoundary(nit, faceQuadInner, l, coeffEn_);

          for(int i=0; i<dimRange; ++i)
          {
            coeffEn_[i].mv(unitNormal,norm[i]);
          }
        }
        else 
        {
          for(int i=0; i<dimRange; ++i)
          {
            for(int j=0; j<dimDomain; ++j)
            {
              norm[i][j] = unitNormal[j];
            }
          }
        }

        double ldt = 0.0 ;
        // overall beta factor (also needed for diffusion time step)
        // facBeta has atually no scaling with h ( 1/|e| integrated over e = 1)
        const RangeType facBeta = factorBeta(intelFactor, coeffEn_, affineIntersection, ldt );
        const bool hasPenaltyTerm = penaltyNotZero_ && ( facBeta.infinity_norm() > 0 );

        wspeedS += ldt * faceQuadInner.weight(l);

        // evaluate base functions on entity 
        evaluateBaseFunctions(en, geo, bsetEn, numDiffDofs, faceQuadInner, l, norm, tau_, phi_ );  

#if USE_COMPACT_LDG 
        // only add compact LDG values on dicrichlet boundary 
        if( hasLifting() && ! spc_.continuous() && bndType.isDirichletType() )
        {
          addCompactDG = true; 
          GradRangeType& tmp = rRets_[0];

          MassVectorType bndRhs, intRhs, intRes, bndRes;

          // get numbre of base functions 
          for(int m=0; m<numGradBase; ++m)
          {  
            // eval base functions 
            enSet.evaluate(m, faceQuadInner[l], tmp );

            // mulitply tmp with unitnormal 
            multiplyVec(unitNormal, tmp, eta_ [ m ]); 
          }

          const double bndIntel = intel;

          // calculate coefficients for r_D( g_D ) and r_D ( u_h )
          for(int k=0; k<numDofs; ++k)
          {
            RangeType phiVal (phi_[k]);
            phiVal *= bndIntel;
            
            RangeType phiBnd ( boundaryValue );
            phiBnd *= bndIntel;
            
            // calculate coefficients 
            for(int m=0; m<numGradBase; ++m)
            {
              intRhs[m] = -(phiVal * eta_[m]); 
              bndRhs[m] = -(phiBnd * eta_[m]);
            }

            if( applyInverseMass_ ) 
            {
              inverseMassEn_.mv( intRhs , intRes );
              inverseMassEn_.mv( bndRhs , bndRes );
            }
            else 
            {
              // apply to vectors 
              for(int m=0; m<numGradBase; ++m)
              {
                bndRes[m] = bndRhs[m] * massInv;
                intRes[m] = intRhs[m] * massInv;
              }
            }
            
            TemporaryLocalFunctionType& r_e = *(r_e_[k]);
            TemporaryLocalFunctionType& r_e_bnd = *(r_e_neigh_[k]);
            for(int m=0; m<numGradBase; ++m)
            {
              r_e[m]     += intRes[m];
              r_e_bnd[m] += bndRes[m];
            }
          }
        }
#endif
               
        // if not Babuska-Zlamal method, add boundary terms 
        {
          // only change right hand side if exists 
          if( bndType.isDirichletNonZero() && rhs_ )
          {
            // fill right hand side  
            for(int k=0; k<numDofs; ++k)
            {  
              // only valid for dim range = 1
              RangeFieldType rhsVal1 = boundaryValue[ k % dimRange ] * tau_[k];

              rhsVal1 *= bilinIntel;
              singleRhs[k] += rhsVal1;
            }
          }

          // only on non Neumann type boundaries
          if( matrixEnPtr && bndType.isDirichletType() )
          {
            // fill matrix entries 
            for(int dk=0; dk<numDiffDofs; ++dk)
            {  
              for (int dj = 0; dj < numDiffDofs; ++dj) 
              {
                int k = dk * dimRange ;
                int j = dj * dimRange ;

                for(int i=0; i<dimRange; ++i, ++k , ++j ) 
                {
                  {
                    // grad w * v 
                    RangeFieldType val = tau_[j] * phi_[k][i];
                    val *= -intel;
                    matrixEn.add( k , j , val );
                  }
                
                  {
                    // w * grad v
                    RangeFieldType val = tau_[k] * phi_[j][i];
                    val *= bilinIntel;
                    matrixEn.add( k , j , val );
                  }
                }
              }
            }
          }
        }
            
        // dirichlet boundary values for u 
        // only change right hand side if exists 
        if(bndType.isNeumannNonZero() && rhs_ )
        {
          // fill right hand side 
          for(int dk=0; dk<numDiffDofs; ++dk)
          {  
            int k = dk * dimRange ;
            for(int i=0; i<dimRange; ++i, ++k ) 
            {
              // only valid for dim range = 1
              RangeFieldType rhsVal = boundaryValue[i] * phi_[k][i];

              rhsVal *= intelFactor;
              singleRhs[k] += rhsVal;
            }
          }
        }

        if( hasPenaltyTerm )
        {
          // stabilization 
          if( matrixEnPtr && bndType.isDirichletType() )
          {
            // fill matrix entries 
            // fill matrix entries 
            for(int dk=0; dk<numDiffDofs; ++dk)
            {  
              for (int dj = 0; dj < numDiffDofs; ++dj) 
              {
                int k = dk * dimRange ;
                int j = dj * dimRange ;

                for(int i=0; i<dimRange; ++i, ++k , ++j ) 
                {
                  // phi_j * phi_k 
                  RangeFieldType phiVal = phi_[j][i] * phi_[k][i]; 
                  phiVal *= facBeta[i];
                  matrixEn.add( k , j , phiVal );
                }
              }
            }
            
            // dirichlet boundary values for u 
            // only change right hand side if exists 
            if(bndType.isDirichletNonZero() && rhs_ )
            {
              // fill right hand side 
              for (int dj = 0; dj < numDiffDofs; ++dj) 
              {
                int j = dj * dimRange ;

                for(int i=0; i<dimRange; ++i, ++j ) 
                {
                  // phi_j * phi_k 
                  const RangeFieldType phiVal = phi_[j][i] * boundaryValue[ i ] * facBeta[i];
                  singleRhs[ j ] += phiVal;
                }
              }
            } 
          }
        }
      }

#if USE_COMPACT_LDG 
      // now add lifting operators if we use compact LDG 
      if( addCompactDG )
      {
        // only change right hand side if exists 
        if( rhs_ )
        {
          addBoundaryLifting(nit,faceQuadInner,numDofs,
              bsetEn,r_e_neigh_, singleRhs );
        }
        
        // only do this if matrix is assembled 
        if( matrixEnPtr ) 
        {
          CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
          addLiftingOperator(coeffCaller,en,
                             geo,volQuad,
                             numDofs,r_e_,r_e_,matrixEn);
        }
      }
#endif

    } // end applyLocalBoundary

    double compBetaK(const FluxRangeType& K) const 
    {
      DomainType eigenvalues; 
      // calculate eigenValues 
      FMatrixHelp :: eigenValues( K, eigenvalues );

      assert( K.cols == K.rows );
      if( K.rows == 1 ) 
      {
        return eigenvalues[0]; 
      }
      else 
      {
        return SQR(eigenvalues[K.rows-1]) / eigenvalues[0];
      }
    }

    // --factorBeta
    RangeType factorBeta(const double intel, 
                         const CoeffMatrixType& K, 
                         const bool affineIntersection,
                         double& wspeedL ) const
    {
      RangeType betaFactor(0);

      if( theoryParams_ && affineIntersection ) 
        return betaFactor;

      for(int i=0; i<dimRange; ++i) 
      {
        // get sqrt( maxLambda^2 / minLambda ) 
        const double betS = compBetaK( K[i] ); 
        //std::cout << betS << "\n";

        //////////////////////////////////////////////
        //  store local diffusion time step 
        //////////////////////////////////////////////
        wspeedL = betS ; 

        // set betaFactor 
        betaFactor[i] = (penalty_ * intel * betS);
        //betaFactor[i] = (globalBeta_ * betS);
      }
      return betaFactor; 
    }

    // --factorBeta
    RangeType factorBeta(const double intel, 
                         const CoeffMatrixType& enKC, 
                         const CoeffMatrixType& nbKC,
                         const bool affineIntersection,
                         double& wspeedL ) const
    {
      RangeType betaFactor(0);
      CoeffMatrixType KC; 
      for(int i=0; i<dimRange; ++i) 
      {
        const FluxRangeType& enK = enKC[i];
        const FluxRangeType& nbK = nbKC[i];

        FluxRangeType& K = KC[i];
        // use average diffusion matrix 
        for (int k=0; k<dimDomain; ++k) 
        {
          for(int j=0; j<dimDomain; ++j) 
          {
            K[k][j] = 0.5 * (enK[k][j] + nbK[k][j]);
          }
        }

        // get ( maxLambda^2 / minLambda ) 
        double betS = compBetaK( K ); 

        //////////////////////////////////////////////
        //  store local diffusion time step 
        //////////////////////////////////////////////
        wspeedL = betS ; 

        const double betEn = compBetaK( enK ); 
        const double betNb = compBetaK( nbK ); 

        const double jump = std::tanh( std::abs( betEn - betNb ) );

        // only for small values of betS apply betS in smooth cases 
        const double betN = ( theoryParams_ && affineIntersection ) ? 
                            ( 0.0 ) : betS ;
                            //(std :: min( betS, 1.0 ));
        //std::cout << betN << "\n";

        // betS becomes 1 if the eigen values of both matrices are the same 
        betS = betS * jump + (1.0 - jump) * betN;

        // set betaFactor 
        betaFactor[i] = (penalty_ * intel * betS);
        //betaFactor[i] = (globalBeta_ * betS);
      }

      return betaFactor; 
    }

    // resize memory for r_e and l_e functions 
    void resizeTemporaryFunctions(const EntityType& en,
                                  TemporaryLocalFunctionArrayType& r_e,
                                  const size_t numDofs, const size_t numGradBase) const 
    {
      if( r_e.size() < numDofs )
      {
        r_e.resize(0);
        r_e.resize( numDofs );
        rRets_.resize( std::max( numDofs , numGradBase ));
        rRetsCoeff_.resize( numDofs );
        for(size_t i=0; i<numDofs; ++i) 
        {  
          std::auto_ptr< TemporaryLocalFunctionType > 
            ptr(  new TemporaryLocalFunctionType ( gradientSpace_ ) );
          r_e[i] = ptr; 
        }
      }

      if( eta_.size() < numGradBase )
      {
        eta_.resize( numGradBase );
      }

      for(size_t i=0; i<numDofs; ++i) 
      {
        TemporaryLocalFunctionType& re = (*r_e[i]);
        re.init ( en );
        // set all dofs to zero 
        re.clear();
      }
    }


    template <class QuadratureImp> 
    bool determineDirection( const IntersectionType & nit,
                             const QuadratureImp& faceQuad,
                             const Geometry& enGeo,
                             const Geometry& nbGeo ) const 
    {
      if( areaSwitch_ ) 
      {
        return nbGeo.volume() > enGeo.volume();
        //return enGeo.volume() > nbGeo.volume();
      }
      else 
      {
        DomainType normal( nit.integrationOuterNormal( faceQuad.localPoint(0) ) );
        assert( std::abs( normal * upwind_ ) > 0 );
        return (normal * upwind_ < 0 ); 
      }
    }


    template <class QuadratureImp> 
    double applyLocalNeighbor(const IntersectionType & nit, 
                              ConstEntityType & en, 
                              const Geometry& enGeo, 
                              ConstEntityType & nb,
                              VolumeQuadratureType & volQuad,
                              const QuadratureImp & faceQuadInner, 
                              const QuadratureImp & faceQuadOuter, 
                              const BaseFunctionSetType & bsetEn, 
                              LocalMatrixType & matrixEn,
                              SingleLFType& singleRhs,
                              const bool interior, 
                              double& wspeedS 
                             ) const
    {
      const int numDofs = bsetEn.numBaseFunctions();

      // make neighbor known to model caller 
      caller_.setNeighbor(nb);
  
      // get neighbors geometry 
      const Geometry& nbGeo = nb.geometry();

      const bool affineIntersection = nit.geometry().affine() ;
      if( ! affineIntersection ) 
      {
        affineGeoms_ = false ;
      }

      ////////////////////////////////////////////////////////////
      RangeFieldType resultLeft(0.0);
      RangeFieldType resultRight(0.0);

      RangeFieldType phi_j;
      RangeFieldType phiNeigh_j;
      RangeFieldType phiEn;
      RangeFieldType phiNeigh;

      // reuse cache mem 
      assert( psi_.size() > 0 );
      assert( coeffPsi_.size() > 0 );
      JacobianRangeType& normEn = psi_[0];
      JacobianRangeType& normNb = coeffPsi_[0];

      // create matrix handles for neighbor 
      LocalMatrixType matrixNb = matrixObj_.localMatrix( en, nb );

#ifndef DG_DOUBLE_FEATURE
      if( ! interior ) 
      {
        // create matrix handles for neighbor 
        LocalMatrixType nbMatrix = matrixObj_.localMatrix( nb, nb ); 
      
        // set matrix to id matrix for ghost cells (otherwise 
        // one has problems with preconditioning)
        for(int k=0; k<numDofs; ++k) 
        {
          nbMatrix.set( k, k, 1);
        }
      }
#endif

      // calculate outer normal 
      const bool useInterior = determineDirection( nit, faceQuadInner, enGeo, nbGeo );

#ifdef DG_DOUBLE_FEATURE
      // create matrix handles for neighbor (when called with ghost do nothing)
      LocalMatrixType enMatrix = matrixObj_.localMatrix( nb, (interior) ? en : nb ); 

      // create matrix handles for neighbor 
      LocalMatrixType nbMatrix = matrixObj_.localMatrix( nb, nb ); 
      
      // set matrix to id matrix 
      if( ! interior ) 
      {
        for(int k=0; k<numDofs; ++k) 
        {
          nbMatrix.set( k, k, 1);
        }
      }

      // only calculate on when useInterior is true 
      if( ! useInterior ) return nbGeo.volume();

#elif USE_COMPACT_LDG 
      const double massInvNb = (applyInverseMass_) ? 
        1.0 : localGradientMass_.getAffineMassFactor( nbGeo );
#endif 

      // get base function set 
      const BaseFunctionSetType bsetNeigh = spc_.baseFunctionSet(nb);

#if USE_COMPACT_LDG 

      // C_12 switch 
      const RangeFieldType C_12 = ( useInterior ) ? -0.5 : 0.5;

      const double massInvEn = (applyInverseMass_) ? 
        1.0 : localGradientMass_.getAffineMassFactor( enGeo );

      typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType BaseFunctionSetType;
      const BaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
      const BaseFunctionSetType nbSet = gradientSpace_.baseFunctionSet( nb );

      // get number of base functions for gradient space 
      const int numGradBase = enSet.numBaseFunctions();

      // if we use compact LDG initialize helper functions 
      if( hasLifting() ) 
      {
        // resize and reset temporary functions 
        resizeTemporaryFunctions( en, r_e_ , numDofs, numGradBase );
        resizeTemporaryFunctions( en, r_e_neigh_ , numDofs, numGradBase );

        resizeTemporaryFunctions( nb, r_e_nb_ , numDofs, numGradBase );
        resizeTemporaryFunctions( nb, r_e_neigh_nb_ , numDofs, numGradBase );
        
        // only need this matrix when useInterior is false 
        if( ! useInterior && applyInverseMass_ )
        {
          // for nb 
          inverseMassNb_ = 0;
          VolumeQuadratureType nbQuad(nb, volumeQuadOrd_);
          getMassMatrix( nbGeo, nbQuad,nbSet,numGradBase,rRets_,inverseMassNb_);

          // invert matrix 
          inverseMassNb_.invert();
        }
      }
#endif

      // number of different dofs 
      const int numDiffDofs = numDofs / dimRange ;

      // loop over all quadrature points 
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // cacluate outer normal 
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
        unitNormal *= 1.0/faceVol; 

        // make sure switch is the same for hole face 
        assert( ( method_ == method_cdg ) ? (unitNormal * upwind_  < 0) == useInterior : true ); 

        // make sure we have the same factors 
        assert( std::abs(faceQuadInner.weight(l) - faceQuadOuter.weight(l)) < 1e-10);
        // integration element factor 
        const double intelFactor = faceQuadInner.weight(l);
        
        const double intel = faceVol * intelFactor; 
#ifdef DG_DOUBLE_FEATURE
        // use opposite signs here
        const double outerIntel = -faceVol * faceQuadOuter.weight(l); 
        // intel switching between bilinear from B_+ and B_-  
        const double outerBilinIntel = (bilinearPlus_) ? outerIntel : -outerIntel;
#endif
        // intel switching between bilinear from B_+ and B_-  
        const double bilinIntel = (bilinearPlus_) ? intel : -intel;

        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        if(problem_.hasCoefficient())
        {
          // call anayltical flux of discrete model 
          caller_.evaluateCoefficientFace(nit,
                  faceQuadInner, faceQuadOuter, l, coeffEn_, coeffNb_);

          for(int i=0; i<dimRange; ++i)
          {
            coeffEn_[i].mv(unitNormal, normEn[i]);
            coeffNb_[i].mv(unitNormal, normNb[i]);
          }
        }
        else 
        {
          for(int i=0; i<dimRange; ++i)
          {
            // set to unit matrix 
            setToUnitMatrix(coeffEn_[i], coeffNb_[i]);

            for(int j=0; j<dimDomain; ++j)
            {
              normEn[i][j] = unitNormal[j];
              normNb[i][j] = unitNormal[j];
            }
          }
        }

        double ldt = 0.0;

        // overall beta factor (also needed for diffusion time step)
        const RangeType facBeta = factorBeta(intelFactor , coeffEn_, coeffNb_, affineIntersection, ldt );
        const bool hasPenaltyTerm = penaltyNotZero_ && ( facBeta.infinity_norm() > 0 );
        //if( hasPenaltyTerm ) 
        //  std::cout << facBeta << " factor ets on interior intersection " << std::endl;

        wspeedS += ldt * faceQuadInner.weight(l) ;

        // evaluate base functions on entity    
        evaluateBaseFunctions(en, enGeo, bsetEn, numDiffDofs, faceQuadInner, l, normEn, tau_, phi_ );  

        // evaluate base functions on neighbor 
        evaluateBaseFunctions(nb, nbGeo, bsetNeigh, numDiffDofs, faceQuadOuter, l, normNb, tauNeigh_, phiNeigh_ );  
           
        // loop over all different dofs 
        for(int dk=0; dk<numDiffDofs; ++dk)
        {
          for(int dj=0; dj<numDiffDofs; ++dj)
          {
            int k = dk * dimRange; 
            int j = dj * dimRange; 

            // loop over dimRange 
            for(int i=0; i<dimRange; ++i, ++k, ++j ) 
            {
              // view from inner entity en 
              // v^+ * (grad w^+  + grad w^-)
              numericalFlux2(phi_[k] , tau_[j] , tauNeigh_[j] , resultLeft, resultRight);

              resultLeft *= -intel;
              matrixEn.add( k , j , resultLeft );

              resultRight *= -intel;
              matrixNb.add( k , j , resultRight );

              // view from inner entity en 
              // grad v^+ * ( w^+  - w^-)
              numericalFlux(tau_[k] , phi_[j] , phiNeigh_[j] , resultLeft, resultRight);

              resultLeft *= bilinIntel;
              matrixEn.add( k , j , resultLeft );

              resultRight *= bilinIntel;
              matrixNb.add( k , j , resultRight );
#ifdef DG_DOUBLE_FEATURE 
              // this part should only be calculated if neighboring
              // entity has partition type interior 
              if( interior ) 
              {
                // view from outer entity nb 
                // v^+ * (grad w^+  + grad w^-)
                numericalFlux2(phiNeigh_[k] , tauNeigh_[j] , tau_[j] , resultLeft, resultRight);

                resultLeft *= -outerIntel;
                nbMatrix.add( k , j , resultLeft );

                resultRight *= -outerIntel;
                enMatrix.add( k , j , resultRight );

                // view from outer entity nb 
                // v^+ * (grad w^+  + grad w^-)
                numericalFlux(tauNeigh_[k] , phiNeigh_[j] , phi_[j] , resultLeft, resultRight);

                resultLeft *= outerBilinIntel;
                nbMatrix.add( k , j , resultLeft );

                resultRight *= outerBilinIntel;
                enMatrix.add( k , j , resultRight );
              } // end interior 
#endif // end DG_DOUBLE_FEATURE 
            } // end for i 
          } // end for dj 
        } // end for dk 

        ////////////////////////////////
        // C_11 (beta) Stabilization 
        ////////////////////////////////
        if( hasPenaltyTerm )
        {
          for(int dk=0; dk<numDiffDofs; ++dk)
          {
            for(int dj=0; dj<numDiffDofs; ++dj)
            {
              int k = dk * dimRange; 
              int j = dj * dimRange; 

              // loop over dimRange 
              for(int i=0; i<dimRange; ++i, ++k, ++j ) 
              {
                // phi_j * phi_k on entity 
                phi_j    = phi_[j][i] * phi_[k][i]; 
                // product with nb 
                phiNeigh = phiNeigh_[j][i] * phi_[k][i]; //bsetNeigh.evaluateSingle(j,faceQuadOuter,l, phi_[k] );      
                
                // phi_j * phi_k on neighbour  
                phiNeigh_j = phiNeigh_[j][i] * phiNeigh_[k][i];//bsetNeigh.evaluateSingle(j,faceQuadOuter,l, phiNeigh_[k] );      
                // product with nb 
                phiEn = phi_[j][i] * phiNeigh_[k][i]; //bsetEn.evaluateSingle(j,faceQuadInner,l, phiNeigh_[k]); 

                // view from inner entity en 
                numericalFluxStab(phi_j, phiNeigh , resultLeft, resultRight);

                RangeFieldType valLeft = facBeta[i];
                valLeft  *= resultLeft;

                matrixEn.add( k , j , valLeft );

                RangeFieldType valRight = facBeta[i];
                valRight *= resultRight;

                matrixNb.add( k , j , valRight );
                
#ifdef DG_DOUBLE_FEATURE
                if( interior )
                {
                  // view from outer entity nb
                  numericalFluxStab(phiNeigh_j, phiEn , resultLeft, resultRight);

                  RangeFieldType valLeft = facBeta[i];
                  valLeft  *= resultLeft;

                  nbMatrix.add( k , j , valLeft );

                  RangeFieldType valRight = facBeta[i];
                  valRight *= resultRight;

                  enMatrix.add( k , j , valRight );
                } // end interior
#endif
              } // end for i
            } // end for dj 
          } // end for dk 
        } // end betaNotZero 

#if USE_COMPACT_LDG 
        // lifting operators 
        if( hasLifting() )
        {
          assert( rRets_.size() > 0 );
          GradRangeType& tmp = rRets_[0];

          // evaluate eta for this entity 
          if( useInterior ) 
          {
            for(int m=0; m<numGradBase; ++m)
            {  
              // eval base functions 
              enSet.evaluate(m, faceQuadInner[l], tmp );

              // apply unit normal 
              multiplyVec(unitNormal, tmp, eta_ [ m ]); 
              //eta_[m] = tmp * unitNormal;
            }
          }
          else 
          {
            // evaluate eta on neighbor 
#ifdef DG_DOUBLE_FEATURE 
            assert( false );
            // we should not get here
            abort();
#endif
            DomainType minusUnitNormal( unitNormal );
            minusUnitNormal *= -1.0;

            for(int m=0; m<numGradBase; ++m)
            {  
              // neighbor stuff 
              nbSet.evaluate(m, faceQuadOuter[l], tmp ); 

              // apply -1 * unit Normal
              multiplyVec(minusUnitNormal, tmp, eta_ [ m ]); 
              // eta_[m] = -(tmp * unitNormal);
            }
          }

          for(int k=0; k<numDofs; ++k)
          {
            // calculate [u] 
            // which is the jump of phi 
            RangeType phiDiff (phi_[k]);
            // scale with integration element 
            phiDiff *= intel;

            // which is the jump of phi 
            RangeType phiDiffNb (phiNeigh_[k]);
            // scale with integration element 
            phiDiffNb *= intel;

            if( useInterior )
            {
              if( applyInverseMass_ ) 
              {
                // calculate coefficients 
                for(int m=0; m<numGradBase; ++m)
                {
                  rhsEn_[m] = - (phiDiff  * eta_[m]);
                  rhsNb_[m] =   phiDiffNb * eta_[m];
                  // this should be rhsEn[m] = -0.5 * value - C_12 * value;
                  // but 0.5 + C_12 = 1 
                  // this should be 
                  // rhsNb[m] = 0.5 * value + C_12 * value; 
                  // but 0.5 + C_12 = 1 
                }

                // multiply with inverse 
                inverseMassEn_.mv( rhsEn_, resEn_ );

                // multiply with inverse
                inverseMassEn_.mv( rhsNb_, resNb_ );

                TemporaryLocalFunctionType& r_e       = *(r_e_)[k];
                TemporaryLocalFunctionType& r_e_neigh = *(r_e_neigh_)[k];

                // calculate coefficients 
                for(int m=0; m<numGradBase; ++m)
                {
                  r_e[m]       += resEn_[m];  
                  r_e_neigh[m] += resNb_[m];  
                }
              }
              else 
              {
                phiDiff   *= massInvEn ;
                phiDiffNb *= massInvEn ;

                TemporaryLocalFunctionType& r_e       = *(r_e_)[k];
                TemporaryLocalFunctionType& r_e_neigh = *(r_e_neigh_)[k];

                // calculate coefficients 
                for(int m=0; m<numGradBase; ++m)
                {
                  r_e[m]       -= ( phiDiff  * eta_[m] );
                  r_e_neigh[m] += ( phiDiffNb * eta_[m] );
                }
              }
#ifndef DG_DOUBLE_FEATURE
            }
            else
            {
              RangeType& phiDiffEn    = phiDiff; 
              RangeType& phiDiffNeigh = phiDiffNb; 
              const MassMatrixType& inverseMassNeigh = inverseMassNb_;
              if( ! applyInverseMass_ ) 
              {
                phiDiffEn    *= massInvNb ;
                phiDiffNeigh *= massInvNb ;
              }
#else 
              // switch values for double feature 
              RangeType& phiDiffEn    = phiDiffNb; 
              RangeType& phiDiffNeigh = phiDiff; 
              const MassMatrixType& inverseMassNeigh = inverseMassEn_;
#endif
              
              if( applyInverseMass_ ) 
              {
                // calculate coefficients 
                for(int m=0; m<numGradBase; ++m)
                {
                  //const RangeFieldType value = phiDiffEn * (eta_[m]);
                  // this should be 
                  // rhsEn[m] = 0.5 * value - C_12_nb * value; 
                  // but the 0.5 - (-C_12) = 1 
                  rhsEn_[m] = phiDiffEn * eta_[m]; 
                  //const RangeFieldType value = phiDiffNeigh * (eta_[m]);
                  // this should be 
                  // rhsNb[m] = -0.5 * value + C_12_nb * value; 
                  // but -0.5 + (-C_12) = -1 
                  rhsNb_[m] = -(phiDiffNeigh * eta_[m]) ; 
                }

                // multiply with inverse (here use mv to initialize res)
                inverseMassNeigh.mv( rhsEn_, resEn_ );

                // add nb part (here use umv to add to res)
                inverseMassNeigh.mv( rhsNb_, resNb_ );

                TemporaryLocalFunctionType& r_e_nb       = *(r_e_nb_)[k];
                TemporaryLocalFunctionType& r_e_neigh_nb = *(r_e_neigh_nb_)[k];
                // calculate coefficients 
                for(int m=0; m<numGradBase; ++m)
                {
                  r_e_nb[m]       += resEn_[m];  
                  r_e_neigh_nb[m] += resNb_[m];  
                }
              }
              else 
              {
                TemporaryLocalFunctionType& r_e_nb       = *(r_e_nb_)[k];
                TemporaryLocalFunctionType& r_e_neigh_nb = *(r_e_neigh_nb_)[k];
                // calculate coefficients 
                for(int m=0; m<numGradBase; ++m)
                {
                  r_e_nb[m]       += (phiDiffEn * eta_[m]);
                  r_e_neigh_nb[m] -= (phiDiffNeigh * eta_[m]);  
                }
              }
            }
          }
        }

        ////////////////////////////////////
        //  C12 stabilization 
        ///////////////////////////////////
        if( method_ == method_cdg )
        {
          const RangeFieldType C_12_bilin = C_12 * bilinIntel;
          // loop over all different dofs 
          for(int dk=0; dk<numDiffDofs; ++dk)
          {
            for(int dj=0; dj<numDiffDofs; ++dj)
            {
              int k = dk * dimRange; 
              int j = dj * dimRange; 

              // loop over dimRange 
              for(int i=0; i<dimRange; ++i, ++k, ++j ) 
              {
                // view from inner entity en 
                numericalFlux_C12(phi_[k], tau_[j] , tauNeigh_[j] , resultLeft, resultRight);

                resultLeft *= C_12_bilin; 
                matrixEn.add( k , j , resultLeft );

                resultRight *= C_12_bilin;
                matrixNb.add( k , j , resultRight );

                // view from inner entity en 
                numericalFlux2_C12(tau_[k] , phi_[j] , phiNeigh_[j] , resultLeft, resultRight);

                resultLeft *= C_12_bilin;
                matrixEn.add( k , j , resultLeft );

                resultRight *= C_12_bilin;
                matrixNb.add( k , j , resultRight );

#ifdef DG_DOUBLE_FEATURE 
                // this part should only be calculated if neighboring
                // entity has partition type interior 
                if( interior ) 
                {
                  // view from outer entity nb 
                  numericalFlux_C12(phiNeigh_[k], tauNeigh_[j] , tau_[j] , resultLeft, resultRight);

                  resultLeft *= C_12_bilin; 
                  nbMatrix.add( k , j , resultLeft );

                  resultRight *= C_12_bilin;
                  enMatrix.add( k , j , resultRight  );

                  // view from outer entity nb 
                  numericalFlux2_C12(tauNeigh_[k] , phiNeigh_[j] , phi_[j] , resultLeft, resultRight);

                  resultLeft *= C_12_bilin; 
                  nbMatrix.add( k , j , resultLeft );

                  resultRight *= C_12_bilin;
                  enMatrix.add( k , j , resultRight );
                } // end interior 
#endif
              } // end for i 
            } // end for dj 
          } // end for dk 
        } // compact LDG 
#endif // end USE_COMPACT_LDG 
      } // end loop quadrature points 

#if USE_COMPACT_LDG 
      if( hasLifting() )
      {
        CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
        if( useInterior )
        {
          addLiftingOperator(coeffCaller,en,
                             enGeo,volQuad,
                             numDofs,r_e_,r_e_,matrixEn,
                             r_e_,r_e_neigh_,matrixNb, true );
          
#ifndef DG_DOUBLE_FEATURE
        } 
        else 
        {
          VolumeQuadratureType nbQuad(nb, volumeQuadOrd_);
          ConstEntityType& neigh = nb;
          LocalMatrixType& matrixNeigh   = matrixEn;
          LocalMatrixType& matrixNeighNb = matrixNb;
          const Geometry& bngeo = nbGeo; 
          const bool evalAgain = true ;
#else
          VolumeQuadratureType& nbQuad = volQuad; 
          ConstEntityType& neigh = en ;
          LocalMatrixType& matrixNeigh   = nbMatrix;
          LocalMatrixType& matrixNeighNb = enMatrix;
          const Geometry& bngeo = enGeo; 
          const bool evalAgain = false ;
#endif
          addLiftingOperator(coeffCaller,neigh,
                             bngeo, nbQuad,
                             numDofs, r_e_nb_,r_e_nb_, matrixNeigh,
                             r_e_nb_,r_e_neigh_nb_, matrixNeighNb,
                             evalAgain );
          
        }
      } // end compactLDG 
#endif

      return nbGeo.volume();
    } // end applyLocalNeighbor 

    template <class QuadratureType>  
    void evaluateBaseFunctions(const EntityType& en,
                               const Geometry& geo,
                               const BaseFunctionSetType& bset,
                               const int numDiffDofs, 
                               const QuadratureType& quad,
                               const int l,
                               const JacobianRangeType& factor,
                               MutableArray<RangeFieldType>& tau,
                               MutableArray<RangeType>& phi) const 
    {
      const JacobianInverseType& inv =
        geo.jacobianInverseTransposed( quad.point(l) );

      JacobianRangeType grad, invGrad;
      for(int dk=0; dk<numDiffDofs; ++dk) 
      {
        int k = dk * dimRange ;

        // eval base functions 
        bset.evaluate( k, quad[l], phi[k] );

        // set other base functions 
        for(int i=1; i<dimRange; ++i) 
        {
          const int j = k + i;
          phi[ j ]    = 0;
          phi[ j ][i] = phi[k][0];
        }

        // evaluate gradient 
        bset.jacobian( k, quad[l], grad );

        // apply jacobian inverse  
        inv.mv( grad[0], invGrad[0] );

        for(int i=0; i<dimRange; ++i , ++k )
        {
          tau[k] = invGrad[ 0 ] * factor[ i ];
        }
      }
    }

    template <class BaseFunctionSet, 
              class LocalStorageType, 
              class MassMatrixType>
    void getMassMatrix(const Geometry& geo,
                       VolumeQuadratureType& volQuad,
                       BaseFunctionSet& set,
                       const int numBase,
                       LocalStorageType& tmp,
                       MassMatrixType& massMatrix) const
    {
      const int volNop = volQuad.nop();
      for(int qp=0; qp<volNop; ++qp) 
      {
        // calculate integration weight 
        const double intel = volQuad.weight(qp)
           * geo.integrationElement(volQuad.point(qp));

        for(int m=0; m<numBase; ++m)
        {  
          // eval base functions 
          set.evaluate(m, volQuad[qp], tmp[m] );
        }

        for(int m=0; m<numBase; ++m)
        {
          {
            double val = intel * (tmp[m] * tmp[m]);
            massMatrix[m][m] += val;
          }
          
          for(int k=m+1; k<numBase; ++k) 
          {
            double val = intel * (tmp[m] * tmp[k]);
            massMatrix[k][m] += val;
            massMatrix[m][k] += val;
          }
        }
      }
    }

#if USE_COMPACT_LDG 
    template <class QuadratureType>
    void addBoundaryLifting(const IntersectionType& nit, 
                            QuadratureType& faceQuadInner,
                            const int numDofs, 
                            const BaseFunctionSetType& bsetEn,
                            TemporaryLocalFunctionArrayType& r_e_array,
                            SingleLFType& singleRhs) const
    {
      assert( nit.boundary() );
      
      const int numDiffDofs = numDofs / dimRange ; 
      // loop over quadrature points 
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calculate normal 
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
        unitNormal *= 1.0/faceVol;

        const double intelFactor = faceQuadInner.weight(l);
        
        // integration element factor 
        const double intel = intelFactor * faceVol * liftFactor_ ;

        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        assert( psi_.size() > 0 );
        JacobianRangeType& coefficient = psi_[0];
        {
          // evaluate coefficient on boundary
          caller_.evaluateCoefficientBoundary(nit, faceQuadInner, l, coeffEn_ );

          // apply multiply unitNormal with coefficient 
          for(int i=0; i<dimRange; ++i)
          {
            coeffEn_[i].mv(unitNormal, coefficient[i]);
          }
        }

        GradRangeType& tmp = rRets_[0];

        for(int dk=0; dk<numDiffDofs; ++dk) 
        {
          int k = dk * dimRange;  
          for( int i=0; i<dimRange; ++i, ++k) 
          {
            // eval r_D 
            TemporaryLocalFunctionType& r_e_bnd = *(r_e_array[k]);
            r_e_bnd.evaluate(faceQuadInner[l], tmp);
        
            // eval phi 
            bsetEn.evaluate(k, faceQuadInner[l] , phi_[k]);

            int idx = i * dimDomain;
            RangeFieldType val = 0; 
            for( int j=0; j<dimDomain; ++j, ++idx) 
              val += tmp[ idx ] * coefficient[i][j];
            //RangeFieldType val = tmp * norm[0]; 

            val *= intel * phi_[k][i];
            singleRhs[k] -= val;
          }
        }
      }
    }

    // --addLiftingOperator 
    template <class CoeffCallerType>
    void evaluateCoefficients(CoeffCallerType& coeffCaller,
                              ConstEntityType& en, 
                              const Geometry& geo,
                              VolumeQuadratureType& volQuad) const 
    {
      const size_t volNop = volQuad.nop();

      if( coeffs_.size() != volNop ) 
      {
        coeffs_.resize( volNop );
        weights_.resize( volNop );
      }

      assert( coeffs_.size() == volNop );
      assert( weights_.size() == volNop );

      for (size_t l = 0; l < volNop ; ++l) 
      {
        // calculate integration weight 
        weights_[l] = volQuad.weight(l)
          * geo.integrationElement(volQuad.point(l));

        // evaluate diffusion coefficient 
        coeffCaller.evaluateCoefficient(caller_, en, volQuad, l, coeffs_[ l ] );
      }
    }

    // --addLiftingOperator 
    template <class CoeffCallerType, class BaseFuncArrayType>
    void addLiftingToMatrix(CoeffCallerType& coeffCaller,
                            const int l,
                            const int numDofs, 
                            const int numGradBase,
                            const BaseFuncArrayType& baseFct, 
                            TemporaryLocalFunctionArrayType& r_e_en,
                            TemporaryLocalFunctionArrayType& r_e_nb,
                            LocalMatrixType& matrix) const
    {
      const int numDiffDofs = numDofs / dimRange ;

      GradRangeType cof ;
      for(int k=0; k<numDofs; ++k) 
      {
        cof = 0;
        TemporaryLocalFunctionType& r_en = *(r_e_en[k]);
        TemporaryLocalFunctionType& r_nb = *(r_e_nb[k]);

        GradRangeType& ret = rRets_[k];
        ret = 0;

        // could be implemented more efficiently for dimRange>1 
        for( int m=0; m<numGradBase; ++m) 
        {
          cof.axpy( r_en[m] , baseFct[m] );
          ret.axpy( r_nb[m] , baseFct[m] );
        }

        // apply diffusion coefficient 
        coeffCaller.applyCoefficient(coeffs_[ l ], cof, rRetsCoeff_[k] );
      }

      // calculate integration weight 
      const double intel = weights_[l] * liftFactor_ ; 

      for(int dk=0; dk<numDiffDofs; ++dk)
      {
        for(int dj=0; dj<numDiffDofs; ++dj)
        {
          int k = dk * dimRange; 
          int j = dj * dimRange; 

          // loop over dimRange 
          for(int i=0; i<dimRange; ++i, ++k, ++j ) 
          {
            RangeFieldType val = intel * multiplyVec( rRetsCoeff_[k], rRets_[j], i);
            matrix.add(k, j, val);
          }
        }
      }
    }

    // --addLiftingOperator 
    template <class CoeffCallerType>
    void addLiftingOperator(CoeffCallerType& coeffCaller,
                            ConstEntityType& en, 
                            const Geometry& geo,
                            VolumeQuadratureType& volQuad,
                            const int numDofs, 
                            TemporaryLocalFunctionArrayType& r_e_en,
                            TemporaryLocalFunctionArrayType& r_e_nb,
                            LocalMatrixType& matrixEn,
                            TemporaryLocalFunctionArrayType& r_e_en_nb,
                            TemporaryLocalFunctionArrayType& r_e_nb_nb,
                            LocalMatrixType& matrixNb,
                            const bool evaluateCoeff) const 
    {
      const int volNop = volQuad.nop();

      if( evaluateCoeff ) 
      {
        // evaluate coefficients for all qaud points 
        evaluateCoefficients(coeffCaller,en,geo,volQuad);
      }

      typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType GradBaseFunctionSetType;
      const GradBaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
      const size_t numGradBase = enSet.numBaseFunctions();

      if( baseFct_.size() != numGradBase )
        baseFct_.resize( numGradBase );

      for (int l = 0; l < volNop ; ++l) 
      {
        // evaluate base functions 
        for(size_t m=0; m<numGradBase; ++m) 
        {
          enSet.evaluate( m, volQuad[ l ], baseFct_[ m ] );
        }

        addLiftingToMatrix(coeffCaller, l, numDofs, numGradBase, 
                           baseFct_, r_e_en, r_e_nb, matrixEn);

        addLiftingToMatrix(coeffCaller, l, numDofs, numGradBase, 
                           baseFct_, r_e_en_nb, r_e_nb_nb, matrixNb);
      }
    }

    // --addLiftingOperator 
    template <class CoeffCallerType>
    void addLiftingOperator(CoeffCallerType& coeffCaller,
                            ConstEntityType& en, 
                            const Geometry& geo,
                            VolumeQuadratureType& volQuad,
                            const int numDofs, 
                            TemporaryLocalFunctionArrayType& r_e_en,
                            TemporaryLocalFunctionArrayType& r_e_nb,
                            LocalMatrixType& matrixEn) const
    {
      const int volNop = volQuad.nop();

      // evaluate coefficients for all qaud points 
      evaluateCoefficients(coeffCaller,en,geo,volQuad);

      typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType GradBaseFunctionSetType;
      const GradBaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
      const size_t numGradBase = enSet.numBaseFunctions();

      if( baseFct_.size() != numGradBase )
        baseFct_.resize( numGradBase );

      for (int l = 0; l < volNop ; ++l) 
      {
        // evaluate base functions 
        for( size_t m=0; m<numGradBase; ++m) 
        {
          enSet.evaluate( m, volQuad[ l ], baseFct_[ m ] );
        }

        addLiftingToMatrix(coeffCaller, l, numDofs, numGradBase, 
                           baseFct_, r_e_en, r_e_nb, matrixEn);
      }
    }
#endif

  private:  
    void numericalFlux(const RangeFieldType & grad, 
                       const RangeType & phiLeft,
                       const RangeType & phiRight, 
                       RangeFieldType & resultLeft,
                       RangeFieldType & resultRight) const
    {
      resultLeft  = phiLeft[0];
      resultRight = phiRight[0];

      for(int i=1; i<dimRange; ++i) 
      {
        resultLeft  += phiLeft[i];
        resultRight += phiRight[i]; 
      }

      resultLeft  *=  0.5 * grad;
      // need negative value of phiRight 
      resultRight *= -0.5 * grad; 
    }
                       
    void numericalFlux2(const RangeType & phi,
                        const RangeFieldType & gradLeft, 
                        const RangeFieldType & gradRight,
                        RangeFieldType & resultLeft,
                        RangeFieldType & resultRight) const
    {
      resultLeft  = phi[0];
      resultRight = phi[0];

      for(int i=1; i<dimRange; ++i) 
      {
        resultLeft  += phi[i];
        resultRight += phi[i]; 
      }

      resultLeft  *= 0.5 * gradLeft; 
      resultRight *= 0.5 * gradRight;
    }
                       
    void numericalFluxStab(const RangeFieldType & phiLeft,
                           const RangeFieldType & phiRight,
                           RangeFieldType & resultLeft,
                           RangeFieldType & resultRight) const
    {
      resultLeft  =  phiLeft; 
      resultRight = -phiRight; 
    }

    void numericalFlux_C12(const RangeType & phi,
                        const RangeFieldType & gradLeft,
                        const RangeFieldType & gradRight,
                        RangeFieldType & resultLeft,
                        RangeFieldType & resultRight) const
    {
      resultLeft  = phi[0];
      resultRight = phi[0];

      for(int i=1; i<dimRange; ++i) 
      {
        resultLeft  += phi[i];
        resultRight += phi[i]; 
      }

      resultLeft  *=  gradLeft;
      resultRight *= -gradRight;
    }

    void numericalFlux2_C12(const RangeFieldType & grad,
                            const RangeType & phiLeft,
                            const RangeType & phiRight,
                            RangeFieldType & resultLeft,
                            RangeFieldType & resultRight) const
    {
      resultLeft  = phiLeft[0];
      resultRight = phiRight[0];

      for(int i=1; i<dimRange; ++i) 
      {
        resultLeft  += phiLeft[i];
        resultRight += phiRight[i]; 
      }

      resultLeft  *=  grad;
      resultRight *= -grad;
    }
    
    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
    //! empty constructor not defined 
    CDGPrimalOperatorImpl();
    //! copy constructor not defined 
    CDGPrimalOperatorImpl(const CDGPrimalOperatorImpl&);

  protected:
    mutable DiscreteModelCallerType caller_;

    DiscreteModelType& problem_; 
         
    
    mutable ArgumentType* arg_;
    mutable DestinationType* rhs_;
    mutable const DestinationType* uh_;

    const DiscreteFunctionSpaceType& spc_;
    GridPartType & gridPart_;
    const DiscreteGradientSpaceType gradientSpace_;

    const GridWidthType& gridWidth_;
    DomainType upwind_ ;
    
    mutable LocalMassMatrixType localGradientMass_;
    const bool applyInverseMass_;
    mutable SingleLFType rhsFunc_;

    const int volumeQuadOrd_;
    const int faceQuadOrd_;

    mutable MatrixObjectType matrixObj_;

    // return type of analyticalFlux 
    mutable CoeffMatrixType coeffEn_;
    mutable CoeffMatrixType coeffNb_;

    mutable MassMatrixType inverseMassEn_;
    mutable MassMatrixType inverseMassNb_;

    mutable MassVectorType rhsEn_;
    mutable MassVectorType rhsNb_;

    mutable MassVectorType resEn_;
    mutable MassVectorType resNb_;


    // caches for base function evaluation 
    mutable MutableArray<RangeFieldType> tau_;
    mutable MutableArray<RangeFieldType> tauNeigh_;
    mutable MutableArray<RangeType> phi_;
    mutable MutableArray<RangeType> phiNeigh_;
    mutable MutableArray<JacobianRangeType> psi_;
    mutable MutableArray<JacobianRangeType> coeffPsi_;

    mutable MutableArray<RangeType> eta_;

    mutable MutableArray<GradRangeType> rRets_;
    mutable MutableArray<GradRangeType> rRetsCoeff_;

    mutable TemporaryLocalFunctionArrayType r_e_;
    mutable TemporaryLocalFunctionArrayType r_e_nb_;

    mutable TemporaryLocalFunctionArrayType r_e_neigh_;
    mutable TemporaryLocalFunctionArrayType r_e_neigh_nb_;


    mutable MutableArray< CoeffMatrixType > coeffs_; 
    mutable MutableArray< double > weights_ ; 
    mutable MutableArray< GradRangeType > baseFct_; 

    mutable JacobianRangeType psitmp_;

    mutable bool matrixAssembled_;
    double power_;
    mutable double tau_1_;
    mutable double theta_;
    mutable double dtMin_;
    const double minLimit_;
    int sequence_ ;
    // if true B_+ is used otherwise B_-
    mutable size_t numberOfElements_;
    const DGDiffusionFluxIdentifier method_;
    const double penalty_;
    const bool penaltyNotZero_;
    const bool bilinearPlus_;
    const bool areaSwitch_;

    const int theoryParams_ ;
    const double liftFactor_; 
    mutable bool affineGeoms_;

    mutable double globalBeta_;
  };

  template <class DiscreteModelImp, 
            class PreviousPassImp, 
            class MatrixObjectTraits, 
            int passId >
  class CDGPrimalOperator
    : public CDGPrimalOperatorImpl< DiscreteModelImp, 
                                    PreviousPassImp, 
                                    MatrixObjectTraits, 
                                    passId > ,
     public OEMSolver :: PreconditionInterface                               
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef CDGPrimalOperatorImpl<DiscreteModelImp,
            PreviousPassImp,MatrixObjectTraits, passId > BaseType;

    typedef CDGPrimalOperator<DiscreteModelImp,
            PreviousPassImp,MatrixObjectTraits, passId > ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;

    typedef PreviousPassImp PreviousPassType;

    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename BaseType :: MatrixObjectType MatrixObjectType;
    typedef typename BaseType :: MatrixObjectType::MatrixType MatrixType;
    typedef typename BaseType :: MatrixObjectType::PreconditionMatrixType PreconditionMatrixType;
    
  protected:  
    using BaseType :: matrixObj_ ;
  public:
    //- Public methods
    /**  \brief Constructor
     \param problem Actual problem definition (see problem.hh)
     \param[in]  gradPass  types and discrete model for the gradient pass
     \param pass Previous pass
     \param spc Space belonging to the discrete function local to this pass
     \param paramFile parameter file to read necessary parameters, if empty 
             default parameters will be applied 
    
     \note Available methods are (chosen by parameters B_{+,-}, beta, and CDG-BZ)
          - Interior Penalty : B_{+,-}: 0 , beta: > 0 (big) , CDG-BZ: 0 
          - Baumann-Oden     : B_{+,-}: 1 , beta: = 0       , CDG-BZ: 0 (needs polOrd > 1) 
          - NIPG             : B_{+,-}: 1 , beta: > 0       , CDG-BZ: 0
          - Compact LDG (CDG): B_{+,-}: 0 , beta: > 0       , CDG-BZ: 1
    */         
    CDGPrimalOperator(DiscreteModelType& problem, 
                      PreviousPassType& pass, 
                      const DiscreteFunctionSpaceType& spc,
                      const std::string paramFile = "")
      : BaseType(problem, pass, spc, paramFile)
    {
    }
  public:
    virtual ~CDGPrimalOperator() {}

    //! return refernence to system matrix, used by Solvers
    const MatrixObjectType & systemMatrix () const { return matrixObj_; }
    
    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType & preconditionMatrix () const { 
      return matrixObj_.preconditionMatrix(); 
    }

    //! returns true if preconditioning matrix has been build 
    bool hasPreconditionMatrix() const  { 
      return matrixObj_.hasPreconditionMatrix(); 
    }
  };
#undef DG_DOUBLE_FEATURE  
#undef USE_COMPACT_LDG 
} // end namespace Dune
#endif
