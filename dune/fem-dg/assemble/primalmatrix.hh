#ifndef PRIMALMATRIXASSEMBLY_HH
#define PRIMALMATRIXASSEMBLY_HH

#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/compatibility.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

#include <dune/fem-dg/algorithm/sub/steadystate.hh>
#include <dune/fem-dg/algorithm/sub/elliptic.hh>
#include <dune/fem-dg/misc/uniquefunctionname.hh>
#include <dune/fem-dg/operator/fluxes/diffusion/dgprimalfluxes.hh>
#include "assemblertraits.hh"

#include <dune/fem/operator/common/spaceoperatorif.hh>

namespace Dune
{
namespace Fem
{

  ////////////////////////////////////////
  //ASSEMBLER
  template< class DomainDiscreteFunction, class RangeDiscreteFunction >
  class AssemblerInterface
    : public Operator< DomainDiscreteFunction, RangeDiscreteFunction >
  {

    typedef Operator< DomainDiscreteFunction, RangeDiscreteFunction > BaseType;

  public:
    typedef typename BaseType::DomainFunctionType             DomainDiscreteFunctionType;
    typedef typename BaseType::RangeFunctionType              RangeDiscreteFunctionType;

    //f           [affine or nonlinear model]
    virtual void operator()( RangeDiscreteFunctionType& rhs ) const = 0;

    //F(u)        [affine or nonlinear model]
    virtual void operator()( const DomainDiscreteFunctionType& u, RangeDiscreteFunctionType& dest ) const = 0;

    //(F(u), f)   [affine or nonlinear model]
    virtual void operator()( const DomainDiscreteFunctionType& u, RangeDiscreteFunctionType& dest, RangeDiscreteFunctionType& rhs ) const = 0;
  };

  template< class DomainDiscreteFunction, class RangeDiscreteFunction >
  class AssemblerDefault
    : public AssemblerInterface< DomainDiscreteFunction, RangeDiscreteFunction >
  {

    typedef AssemblerInterface< DomainDiscreteFunction, RangeDiscreteFunction > BaseType;

    using BaseType::operator();
  public:
    typedef typename BaseType::DomainDiscreteFunctionType DomainDiscreteFunctionType;
    typedef typename BaseType::RangeDiscreteFunctionType  RangeDiscreteFunctionType;

    //(F(u), f)   [affine or nonlinear model]
    virtual void operator()( const DomainDiscreteFunctionType& u, RangeDiscreteFunctionType& dest, RangeDiscreteFunctionType& rhs ) const
    {
      //compute rhs
      operator()( rhs );
      //compute lhs
      operator()( u, dest );
    }
  };

  template< class JacobianAssemblerImp >
  class JacobianAssemblerToAssembler
    : public AssemblerDefault< typename JacobianAssemblerImp::DomainDiscreteFunction,
                               typename JacobianAssemblerImp::RangeDiscreteFunction >
  {

    typedef AssemblerDefault< typename JacobianAssemblerImp::DomainDiscreteFunction,
                              typename JacobianAssemblerImp::RangeDiscreteFunction > BaseType;

  public:
    typedef JacobianAssemblerImp JacobianAssemblerType;

    typedef typename BaseType::DomainDiscreteFunctionType        DomainDiscreteFunctionType;
    typedef typename BaseType::RangeDiscreteFunctionType         RangeDiscreteFunctionType;
    typedef typename JacobianAssemblerType::JacobianOperatorType JacobianOperatorType;


    JacobianAssemblerToAssembler( const JacobianAssemblerType& jacAssembler, JacobianOperatorType& jac )
      : jacAssembler_( jacAssembler ),
        jac_( jac )
    {}

    //f           [affine or nonlinear model]
    virtual void operator()( RangeDiscreteFunctionType& dest ) const
    {
      //DomainDiscreteFunctionType nullFunc( "nullfunc", op_.space() );
      //nullFunc.clear();

      //jacAssembler_.jacobian( jac_ );
      //jac_.apply( nullFunc, dest );
    }

    //F(u)        [affine or nonlinear model]
    virtual void operator()( const DomainDiscreteFunctionType& u, RangeDiscreteFunctionType& dest ) const
    {
      jacAssembler_.jacobian( jac_ );
      jac_.apply( u, dest );
    }


  private:
    const JacobianAssemblerType& jacAssembler_;
    JacobianOperatorType& jac_;
  };




  ///////////////////////////////////////////////////
  //JACOBIANASSEMBLER
  template< class JacobianOperatorImp >
  class JacobianAssemblerInterface
  {

  public:
    typedef DifferentiableOperator< JacobianOperatorImp > BaseType;

    typedef typename BaseType::DomainFunctionType   DomainDiscreteFunctionType;
    typedef typename BaseType::RangeFunctionType    RangeDiscreteFunctionType;

    typedef typename BaseType::JacobianOperatorType JacobianOperatorType;


    virtual void prepare( JacobianOperatorType& jac ) const = 0;

    //A           [affine model]
    virtual void jacobian( JacobianOperatorType& jac ) const = 0;

    //F_u*        [nonlinear model]
    virtual void jacobian( const DomainDiscreteFunctionType& u, JacobianOperatorType& jac ) const = 0;

    //(A, f)      [affine model]
    virtual void jacobian( JacobianOperatorType& jac, RangeDiscreteFunctionType& affine ) const = 0;

    //(F_u*, f)   [nonlinear model]
    virtual void jacobian( const DomainDiscreteFunctionType& u, JacobianOperatorType& jac, RangeDiscreteFunctionType& affine ) const = 0;

    virtual void finalize( JacobianOperatorType& jac ) const = 0;
  };

  template< class JacobianOperatorImp >
  class JacobianAssemblerDefault
    : public JacobianAssemblerInterface< JacobianOperatorImp >
  {

    typedef JacobianAssemblerInterface< JacobianOperatorImp > BaseType;

    using BaseType::jacobian;
  public:

    typedef typename BaseType::DomainDiscreteFunctionType DomainDiscreteFunctionType;
    typedef typename BaseType::RangeDiscreteFunctionType  RangeDiscreteFunctionType;

    typedef typename BaseType::JacobianOperatorType       JacobianOperatorType;

    virtual void prepare( JacobianOperatorType& jac ) const
    {}

    //F_u*        [nonlinear model]
    virtual void jacobian( const DomainDiscreteFunctionType& u, JacobianOperatorType& jac ) const
    {
      jacobian( jac );
    }

    //(A, f)      [affine model]
    virtual void jacobian( JacobianOperatorType& jac, RangeDiscreteFunctionType& affine ) const
    {
      jacobian( jac );
    }

    //(F_u*, f)   [nonlinear model]
    virtual void jacobian( const DomainDiscreteFunctionType& u, JacobianOperatorType& jac, RangeDiscreteFunctionType& affine ) const
    {
      jacobian( u, jac );
    }

    virtual void finalize( JacobianOperatorType& jac ) const
    {}
  };

  template< class DomainDiscreteFunction, class RangeDiscreteFunction,
            class OperatorImp = AutomaticDifferenceOperator< DomainDiscreteFunction, RangeDiscreteFunction > >
  class JacobianAutoDiffAssembler
    : public JacobianAssemblerDefault< OperatorImp >
  {
    typedef JacobianAssemblerDefault< OperatorImp > BaseType;

    typedef typename BaseType::DomainDiscreteFunctionType DomainDiscreteFunctionType;
    typedef typename BaseType::RangeDiscreteFunctionType  RangeDiscreteFunctionType;
    typedef typename BaseType::OperatorType               OperatorType;

    typedef typename BaseType::JacobianOperatorType       JacobianOperatorType;

  public:
    JacobianAutoDiffAssembler( OperatorType& op )
      : op_ ( op )
    {}

    void prepare( JacobianOperatorType& jac ) const
    {}

    //A           [affine model]
    virtual void jacobian( JacobianOperatorType& jac ) const
    {
      DomainDiscreteFunctionType nullFunc( "nullfunc", op_.space() );
      nullFunc.clear();
      op_.jacobian( nullFunc, jac );
    }

    //F_u*        [nonlinear model]
    virtual void jacobian( const DomainDiscreteFunctionType& u, JacobianOperatorType& jac ) const
    {
      op_.jacobian( u, jac );
    }

    void finalize( JacobianOperatorType& jac ) const
    {}

  private:
    OperatorType& op_;
  };





  //forward declaration
  template <class Traits>
  class DGPrimalMatrixAssembly;

  template< class Traits >
  class DGPrimalMatrixDiscreteModel;

  template< class Traits >
  class PoissonAssembledLinearOperator
    : public Traits::MatrixContainerType, //this is the linear operator with access to system Matrix....
      public HasMass
  {
    typedef typename Traits::MatrixContainerType                                   BaseType;
  public:

    typedef typename Traits::DomainDiscreteFunctionType                            DomainDiscreteFunctionType;
    typedef typename Traits::RangeDiscreteFunctionType                             RangeDiscreteFunctionType;
    typedef typename Traits::DomainDiscreteFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
    typedef typename Traits::RangeDiscreteFunctionType::DiscreteFunctionSpaceType  RangeSpaceType;

    static constexpr bool assembled = true;

    using BaseType::apply;


    PoissonAssembledLinearOperator( const std::string& name,
                                    const DomainSpaceType &domainSpace,
                                    const RangeSpaceType &rangeSpace,
                                    const MatrixParameter& param = SparseRowMatrixParameter() )
      : BaseType( name, domainSpace, rangeSpace ),
        lambda_( 0 )
    {
    }

    template< class Splitting >
    bool split( const Splitting& split )
    {
      return false;
    }

    const double &lambda () const
    {
      return lambda_;
    }

    void setLambda ( double lambda )
    {
      lambda_ = lambda;
    }

    virtual void operator()( const DomainDiscreteFunctionType &arg, RangeDiscreteFunctionType &dest ) const
    {
      apply( arg, dest );
      //dest *= lambda_; //and other manipulations...?!
    }

    const BaseType &systemMatrix() const
    {
      return *this;
    }

    BaseType &systemMatrix()
    {
      return *this;
    }

    void communicate()
    {}

  private:
    double lambda_;
  };


  //----------------------------
  //for SubEvolutionAlgorithm
  //----------------------------
  template< class Traits, class Assembler = DGPrimalMatrixAssembly<Traits>, class JacobianAssembler = DGPrimalMatrixAssembly<Traits> >
  class AssembledPoissonSpaceOperator
    : public SpaceOperatorInterface< typename Traits::DomainDiscreteFunctionType, PoissonAssembledLinearOperator< Traits > >
  {
    typedef AssembledPoissonSpaceOperator< Traits, JacobianAssembler > ThisType;

    typedef SpaceOperatorInterface< typename Traits::DomainDiscreteFunctionType, PoissonAssembledLinearOperator< Traits > > BaseType;

  public:
    typedef Traits                                                         TraitsType;

    typedef typename Traits::GridPartType                                  GridPartType;
    typedef typename Traits::ModelType                                     ModelType;
    typedef typename ModelType::ProblemType                                ProblemType;

    typedef Assembler                                                      AssemblerType;
    typedef JacobianAssembler                                              JacobianAssemblerType;

    typedef PoissonAssembledLinearOperator< Traits>                        JacobianOperatorType;

    //typedef DGHelmholtzJacobianOperator< JacobianOperatorType >          NewJacobianOperatorType;

    typedef typename Traits::DomainDiscreteFunctionType                    DomainDiscreteFunctionType;
    typedef typename Traits::RangeDiscreteFunctionType                     RangeDiscreteFunctionType;

    typedef typename DomainDiscreteFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
    typedef typename RangeDiscreteFunctionType::DiscreteFunctionSpaceType  RangeDiscreteFunctionSpaceType;

    typedef RangeDiscreteFunctionType                                      RangeFunctionType;
    typedef DomainDiscreteFunctionType                                     DomainFunctionType;

    typedef DomainDiscreteFunctionSpaceType                                DiscreteFunctionSpaceType;

    typedef std::tuple<>                                                   ExtraParameterTupleType;



    AssembledPoissonSpaceOperator( GridPartType& gridPart, ProblemType& problem,
                                   ExtraParameterTupleType tuple = ExtraParameterTupleType(),
                                   const std::string name = "" )
      : gridPart_( gridPart ),
        problem_( problem ),
        model_( problem_ ),
        space_( gridPart_ ),
        assembler_( gridPart_, model_ ),
        jacobianAssembler_( gridPart_, model_ )
    {}

    template< class Splitting >
    bool split( const Splitting& split )
    {
      //setAlpha( lhs.getAlpha(), rhs.getAlpha() );
      //setBeta( lhs.getBeta(), rhs.getBeta() );
      //setLambda( lhs.getLambda(), rhs.getLambda() );
      //setDeltaT( lhs.getDeltaT(), rhs.getDeltaT() );
      //setScaling( lhs.getScaling(), rhs.getScaling() );

      //lhs.checkSpaceOp( this );
      //rhs.checkSpaceOp( this );
      return false;
    }

    virtual bool split( const LambdaSplitting& split )
    {
      //do lambda splitting
      return true;
    }

    //virtual bool split( const DeltaTSplitting& split )
    //{
    //  // do DeltaTSplitting
    //  return true;
    //}

    template< class ScaleImp >
    bool scale( const DomainDiscreteFunctionType& df, const ScaleImp& scal )
    {
      return false;
    }

    //bool scale( const DomainFunctionType& df, const VelocityPressureSplitting& split )
    //{
    //  scal.apply( df );
    //}


    void setTime(const double time)
    {
      model_.setTime( time );
    }

    double timeStepEstimate() const
    {
      //TODO implementation needed!
      return 0.0001;
    }

    //! evaluate the spatial operator
    void operator()( const DomainDiscreteFunctionType& u, RangeDiscreteFunctionType& dest ) const
    {
      assembler_( u, dest );
    }

    //! evaluate rhs of the spatial operator
    void operator()( RangeDiscreteFunctionType& dest ) const
    {
      assembler_( dest );
    }

    inline const DiscreteFunctionSpaceType& space() const
    {
	    return space_;
    }

    inline DiscreteFunctionSpaceType& space()
    {
	    return space_;
    }

    void jacobian ( const DomainDiscreteFunctionType &u, JacobianOperatorType &jac ) const
    {
      jacobianAssembler_.prepare( jac );

      const bool isAffine = true;

      if( isAffine )
      {
        jacobianAssembler_.jacobian( jac );
      }
      else
      {
        //nonlinear assembler
        jacobianAssembler_.jacobian( u, jac );
      }

      jacobianAssembler_.finalize( jac );
    }

    //for RK-Solver
    //========================
    inline void switchupwind()
    {
      //TODO implement if needed
    }

    inline double computeTime() const
    {
    }

    inline size_t numberOfElements () const
    {
      return numberOfElements_;
    }

    //for algorithm
    //==========================
    const ModelType& model() const { return model_; }


  private:
    GridPartType& gridPart_;
    ProblemType& problem_;
    ModelType model_;
    DiscreteFunctionSpaceType space_;
    AssemblerType assembler_;
    JacobianAssemblerType jacobianAssembler_;
    mutable size_t numberOfElements_;

  };



  //----------------------------
  //for SubSteadyStateAlgorithm
  //----------------------------
  template< class Traits, class JacobianAssemblerImp >
  class AssembledPoissonOperator
  {
    typedef JacobianAssemblerImp JacobianAssemblerType;

    typedef typename JacobianAssemblerType::ContainerType ContainerType;

    typedef typename JacobianAssemblerType::ModelType ModelType;

  public:
                                      //GridPartType& gridPart, ProblemType& problem,
                                      //ExtraParameterTupleType& tuple,
                                      //const std::string name = "" )

    AssembledPoissonOperator( ContainerType& container,
                              const ModelType& model,
                              const bool calculateFluxes = true,
                              const bool strongBC = false )
      : container_( container ),
        jacAssembler_( container_, model, calculateFluxes, strongBC )
    {}

    void assemble() const
    {
      jacAssembler_.jacobian( container_.matrix(), container_.rhs() );
    }


    JacobianAssemblerType& assembler()
    {
      return jacAssembler_;
    }

  private:
    ContainerType& container_;
    JacobianAssemblerType jacAssembler_;
  };














  template <class Entity,
            class Quadrature,
            class Domain >
  class EntityStorage
  {
  public:
    typedef Entity                            EntityType;
    typedef Quadrature                        QuadratureType;
    typedef Domain                                         DomainType;
    typedef typename EntityType::Geometry::LocalCoordinate LocalDomainType;


    EntityStorage( const EntityType& entity, const double volume,
                   const DomainType& position, const LocalDomainType& local )
      : en_(entity),
        volume_(volume),
        position_( position ),
        localPosition_( local )
    {}

    const EntityType &entity() const { return en_; }
    const double volume() const { return volume_; }
    const LocalDomainType& localPosition() const { return localPosition_; }
    const DomainType& position() const { return position_; }

  private:
    const EntityType &en_;
    const double volume_;
    const DomainType& position_;
    const LocalDomainType& localPosition_;
  };

  template <class Entity, class Intersection, class Quadrature >
  class IntersectionStorage
  {
  public:
    typedef Entity                            EntityType;
    typedef Intersection                      IntersectionType;
    typedef Quadrature                        QuadratureType;

    IntersectionStorage( const IntersectionType& intersection,
                         const EntityType &entity,
                         const QuadratureType& quad,
                         const double volume )
      : intersection_( intersection ),
        entity_(entity),
        quad_( quad ),
        volume_(volume)
    {}

    const IntersectionType& intersection() const { return intersection_; }

    const EntityType &entity() const { return entity_; }

    const double volume() const { return volume_; }

    const Quadrature& quadrature() const { return quad_; }


    private:
    const IntersectionType& intersection_;
    const EntityType &entity_;
    const QuadratureType& quad_;
    const double volume_;
  };




  /**
   * \brief Assembles the primal DG matrix.
   *
   * \ingroup AssemblyOperator
   */
  template <class Traits>
  class DGPrimalMatrixAssembly
    : public AssemblerDefault< typename Traits::DomainDiscreteFunctionType, typename Traits::RangeDiscreteFunctionType >,
      public JacobianAssemblerDefault< typename Traits::MatrixContainerType >
  {
    public:
    typedef typename Traits::DomainDiscreteFunctionType           DomainDiscreteFunctionType;
    typedef typename Traits::RangeDiscreteFunctionType            RangeDiscreteFunctionType;
    typedef typename Traits::MatrixContainerType                  LinearOperatorType;
    typedef LinearOperatorType                                    MatrixType;
    typedef typename Traits::DomainDiscreteFunctionType           DestinationType;

    typedef typename Traits::ModelType                            ModelType;
    static const bool hasDiffusion = ModelType::hasDiffusion;



    //static const bool hasAdvection = ModelType::hasAdvection;
    //static const bool hasSource = ModelType::hasSource;
    //static const bool hasNonLinearDiffusion = ModelType::hasNonLinearDiffusion;
    //static const bool hasNonLinearAdvection = ModelType::hasNonLinearAdvection;
    //static const bool hasNonLinearSource = ModelType::hasNonLinearSource;

    typedef typename DestinationType::DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::IteratorType      IteratorType;
    typedef typename IteratorType::Entity                         EntityType;
    typedef typename EntityType::Geometry                         GeometryType;

    typedef typename DiscreteFunctionSpaceType::GridPartType      GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType        DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType         RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType   DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType    RangeFieldType;

    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    typedef typename DestinationType::LocalFunctionType           LocalFunctionType;

    typedef typename GridPartType::IntersectionIteratorType       IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection       IntersectionType;
    typedef typename IntersectionType::Geometry                   IntersectionGeometryType;

    typedef Fem::Parameter                                        ParameterType;

    // need treatment of non conforming grids
    typedef typename Traits::FaceQuadratureType                   FaceQuadratureType;
    typedef typename Traits::VolumeQuadratureType                 VolumeQuadratureType;


    typedef IntersectionStorage< EntityType, IntersectionType, FaceQuadratureType >
                                                                  IntersectionStorageType;
    typedef EntityStorage< EntityType, VolumeQuadratureType, DomainType >
                                                                  EntityStorageType;

    typedef typename MatrixType::LocalMatrixType                  LocalMatrixType;

    typedef typename Traits::AdvectionFluxType                    AdvectionFluxType;
    typedef ExtendedDGPrimalDiffusionFlux<DiscreteFunctionSpaceType, ModelType,
                                          typename Traits::DiffusionFluxType::ParameterType >
                                                                  DiffusionFluxType;

    typedef DefaultAssemblerTraits< Traits::polynomialOrder,
                                    typename Traits::AnalyticalTraitsType,
                                    typename Traits::MatrixContainerType,
                                    typename Traits::AdvectionFluxType,
                                    DiffusionFluxType,
                                    typename Traits::DomainDiscreteFunctionType,
                                    typename Traits::RangeDiscreteFunctionType > NewAssemblerTraits;

    typedef DGPrimalMatrixDiscreteModel< NewAssemblerTraits >     DiscreteModelType;


    typedef Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType >
                                                                  TemporaryLocalFunctionType;
    typedef Dune::Fem::TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType >
                                                                  TemporaryLocalMatrixType;

    struct RangeValues
    {
      typedef std::vector< std::vector< RangeType > > VectorType;

      RangeValues(const VectorType &vec, int col = -1 ) : zero_(0), col_(col), vec_(vec)
      {}

      RangeValues() : zero_(0), col_(-1), vec_(VectorType())
      {}

      const RangeType& at( int i ) const { return this->operator[] ( i ); }
      const RangeType &operator[](int row) const
      {
        if ( col_ == -1 )
          return zero_;
        else
          return vec_[row][col_];
      }
    private:
      const RangeType zero_;
      const int col_;
      const VectorType& vec_;
    };

    struct JacobianRangeValues
    {
      typedef std::vector< std::vector< JacobianRangeType > > VectorType;

      JacobianRangeValues(const VectorType &vec, int col = -1 ) : zero_(0), col_(col), vec_(vec)
      {}

      JacobianRangeValues() : zero_(0), col_(-1), vec_(VectorType())
      {}

      const JacobianRangeType& at( int i ) const { return this->operator[] ( i ); }
      const JacobianRangeType &operator[](int row) const
      {
        if ( col_ == -1 )
          return zero_;
        else
          return vec_[row][col_];
      }
    private:
      const JacobianRangeType zero_;
      const int col_;
      const VectorType& vec_;
    };


    template< class GlobalMatrixImp >
    struct LocalMatrixStorage
    {
      template< class DomainSpaceImp, class RangeSpaceImp >
      LocalMatrixStorage( GlobalMatrixImp* const& matrix,
                          DomainSpaceImp& domainSpace,
                          RangeSpaceImp& rangeSpace  )
        : localMatrix_( domainSpace, rangeSpace ),
          matrix_( matrix ),
          locked_( 0 )
      {
      }

      template< class LocalMatrixImp >
      void init( const EntityType& domainEntity, const EntityType& rangeEntity, LocalMatrixImp& localMatrix )
      {
        if( matrix_ )
        {
          if( locked_ <= 0 )
          {
            localMatrix.init( domainEntity, rangeEntity );
            localMatrix.clear();
          }
          locked_++;
        }
      }

      template< class LocalMatrixImp >
      void finalize( LocalMatrixImp& localMatrix  )
      {
        if( matrix_ )
        {
          assert( locked_ > 0 );
          if( locked_ <= 0 )
          {
            const EntityType &domainEntity = localMatrix.domainEntity();
            const EntityType &rangeEntity = localMatrix.rangeEntity();
            (matrix_)->addLocalMatrix( domainEntity, rangeEntity, localMatrix );
          }
          locked_--;
        }
      }

      void init( const EntityType& domainEntity, const EntityType& rangeEntity )
      {
        init( domainEntity, rangeEntity, localMatrix_ );
      }

      void finalize()
      {
        finalize( localMatrix_ );
      }

      const TemporaryLocalMatrixType& operator()() const
      {
        return localMatrix_;
      }

      TemporaryLocalMatrixType& operator()()
      {
        return localMatrix_;
      }

    private:
      TemporaryLocalMatrixType localMatrix_;
      GlobalMatrixImp* const& matrix_;
      int locked_;
    };

    template< class DestinationImp >
    struct LocalFunctionStorage
    {
      template< class DomainSpaceImp >
      LocalFunctionStorage( DestinationImp* const& rhs, DomainSpaceImp& space )
        : localFunction_( space ),
          rhs_( rhs ),
          locked_( 0 )
      {
      }

      template< class LocalFunctionImp >
      void init( const EntityType& entity, LocalFunctionImp& localFunction )
      {
        if( rhs_ )
        {
          if( locked_ <= 0 )
          {
            localFunction.init( entity );
            localFunction.clear();
          }
          locked_++;
        }
      }

      template< class LocalFunctionImp >
      void finalize( LocalFunctionImp& localFunction  )
      {
        if( rhs_ )
        {
          assert( locked_ > 0 );
          if( locked_ <= 0 )
          {
            const EntityType &entity = localFunction.entity();
            rhs_->addLocalDofs( entity, localFunction.localDofVector() );
          }
          locked_--;
        }
      }

      void init( const EntityType& entity )
      {
        init( entity, localFunction_ );
      }

      void finalize()
      {
        finalize( localFunction_ );
      }

      const TemporaryLocalFunctionType& operator()() const
      {
        return localFunction_;
      }

      TemporaryLocalFunctionType& operator()()
      {
        return localFunction_;
      }

    private:
      TemporaryLocalFunctionType localFunction_;
      DestinationImp* const& rhs_;
      int locked_;
    };


    public:

    typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;


    //! constructor for DG matrix assembly
    DGPrimalMatrixAssembly( GridPartType& gridPart,
                            const ModelType& model,
                            const bool calculateFluxes = true,
                            const bool strongBC = false )
      : gridPart_( gridPart ),
        model_( model ),
        space_( gridPart_ ),
        stencil_( space_, space_ ),
        time_( 0 ),
        matrix_( nullptr ),
        rhs_( nullptr ),
        argU_( nullptr ),
        localMatrixStorageEn_( matrix_, space_, space_ ),
        localMatrixStorageNb_( matrix_, space_, space_ ),
        localMatrixStorageNbEn_( matrix_, space_, space_ ),
        localMatrixStorageNbNb_( matrix_, space_, space_ ),
        localFunctionStorage_( rhs_, space_ ),
        discreteModel_( model_, space_ ),
        calculateFluxes_( calculateFluxes ),
        useStrongBoundaryCondition_( strongBC ),
        uZero_( 0 ),
        uJacZero_( 0 ),
        null_(),
        jacNull_()
    {
    }

    const DiscreteModelType &discreteModel () const { return discreteModel_; }
    DiscreteModelType &discreteModel () { return discreteModel_; }

    const ModelType &model () const { return model_; }
    ModelType &model () { return model_; }

    const StencilType& stencil() const
    {
      return stencil_;
    }


    const DiscreteFunctionSpaceType &space() const
    {
      return space_;
    }

    //const typename DiffusionFluxType::DiscreteGradientSpaceType &gradientSpace() const
    //{
    //  return diffusionFlux_.gradientSpace();
    //}

    size_t maxNumScalarBasisFunctions( const DiscreteFunctionSpaceType& space ) const
    {
      return space.blockMapper().maxNumDofs() * DiscreteFunctionSpaceType :: localBlockSize ;
    }

    void setTime ( double time )
    {
      time_ = time;
    }

    void prepare( MatrixType& matrix ) const
    {
      matrix.reserve( stencil() );
      matrix.clear();
    }

    void finalize( MatrixType& matrix ) const
    {
      // finish matrix build process
      matrix.communicate();
    }

    //f           [affine or nonlinear model]
    void operator()( RangeDiscreteFunctionType& rhs ) const
    {
      setArg();
      setMatrix();
      setRhs( rhs );

      assembleRhs();
    }

    //F(u)        [affine or nonlinear model]
    void operator()( const DomainDiscreteFunctionType& u, RangeDiscreteFunctionType& dest ) const
    {
      //MatrixType matrix;
    }

    /*
     * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
     * implementation.
     *
     * Linearize around u?
     */
    void jacobian( const DomainDiscreteFunctionType& u, MatrixType& matrix, DestinationType& rhs ) const
    {
      setArg( u );
      setMatrix( matrix );
      setRhs( rhs );
      assemble();
    }

    void jacobian( const DomainDiscreteFunctionType& u, MatrixType& matrix ) const
    {
      setArg( u );
      setMatrix( matrix );
      setRhs();
      assemble();
    }


    void jacobian( const DomainDiscreteFunctionType& u, DestinationType& rhs ) const
    {
      setArg( u );
      setMatrix();
      setRhs( rhs );
      assemble();
    }

    void jacobian( MatrixType& matrix, DestinationType& rhs ) const
    {
      setArg();
      setMatrix( matrix );
      setRhs( rhs );
      assemble();
    }

    void jacobian( MatrixType& matrix ) const
    {
      setArg();
      setMatrix( matrix );
      setRhs();
      assemble();
    }



    /*
     * Assemble Matrix for Elliptic Problem using the DGPrimalDIffusionFlux
     * implementation.
     */
    void assemble() const
    {
      Dune::Timer timer ;

      if( matrix_ )
        matrix_->clear();
      if( rhs_ )
        rhs_->clear();

      discreteModel().initialize();

      for( const EntityType &entity : space_ )
        assembleLocal( entity  );

      if( Dune::Fem::Parameter::verbose() )
        std::cout << "DG( " << space_.grid().size( 0 ) << " ) matrix assemble took " << timer.elapsed() << " sec." << std::endl;
    }

    void assembleLocal( const EntityType& en ) const
    {
      //init
      localFunctionStorage_.init( en );
      localMatrixStorageEn_.init( en, en );

      volumeIntegral( en, localMatrixStorageEn_(), localFunctionStorage_() );

      for( const IntersectionType &intersection : intersections( space().gridPart(), en ) )
      {
        assembleLocal( intersection );
      }

      //update
      localMatrixStorageEn_.finalize();
      localFunctionStorage_.finalize();
    }

    void assembleLocal( const IntersectionType& intersection ) const
    {
      const EntityType& entity = intersection.inside();

      //init
      localFunctionStorage_.init( entity );
      localMatrixStorageEn_.init( entity, entity );

      for( const IntersectionType &intersection : intersections( space().gridPart(), entity ) )
      {
        if( isSurfaceIntegral( intersection ) )
        {
          const EntityType& neighbor = intersection.outside();

          //init
          localMatrixStorageNb_.init( entity, neighbor );
          localMatrixStorageNbEn_.init( neighbor, entity );
          localMatrixStorageNbNb_.init( neighbor, neighbor );

          if( intersection.conforming() )
            surfaceIntegral< true >( intersection, localMatrixStorageEn_(), localMatrixStorageNb_(), localMatrixStorageNbEn_(), localMatrixStorageNbNb_(), localFunctionStorage_() );
          else
            surfaceIntegral< false >( intersection, localMatrixStorageEn_(), localMatrixStorageNb_(), localMatrixStorageNbEn_(), localMatrixStorageNbNb_(), localFunctionStorage_() );

          //update
          localMatrixStorageNb_.finalize();
          localMatrixStorageNbEn_.finalize();
          localMatrixStorageNbNb_.finalize();
        }
        else if( isBoundaryIntegral( intersection) )
          boundaryIntegral( intersection, localMatrixStorageEn_(), localFunctionStorage_() );

      }

      //update
      localMatrixStorageEn_.finalize();
      localFunctionStorage_.finalize();
    }

    // assemble vector containing boundary fluxes for right hand side
    void assembleRhs() const
    {
      if( rhs_ )
        rhs_->clear();

      if( hasDiffusion )
        discreteModel().initialize();

      for( const EntityType &entity : space() )
      {
        //init
        localFunctionStorage_.init( entity );
        localMatrixStorageEn_.init( entity, entity );

        for( const IntersectionType &intersection : intersections( space().gridPart(), entity ) )
          if( !isSurfaceIntegral( intersection ) && !isBoundaryIntegral( intersection) )
            boundaryIntegral( intersection, localMatrixStorageEn_(), localFunctionStorage_() );

        //update
        localMatrixStorageEn_.finalize();
        localFunctionStorage_.finalize();
      }
    }

    template <class LocalMatrix, class LocalFunction>
    void volumeIntegral( const EntityType& entity,
                         LocalMatrix& localMatrix,
                         LocalFunction& rhsLocal ) const
    {
      typedef RangeType           RangeTuple;
      typedef JacobianRangeType   JacobianTuple;
      typedef ElementQuadraturePointContext< EntityType, VolumeQuadratureType, RangeTuple, JacobianTuple >      LocalEvaluationType;

      const GeometryType &geometry = entity.geometry();

      VolumeQuadratureType quadrature( entity, elementQuadOrder( space_.order( entity ) ) );
      const size_t numQuadPoints = quadrature.nop();

      const BasisFunctionSetType &baseSet = localMatrix.domainBasisFunctionSet();
      const unsigned int numBasisFunctionsEn = baseSet.size();

      const double volume = geometry.volume();

      RangeType arhs;
      JacobianRangeType adrhs;

      resize( numQuadPoints );
      for( size_t pt = 0; pt < numQuadPoints; ++pt )
      {
        baseSet.evaluateAll( quadrature[ pt ], phiEn_ );
        baseSet.jacobianAll( quadrature[ pt ], dphiEn_ );

        LocalEvaluationType local( entity, quadrature, uZero_, uJacZero_, pt, time_, volume );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( local.position() );

        // compute rhs
        discreteModel().rhs( local, arhs, adrhs );

        if( matrix_ )
        {
          for( unsigned int i = 0; i < numBasisFunctionsEn; ++i )
          {
            RangeType aphi(0);
            JacobianRangeType adphi(0);
            //TODO save time step restriction?
            const double dtEst = discreteModel().analyticalFluxAndSource( local, phiEn_[i], dphiEn_[i], aphi, adphi);

            // subtract affine part and move to left hand side
            aphi -= arhs;
            aphi *= -1;
            adphi -= adrhs;

            // get column object and call axpy method
            localMatrix.column( i ).axpy( phiEn_, dphiEn_, aphi, adphi, weight );
          }
        }

        if( rhs_ )
        {
          // assemble rhs
          arhs  *=  weight;
          adrhs *= -weight;
          rhsLocal.axpy( quadrature[ pt ], arhs, adrhs );
        }

      }
    }

    template< bool conforming, class LocalMatrix, class LocalRHS >
    void surfaceIntegral( const IntersectionType& intersection,
                          LocalMatrix& localMatrixEn,
                          LocalMatrix& localMatrixNb,
                          LocalMatrix& localMatrixNbEn,
                          LocalMatrix& localMatrixNbNb,
                          LocalRHS& rhsLocal ) const
    {
      if ( space().continuous(intersection) ) return;

      // make sure we got the right conforming statement
      assert( intersection.conforming() == conforming );

      // use IntersectionQuadrature to create appropriate face quadratures
      typedef Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
      typedef typename IntersectionQuadratureType::FaceQuadratureType QuadratureImp;

      const EntityType& entity = intersection.inside();
      const EntityType& neighbor = intersection.outside();

      const GeometryType &geometry = entity.geometry();

      const int entityOrder   = space().order( entity );
      const int neighborOrder = space().order( neighbor );
      const int polOrdOnFace = std::max( entityOrder, neighborOrder );
      // create intersection quadrature (without neighbor check)
      IntersectionQuadratureType interQuad( space().gridPart(), intersection, faceQuadOrder( polOrdOnFace ), true);
      // get appropriate references
      const QuadratureImp &faceQuadInside  = interQuad.inside();
      const QuadratureImp &faceQuadOutside = interQuad.outside();

      const size_t numFaceQuadPoints = faceQuadInside.nop();

      const double entityVolume = entity.geometry().volume();
      const double neighborVolume = neighbor.geometry().volume();

      typedef IntersectionStorage< EntityType, IntersectionType, QuadratureImp > LocalIntersectionStorageType;
      LocalIntersectionStorageType leftStorage( intersection, entity, faceQuadInside, entityVolume );
      LocalIntersectionStorageType rightStorage( intersection, neighbor, faceQuadOutside, neighborVolume );

      typedef RangeType           RangeTuple;
      typedef JacobianRangeType   JacobianTuple;
      typedef IntersectionQuadraturePointContext<IntersectionType, EntityType, QuadratureImp, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;

      RangeType null(0);
      JacobianRangeType jacNull(0);

      //initialize intersection
      discreteModel().initializeIntersection( intersection, model_.time(),faceQuadInside, faceQuadOutside, null, null );

      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        IntersectionLocalEvaluationType left( intersection, entity, faceQuadInside, null, jacNull, pt, model_.time(), entityVolume );

        discreteModel().numericalFlux(left, left,
                                      null, null, jacNull, jacNull,
                                      rhsValueEn_[pt], rhsValueNb_[pt], rhsDValueEn_[pt], rhsDValueNb_[pt] );
      }

      if( matrix_ )
      {
        // get local basis function sets for face entries
        const BasisFunctionSetType &baseSet = localMatrixEn.domainBasisFunctionSet();
        const BasisFunctionSetType &baseSetNb = localMatrixNb.rangeBasisFunctionSet();

        // get neighbor's base function set
        const unsigned int numBasisFunctionsEn = baseSet.size();
        const unsigned int numBasisFunctionsNb = baseSetNb.size();

        // only do one sided evaluation if the polynomial orders match
        const bool oneSidedEvaluation = ( numBasisFunctionsEn == numBasisFunctionsNb );

        const bool updateOnNeighbor = space().gridPart().indexSet().index(entity) > space().gridPart().indexSet().index(neighbor) ;

        // only do one sided evaluation if the polynomial orders match
        if( updateOnNeighbor && oneSidedEvaluation )
          return;

        resize( numFaceQuadPoints );
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const double weight = faceQuadInside.weight( pt );

          baseSet.evaluateAll( faceQuadInside[ pt ], phiEn_ );
          baseSet.jacobianAll( faceQuadInside[ pt ], dphiEn_ );
          baseSetNb.evaluateAll( faceQuadOutside[ pt ], phiNb_ );
          baseSetNb.jacobianAll( faceQuadOutside[ pt ], dphiNb_ );

          for( unsigned int i = 0; i < numBasisFunctionsEn; ++i )
          {
            IntersectionLocalEvaluationType right( intersection, neighbor, faceQuadOutside, null, jacNull, pt, model_.time(), neighborVolume );
            {
              IntersectionLocalEvaluationType left( intersection, entity, faceQuadInside, phiEn_[i], dphiEn_[i], pt, model_.time(), entityVolume );
              discreteModel().numericalFlux(left, right,
                                            phiEn_[i], null, dphiEn_[i], jacNull,
                                            valueEn_, valueNb_, dvalueEn_, dvalueNb_ );

              valueEn_ -= rhsValueEn_[pt];
              dvalueEn_ -= rhsDValueEn_[pt];
              valueNb_ -= rhsValueNb_[pt];
              dvalueNb_ -= rhsDValueNb_[pt];

              localMatrixEn.column( i ).axpy( phiEn_, dphiEn_,
                                              valueEn_, dvalueEn_, weight );
              localMatrixNb.column( i ).axpy( phiNb_, dphiNb_,
                                              valueNb_, dvalueNb_, -weight );
            }
            if( oneSidedEvaluation ) //i.e. numBasisFunctionsEn == numBasisFunctionsNb
            {
              IntersectionLocalEvaluationType left( intersection, entity, faceQuadInside, phiNb_[i], dphiNb_[i], pt, model_.time(), entityVolume );
              discreteModel().numericalFlux(left, right,
                                            phiNb_[i], null, dphiNb_[i], jacNull,
                                            valueEn_, valueNb_, dvalueEn_, dvalueNb_ );

              valueEn_ -= rhsValueEn_[pt];
              dvalueEn_ -= rhsDValueEn_[pt];
              valueNb_ -= rhsValueNb_[pt];
              dvalueNb_ -= rhsDValueNb_[pt];
              localMatrixNbNb.column( i ).axpy( phiNb_, dphiNb_,  // +
                                                valueNb_, dvalueNb_, -weight );
              localMatrixNbEn.column( i ).axpy( phiEn_, dphiEn_,  // -
                                                valueEn_, dvalueEn_, weight );
            }
          }
        }
       // // compute fluxes and assemble matrix
       // for( unsigned int i = 0; i < numBasisFunctionsEn; ++i )
       // {
       //   // compute flux for one base function, i.e.:  uLeft=phiEn_[.][i]; uRight=0
       //   discreteModel().numericalFlux(leftStorage, rightStorage,
       //                                 RangeValues(phiEn_,i), null_, JacobianRangeValues(dphiEn_,i),jacNull_,
       //                                 valueEn_, valueNb_, dvalueEn_, dvalueNb_ );

       //   for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
       //   {
       //     const double weight = faceQuadInside.weight( pt );
       //     valueEn_[pt] -= rhsValueEn_[pt];
       //     dvalueEn_[pt] -= rhsDValueEn_[pt];
       //     valueNb_[pt] -= rhsValueNb_[pt];
       //     dvalueNb_[pt] -= rhsDValueNb_[pt];
       //     localMatrixEn.column( i ).axpy( phiEn_[pt], dphiEn_[pt],
       //                                     valueEn_[pt], dvalueEn_[pt], weight );
       //     localMatrixNb.column( i ).axpy( phiNb_[pt], dphiNb_[pt],
       //                                     valueNb_[pt], dvalueNb_[pt], -weight );
       //   }
       // }
      }

      if( rhs_ )
      {
        // now move affine part to right hand side
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const double weight = faceQuadInside.weight( pt );
          rhsValueEn_[pt] *= -weight;
          rhsDValueEn_[pt] *= -weight;
        }
        rhsLocal.axpyQuadrature( faceQuadInside, rhsValueEn_, rhsDValueEn_ );
      }
    }


    template< class LocalMatrix, class LocalRHS >
    void boundaryIntegral( const IntersectionType& intersection,
                           LocalMatrix& localMatrixEn,
                           LocalRHS& rhsLocal ) const
    {
      const EntityType& entity = intersection.inside();
      const GeometryType &geometry = entity.geometry();

      FaceQuadratureType faceQuadInside(space_.gridPart(), intersection,
                                        faceQuadOrder( space_.order( entity ) ),
                                        FaceQuadratureType::INSIDE);

      const size_t numFaceQuadPoints = faceQuadInside.nop();


      const double entityVolume = entity.geometry().volume();

      typedef RangeType           RangeTuple;
      typedef JacobianRangeType   JacobianTuple;
      typedef IntersectionQuadraturePointContext<IntersectionType, EntityType, FaceQuadratureType, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;

      //assemble rhs
      RangeType null(0);
      JacobianRangeType jacNull(0);

      //initialize boundary
      discreteModel().initializeBoundary( intersection, model_.time(), faceQuadInside, null );

      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        IntersectionLocalEvaluationType left( intersection, entity, faceQuadInside, null, jacNull, pt, model_.time(), entityVolume );

        discreteModel().boundaryValues(left, bndValues_[pt]); //inside initIntersection()?

        // first compute affine part of boundary flux
        discreteModel().boundaryFlux(left, null, bndValues_[pt], jacNull,
                                     rhsValueNb_[pt],rhsDValueNb_[pt]);
      }

      if( matrix_ )
      {
        const BasisFunctionSetType &baseSet = localMatrixEn.domainBasisFunctionSet();
        const unsigned int numBasisFunctionsEn = baseSet.size();

        resize( numFaceQuadPoints );
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const double weight = faceQuadInside.weight( pt );

          baseSet.evaluateAll( faceQuadInside[ pt ], phiEn_ );
          baseSet.jacobianAll( faceQuadInside[ pt ], dphiEn_ );

          //TODO initialize boundary
          //TODO discreteModel().initializeBoundary( intersection, model_.time(), faceQuadInside, phiEn_ );

          for( unsigned int i = 0; i < numBasisFunctionsEn; ++i )
          {
            IntersectionLocalEvaluationType left( intersection, entity, faceQuadInside, phiEn_[i], dphiEn_[i], pt, model_.time(), entityVolume );
            discreteModel().boundaryFlux(left,
                                         phiEn_[i], bndValues_[pt], dphiEn_[i],
                                         valueEn_, dvalueEn_ );

            valueEn_ -= rhsValueNb_[pt];
            dvalueEn_ -= rhsDValueNb_[pt];

            localMatrixEn.column( i ).axpy( phiEn_, dphiEn_,
                                            valueEn_, dvalueEn_, weight );

          }
        }
      }

      if( rhs_ )
      {
        // now move affine part to right hand side
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const double weight = faceQuadInside.weight( pt );
          rhsValueNb_[ pt ] *= -weight;
          rhsDValueNb_[ pt ] *= -weight;
        }
        rhsLocal.axpyQuadrature( faceQuadInside, rhsValueNb_, rhsDValueNb_ );
      }
    }

    template< class Quadrature, class BasisFunctionSet, class RetVector, class RetJacobianVector >
    int evaluateBasisFunctions( const Quadrature& quad, const BasisFunctionSet basisSet, RetVector& ret, RetJacobianVector& dret ) const
    {
      const size_t numQuadPoints = quad.nop();
      resize( numQuadPoints );

      // evaluate base functions
      for( size_t pt = 0; pt < numQuadPoints; ++pt )
      {
        basisSet.evaluateAll( quad[ pt ], ret[ pt ] );
        basisSet.jacobianAll( quad[ pt ], dret[ pt ] );
      }
      return numQuadPoints;
    }

    bool isSurfaceIntegral( const IntersectionType& intersection ) const
    {
      return intersection.neighbor() && calculateFluxes_;
    }

    bool isBoundaryIntegral( const IntersectionType& intersection ) const
    {
      return intersection.boundary() && ! useStrongBoundaryCondition_;
    }


    void setMatrix( const MatrixType& matrix ) const
    {
      //std::cout << "setMatrix: " << &matrix << std::endl;
      matrix_ = const_cast< MatrixType* >(&matrix);
    }

    void setRhs( const DestinationType& rhs ) const
    {
      //std::cout << "setRhs: " << &rhs << std::endl;
      rhs_ = const_cast< DestinationType* >(&rhs);
    }

    void setArg( const DestinationType& argU ) const
    {
      //std::cout << "setArg: " << &argU << std::endl;
      argU_ = const_cast< DestinationType* >(&argU);
    }

    void setMatrix() const
    {
      //std::cout << "setMatrix: nullptr " << std::endl;
      matrix_ = nullptr;
    }

    void setRhs() const
    {
      //std::cout << "setRhs: nullptr " << std::endl;
      rhs_ = nullptr;
    }

    void setArg() const
    {
      //std::cout << "setArg: nullptr " << std::endl;
      argU_ = nullptr;
    }

    //MatrixType* matrix()
    //{
    //  assert( matrix_ );
    //  return matrix_;
    //}

    //DestinationType* rhs()
    //{
    //  assert( rhs_ );
    //  return rhs_;
    //}

    //DestinationType* arg()
    //{
    //  assert( argU_ );
    //  return argU_;
    //}

    void testSymmetry() const
    {
      // not check at the moment
      return;
    }

  private:
    int elementQuadOrder( int polOrder ) const
    {
      return 2*polOrder;
    }
    int faceQuadOrder( int polOrder ) const
    {
      return 2*polOrder;
    }

    void resize( unsigned int numQuadPoints, unsigned int maxNumBasisFunctions = -1 ) const
    {
      if( maxNumBasisFunctions == -1 )
        maxNumBasisFunctions = maxNumScalarBasisFunctions( space() );

      //if (phiEn_.size() >= numQuadPoints)
      //  return;

      //phiEn_.resize( numQuadPoints );
      //for (unsigned int i=0;i<numQuadPoints;++i) phiEn_[i].resize(maxNumBasisFunctions);
      //dphiEn_.resize( numQuadPoints );
      //for (unsigned int i=0;i<numQuadPoints;++i) dphiEn_[i].resize(maxNumBasisFunctions);
      //phiNb_.resize( numQuadPoints );
      //for (unsigned int i=0;i<numQuadPoints;++i) phiNb_[i].resize(maxNumBasisFunctions);
      //dphiNb_.resize( numQuadPoints );
      //for (unsigned int i=0;i<numQuadPoints;++i) dphiNb_[i].resize(maxNumBasisFunctions);
      //valueEn_.resize( numQuadPoints );
      //dvalueEn_.resize( numQuadPoints );
      //valueNb_.resize( numQuadPoints );
      //dvalueNb_.resize( numQuadPoints );
      //

      if( phiEn_.size() < maxNumBasisFunctions )
      {
        phiEn_.resize( maxNumBasisFunctions );
        dphiEn_.resize( maxNumBasisFunctions );
        phiNb_.resize( maxNumBasisFunctions );
        dphiNb_.resize( maxNumBasisFunctions );
      }

      if( rhsValueEn_.size() >= numQuadPoints )
        return;

      rhsValueEn_.resize( numQuadPoints );
      rhsDValueEn_.resize( numQuadPoints );
      rhsValueNb_.resize( numQuadPoints );
      rhsDValueNb_.resize( numQuadPoints );
      bndValues_.resize( numQuadPoints );
    }

    GridPartType&                      gridPart_;
    const ModelType&                   model_;
    DiscreteFunctionSpaceType          space_;
    StencilType                        stencil_;
    double                             time_;

    mutable MatrixType*                       matrix_;
    mutable DestinationType*                  rhs_;
    mutable DestinationType*                  argU_;

    mutable LocalMatrixStorage<MatrixType>                 localMatrixStorageEn_;
    mutable LocalMatrixStorage<MatrixType>                 localMatrixStorageNb_;
    mutable LocalMatrixStorage<MatrixType>                 localMatrixStorageNbEn_;
    mutable LocalMatrixStorage<MatrixType>                 localMatrixStorageNbNb_;
    mutable LocalFunctionStorage<DestinationType>          localFunctionStorage_;

    DiscreteModelType                  discreteModel_;
    const bool                         calculateFluxes_;
    const bool                         useStrongBoundaryCondition_ ;

    const RangeType                    uZero_;
    const JacobianRangeType            uJacZero_;

    RangeValues                        null_;
    JacobianRangeValues                jacNull_;

    // storage for all flux values
    //mutable std::vector< RangeType >         valueEn_;
    //mutable std::vector< JacobianRangeType > dvalueEn_;
    //mutable std::vector< RangeType >         valueNb_;
    //mutable std::vector< JacobianRangeType > dvalueNb_;
    mutable RangeType                        valueEn_;
    mutable JacobianRangeType                dvalueEn_;
    mutable RangeType                        valueNb_;
    mutable JacobianRangeType                dvalueNb_;
    mutable std::vector< RangeType >         rhsValueEn_;
    mutable std::vector< JacobianRangeType > rhsDValueEn_;
    mutable std::vector< RangeType >         rhsValueNb_;
    mutable std::vector< JacobianRangeType > rhsDValueNb_;
    mutable std::vector< RangeType >         bndValues_;

    // store all basis functions
    //mutable std::vector< std::vector< RangeType > >         phiEn_;
    //mutable std::vector< std::vector< JacobianRangeType > > dphiEn_;
    //mutable std::vector< std::vector< RangeType > >         phiNb_;
    //mutable std::vector< std::vector< JacobianRangeType > > dphiNb_;
    mutable std::vector< RangeType >         phiEn_;
    mutable std::vector< JacobianRangeType > dphiEn_;
    mutable std::vector< RangeType >         phiNb_;
    mutable std::vector< JacobianRangeType > dphiNb_;
  };






  template< class Traits >
  class DGPrimalMatrixDiscreteModel
  {

    typedef typename Traits::ModelType                            ModelType;

    typedef typename Traits::AdvectionFluxType                    AdvectionFluxType;
    typedef typename Traits::DiffusionFluxType                    DiffusionFluxType;

    typedef typename Traits::DomainDiscreteFunctionType           DestinationType;

    static const bool hasDiffusion = ModelType::hasDiffusion;
    //static const bool hasAdvection = ModelType::hasAdvection;
    //static const bool hasSource = ModelType::hasSource;
    //static const bool hasNonLinearDiffusion = ModelType::hasNonLinearDiffusion;
    //static const bool hasNonLinearAdvection = ModelType::hasNonLinearAdvection;
    //static const bool hasNonLinearSource = ModelType::hasNonLinearSource;

    typedef typename DestinationType::DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::IteratorType      IteratorType;
    typedef typename IteratorType::Entity                         EntityType;
    typedef typename EntityType::Geometry                         GeometryType;

    typedef typename DiscreteFunctionSpaceType::GridPartType      GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType        DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType         RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType   DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType    RangeFieldType;

    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    typedef typename DestinationType::LocalFunctionType           LocalFunctionType;

    typedef typename GridPartType::IntersectionIteratorType       IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection       IntersectionType;
    typedef typename IntersectionType::Geometry                   IntersectionGeometryType;

    struct VectorToTupleVector
    {
      typedef std::vector< RangeType > VectorType;
      VectorToTupleVector(const VectorType &vec) : vec_(vec)
      {}

      const RangeType& at( int i ) const { return this->operator[] ( i ); }
      const RangeType &operator[](int i) const
      {
        return vec_[ i ];
      }
    private:
      const VectorType& vec_;
    };



    class LinearDiscreteModel
    {


    public:

      //nonlinear model without affine shift
      template <class LocalEvaluation>
      double numericalFlux(const LocalEvaluation& left,
                           const LocalEvaluation& right,
                           const RangeType& uLeft,
                           const RangeType& uRight,
                           const JacobianRangeType& jacLeft,
                           const JacobianRangeType& jacRight,
                           RangeType& gLeft,
                           RangeType& gRight,
                           JacobianRangeType& gDiffLeft,
                           JacobianRangeType& gDiffRight ) const
      {}

      //linear model without affine shift
      template <class LocalEvaluation>
      double numericalFlux(const LocalEvaluation& left,
                           const LocalEvaluation& right,
                           const RangeType& uLinLeft,
                           const RangeType& uLinRight,
                           const JacobianRangeType& jacLinLeft,
                           const JacobianRangeType& jacLinRight,
                           const RangeType& uLeft,
                           const RangeType& uRight,
                           const JacobianRangeType& jacLeft,
                           const JacobianRangeType& jacRight,
                           RangeType& gLeft,
                           RangeType& gRight,
                           JacobianRangeType& gDiffLeft,
                           JacobianRangeType& gDiffRight ) const
      {}

    };



    typedef LinearDiscreteModel LinearDiscreteModelType;

  public:


    DGPrimalMatrixDiscreteModel( const ModelType& model, const DiscreteFunctionSpaceType& space )
      : model_( model ),
        space_( space ),
        advFlux_( model_ ),
        diffusionFlux_( space_.gridPart(), model_ /*,typename Traits::DiffusionFluxType::ParameterType( ParameterKey::generate( "", "dgdiffusionflux." ) ) */ )
    {}

    void initialize() const
    {
      diffusionFlux_.initialize( space_ );
    }


    template <class IntersectionStorageImp,class Value,class LiftingFunction>
    void lifting(const IntersectionStorageImp& left,
                 const IntersectionStorageImp& right,
                 const Value &valueEn,   /* or: localfunction */
                 const Value &valueNb,  /* or: localfunction */
                 LiftingFunction &lifting) const
    {
      VectorToTupleVector valEn( valueEn );
      VectorToTupleVector valNb( valueNb );
      if( hasDiffusion )
      {
        diffusionFlux_.initializeIntersection( left.intersection(), left.entity(), right.entity(), model_.time(), left.quadrature(), right.quadrature(), valEn, valNb, true );
        lifting += diffusionFlux_.getInsideLifting();
      }
    }

    //template <class IntersectionStorageImp, class RangeVector>
    //void boundaryValues( const IntersectionStorageImp local,
    //                     RangeVector &bndValues) const
    //{
    //  const RangeType uZero(0);
    //  const JacobianRangeType uJacZero( 0 );

    //  typedef RangeType           RangeTuple;
    //  typedef JacobianRangeType   JacobianTuple;
    //  typedef IntersectionQuadraturePointContext<typename IntersectionStorageImp::IntersectionType, typename IntersectionStorageImp::EntityType,
    //                                             typename IntersectionStorageImp::QuadratureType, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;

    //  const size_t numFaceQuadPoints = local.quadrature().nop();
    //  for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
    //  {
    //    IntersectionLocalEvaluationType local( local.intersection(), local.entity(), local.quadrature(), uZero, uJacZero, pt, model_.time(), local.volume() );
    //    model_.boundaryValue(local, uZero, bndValues[ pt ]);
    //  }
    //}

    template <class IntersectionStorageImp>
    void boundaryValues( const IntersectionStorageImp local,
                         RangeType &bndValue) const
    {
      const RangeType uZero(0);
      model_.boundaryValue(local, uZero, bndValue);
    }


    template <class QuadratureImp, class ArgumentVector >
    void initializeIntersection(const IntersectionType& it,
                                const double time,
                                const QuadratureImp& quadInner,
                                const QuadratureImp& quadOuter,
                                const ArgumentVector& uLeft,
                                const ArgumentVector& uRight) const
    {
      if( hasDiffusion )
      {
        // call diffusion flux
        diffusionFlux_.initializeIntersection( it, it.inside(), it.outside(),
                                               time, quadInner, quadOuter,
                                               uLeft, uRight );
      }
    }

    template <class QuadratureImp, class ArgumentVector >
    void initializeBoundary(const IntersectionType& it,
                            const double time,
                            const QuadratureImp& quadInner,
                            const ArgumentVector& uLeft) const
    {
      if( hasDiffusion )
      {
        if( diffusionFlux_.hasLifting() )
        {
          typedef IntersectionQuadraturePointContext< IntersectionType, EntityType, QuadratureImp, ArgumentVector, ArgumentVector > IntersectionQuadraturePointContextType;
          IntersectionQuadraturePointContextType local0( it, it.inside(), quadInner, uLeft[ 0 ], uLeft[ 0 ], 0, time, it.inside().geometry().volume() );
          const bool hasBoundaryValue = model_.hasBoundaryValue( local0 );

          const size_t quadNop = quadInner.nop();
          if( uBndVec_.size() < quadNop ) uBndVec_.resize( quadNop );

          for(size_t qp = 0; qp < quadNop; ++qp)
          {
            IntersectionQuadraturePointContextType
              local( it, it.inside(), quadInner, uLeft[ qp ], uLeft[ qp ], qp, time, it.inside().geometry().volume() );

            assert( hasBoundaryValue == model_.hasBoundaryValue( local ) );

            if( hasBoundaryValue )
              model_.boundaryValue(local, uLeft[ qp ], uBndVec_[ qp ] );
            else
              // do something bad to uBndVec as it shouldn't be used
              uBndVec_[ qp ] = std::numeric_limits< double >::quiet_NaN();
          }

          // call diffusion flux
          diffusionFlux_.initializeBoundary( it, it.inside(),
                                             time, quadInner,
                                             uLeft, uBndVec_ );
        }
      }
    }


    template <class IntersectionStorageImp, class RangeInVector, class JacobianRangeInVector, class RangeVector, class JacobianRangeVector>
    void numericalFlux( const IntersectionStorageImp& left,
                        const IntersectionStorageImp& right,
                        const RangeInVector& uLeft,
                        const RangeInVector& uRight,
                        const JacobianRangeInVector& jacLeft,
                        const JacobianRangeInVector& jacRight,
                        RangeVector& gLeft,
                        RangeVector& gRight,
                        JacobianRangeVector& gDiffLeft,
                        JacobianRangeVector& gDiffRight,
                        const bool initializeIntersection = true ) const
    {
      if( initializeIntersection )
      {
        initializeIntersection( left.intersection(), model_.time(),left.quadrature(), right.quadrature(), uLeft, uRight );
      }

      const size_t numFaceQuadPoints = left.quadrature().nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        typedef RangeType           RangeTuple;
        typedef JacobianRangeType   JacobianTuple;
        typedef IntersectionQuadraturePointContext<typename IntersectionStorageImp::IntersectionType, typename IntersectionStorageImp::EntityType,
                                                   typename IntersectionStorageImp::QuadratureType, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;
        IntersectionLocalEvaluationType left( left.intersection(), left.entity(),
                                              left.quadrature(), uLeft[ pt ], jacLeft[ pt ], pt, model_.time(), left.volume() );
        IntersectionLocalEvaluationType right( right.intersection(), right.entity(),
                                               right.quadrature(), uRight[ pt ], jacRight[ pt ], pt, model_.time(), right.volume() );

        numericalFlux( left, right, uLeft[ pt ], uRight[ pt ], jacLeft[ pt ], jacRight[ pt ],
                       gLeft[ pt ], gRight[ pt ], gDiffLeft[ pt ], gDiffRight[ pt ] );
      }
    }

    template <class IntersectionStorageImp, class RangeInVector, class RangeInBndVector, class JacobianRangeInVector, class RangeVector, class JacobianRangeVector>
    void boundaryFlux( const IntersectionStorageImp& local,
                       const RangeInVector& uLeft,
                       const RangeInBndVector& uRight,
                       const JacobianRangeInVector& jacLeft,
                       RangeVector& gLeft,
                       JacobianRangeVector& gDiffLeft ) const
    {
      if( hasDiffusion )
      {
        diffusionFlux_.initializeBoundary( local.intersection(), local.entity(), model_.time(), local.quadrature(), uLeft, uRight );
      }

      const size_t numFaceQuadPoints = local.quadrature().nop();
      for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
      {
        typedef RangeType           RangeTuple;
        typedef JacobianRangeType   JacobianTuple;
        typedef IntersectionQuadraturePointContext<typename IntersectionStorageImp::IntersectionType, typename IntersectionStorageImp::EntityType,
                                                   typename IntersectionStorageImp::QuadratureType, RangeTuple, JacobianTuple > IntersectionLocalEvaluationType;
        IntersectionLocalEvaluationType left( local.intersection(), local.entity(), local.quadrature(), uLeft[ pt ], jacLeft[ pt ], pt, model_.time(), local.volume() );

        boundaryFlux( left, uLeft[ pt ], uRight[ pt ], jacLeft[ pt ], gLeft[ pt ], gDiffLeft[ pt ] );
      }

    }


    //////////nonlinear case (including affine shift)
    ////////template< class ArgLeft, class ArgRight, class LinLeft, class LinRight, class Left, class Right, class RetLeft, class RetRight >
    ////////double numericalFlux( const ArgLeft& left, const ArgRight& right,
    ////////                      const LinLeft& uLinLeft, const LinRight& uLinRight,
    ////////                      const Left& uLeft, const Right& uRight,
    ////////                      RetLeft& retLeft, RetRight& retRight )
    ////////{}

    //////////linear case (including affine shift)
    ////////template< class ArgLeft, class ArgRight, class Left, class Right, class RetLeft, class RetRight >
    ////////double numericalFlux( const ArgLeft& left, const ArgRight& right,
    ////////                      const Left& uLeft, const Right& uRight,
    ////////                      RetLeft& retLeft, RetRight& retRight )
    ////////{}

    //////////affine
    ////////template< class ArgLeft, class ArgRight, class Left, class Right, class RetLeft, class RetRight >
    ////////double numericalFlux( const ArgLeft& left, const ArgRight& right,
    ////////                      const Left& uLeft, const Right& uRight,
    ////////                      RetLeft& retLeft, RetRight& retRight )
    ////////{}

    //nonlinear case (including affine shift)
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         const RangeType& uLinLeft,
                         const RangeType& uLinRight,
                         const JacobianRangeType& jacLinLeft,
                         const JacobianRangeType& jacLinRight,
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const JacobianRangeType& jacLeft,
                         const JacobianRangeType& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {}


    //linear case (including affine shift)
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const JacobianRangeType& jacLeft,
                         const JacobianRangeType& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      if( hasDiffusion )
      {
        diffusionFlux_.numericalFlux( left, right,
                                      uLeft, uRight, jacLeft, jacRight,
                                      gLeft, gRight, gDiffLeft, gDiffRight );
      }
      else
      {
        gLeft  = RangeType(0);
        gRight = RangeType(0);
        gDiffLeft = JacobianRangeType(0);
        gDiffRight = JacobianRangeType(0);
      }

      RangeType gLeftTmp;
      RangeType gRightTmp;

      advFlux_.numericalFlux( left, right,
                              uLeft, uRight, jacLeft, jacRight,
                              gLeftTmp, gRightTmp);
      gLeft += gLeftTmp;
      gRight += gRightTmp;
    }

    //affine shift
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
    }


    template <class LocalEvaluation>
    double boundaryFlux (const LocalEvaluation& local,
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const JacobianRangeType& jacLeft,
                         /*const JacobianRangeType& jacRight,*/
                         RangeType& gLeft,
                         JacobianRangeType& gDiffLeft ) const
    {
      if ( model_.hasBoundaryValue( local ) )
      {
        if( hasDiffusion )
        {
          diffusionFlux_.boundaryFlux( local,
                                       uLeft, uRight, jacLeft,
                                       gLeft, gDiffLeft );
        }
        else
        {
          gLeft  = RangeType(0);
          gDiffLeft = JacobianRangeType(0);
        }

        RangeType gLeftTmp;
        RangeType gRightTmp;
        advFlux_.numericalFlux(local, local,
                               uLeft, uRight, jacLeft, jacLeft,
                               gLeftTmp, gRightTmp );
        gLeft += gLeftTmp;
      }
      else
      {
        model_.boundaryFlux( local,
                             uLeft, jacLeft,
                             gLeft );
        gDiffLeft = 0;
      }
    }

    template <class LocalEvaluation>
    double source( const LocalEvaluation& local,
                   const RangeType& uLeft,
                   const JacobianRangeType& jacLeft,
                   RangeType& gLeft,
                   JacobianRangeType& gDiffLeft ) const
    {
      double dtEst = 0;
      // now assemble source part depending on u (mass part)
      if ( model_.hasStiffSource() )
      {
        const double dtStiff = model_.stiffSource( local, uLeft, jacLeft, gLeft );

        dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
      }
      if ( model_.hasNonStiffSource() )
      {
        RangeType sNonStiff (0);
        const double dtNon = model_.nonStiffSource( local, uLeft, jacLeft, sNonStiff );
        gLeft += sNonStiff;

        dtEst = ( dtNon > 0 ) ? std::min( dtEst, dtNon ) : dtEst;
      }
      // return the fastest wave from source terms
      return dtEst;
    }

    template <class LocalEvaluation>
    void analyticalFlux( const LocalEvaluation& local,
                         const RangeType& uLeft,
                         const JacobianRangeType& jacLeft,
                         JacobianRangeType& gDiffLeft ) const
    {
      model_.diffusion( local, uLeft, jacLeft, gDiffLeft);
      JacobianRangeType brhsphi;
      //if( model_.hasFlux() )
      {
        model_.advection( local, uLeft, jacLeft, brhsphi);
        gDiffLeft -= brhsphi;
      }
    }

    template <class LocalEvaluation>
    double analyticalFluxAndSource( const LocalEvaluation& local,
                                    const RangeType& uLeft,
                                    const JacobianRangeType& jacLeft,
                                    RangeType& gLeft,
                                    JacobianRangeType& gDiffLeft ) const
    {
      gLeft = 0;
      gDiffLeft = 0;
      JacobianRangeType jac(0);
      analyticalFlux( local, uLeft, jacLeft, jac );
      const double dtEst = source( local, uLeft, jacLeft, gLeft, gDiffLeft );
      gDiffLeft += jac;
      return dtEst;
    }

    template <class LocalEvaluation>
    void rhs( const LocalEvaluation& local,
              RangeType& gLeft,
              JacobianRangeType& gDiffLeft ) const
    {
      const RangeType uZero(0);
      const JacobianRangeType uJacZero(0);
      analyticalFluxAndSource( local, uZero, uJacZero, gLeft, gDiffLeft );
    }


    const LinearDiscreteModelType& linear() const
    {
      return linearModel_;
    }

  private:

    const ModelType&                 model_;
    const DiscreteFunctionSpaceType& space_;
    AdvectionFluxType                advFlux_;
    mutable DiffusionFluxType        diffusionFlux_;
    LinearDiscreteModelType          linearModel_;
    mutable std::vector< RangeType > uBndVec_;


  };


}
}

#endif
