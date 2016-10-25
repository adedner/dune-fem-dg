#ifndef DUNE_STOKESASSEMBLER_HH
#define DUNE_STOKESASSEMBLER_HH

#include <dune/fem/operator/linear/spoperator.hh>
#include "tensorhelper.hh"
#include <dune/fem/solver/timeprovider.hh>

//#include "matrixoperator.hh"
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem-dg/operator/dg/passtraits.hh>
#include <dune/fem-dg/assemble/assemblertraits.hh>

#include <dune/fem-dg/algorithm/sub/steadystate.hh>
#include <dune/fem-dg/algorithm/sub/elliptic.hh>

#include <dune/fem-dg/algorithm/sub/container.hh>


namespace Dune
{
namespace Fem
{
#define PRESSURESTABILIZATION 0


  struct MatrixParameterNoPreconditioner
      : public Dune::Fem::MatrixParameter
  {
    typedef Dune::Fem::MatrixParameter BaseType;

    MatrixParameterNoPreconditioner( const std::string keyPrefix = "istl." )
      : keyPrefix_( keyPrefix )
    {}

    virtual double overflowFraction () const
    {
      return Fem::Parameter::getValue< double >( keyPrefix_ + "matrix.overflowfraction", 1.0 );
    }

    virtual int numIterations () const
    {
      return Fem::Parameter::getValue< int >( keyPrefix_ + "preconditioning.iterations", 5 );
    }

    virtual double relaxation () const
    {
      return Fem::Parameter::getValue< double >( keyPrefix_ + "preconditioning.relaxation", 1.1 );
    }

    virtual int method () const
    {
      return 0;
    }

    virtual std::string preconditionName() const
    {
      return "None";
    }

   private:
    std::string keyPrefix_;

  };

  template< int row, int col >
  struct OperatorPart;

  template<>
  struct OperatorPart<0,0>
  {
    template< class R, class C >
    using type = SparseRowLinearOperator<R,C>;
  };
  template<>
  struct OperatorPart<0,1>
  {
    template< class R, class C >
    using type = SparseRowLinearOperator<R,C>;
  };
  template<>
  struct OperatorPart<1,0>
  {
    template< class R, class C >
    using type = SparseRowLinearOperator<R,C>;
  };
  template<>
  struct OperatorPart<1,1>
  {
    template< class R, class C >
    using type = SparseRowLinearOperator<R,C>;
  };



  //template< template<class,class> class MatrixImp, class... DiscreteFunctions >
  //struct SubEllipticContainer
  //  : public TwoArgContainer< LinOperatorSelect< SubEllipticContainerItem, MatrixImp >::template Object,
  //                            SubSteadyStateContainerItem, std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > >
  //{
  //  typedef TwoArgContainer< LinOperatorSelect< SubEllipticContainerItem, MatrixImp >::template Object,
  //                           SubSteadyStateContainerItem, std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > > BaseType;

  //public:
  //  using BaseType::operator();

  //  // constructor: do not touch/delegate everything
  //  template< class ... Args>
  //  SubEllipticContainer( Args&&... args )
  //  : BaseType( args... )
  //  {}

  //};


  //template< template<class,class> class MatrixImp, class UDiscreteFunctionImp, class PDiscreteFunctionImp >
  template< template<class,class> class MatrixImp, class UDiscreteFunctionImp,class PDiscreteFunctionImp >
  class UzawaContainer
    //: public TwoArgContainer< LinOperatorSelect< SubEllipticContainerItem, MatrixImp >::template Object,
    //                          SubSteadyStateContainerItem, std::tuple< DiscreteFunctions... >, std::tuple< DiscreteFunctions... > >
    : public SubEllipticContainer< MatrixImp, UDiscreteFunctionImp, PDiscreteFunctionImp >
  {

    typedef SubEllipticContainer< MatrixImp, UDiscreteFunctionImp, PDiscreteFunctionImp > BaseType;
  public:

    using BaseType::operator();

    // constructor: do not touch/delegate everything
    template< class ... Args>
    UzawaContainer( Args&&... args )
    : BaseType( args... )
    {}

  };

  //template <template< class, class > class MatrixContainerImp, class UDiscreteFunctionImp,class PDiscreteFunctionImp >
  //class UzawaContainer
  //  : public SubEllipticContainer< MatrixContainerImp, UDiscreteFunctionImp, PDiscreteFunctionImp >
  //{

  //  typedef SubEllipticContainer< MatrixContainerImp, UDiscreteFunctionImp, PDiscreteFunctionImp > BaseType;
  //public:

  //  using BaseType::operator();

  //  template< class SameObject >
  //  UzawaContainer( SameObject& obj, const std::string name = "" )
  //  : BaseType( obj, name )
  //  {}

  //  UzawaContainer( const std::string name = "" )
  //  : BaseType( name )
  //  {}

  //};

  //template <class UDiscreteFunctionImp,class PDiscreteFunctionImp/*, class EllipticContainerImp*/  >
  //class UzawaContainer
  //{
  //  typedef UDiscreteFunctionImp                                            UDiscreteFunctionType;
  //  typedef PDiscreteFunctionImp                                            PDiscreteFunctionType;

  //  typedef typename UDiscreteFunctionType::DiscreteFunctionSpaceType       UDiscreteFunctionSpaceType;
  //  typedef typename PDiscreteFunctionType::DiscreteFunctionSpaceType       PDiscreteFunctionSpaceType;

  //  template< class DomainDF, class RangeDF >
  //  using StokesMatrixContainer = Dune::Fem::SparseRowLinearOperator< DomainDF, RangeDF >;

  //  typedef StokesMatrixContainer<UDiscreteFunctionType,UDiscreteFunctionType> UUMatrixType;
  //  typedef StokesMatrixContainer<UDiscreteFunctionType,PDiscreteFunctionType> UPMatrixType;
  //  typedef StokesMatrixContainer<PDiscreteFunctionType,UDiscreteFunctionType> PUMatrixType;
  //  typedef StokesMatrixContainer<PDiscreteFunctionType,PDiscreteFunctionType> PPMatrixType;

  //  //short cut for tuple containing shared_ptr
  //  template< class... Args >
  //  using shared_tuple = std::tuple< std::shared_ptr<Args>... >;

  //  typedef MatrixParameterNoPreconditioner                                 MatrixParameterType;
  //  //typedef std::tuple< shared_tuple< UUMatrixType, UPMatrixType >,
  //  //                    shared_tuple< PUMatrixType, PPMatrixType > >        MatrixTupleType;


  //  //static_assert( std::is_same< typename UDiscreteFunctionSpaceType::GridType,
  //  //                             typename PDiscreteFunctionSpaceType::GridType >::value,
  //  //               "GridTypes must coincide for both discrete functions!" );

  //  //static_assert( std::is_same< typename UDiscreteFunctionSpaceType::GridPartType,
  //  //                             typename PDiscreteFunctionSpaceType::GridPartType >::value,
  //  //               "GridPartTypes must coincide for both discrete functions!" );

  //  //typedef shared_tuple< UDiscreteFunctionType, PDiscreteFunctionType >    DiscreteFunctionTupleType;

  //  //typedef EllipticContainerImp                                            UContainerType;
  //  //typedef SubSteadyStateContainerItem< UDiscreteFunctionType >       UContainerType;
  //  //typedef SubSteadyStateContainerItem< PDiscreteFunctionType >       PContainerType;

  //  template< class Arg >
  //  using Mat = SubEllipticContainerItem< Arg >;

  //  typedef SubEllipticContainer< UDiscreteFunctionType, Mat<UUMatrixType> >       UContainerType;
  //  typedef SubEllipticContainer< PDiscreteFunctionType, Mat<PPMatrixType> >       PContainerType;

  //  typedef std::tuple< shared_tuple< Mat< UUMatrixType >, Mat< UPMatrixType > >,
  //                      shared_tuple< Mat< PUMatrixType >, Mat< PPMatrixType > > > MatrixContainerType;

  //  typedef shared_tuple< UContainerType, PContainerType >                  ContainerList;


  //public:
  //  template< int i, int j >
  //  using Matrix = typename std::tuple_element< i, typename std::tuple_element< j, MatrixContainerType >::type >::type::element_type;
  //  template< int i >
  //  using ContainerItem = typename std::tuple_element< i, ContainerList >::type::element_type;

  //  template< class GridImp >
  //  UzawaContainer( GridImp& grid, const std::string name = "" )
  //  : item_( std::make_shared< UContainerType >( grid ), std::make_shared< PContainerType >( grid ) ),
  //    matrix_( std::make_tuple( std::make_shared< Matrix<0,0> >( (*this)(_0), (*this)(_0), name + "uu" ),
  //                              std::make_shared< Matrix<1,0> >( (*this)(_0), (*this)(_1), name + "up"/*, MatrixParameterType()*/ ) ),
  //             std::make_tuple( std::make_shared< Matrix<0,1> >( (*this)(_1), (*this)(_0), name + "pu"/*, MatrixParameterType()*/ ),
  //                              std::make_shared< Matrix<1,1> >( (*this)(_1), (*this)(_1), name + "pp"/*, MatrixParameterType()*/ ) ) )
  //  {}

  //  template< unsigned long int i >
  //  std::shared_ptr< ContainerItem<i> > operator() ( std::integral_constant<unsigned long int, i> index )
  //  {
  //    return (*std::get<i>( item_ ))(index);
  //  }

  //  template< unsigned long int i, unsigned long int j >
  //  std::shared_ptr< Matrix<i,j> > operator() ( std::integral_constant<unsigned long int, i> row,
  //                                              std::integral_constant<unsigned long int, j> col  )
  //  {
  //    return (*std::get<i>( std::get<j>( matrix_ ) ))(row,col);
  //  }

  //  template< unsigned long int i >
  //  std::shared_ptr< ContainerItem<i> > operator[] ( std::integral_constant<unsigned long int, i> index )
  //  {
  //    return std::get<i>( item_ );
  //  }


  //  // sub Container
  //  template< unsigned long int i... >
  //  std::shared_ptr< ContainerItem<i>... > operator() ( std::tuple< std::integral_constant<unsigned long int, i>... > index )
  //  {
  //    return (*std::get<i>( item_ ))(index);
  //  }

  //  //// sub Container
  //  //template< unsigned long int i..., unsigned long int j... >
  //  //std::shared_ptr< Matrix<i,j>... > operator() ( std::tuple< std::integral_constant<unsigned long int, i>... > row,
  //  //                                               std::tuple< std::integral_constant<unsigned long int, j>... > col  )
  //  //{
  //  //  return (*std::get<i>( std::get<j>( matrix_ ) ))(row,col);
  //  //}



  //  //template< unsigned long int i, unsigned long int j >
  //  //std::shared_ptr< Matrix<i,j> > sub() ( std::integral_constant<unsigned long int, i> row,
  //  //                                       std::integral_constant<unsigned long int, j> col  )
  //  //{
  //  //  return std::get<i>( std::get<j>( matrix_ ) );
  //  //}

  //private:
  //  ContainerList             item_;
  //  MatrixContainerType       matrix_;

  //};




  /**
   * \brief Assembles the Stokes problem.
   *
   * \ingroup AssemblyOperator
   */
  template <class CombAssTraits, class OpTraits >
  class StokesAssembler
  {
  public:
    // ( ------------- )
    // | (0,0) | (0,1) |
    // | ------------- |
    // | (1,0) | (1,1) |
    // ( ------------- )
    //
    typedef typename CombAssTraits::template DomainDiscreteFunction<0,0> DiscreteFunctionType;
    typedef typename CombAssTraits::template RangeDiscreteFunction<1,1>  DiscretePressureFunctionType;

    // some basic tests...
    static_assert( std::is_same< typename CombAssTraits::template RangeDiscreteFunction<1,0>,
                                 typename CombAssTraits::template DomainDiscreteFunction<0,1> >::value,
                   "discrete function spaces does not fit." );
    static_assert( std::is_same< typename CombAssTraits::template RangeDiscreteFunction<0,1>,
                                 typename CombAssTraits::template DomainDiscreteFunction<1,0> >::value,
                   "discrete function spaces does not fit." );


    typedef typename CombAssTraits::template MatrixContainer< DiscreteFunctionType, DiscreteFunctionType >  EllMatrixContainerType;

    template< class DomainDF, class RangeDF >
    using Mat = typename CombAssTraits::template MatrixContainer< DomainDF, RangeDF >;

    typedef UzawaContainer< Mat, DiscreteFunctionType, DiscretePressureFunctionType >
                                                                         ContainerType;

    typedef OpTraits                                                     OperatorTraits;
    typedef typename OperatorTraits::ModelType                           ModelType;
    typedef typename ModelType::ProblemType                              ProblemType;

    //! type of discrete function spaceS
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType     DiscreteFunctionSpaceType;

    typedef typename DiscretePressureFunctionType::DiscreteFunctionSpaceType DiscretePressureSpaceType;
    //! type of problem

    typedef typename DiscreteFunctionSpaceType::RangeFieldType           RangeFieldType;


  protected:

    //! type of the base function sets
    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType     BaseFunctionSetType;
    typedef typename DiscretePressureSpaceType::BasisFunctionSetType     PressureBaseFunctionSetType;

    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType::GridPartType             GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType::GridType                 GridType;

    //! polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder };

    //! The grid's dimension
    enum { dimension = GridType::dimension };

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::DomainType                DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType                 RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType         JacobianRangeType;

    typedef typename DiscretePressureSpaceType::RangeType                 PressureRangeType;
    typedef typename DiscretePressureSpaceType::JacobianRangeType         PressureJacobianRangeType;


    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionType                       IntersectionType;
    typedef typename GridType::template Codim<0>::Entity                  EntityType ;

  public:
    typedef typename OperatorTraits::VolumeQuadratureType                 VolumeQuadratureType;
    typedef typename OperatorTraits::FaceQuadratureType                   FaceQuadratureType;
  protected:

    typedef typename ContainerType::template Item2<0,1>::MatrixType                  PressureGradMatType;
    typedef typename ContainerType::template Item2<1,0>::MatrixType                  PressureDivMatType;
    typedef typename ContainerType::template Item2<1,1>::MatrixType                  PressureStabMatType;

    typedef typename PressureGradMatType::LocalMatrixType                 LocalPressureGradMatType;
    typedef typename PressureDivMatType::LocalMatrixType                  LocalPressureDivMatType;
    typedef typename PressureStabMatType::LocalMatrixType                 LocalPressureStabMatType;

    typedef FieldMatrix<double,dimension,dimension>                       JacobianInverseType;

  public:

    //Constructor
    template< class ContainerImp >
    StokesAssembler( std::shared_ptr< ContainerImp > cont,
                     const ModelType& model,
                     double d11=1.,
                     double d12=1.) :
      //container_( cont ),
      spc_( (*cont)(_0)->solution()->space() ),
      pressurespc_( (*cont)(_1)->solution()->space() ),
      problem_( model.problem() ),
      veloRhs_(  (*cont)(_0)->rhs() ),
      pressureRhs_( (*cont)(_1)->rhs() ),
      volumeQuadOrd_( 2*spc_.order()+1),
      faceQuadOrd_( 2*spc_.order()+1 ),
      pressureGradMatrix_( (*cont)(_0,_1)->matrix() ),//PGM
      pressureDivMatrix_( (*cont)(_1,_0)->matrix() ),//PDM
      pressureStabMatrix_( (*cont)(_1,_1)->matrix() ),//PSM
      d11_(d11),
      d12_(d12),
      time_(0)
    {}

    void divergence(const JacobianRangeType& du, PressureRangeType& divu) const
    {
      divu=0.;
      for( int i=0;i<dimension;++i)
        divu+=du[i][i];
    }

    void setTime ( double time )
    {
      time_ = time;
    }

    void assemble()
    {

      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscretePressureSpaceType,DiscreteFunctionSpaceType> PgStencilType;
      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscretePressureSpaceType> PdStencilType;
      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscretePressureSpaceType,DiscretePressureSpaceType> PsStencilType;

      PgStencilType pgstencil( pressurespc_, spc_ );
      PdStencilType pdstencil( spc_, pressurespc_ );
      PsStencilType psstencil( pressurespc_, pressurespc_ );

      pressureGradMatrix_->reserve( pgstencil );
      pressureGradMatrix_->clear();
      pressureDivMatrix_->reserve( pdstencil );
      pressureDivMatrix_->clear();
      pressureStabMatrix_->reserve( psstencil );
      pressureStabMatrix_->clear();
      pressureRhs_->clear();

      for( const auto& en : elements( spc_.gridPart() ) )
        assembleLocal( en );

      //important: finish build process!
      pressureGradMatrix_->communicate();
      pressureDivMatrix_->communicate();
      pressureStabMatrix_->communicate();
#define SYMMCHECK 0
#if SYMMCHECK
      std::cout << "symcheck" << std::endl;
      int size=spc_.size();
      int pressuresize=pressurespc_.size();

      for(int i=0; i<size; ++i)
      {
        for(int j=0;j<pressuresize; ++j)
        {
          double error=0.0;
          error=fabs(pressureDivMatrix_->matrix()(j,i) - pressureGradMatrix_->matrix()(i,j));
          if(error > 1e-10)
          {
            std::cout<<"Wrong at:"<<i<<" , "<< j<<"error="<<error<<"\n";
            //abort();
          }
        }
      }
      std::cout << std::endl;

      //for(int i=0; i<size; ++i)
      //{
      //  for(int j=0;j<pressuresize; ++j)
      //  {
      //    std::cout.width(5);
      //    std::cout << std::setprecision(2) << pressureDivMatrix_->matrix()(j,i) << " ";
      //  }
      //  std::cout << std::endl;
      //}
      //std::cout << std::endl;
      //std::cout << std::endl;

      //for(int i=0; i<size; ++i)
      //{
      //  for(int j=0;j<pressuresize; ++j)
      //  {
      //    std::cout.width(5);
      //    std::cout << std::setprecision(2) << pressureGradMatrix_->matrix()(i,j) << " ";
      //  }
      //  std::cout << std::endl;
      //}

#endif
      //pressureDivMatrix_->matrix().print(std::cout);
      //std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n";
      //pressureGradMatrix_->matrix().print(std::cout);
      //std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n";
      //pressureStabMatrix_->matrix().print(std::cout);
    }


    template<class EntityType>
    void assembleLocal (const EntityType& en) const
    {

      GridPartType& gridPart = spc_.gridPart();

      const auto& bsetEn         = spc_.basisFunctionSet(en);
      const auto& pressurebsetEn = pressurespc_.basisFunctionSet(en);
      const size_t numBaseFunctions=bsetEn.size();
      const size_t numPBaseFunctions=pressurebsetEn.size();

      const auto& geo = en.geometry();

      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      std::vector<RangeType> u(numBaseFunctions);
      std::vector<JacobianRangeType> du(numBaseFunctions);

      std::vector<PressureRangeType> p(numPBaseFunctions);
      std::vector<PressureJacobianRangeType> dp(numPBaseFunctions);


      LocalPressureGradMatType enPGrad = pressureGradMatrix_->localMatrix(en,en);
      LocalPressureDivMatType  enPDiv  = pressureDivMatrix_->localMatrix(en,en);
#if PRESSURESTABILIZATION
      LocalPressureStabMatType enPStab = pressureStabMatrix_->localMatrix(en,en);
#endif
      JacobianInverseType inv;

      for( const auto qp : volQuad )
      {
        const auto& x = qp.position();
        inv = geo.jacobianInverseTransposed(x);

        bsetEn.evaluateAll(qp,u);
        bsetEn.jacobianAll(qp,du);

        pressurebsetEn.evaluateAll(qp,p);
        pressurebsetEn.jacobianAll(qp,dp);

        double weight = qp.weight()*geo.integrationElement(x);

        for(unsigned int k = 0;k < numBaseFunctions ;++k)
        {
          PressureRangeType divu(0.0);
          divergence(du[k],divu);

          // u = u, v
          // p = p, q
          //pressureGradientMatrix , pressureDivMatrix
          for(unsigned int n = 0; n < numPBaseFunctions ; ++n)
          {
            //eval (-p_n*div u_j)
            double PGM=p[n]*divu*weight;
            PGM*=-1;
            enPGrad.add(k,n,PGM);

            //eval -(-u_j*grad p_n)
            double PDM =(u[k]*dp[n][0])*weight;
            enPDiv.add(n,k,PDM);
          }
        }
      }

      typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingFaceQuadratureType;

      for (const auto& intersection : intersections(gridPart, en) )
      {

        if(intersection.neighbor())
        {
          const EntityType& nb = intersection.outside();
          if(intersection.conforming() )
          {
            typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, true > IntersectionQuadratureType;

            IntersectionQuadratureType faceQuad(gridPart,intersection, faceQuadOrd_);
            applyNeighbor(intersection, en, nb, faceQuad, enPGrad, enPDiv
#if PRESSURESTABILIZATION
                , enPStab
#endif
                );
          }
          else
          {
            typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, false > IntersectionQuadratureType;
            IntersectionQuadratureType faceQuad(gridPart,intersection, faceQuadOrd_);
            applyNeighbor(intersection, en, nb, faceQuad, enPGrad, enPDiv
#if PRESSURESTABILIZATION
                , enPStab
#endif
                );
          }
        }
        else if(intersection.boundary())
        {
          if(intersection.conforming() )
          {
            FaceQuadratureType faceQuadInner(gridPart,intersection, faceQuadOrd_, FaceQuadratureType::INSIDE);
            applyBoundary(intersection, en, faceQuadInner, enPGrad, enPDiv
#if PRESSURESTABILIZATION
                , enPStab
#endif
                );
          }
          else
          {
            NonConformingFaceQuadratureType   faceQuadInner(gridPart, intersection, faceQuadOrd_, NonConformingFaceQuadratureType::INSIDE);
            applyBoundary(intersection, en, faceQuadInner, enPGrad, enPDiv
#if PRESSURESTABILIZATION
                , enPStab
#endif
                );

          }
        }
      }
    }

  private:
    template<class Quad,class Entity>
    void applyNeighbor(const IntersectionType &edge,
                       const Entity &en,
                       const Entity &nb,
                       const Quad &faceQuad,
                       LocalPressureGradMatType&  enPGrad,
                       LocalPressureDivMatType&   enPDiv
#if PRESSURESTABILIZATION
                       , LocalPressureStabMatType&  enPStab
#endif
                       ) const
    {
      const auto& faceQuadInner  = faceQuad.inside();
      const auto& faceQuadOuter  = faceQuad.outside();

      const auto& bsetEn         = spc_.basisFunctionSet(en);
      const auto& bsetNb         = spc_.basisFunctionSet(nb);
      const auto& pressurebsetEn = pressurespc_.basisFunctionSet(en);
      const auto& pressurebsetNb = pressurespc_.basisFunctionSet(nb);
      const size_t numBaseFunctionsEn=bsetEn.size();
      const size_t numBaseFunctionsNb=bsetNb.size();
      const size_t numPBaseFunctionsEn=pressurebsetEn.size();
      const size_t numPBaseFunctionsNb=pressurebsetNb.size();

      std::vector<RangeType> uEn(numBaseFunctionsEn);
      std::vector<JacobianRangeType> duEn(numBaseFunctionsEn);

      std::vector<PressureRangeType> pEn(numPBaseFunctionsEn);
      std::vector<PressureJacobianRangeType> dpEn(numPBaseFunctionsEn);

      std::vector<RangeType> uNb(numBaseFunctionsNb);
      std::vector<JacobianRangeType> duNb(numBaseFunctionsNb);

      std::vector<PressureRangeType> pNb(numPBaseFunctionsNb);
      std::vector<PressureJacobianRangeType> dpNb(numPBaseFunctionsNb);

      PressureRangeType uNormal(0.);
      RangeType pNormal(0.);

      LocalPressureGradMatType  nbPressGrad=pressureGradMatrix_->localMatrix(nb,en);
      LocalPressureDivMatType   nbPressDiv=pressureDivMatrix_->localMatrix(nb,en);

#if PRESSURESTABILIZATION
      LocalPressureStabMatType  nbPressStab=pressureStabMatrix_->localMatrix(nb,en);
#endif
      //for( const auto qpp : faceQuad )
      //{
      //  const auto& qp = qpp.first;
      //  const auto& qpOut = qpp.second;
      for( const auto qp : faceQuadInner )
      {
        bsetEn.evaluateAll(qp,uEn);
        bsetNb.evaluateAll(faceQuadOuter[qp.index()],uNb);
        //bsetNb.evaluateAll(qpOut,uNb);
        pressurebsetEn.evaluateAll(qp,pEn);
        pressurebsetNb.evaluateAll(faceQuadOuter[qp.index()],pNb);
        //pressurebsetNb.evaluateAll(qpOut,pNb);

        DomainType normal=edge.integrationOuterNormal(qp.localPosition());

        double weight=qp.weight();

        for(unsigned int j=0; j< numBaseFunctionsEn; ++j)
        {
          for(unsigned int n = 0; n < numPBaseFunctionsEn ; ++n)
          {
            // ******************************************************************************
            // u+.p+*n+
            // ******************************************************************************
            pNormal=normal;
            pNormal*=pEn[n][0];
            double PGM_en=(uEn[j]*pNormal*weight);
            PGM_en*=0.5;
            enPGrad.add(j,n,PGM_en);

            // ******************************************************************************
            // -p+*u+.n+
            // ******************************************************************************
            uNormal[0]=uEn[j]*normal;
            double PDM_en=(pEn[n]*uNormal)*weight;
            PDM_en*=-0.5;
            enPDiv.add(n,j,PDM_en);

#if PRESSURESTABILIZATION
            if(j==0)
            {
              // ******************************************************************************
              // -p+*[p]*n+ * d11
              // ******************************************************************************
              for(unsigned int m=0;m<numPBaseFunctionsEn;++m)
              {
                double PSM_nb=pEn[m]*pEn[n]*weight*d11_;
                nbPressStab.add(m,n,PSM_nb);
              }
              for(unsigned int m=0;m<numPBaseFunctionsNb;++m)
              {
                double PSM_en=pEn[n]*pNb[m]*weight*d11_;
                enPStab.add(m,n,PSM_en);
              }

            }
#endif
          }
        }

        for(unsigned int j=0; j< numBaseFunctionsEn; ++j)
        {
          for(unsigned int n = 0; n < numPBaseFunctionsNb ; ++n)
          {
            // ******************************************************************************
            // -*u+.p-*n+
            // ******************************************************************************
            pNormal=normal;
            pNormal*=pNb[n][0];
            double PGM_nb=(uEn[j]*pNormal*weight);
            PGM_nb*=0.5;
            nbPressGrad.add(j,n,PGM_nb);
          }
        }

        for(unsigned int j=0; j< numBaseFunctionsNb; ++j)
        {
          for(unsigned int n = 0; n < numPBaseFunctionsEn ; ++n)
          {
            // ******************************************************************************
            // -p+*u-.n+
            // ******************************************************************************

            uNormal[0]=uNb[j]*normal;
            double PDM_nb =( pEn[n]*uNormal[0])*weight;
            PDM_nb*=-0.5;
            nbPressDiv.add(n,j,PDM_nb);

#if PRESSURESTABILIZATION
            if(j==0)
            {
              // ******************************************************************************
              // -p+*[p]*n+ * d11
              // ******************************************************************************
              for(unsigned int m=0;m<numPBaseFunctionsEn;++m)
              {
                double PSM_nb=pEn[m]*pEn[n]*weight*d11_;
                nbPressStab.add(m,n,PSM_nb);
              }
              for(unsigned int m=0;m<numPBaseFunctionsNb;++m)
              {
                double PSM_en=pEn[n]*pNb[m]*weight*d11_;
                enPStab.add(m,n,PSM_en);
              }
            }
#endif
          }
        }
      }
    }



    template<class Quad,class Entity>
    void applyBoundary(const IntersectionType &edge,
                       const Entity &en,
                       const Quad &faceQuadInner,
                       LocalPressureGradMatType&  enPGrad,
                       LocalPressureDivMatType&   enPDiv
#if PRESSURESTABILIZATION
                       , LocalPressureStabMatType&  enPStab
#endif
                       ) const
    {
      const auto& bsetEn             = spc_.basisFunctionSet(en);
      const auto& pressurebsetEn     = pressurespc_.basisFunctionSet(en);
      const size_t numBaseFunctions=bsetEn.size();
      const size_t numPBaseFunctions=pressurebsetEn.size();
      std::vector<RangeType> uEn(numBaseFunctions);
      std::vector<JacobianRangeType> duEn(numBaseFunctions);

      std::vector<PressureRangeType> pEn(numPBaseFunctions);
      std::vector<PressureJacobianRangeType> dpEn(numPBaseFunctions);


      auto localPressure = pressureRhs_->localFunction(en);
      RangeType dirichletValue(0.0);
      RangeType pNormal(0.0);

      for( const auto qp : faceQuadInner )
      {
        bsetEn.evaluateAll(qp,uEn);
        pressurebsetEn.evaluateAll(qp,pEn);

        auto normal = edge.integrationOuterNormal(qp.localPosition());

        double weight=qp.weight();
        DomainType quadInEn=edge.geometry().global(qp.localPosition());

        problem_.template get<0>().g(quadInEn,dirichletValue);
        double pressureDirichlet;

        for(unsigned int n = 0; n <(unsigned int)localPressure.numDofs() ; ++n)
        {
          // ******************************************************************************
          // p+.u_D+*n+
          // ******************************************************************************
          pressureDirichlet=normal*dirichletValue;
          pressureDirichlet*=pEn[n][0];
          pressureDirichlet*=weight;
          //pressureDirichlet*=-1.;
          localPressure[n]+=pressureDirichlet;
        }

        for(unsigned int j=0; j< numBaseFunctions; ++j)
        {
          for(unsigned int n = 0; n <numPBaseFunctions ; ++n)
          {
            // ******************************************************************************
            // u+.p+*n+
            // ******************************************************************************
            pNormal=normal;
            pNormal*=pEn[n][0];
            double PGM_en=(pNormal*uEn[j]) *weight;
            enPGrad.add(j,n,PGM_en);
          }
        }
      }
    }


  private:
    //std::shared_ptr< ContainerType >             container_;
    const DiscreteFunctionSpaceType&             spc_;
    const DiscretePressureSpaceType&             pressurespc_;
    const ProblemType&                           problem_;

    std::shared_ptr< DiscreteFunctionType >      veloRhs_;
    std::shared_ptr< DiscretePressureFunctionType > pressureRhs_;
    int volumeQuadOrd_;
    int faceQuadOrd_;

    std::shared_ptr< PressureGradMatType >       pressureGradMatrix_;
    std::shared_ptr< PressureDivMatType >        pressureDivMatrix_;
    std::shared_ptr< PressureStabMatType >       pressureStabMatrix_;

    double                                       d11_;
    double                                       d12_;
    double                                       time_;
    DomainType                                   direction_;
  };
}
}
#endif
