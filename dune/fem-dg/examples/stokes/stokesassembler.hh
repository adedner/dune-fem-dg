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

  //template <class UDiscreteFunctionImp,class PDiscreteFunctionImp, EllipticContainerImp, template< class, class> class MatrixContainerImp  >

  template <class UDiscreteFunctionImp,class PDiscreteFunctionImp, class EllipticContainerImp  >
  class UzawaContainer
  {
    typedef UDiscreteFunctionImp                                            UDiscreteFunctionType;
    typedef PDiscreteFunctionImp                                            PDiscreteFunctionType;

    typedef typename UDiscreteFunctionType::DiscreteFunctionSpaceType       UDiscreteFunctionSpaceType;
    typedef typename PDiscreteFunctionType::DiscreteFunctionSpaceType       PDiscreteFunctionSpaceType;

    template< class DomainDF, class RangeDF >
    using StokesMatrixContainer = Dune::Fem::SparseRowLinearOperator< DomainDF, RangeDF >;

    typedef StokesMatrixContainer<UDiscreteFunctionType,UDiscreteFunctionType> UUMatrixType;
    typedef StokesMatrixContainer<UDiscreteFunctionType,PDiscreteFunctionType> UPMatrixType;
    typedef StokesMatrixContainer<PDiscreteFunctionType,UDiscreteFunctionType> PUMatrixType;
    typedef StokesMatrixContainer<PDiscreteFunctionType,PDiscreteFunctionType> PPMatrixType;

    typedef EllipticContainerImp                                            UContainerType;

    typedef Fem::SubSteadyStateContainer< PDiscreteFunctionType >           PContainerType;

    //short cut for tuple containing shared_ptr
    template< class... Args >
    using shared_tuple = std::tuple< std::shared_ptr<Args>... >;

    static_assert( std::is_same< typename UDiscreteFunctionSpaceType::GridType,
                                 typename PDiscreteFunctionSpaceType::GridType >::value,
                   "GridTypes must coincide for both discrete functions!" );

    static_assert( std::is_same< typename UDiscreteFunctionSpaceType::GridPartType,
                                 typename PDiscreteFunctionSpaceType::GridPartType >::value,
                   "GridPartTypes must coincide for both discrete functions!" );


    typedef std::tuple< UDiscreteFunctionSpaceType, PDiscreteFunctionSpaceType >
                                                                            DiscreteFunctionSpaceTupleType;

    typedef shared_tuple< UDiscreteFunctionType, PDiscreteFunctionType >    DiscreteFunctionTupleType;

    typedef std::tuple< shared_tuple< UUMatrixType, UPMatrixType >,
                        shared_tuple< PUMatrixType, PPMatrixType > >        MatrixTupleType;

    typedef std::tuple< std::shared_ptr< UContainerType >, std::shared_ptr< PContainerType > >
                                                                            ContainerList;

    typedef MatrixParameterNoPreconditioner                                 MatrixParameterType;
  public:
    typedef typename UDiscreteFunctionSpaceType::GridType                   GridType;
    typedef typename UDiscreteFunctionSpaceType::GridPartType               GridPartType;

    template< int i >
    using DiscreteFunctionSpace = typename std::tuple_element< i, DiscreteFunctionSpaceTupleType >::type;
    template< int i >
    using DiscreteFunction = typename std::tuple_element< i, DiscreteFunctionTupleType >::type::element_type;
    template< int i, int j >
    using Matrix = typename std::tuple_element< i, typename std::tuple_element< j, MatrixTupleType >::type >::type::element_type;
    template< int i >
    using ContainerElement = typename std::tuple_element< i, ContainerList >::type::element_type;

    UzawaContainer( GridType& grid, const std::string name = "" )
    : gridPart_( grid ),
      list_( std::make_shared< UContainerType >( grid ), std::make_shared< PContainerType >( grid ) ),
      matrix_( std::make_tuple( std::shared_ptr< UUMatrixType >( new UUMatrixType( name + "uu",space<0>(), space<0>() ) ),
                                std::shared_ptr< UPMatrixType >( new UPMatrixType( name + "up",space<0>(), space<1>(), MatrixParameterType() ) ) ),
               std::make_tuple( std::shared_ptr< PUMatrixType >( new PUMatrixType( name + "pu",space<1>(), space<0>(), MatrixParameterType() ) ),
                                std::shared_ptr< PPMatrixType >( new PPMatrixType( name + "pp",space<1>(), space<1>(), MatrixParameterType() ) ) ) )
    {}

    //container list
    template< int i >
    const ContainerElement<i>& adapter() const
    {
      assert( std::get<i>( list_ ) );
      return *std::get<i>( list_ );
    }
    template< int i >
    ContainerElement<i>& adapter()
    {
      assert( std::get<i>( list_ ) );
      return *std::get<i>( list_ );
    }

    //spaces
    template< int i >
    const DiscreteFunctionSpace<i>& space() const
    {
      return adapter<i>().space();
    }
    template< int i >
    DiscreteFunctionSpace<i>& space()
    {
      return adapter<i>().space();
    }

    //solution
    template< int i >
    std::shared_ptr< DiscreteFunction<i> > solution() const
    {
      return adapter<i>().solution();
    }
    template< int i >
    void setSolution( std::shared_ptr< DiscreteFunction<i> > solution )
    {
      adapter<i>().solution() = solution;
    }

    //rhs
    template< int i >
    std::shared_ptr< DiscreteFunction<i> > rhs() const
    {
      return adapter<i>().rhs();
    }
    template< int i >
    void setRhs( std::shared_ptr< DiscreteFunction<i> > rhs )
    {
      adapter<i>().rhs() = rhs;
    }

    //matrix
    template< int i, int j >
    std::shared_ptr< Matrix<i,j> > matrix() const
    {
      return std::get<i>( std::get<j>( matrix_ ) );
    }
    template< int i, int j >
    void setMatrix( std::shared_ptr< Matrix<i,j> > matrix )
    {
      std::get<i>( std::get<j>( matrix_ ) ) = matrix;
    }

  private:
    GridPartType              gridPart_;

    ContainerList             list_;
    MatrixTupleType           matrix_;

  };


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

    typedef Fem::SubEllipticContainer< DiscreteFunctionType, EllMatrixContainerType >
    //typedef Fem::SubEllipticContainer< DiscreteFunctionType, Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > >
                                                                         EllContainerType;

    typedef UzawaContainer< DiscreteFunctionType, DiscretePressureFunctionType, EllContainerType >
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
    typedef typename DiscreteFunctionSpaceType::IteratorType              IteratorType;

    typedef typename DiscretePressureSpaceType::RangeType                 PressureRangeType;
    typedef typename DiscretePressureSpaceType::JacobianRangeType         PressureJacobianRangeType;


    typedef typename DiscreteFunctionType::LocalFunctionType              LocalFuncType;
    typedef typename DiscretePressureFunctionType::LocalFunctionType      LocalPressureType;

    // Types extracted from the underlying grid
    typedef typename GridType::Traits::LeafIntersectionIterator           IntersectionIterator;

    typedef typename GridPartType::Traits::IndexSetType                   IndexSetType;
    typedef typename IntersectionIterator::Intersection                   IntersectionType;
    typedef typename GridType::template Codim<0>::Entity                  EntityType ;
    typedef typename GridType::template Codim<0>::Geometry                GeometryType;
    //! type of quadrature to be used

  public:
    typedef typename OperatorTraits::VolumeQuadratureType                 VolumeQuadratureType;
    typedef typename OperatorTraits::FaceQuadratureType                   FaceQuadratureType;
  protected:

    typedef typename ContainerType::template Matrix<0,1>                  PressureGradMatType;
    typedef typename ContainerType::template Matrix<1,0>                  PressureDivMatType;
    typedef typename ContainerType::template Matrix<1,1>                  PressureStabMatType;

    //typedef typename CombAssTraits::template Matrix<0,1>                  PressureGradMatType;
    //typedef typename CombAssTraits::template Matrix<1,0>                  PressureDivMatType;
    //typedef typename CombAssTraits::template Matrix<1,1>                  PressureStabMatType;

    typedef typename PressureGradMatType::LocalMatrixType                 LocalPressureGradMatType;
    typedef typename PressureDivMatType::LocalMatrixType                  LocalPressureDivMatType;
    typedef typename PressureStabMatType::LocalMatrixType                 LocalPressureStabMatType;

    typedef FieldMatrix<double,dimension,dimension>                       JacobianInverseType;

  public:

    //Constructor
    StokesAssembler( ContainerType& container,
                     const ModelType& model,
                     double d11=1.,
                     double d12=1.) :
      container_( container ),
      spc_( container_.template space<0>() ),
      pressurespc_( container_.template space<1>() ),
      problem_( model.problem() ),
      veloRhs_( container_.template rhs<0>() ),
      pressureRhs_( container_.template rhs<1>() ),
      volumeQuadOrd_( 2*spc_.order()+1),
      faceQuadOrd_( 2*spc_.order()+1 ),
      pressureGradMatrix_( container_.template matrix<0,1>() ),//PGM
      pressureDivMatrix_( container_.template matrix<1,0>() ),//PDM
      pressureStabMatrix_( container_.template matrix<1,1>() ),//PSM
      d11_(d11),
      d12_(d12)
    {}

    void divergence(const JacobianRangeType& du, PressureRangeType& divu) const
    {
      divu=0.;
      for( int i=0;i<dimension;++i)
        divu+=du[i][i];
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

      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      IteratorType end = spc_.end();
      for(IteratorType it = spc_.begin(); it != end; ++it)
      {
        assembleLocal( *it );
      }

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
      //pressureGradMatrix_->.matrix().print(std::cout);
      //std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n";
      //pressureStabMatrix_->matrix().print(std::cout);
    }


    template<class EntityType>
    void assembleLocal (const EntityType& en) const
    {

      GridPartType& gridPart = spc_.gridPart();

      const BaseFunctionSetType& bsetEn = spc_.basisFunctionSet(en);
      const PressureBaseFunctionSetType& pressurebsetEn = pressurespc_.basisFunctionSet(en);
      const size_t numBaseFunctions=bsetEn.size();
      const size_t numPBaseFunctions=pressurebsetEn.size();

      //- typedefs
      typedef typename EntityType :: Geometry Geometry;
      const Geometry& geo=en.geometry();

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

      const size_t quadNop = volQuad.nop();
      for (size_t l = 0; l < quadNop ; ++l)
      {
        const typename  VolumeQuadratureType::CoordinateType& x=volQuad.point(l);
        inv = geo.jacobianInverseTransposed(x);
        bsetEn.evaluateAll(volQuad[l],u);
        bsetEn.jacobianAll(volQuad[l],du);

        pressurebsetEn.evaluateAll(volQuad[l],p);
        pressurebsetEn.jacobianAll(volQuad[l],dp);

        double quadWeight = volQuad.weight(l)*geo.integrationElement(x);

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
            double PGM=p[n]*divu*quadWeight;
            PGM*=-1;
            enPGrad.add(k,n,PGM);

            //eval -(-u_j*grad p_n)
            double PDM =(u[k]*dp[n][0])*quadWeight;
            enPDiv.add(n,k,PDM);
          }
        }
      }

      IntersectionIterator endnit = gridPart.iend(en);
      typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingFaceQuadratureType;

      for (IntersectionIterator nit = gridPart.ibegin(en); nit != endnit; ++nit)
      {
        const IntersectionType& edge=*nit;

        if(edge.neighbor())
        {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
          const EntityType& nb= edge.outside();
#else
         //Access to the neighbor Element
          EntityPointerType neighEp=edge.outside();
          const EntityType& nb=  *neighEp;
#endif

          if(edge.conforming() )
          {
            FaceQuadratureType faceQuadInner(gridPart,edge, faceQuadOrd_, FaceQuadratureType::INSIDE);
            FaceQuadratureType faceQuadOuter(gridPart,edge, faceQuadOrd_, FaceQuadratureType::OUTSIDE);
            applyNeighbor(edge, en, nb, faceQuadInner, faceQuadOuter, enPGrad, enPDiv
#if PRESSURESTABILIZATION
                , enPStab
#endif
                );
          }
          else
          {
            NonConformingFaceQuadratureType  faceQuadOuter(gridPart, edge, faceQuadOrd_, NonConformingFaceQuadratureType::OUTSIDE);
            NonConformingFaceQuadratureType  faceQuadInner(gridPart, edge, faceQuadOrd_, NonConformingFaceQuadratureType::INSIDE);
            applyNeighbor(edge, en, nb, faceQuadInner, faceQuadOuter, enPGrad, enPDiv
#if PRESSURESTABILIZATION
                , enPStab
#endif
                );
          }
        }
        else if(edge.boundary())
        {
          if(edge.conforming() )
          {
            FaceQuadratureType faceQuadInner(gridPart,edge, faceQuadOrd_, FaceQuadratureType::INSIDE);
            applyBoundary(edge, en, faceQuadInner, enPGrad, enPDiv
#if PRESSURESTABILIZATION
                , enPStab
#endif
                );
          }
          else
          {
            NonConformingFaceQuadratureType   faceQuadInner(gridPart, edge, faceQuadOrd_, NonConformingFaceQuadratureType::INSIDE);
            applyBoundary(edge, en, faceQuadInner, enPGrad, enPDiv
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
                       const Quad &faceQuadInner,
                       const Quad &faceQuadOuter,
                       LocalPressureGradMatType&  enPGrad,
                       LocalPressureDivMatType&   enPDiv
#if PRESSURESTABILIZATION
                       , LocalPressureStabMatType&  enPStab
#endif
                       ) const
    {
      const BaseFunctionSetType& bsetEn = spc_.basisFunctionSet(en);
      const BaseFunctionSetType& bsetNb = spc_.basisFunctionSet(nb);
      const PressureBaseFunctionSetType& pressurebsetEn=pressurespc_.basisFunctionSet(en);
      const PressureBaseFunctionSetType& pressurebsetNb=pressurespc_.basisFunctionSet(nb);
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

      const int quadNop = faceQuadInner.nop();

      LocalPressureGradMatType  nbPressGrad=pressureGradMatrix_->localMatrix(nb,en);
      LocalPressureDivMatType   nbPressDiv=pressureDivMatrix_->localMatrix(nb,en);

#if PRESSURESTABILIZATION
      LocalPressureStabMatType  nbPressStab=pressureStabMatrix_->localMatrix(nb,en);
#endif
      for (int l = 0; l < quadNop ; ++l)
      {
        bsetEn.evaluateAll(faceQuadInner[l],uEn);
        bsetNb.evaluateAll(faceQuadOuter[l],uNb);
        pressurebsetEn.evaluateAll(faceQuadInner[l],pEn);
        pressurebsetNb.evaluateAll(faceQuadOuter[l],pNb);

        DomainType normal=edge.integrationOuterNormal(faceQuadInner.localPoint(l));

        double intWeight=faceQuadInner.weight(l);

        for(unsigned int j=0; j< numBaseFunctionsEn; ++j)
        {
          for(unsigned int n = 0; n < numPBaseFunctionsEn ; ++n)
          {
            // ******************************************************************************
            // u+.p+*n+
            // ******************************************************************************
            pNormal=normal;
            pNormal*=pEn[n][0];
            double PGM_en=(uEn[j]*pNormal*intWeight);
            PGM_en*=0.5;
            enPGrad.add(j,n,PGM_en);

            // ******************************************************************************
            // -p+*u+.n+
            // ******************************************************************************
            uNormal[0]=uEn[j]*normal;
            double PDM_en=(pEn[n]*uNormal)*intWeight;
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
                double PSM_nb=pEn[m]*pEn[n]*intWeight*d11_;
                nbPressStab.add(m,n,PSM_nb);
              }
              for(unsigned int m=0;m<numPBaseFunctionsNb;++m)
              {
                double PSM_en=pEn[n]*pNb[m]*intWeight*d11_;
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
            double PGM_nb=(uEn[j]*pNormal*intWeight);
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
            double PDM_nb =( pEn[n]*uNormal[0])*intWeight;
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
                double PSM_nb=pEn[m]*pEn[n]*intWeight*d11_;
                nbPressStab.add(m,n,PSM_nb);
              }
              for(unsigned int m=0;m<numPBaseFunctionsNb;++m)
              {
                double PSM_en=pEn[n]*pNb[m]*intWeight*d11_;
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
      const BaseFunctionSetType& bsetEn = spc_.basisFunctionSet(en);
      const PressureBaseFunctionSetType& pressurebsetEn=pressurespc_.basisFunctionSet(en);
      const size_t numBaseFunctions=bsetEn.size();
      const size_t numPBaseFunctions=pressurebsetEn.size();
      std::vector<RangeType> uEn(numBaseFunctions);
      std::vector<JacobianRangeType> duEn(numBaseFunctions);

      std::vector<PressureRangeType> pEn(numPBaseFunctions);
      std::vector<PressureJacobianRangeType> dpEn(numPBaseFunctions);


      LocalPressureType localPressure=pressureRhs_->localFunction(en);
      RangeType dirichletValue(0.0);
      RangeType pNormal(0.0);

      int quadNop=faceQuadInner.nop();

      for (int l = 0; l < quadNop ; ++l)
      {
        bsetEn.evaluateAll(faceQuadInner[l],uEn);
        pressurebsetEn.evaluateAll(faceQuadInner[l],pEn);

        DomainType normal=edge.integrationOuterNormal(faceQuadInner.localPoint(l));

        double intWeight=faceQuadInner.weight(l);
        DomainType quadInEn=edge.geometry().global(faceQuadInner.localPoint(l));

        problem_.template get<0>().g(quadInEn,dirichletValue);
        double pressureDirichlet;

        for(unsigned int n = 0; n <(unsigned int)localPressure.numDofs() ; ++n)
        {
          // ******************************************************************************
          // p+.u_D+*n+
          // ******************************************************************************
          pressureDirichlet=normal*dirichletValue;
          pressureDirichlet*=pEn[n][0];
          pressureDirichlet*=intWeight;
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
            double PGM_en=(pNormal*uEn[j]) *intWeight;
            enPGrad.add(j,n,PGM_en);
          }
        }
      }
    }


  private:
    ContainerType&                               container_;
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
    DomainType                                   direction_;
  };
}
}
#endif
