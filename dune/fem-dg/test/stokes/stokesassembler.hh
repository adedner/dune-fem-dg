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

namespace Dune {
#define MATRIXBUG 1

  //! implementation of the operator
  template <class DiscreteFunction,class DiscretePressureFunction, class Traits>
  class StokesAssembler
  {
  public:
    typedef typename Traits :: InitialDataType  ProblemType;
    //! type of discrete functions
    typedef DiscreteFunction DiscreteFunctionType;
    typedef DiscretePressureFunction DiscretePressureFunctionType;

    //! type of discrete function spaceS
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteSpaceType;

    typedef typename DiscretePressureFunction :: DiscreteFunctionSpaceType
    DiscretePressureSpaceType;
    //! type of problem
    //  typedef StokesProblemInterface< typename DiscreteFunctionSpaceType :: FunctionSpaceType,
    //            typename DiscretePressureSpaceType :: FunctionSpaceType >

    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
    RangeFieldType;

  protected:

    //! type of the base function sets
    typedef typename DiscreteFunctionSpaceType :: BasisFunctionSetType
    BaseFunctionSetType;
    typedef typename DiscretePressureSpaceType:: BasisFunctionSetType
    PressureBaseFunctionSetType;

    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

    //! The grid's dimension
    enum { dimension = GridType :: dimension };

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    typedef typename DiscretePressureSpaceType::RangeType PressureRangeType;
    typedef typename DiscretePressureSpaceType::JacobianRangeType PressureJacobianRangeType;


    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    typedef typename DiscretePressureFunctionType::LocalFunctionType LocalPressureType;

    // Types extracted from the underlying grid
    typedef typename GridType::Traits::LeafIntersectionIterator IntersectionIterator;

    typedef typename GridPartType::Traits::IndexSetType IndexSetType;
    typedef typename IntersectionIterator::Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Entity   EntityType ;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;
    //! type of quadrature to be used

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
#if 0
#warning USING ISTL
#if DGSCHEME // for all dg schemes including pdg (later not working)
  typedef Dune::Fem::DGMatrixTraits< DiscreteFunctionSpaceType > MyMatrixTraits;
#else
  typedef Dune::Fem::LagrangeMatrixTraits< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType, true > MyMatrixTraits;
#endif
//#else
    struct MyMatrixTraits
      : public MatrixTraits,
        public Fem :: OverloadedSparseRowMatrix
    {
      // define new type of SparseRowMatrix to be used in SparseRowMatrixObject
      typedef Dune::SparseRowMatrixExtra< double > MatrixType ;
    };
#endif

  public:
#if 0
    typedef Dune::Fem::ISTLMatrixOperator< DiscretePressurteFunctionType, DiscreteFunctionType, MyMatrixTraits >PressureGradMatType;
    typedef Dune::Fem::ISTLMatrixOperator< DiscreteFunctionType, DiscretePressureFunctionType,MyMatrixTraits > PressureDivMatType;
    typedef Dune::Fem::SparseRowMatrixOperator< DiscretePressureFunctionType,  DiscretePressureFunctionType, MyMatrixTraits >PressureStabMatType;
#else
    typedef Dune::Fem::SparseRowLinearOperator< DiscretePressureFunctionType, DiscreteFunctionType >PressureGradMatType;
    typedef Dune::Fem::SparseRowLinearOperator<DiscreteFunctionType, DiscretePressureFunctionType > PressureDivMatType;
    typedef Dune::Fem::SparseRowLinearOperator< DiscretePressureFunctionType,  DiscretePressureFunctionType >PressureStabMatType;
#endif
    typedef typename PressureGradMatType::LocalMatrixType LocalPressureGradMatType;
    typedef typename PressureDivMatType::LocalMatrixType LocalPressureDivMatType;
    typedef typename PressureStabMatType::LocalMatrixType LocalPressureStabMatType;
    //typedef typename PressureGradMatType::MatrixType  PGMType;
    //typedef typename PressureDivMatType::MatrixType  PDMType;:
    //typedef typename PressureStabMatType::MatrixType PSMType;
    typedef PressureGradMatType PGMType;
    typedef PressureDivMatType PDMType;
    typedef PressureStabMatType PSMType;

    typedef FieldMatrix<double,dimension,dimension> JacobianInverseType;
  private:
    //Class Members
    const DiscreteFunctionSpaceType& spc_;
    const ProblemType& problem_;
    const DiscretePressureSpaceType& pressurespc_;
    mutable DiscreteFunctionType veloRhs_;
    mutable DiscretePressureFunctionType pressureRhs_;
    int volumeQuadOrd_,faceQuadOrd_;

    PressureGradMatType pressureGradMatrix_;
    PressureDivMatType  pressureDivMatrix_;
    PressureStabMatType pressureStabMatrix_;

    double d11_;
    double d12_;
    DomainType direction_;
  public:

    //Constructor
    StokesAssembler( const DiscreteFunctionSpaceType& spc,
                     const DiscretePressureSpaceType& pressurespc,
                     const ProblemType& problem,
                     double d11=1.,
                     double d12=1.) :
      spc_(spc),
      problem_(problem),
      pressurespc_( pressurespc ),
      veloRhs_("VelocityRhs",spc_),
      pressureRhs_("PressureRhs",pressurespc_),
      volumeQuadOrd_( 2*spc_.order()+1),
      faceQuadOrd_( 2*spc_.order()+1 ),
      pressureGradMatrix_("pgm",pressurespc_,spc_),//PGM
      pressureDivMatrix_("pdm",spc_,pressurespc_),//PDM
      pressureStabMatrix_("psm",pressurespc_,pressurespc_),//PSM
      d11_(d11),
      d12_(d12)
    {}

    const DiscretePressureSpaceType& pressurespc()const
    {return pressurespc_;}
    const DiscreteFunctionSpaceType& spc() const
    {return spc_;}

    const PressureGradMatType& getBOP() const  { return pressureGradMatrix_; }
    const PressureDivMatType&  getBTOP() const { return pressureDivMatrix_; }
    const PressureStabMatType& getCOP() const  { return pressureStabMatrix_; }

    DiscreteFunctionType& veloRhs() const { return veloRhs_; }
    DiscretePressureFunctionType& pressureRhs() const { return pressureRhs_; }

    void divergence(const JacobianRangeType& du, PressureRangeType& divu) const
    {
      divu=0.;
      for( int i=0;i<dimension;++i)
        divu+=du[i][i];
    }

    void assemble(const ProblemType& prob )
    {

      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscretePressureSpaceType,DiscreteFunctionSpaceType> PgStencilType;
      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscretePressureSpaceType> PdStencilType;
      typedef Dune::Fem::DiagonalAndNeighborStencil<DiscretePressureSpaceType,DiscretePressureSpaceType> PsStencilType;

      PgStencilType pgstencil( pressurespc_, spc_ );
      PdStencilType pdstencil( spc_, pressurespc_ );
      PsStencilType psstencil( pressurespc_, pressurespc_ );

      pressureGradMatrix_.reserve( pgstencil );
      pressureGradMatrix_.clear();
      pressureDivMatrix_.reserve( pdstencil );
      pressureDivMatrix_.clear();
#if PRESSURESTABILIZATION
      pressureStabMatrix_.reserve( psstencil );
      pressureStabMatrix_.clear();
#endif
      pressureRhs_.clear();

      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      IteratorType end = spc_.end();
      for(IteratorType it = spc_.begin(); it != end; ++it)
      {
        assembleLocal( *it,prob );
      }

#define SYMMCHECK 0
#if SYMMCHECK
      int size=spc_.size();
      int pressuresize=pressurespc_.size();

      for(int i=0; i<size; ++i)
      {
        for(int j=0;j<pressuresize; ++j)
        {
          double error=0.0;
          error=fabs(pressureDivMatrix_.matrix()(j,i) - pressureGradMatrix_.matrix()(i,j));
          if(error > 1e-10)
          {
            std::cout<<"Wrong at:"<<i<<" , "<< j<<"error="<<error<<"\n";
            abort();
          }
        }
      }
      //      pressureDivMatrix_.matrix().print(std::cout);
      //      std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n";
      //      pressureGradMatrix_.matrix().print(std::cout);
#endif
    }


    template<class EntityType>
    void assembleLocal(const EntityType& en/*,DiscreteFunctionType& rhs*/,const ProblemType& problem) const
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


      LocalPressureGradMatType enPGrad = pressureGradMatrix_.localMatrix(en,en);
      LocalPressureDivMatType  enPDiv  = pressureDivMatrix_.localMatrix(en,en);
      LocalPressureStabMatType enPStab = pressureStabMatrix_.localMatrix(en,en);
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
            applyNeighbor(edge, en, nb, faceQuadInner, faceQuadOuter, enPGrad, enPDiv, enPStab);
          }
          else
          {
            NonConformingFaceQuadratureType  faceQuadOuter(gridPart, edge, faceQuadOrd_, NonConformingFaceQuadratureType::OUTSIDE);
            NonConformingFaceQuadratureType  faceQuadInner(gridPart, edge, faceQuadOrd_, NonConformingFaceQuadratureType::INSIDE);
            applyNeighbor(edge, en, nb, faceQuadInner, faceQuadOuter, enPGrad, enPDiv, enPStab);
          }
        }
        else if(edge.boundary())
        {
          if(edge.conforming() )
          {
            FaceQuadratureType faceQuadInner(gridPart,edge, faceQuadOrd_, FaceQuadratureType::INSIDE);
            applyBoundary(edge, en, faceQuadInner, enPGrad, enPDiv, enPStab);
          }
          else
          {
            NonConformingFaceQuadratureType   faceQuadInner(gridPart, edge, faceQuadOrd_, NonConformingFaceQuadratureType::INSIDE);
            applyBoundary(edge, en, faceQuadInner, enPGrad, enPDiv, enPStab);

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
                       LocalPressureDivMatType&   enPDiv,
                       LocalPressureStabMatType&  enPStab
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

      LocalPressureGradMatType  nbPressGrad=pressureGradMatrix_.localMatrix(nb,en);
      LocalPressureDivMatType   nbPressDiv=pressureDivMatrix_.localMatrix(nb,en);

#if PRESSURESTABILIZATION
      LocalPressureStabMatType  nbPressStab=pressureStabMatrix_.localMatrix(nb,en);
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
                       LocalPressureDivMatType&   enPDiv,
                       LocalPressureStabMatType&  enPStab
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


      LocalPressureType localPressure=pressureRhs_.localFunction(en);
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

        problem_.g(quadInEn,dirichletValue);
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
  };
}
#endif
