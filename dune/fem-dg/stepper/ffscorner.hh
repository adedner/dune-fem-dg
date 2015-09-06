#ifndef DUNE_FFSFILTER_HH
#define DUNE_FFSFILTER_HH

//- System includes
#include <vector>
#include <cassert>

//- Dune includes
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune {
  template <class GridPartType, int dim > 
  class CornerFinderBase
  {
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
    typedef typename GridPartType :: GridType GridType; 
    typedef typename GridPartType :: IndexSetType IndexSetType; 

    typedef typename GridType :: template Codim<0> :: Entity EntityType; 
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType; 

    typedef typename GridType :: ctype ctype ;
    enum { dimworld = GridType :: dimensionworld };
    typedef FieldVector<ctype, dimworld> DomainType;

    class ReferenceAndCopies 
    {
    protected:  
      mutable EntityPointerType* refCell_;
      mutable std::vector< EntityPointerType* > pointers_;
    public:  
      ReferenceAndCopies() 
        : refCell_ (0), pointers_ () 
      {}

      template <class EntityPtrType> 
      void addReferenceCell( const EntityPtrType& refCell ) 
      {
        if ( refCell_ != 0 ) 
        {
          std::cerr << "ERROR: duplicate ref cell! \n";
          abort();
        }

        //std::cout << "Pushing reference cell! \n";
        refCell_ = new EntityPointerType(refCell);
      }
      
      template <class EntityPtrType>
      void addCopyCell( const EntityPtrType& cell ) 
      {
        pointers_.push_back( new EntityPointerType(cell) );
      }
        
      void clear() 
      {
        for(size_t i=0; i<pointers_.size(); ++i) 
        {
          delete pointers_[i];
        }
        pointers_.clear();
      }

      size_t size() const { return pointers_.size(); }

      bool foundRefCell () const { return refCell_ != 0; }

      EntityPointerType& referenceCell() const 
      {
        assert( refCell_ );
        return *(refCell_);
      }
      
      EntityPointerType& copyCell(const size_t i) const 
      {
        assert( i < pointers_.size() );
        return *(pointers_[i]);
      }
    };

    typedef ReferenceAndCopies ReferenceAndCopiesType;

    const GridPartType& gridPart_;
    Dune::Fem::GridWidthProvider< GridType > gridWidth_;

    std::vector<ReferenceAndCopiesType> cellPtrVec_;
    int refRank_ ;

    double hx_;
    DomainType refPoint_ ;
    
    DomainType upperBound_ ;
    DomainType lowerBound_ ;

  public:  
    CornerFinderBase(const GridPartType& gridPart) 
    : gridPart_( gridPart ),
      gridWidth_( & gridPart_.grid() ),
      cellPtrVec_(),
      refRank_( -1 ),
      hx_(0.0),
      refPoint_(0),
      upperBound_(0.2),
      lowerBound_(0.0)
    {
      upperBound_[0] = 0.6;

      // search reference and copy cells 
      findElements(); 
    }

    void findElements() 
    {
      hx_ = 0.;
      double width = gridWidth_.gridWidth();

      DomainType upperBound = upperBound_;
      DomainType lowerBound = lowerBound_;

      // first in y direction 
      size_t dir = dimworld-1;
      size_t normDir = 0;

      bool finished = false;
      double dx = 0.0;
      size_t count = 0;
      
      while( ! finished ) 
      {
        upperBound = upperBound_;
        upperBound[dir] -= dx;
        
        for(int i=0; i<dimworld; ++i) 
          lowerBound[i] = upperBound[i] - width;

        //std::cout << "Upper bound is " << upperBound << "\n";
        //std::cout << "Lower bound is " << lowerBound << "\n";

        cellPtrVec_.push_back( ReferenceAndCopiesType() );

        ReferenceAndCopiesType& cellPtr = cellPtrVec_[count];
      
        bool foundRefCell = false;

        const IteratorType endit = gridPart_.template end<0> ();
        for(IteratorType it = gridPart_.template begin<0> ();
            it != endit; ++it)
        {
          const EntityType& en = *it;
          assert( (foundRefCell) ? ( ! isReferenceCell( en, upperBound , lowerBound ) ) : true );
          if( isReferenceCell( en, upperBound, lowerBound ) && ! foundRefCell ) 
          {
            // only interior cells can be reference cells 
            if( it->partitionType() == InteriorEntity ) 
            {
              cellPtr.addReferenceCell( it );
              foundRefCell = true; 
              refRank_ = gridPart_.grid().comm().rank();
            }
          }
        }

        // in serial runs reference cell should have been found 
        assert( ( gridPart_.grid().comm().size() == 1 ) ? foundRefCell : true ); 

        // search for copy cells 
        {
          cellPtr.clear();
          
          for(IteratorType it = gridPart_.template begin<0> ();
              it != endit; ++it)
          {
            const EntityType& en = *it;
            if( isCopyCell( en, width, dir, normDir, upperBound, lowerBound ) ) 
            {
              cellPtr.addCopyCell( it );
            }
            // if not already found a reference cell use this one 
            else if( ! foundRefCell ) 
            {
              cellPtr.addReferenceCell( it );
              foundRefCell = true; 
            }
          }
        }

        // if (!foundRefCell || pointers_.size()<=5) 
        // if( Parameter :: verbose() ) 
        {
          std::cout << " Reference cell: " << foundRefCell << std::endl
                    << " nof copy cells: " << cellPtr.size()  
                    << std::endl;
          // abort();
        }

        // make sure only one rank was found 
#ifndef NDEBUG 
        {
          int val = ( refRank_ != -1 ) ? 1 : 0; 
          val = gridPart_.grid().comm().sum( val );
          assert( val == 1 );
        }
#endif 

        {
          refRank_ = gridPart_.grid().comm().max( refRank_ );
        }

        // in 2d only one loop needed 
        if( dimworld < 3 ) break ; 

        dx += width; 
        if( std::abs(dx - upperBound_[dir]) < 1e-8 ) 
        {
          --dir;
          dx = 0.0;
          ++normDir;
        }
        
        if( dir == 0 ) break; 

        ++count;
      }
    }

    bool isReferenceCell( const EntityType& en ,
                          const DomainType& upperBound,
                          const DomainType& lowerBound) 
    {
      typedef typename EntityType :: Geometry Geometry;
      const Geometry& geo = en.geometry();
      const DomainType baryCenter = geo.center();

      FieldVector<ctype, dimworld-1> mid(0.5);

      // first check barycenter is above upper bound
      for(int i=0; i<dimworld; ++i) 
      {
        if( baryCenter[i] >= upperBound_[i] ) return false;
      }
      
      // first check barycenter is above upper bound
      for(int i=0; i<dimworld; ++i) 
      {
        if( baryCenter[i] >= upperBound[i] ) return false;
      }
      
      if( geo.type().isCube() ) 
      {
        // then check barycenter is below lower bound
        for(int i=0; i<dimworld; ++i) 
        {
          if( baryCenter[i] <= lowerBound[i] ) return false;
        }
      }
      
      bool foundId4 = false;
      bool rightBound = false;
      bool nbBaryUpper  = false;

      typedef typename GridPartType :: IntersectionIteratorType  IntersectionIteratorType;
      IntersectionIteratorType endnit = gridPart_.iend( en );
      for(IntersectionIteratorType nit = gridPart_.ibegin( en ); 
          nit != endnit; ++nit) 
      {
        typedef typename IntersectionIteratorType::Intersection IntersectionType;
        const IntersectionType &inter = *nit;
        if( inter.boundary () )
        {
          if( inter.boundaryId() >= 3 ) foundId4 = true ;
          if( foundId4 ) 
          {
            DomainType normal = inter.outerNormal( mid );
            if( normal[0] > 0 ) 
            {
              typedef typename IntersectionType :: Geometry LocalGeometryType;
              const LocalGeometryType& interGeo = inter.geometry();   
              std::cout << "ELEMENT FOUND: "
                << interGeo.corners() << std::endl;
              for (int p=0; p<interGeo.corners(); ++p) 
              {
                std::cout << p << " ";
                for(int i=0; i<dimworld; ++i) 
                {
                  std::cout << std::abs(interGeo.corner(p)[i]-upperBound[i]) << " ";
                }
                std::cout << std::endl;

                if ( (interGeo.corner(p) - upperBound).two_norm()  < 1e-8 )
                {
                  rightBound = true;
                  // get grid width
                  hx_ = std::pow(geo.volume(), ((double) 1./dimworld) );
                }
              }
            }
          }
        }
        nbBaryUpper = true;
      }
      return (foundId4 && rightBound && nbBaryUpper);
    }

    bool isCopyCell(const EntityType& en, 
                    const double width,
                    const int dir,
                    const int normDir,
                    const DomainType& upperBound,
                    const DomainType& lowerBound) 
    { 
      const DomainType baryCenter = en.geometry().center();

      if( dimworld > 2 ) 
      {
        if( baryCenter[dir] >= upperBound[dir] ) return false;
        if( baryCenter[dir] <= lowerBound[dir] ) return false;
      }

      DomainType x(baryCenter);
      x -= upperBound; 

      if( x[0] >= 0.0 ) 
      {
        if( dimworld > 2 ) x[dir] = 0;
        x[0] -= hx_;

        const double rad = 3.5*hx_; // Struct. < sqrt(13/2)=2.55
        const double r = x.two_norm();
        if (r<rad) return true;
      }
      return false;
    }

    template <class RangeType>
    void evalRef(const RangeType& valRef,
                 const RangeType& valU,
                 RangeType& ret) 
    {
      const int e = dimworld+1;
      const double gamma = 1.4;
      RangeType primU(valU),primRet;
      // convert to prim and compute eval |v|^2
      double kinU = 0.0;
      double rhokinRef = 0.0;
      for (int i=0;i<dimworld;++i)
      {
        primU[i+1] /= valU[0];
        kinU += primU[i+1]*primU[i+1];
        rhokinRef += (valRef[i+1]*valRef[i+1]);
      }
      kinU *= 0.5;
      rhokinRef *= 0.5;
      primU[e] = (valU[e]-kinU*valU[0])*(gamma-1.);

      double pRef = (valRef[e]-rhokinRef/valRef[0])*(gamma-1.);
      double SRef = log(pRef*pow(valRef[0],-gamma));
      double HRef = (valRef[e]+pRef)/valRef[0];
      
      primRet[e] = primU[e];
      primRet[0] = exp((log(primRet[e])-SRef)/gamma);
      double a =  sqrt(
                   (HRef-gamma/(gamma-1.)*primRet[e]/primRet[0])/kinU 
                  );
      assert(a==a);

      kinU = 0.0;
      for(int i=0; i<dimworld; ++i) 
      {
        primRet[i+1] = primU[i+1]*a;
        kinU += primRet[i+1]*primRet[i+1];
      }
      kinU *= 0.5;
      
      ret = primRet;
      for (int i=0;i<dimworld;++i)
        ret[i+1] *= primRet[0];
      ret[e] = kinU*ret[0] + primRet[e]/(gamma-1.);
      assert(ret[0]==ret[0] && ret[1]==ret[1]);
    }

    template <class DiscreteFunctionType>
    void setFunction(DiscreteFunctionType& func)
    {
      std::cout << "ffscorner.hh: setFunction" << std::endl;
      typedef Dune::Fem::CachingQuadrature<GridPartType,0> QuadratureType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BaseFunctionSetType ; 
      typedef typename DiscreteFunctionType :: LocalFunctionType
        LocalFunctionType;
      typedef typename LocalFunctionType :: RangeType RangeType;

      RangeType ret,phi;

      // evaluate reference cell 
      if( cellPtrVec_.size() != 1 )
      {
        std::cerr << "More than one cell pointer vector in CornerFinder! \n";
        abort();
      }

      //for(size_t j=0; j<cellPtrVec_.size(); ++j) 
      size_t j = 0;

      typedef typename GridType :: Traits :: CollectiveCommunication  Communicator; 
      const Communicator& comm = gridPart_.grid().comm();

      typedef Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > TemporaryLocalFunctionType;
      TemporaryLocalFunctionType referenceSolution( func.space() );

      // on the rank with the reference cell 
      // if( refRank_ == comm.rank() )
      ReferenceAndCopiesType& cellPtr = cellPtrVec_[j];
      {

        assert( cellPtr.foundRefCell () );

        EntityPointerType& refCell = cellPtr.referenceCell();
        assert( refCell->isLeaf() );
        
        // get local function 
        LocalFunctionType lfRef = func.localFunction( *refCell );
        referenceSolution.init( *refCell );

        // copy to vector 
        const int numDofs = lfRef.numDofs();
        std::vector< double > values(numDofs, 0.0); 

        // on reference rank get values of function 
        if( refRank_ == comm.rank() )
        {
          assert(refCell->isLeaf());
          // unmark element from anything 
          func.space().gridPart().grid().mark(0, *refCell );

          referenceSolution.assign( lfRef );
          for(int i=0; i<numDofs; ++i) 
            values[i] = lfRef[i];
        }

        // copy values from ref rank to others 
        comm.broadcast( &values[0], numDofs, refRank_ );

        // copy values to reference solution 
        if( refRank_ != comm.rank() )
        {
          for(int i=0; i<numDofs; ++i) 
            referenceSolution[i] = values[i];
        }
      }

      if( cellPtr.size() > 0 ) 
      {
        assert( cellPtr.foundRefCell () );

        EntityPointerType& refCell = cellPtr.referenceCell();
        assert(refCell->isLeaf());
        
        const int quadOrd = 2 * func.space().order();
        QuadratureType quad( *refCell , quadOrd);
        
        const int quadNop = quad.nop();
        RangeType uRef[quadNop], uOld[quadNop];

        // evaluate reference solution 
        referenceSolution.init( *refCell );
        for(int qP = 0; qP < quadNop ; ++qP) 
          referenceSolution.evaluate(quad[qP],uRef[qP]);
        
        for(size_t i=0; i< cellPtr.size(); ++i) 
        {
          // get copy cell 
          EntityPointerType& cell = cellPtr.copyCell(i);
          assert(cell->isLeaf());
          
          // unmark element from anything 
          func.space().gridPart().grid().mark(0, *cell );
          
          // get entity 
          const EntityType& en = *cell; 
          
          // get quadrature 
          QuadratureType quad(en, quadOrd);

          LocalFunctionType lf = func.localFunction( en );
          for(int qP = 0; qP < quadNop ; ++qP) 
            lf.evaluate(quad[qP],uOld[qP]);

          lf.clear();

          // only to check 
          //lf[0] = (j+1) * 2;
          
          // get base function set 
          const BaseFunctionSetType & baseset = lf.baseFunctionSet();
          const int numDofs = lf.numDofs();
          for(int qP = 0; qP < quadNop ; ++qP) 
          {
            // evaluate function 
            evalRef(uRef[qP],uOld[qP],ret);  
            const double intel =  quad.weight(qP); // affine case 
            ret *= intel;
            // do projection 
            lf.axpy(quad[qP], ret);
            /*
            for(int i=0; i<numDofs; ++i) 
            {
              baseset.evaluate(i, quad[qP], phi);
              lf[i] += intel * (ret * phi) ;
            }
            */
          }
        }
      }

      // communicate new changed solution  
      func.space().communicate( func );

    } // end setFunction 
  };

  // do nothing for 3d 
  template <class GridPartType> 
  class CornerFinderBase<GridPartType,3> 
  { 
  public:
    CornerFinderBase(const GridPartType& gp) {}
    
    template <class DiscreteFunctionType>
    void setFunction(DiscreteFunctionType& func)
    {
    }
  };

  template <class GridPartType> 
  class CornerFinder : public CornerFinderBase<GridPartType,
                         GridPartType :: GridType :: dimension > 
  { 
    typedef CornerFinderBase<GridPartType,
            GridPartType :: GridType :: dimension > BaseType ;
  public:  
    CornerFinder(const GridPartType& gp) : BaseType(gp) {}
  };

}  // end namespace Dune
#endif
