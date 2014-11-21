#ifndef DUNE_ROWWISEMATRIX_HH
#define DUNE_ROWWISEMATRIX_HH

//- system includes
#include <stack>
#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>

//- local includes
// #include "blockmatrix.hh"

namespace Dune
{

//*****************************************************************
//
//  --RowwiseMatrixObject
//
//*****************************************************************
template <class RowSpaceImp, class ColumnSpaceImp, class TraitsImp>
class RowwiseMatrixObject;

template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
struct RowwiseMatrixTraits
{
  typedef RowSpaceImp RowSpaceType;
  typedef ColSpaceImp ColumnSpaceType;
  typedef RowwiseMatrixTraits<RowSpaceType,ColumnSpaceType> ThisType;

  template <class OperatorTraits>
  struct MatrixObject
  {
    typedef RowwiseMatrixObject<RowSpaceType,ColumnSpaceType, OperatorTraits> MatrixObjectType;
  };
};



//! matrix object holding a blockamtrix
template <class RowSpaceImp, class ColumnSpaceImp, class TraitsImp>
class RowwiseMatrixObject
{
public:
  typedef RowSpaceImp RowSpaceType;
  typedef ColumnSpaceImp ColumnSpaceType ;

  //! number of rows of blocks
  enum { littleRows = RowSpaceType :: localBlockSize };
  //! number of columns of blocks
  enum { littleCols = ColumnSpaceType :: localBlockSize };

  typedef typename RowSpaceType::GridType::template Codim<0>::Entity EntityType;

  typedef RowwiseMatrixObject<RowSpaceType,ColumnSpaceType,TraitsImp> ThisType;

  //typedef RowwiseMatrix<double> MatrixType;
  //typedef MatrixType PreconditionMatrixType;

  template <class MatrixObjectImp>
  class LocalMatrix;

  struct LocalMatrixTraits
  {
    typedef RowSpaceImp DomainSpaceType ;
    typedef ColumnSpaceImp RangeSpaceType;
    typedef typename RowSpaceImp :: RangeFieldType RangeFieldType;
    typedef LocalMatrix<ThisType>  LocalMatrixType;
    //! number of rows of blocks
    enum { littleRows = RowSpaceType :: localBlockSize };
    //! number of columns of blocks
    enum { littleCols = ColumnSpaceType :: localBlockSize };

    typedef FieldMatrix< typename RowSpaceType::RangeFieldType,
                         littleRows, littleCols > LittleMatrixType ;

    typedef LittleMatrixType LittleBlockType ;

  };

  typedef typename LocalMatrixTraits :: LittleMatrixType LittleMatrixType ;
  //////////////////////////////////////////////////
  //
  //  LocalMatrix implementation
  //
  //////////////////////////////////////////////////
  template <class MatrixObjectImp>
  class LocalMatrix : public Fem::LocalMatrixDefault<LocalMatrixTraits>
  {
  public:
    //! type of base class
    typedef Fem::LocalMatrixDefault<LocalMatrixTraits> BaseType;

    //! type of matrix object
    typedef MatrixObjectImp MatrixObjectType;

    //! type of entries of little blocks
    typedef typename RowSpaceType :: RangeFieldType DofType;

    //! type of little blocks
    typedef typename MatrixObjectImp :: LittleMatrixType LittleMatrixType;

    //! type of row mapper
    typedef typename MatrixObjectType :: RowMapperType RowMapperType;
    //! type of col mapper
    typedef typename MatrixObjectType :: ColMapperType ColMapperType;

  private:
    // matrix to build
    const MatrixObjectType & matrixObj_;

    // point to local matrix
    mutable LittleMatrixType* matrix_;
  public:
    LocalMatrix(const MatrixObjectType & mObj,
                const RowSpaceType & rowSpace,
                const ColumnSpaceType & colSpace)
      : BaseType(rowSpace,colSpace)
      , matrixObj_(mObj)
      , matrix_(0)
    {
    }

    LocalMatrix(const LocalMatrix& org)
      : BaseType( org )
      , matrixObj_(org.matrixObj_)
      , matrix_(org.matrix_)
    {
    }

    //! initialize local matrix to entities
    void init(const EntityType& rowEntity,
              const EntityType& colEntity)
    {
      // initialize base functions sets
      BaseType :: init ( rowEntity , colEntity );
      matrix_ = & matrixObj_.getLocalMatrix( rowEntity, colEntity );
    }

  private:
    //! return matrix
    LittleMatrixType& matrix() const
    {
      assert( matrix_ );
      return *matrix_;
    }

    // check whether given (row,col) pair is valid
    void check(int localRow, int localCol) const
    {
      assert( localRow >= 0 && localRow < matrix().rows );
      assert( localCol >= 0 && localCol < matrix().cols );
    }

  public:
    //! add value to matrix
    void add(int localRow, int localCol , const DofType value)
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      matrix()[localRow][localCol] += value;
    }

    //! return matrix entry
    const DofType get(int localRow, int localCol ) const
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      return matrix()[localRow][localCol];
    }

    //! set matrix enrty to value
    void set(int localRow, int localCol, const DofType value)
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      matrix()[localRow][localCol] = value;
    }

    //! set matrix enrty to value
    void unitRow(const int localRow)
    {
      const int col = this->columns();
      for(int localCol=0; localCol<col; ++localCol)
      {
        matrix()[localRow][localCol] = 0;
      }
      matrix()[localRow][localRow] = 1;
    }

    //! clear all entries belonging to local matrix
    void clear ()
    {
      matrix() = 0.0;
    }

    //! scale local matrix
    void scale ( const DofType& value )
    {
      matrix() *= value ;
    }

    //! resort does nothing because already sorted in good way
    void resort ()
    {
    }
  };
  //////////////////////////////////////////////////////////
  //  end LocalMatrix
  //////////////////////////////////////////////////////////

public:
  //! type of local matrix
  typedef LocalMatrix<ThisType> ObjectType;
  typedef ThisType LocalMatrixFactoryType;
  typedef Fem :: ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
  //! type of local matrix
  typedef Fem :: LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

  typedef typename RowSpaceType :: BlockMapperType RowMapperType;
  typedef typename ColumnSpaceType :: BlockMapperType ColMapperType;

  typedef ThisType PreconditionMatrixType;
  typedef ThisType MatrixType ;

  typedef Fem :: AdaptiveDiscreteFunction<RowSpaceType> DestinationType;

  typedef Fem :: CommunicationManager<RowSpaceType> CommunicationManagerType;

  // row space
  const RowSpaceType & rowSpace_;
  // column space
  const ColumnSpaceType & colSpace_;

  // sepcial row mapper
  const RowMapperType& rowMapper_;
  // special col mapper
  const ColMapperType& colMapper_;

  // vector with local matrices
  mutable std::vector< LittleMatrixType > matrices_ ;
  typedef std::map< int, size_t >  ElNumberType ;
  //typedef std::map< int, ElNumberType > ElNumberMapType;

  //mutable ElNumberMapType elNumbers_ ;
  mutable ElNumberType elNumbers_ ;

  enum { dimension = RowSpaceType :: dimDomain };

  // local matrix stack
  mutable LocalMatrixStackType localMatrixStack_;

  //! constructor
  //! \param rowSpace space defining row structure
  //! \param colSpace space defining column structure
  //! \param paramfile parameter file to read variables
  //!         - Preconditioning: {0 == no, 1 == yes} used is SSOR
  RowwiseMatrixObject(const RowSpaceType & rowSpace,
                       const ColumnSpaceType & colSpace,
                       const std::string paramfile = "")
    : rowSpace_(rowSpace)
    , colSpace_(colSpace)
    , rowMapper_( rowSpace.blockMapper() )
    , colMapper_( colSpace.blockMapper() )
    , matrices_( 25 ) // reserve elements ( element + max neighboring elements )
    , elNumbers_ ()
    , localMatrixStack_(*this)
  {
    assert( rowSpace_.indexSet().size(0) ==
            colSpace_.indexSet().size(0) );
  }

  //! get local matrix for column of entity
  LittleMatrixType& getLocalMatrix(const EntityType& rowEntity,
                                   const EntityType& colEntity) const
  {
    assert( rowMapper_.maxNumDofs() == 1 );
    assert( colMapper_.maxNumDofs() == 1 );

    //const int row = rowMapper_.mapToGlobal(rowEntity, 0);
    //ElNumberType& colNumbers = elNumbers_[ row ] ;
    ElNumberType& colNumbers = elNumbers_ ;

    {
      const int idx = colMapper_.mapToGlobal(colEntity, 0);
      // find number of entity
      typename ElNumberType :: iterator it = colNumbers.find( idx );
      // if number not exists create new
      if( it == colNumbers.end() )
      {
        assert( colNumbers.size() < matrices_.size() );
        // next number is current size
        const size_t next = colNumbers.size() ;
        // insert to map
        colNumbers[ idx ] = next;

        // check size
        const size_t nextSize = next + 1;
        if( matrices_.size() < nextSize )
          matrices_.resize( nextSize );

        LittleMatrixType& matrix = matrices_[ next ];

        // initialize matrix
        matrix = 0;
        return matrix ;
      }
      else
      {
        assert( (*it).second < matrices_.size() );
        return matrices_[ (*it).second ];
      }
    }
  }

  //! apply matrix to discrete function
  template< class DomainFunction, class RangeFunction >
  void applyLocal( const EntityType& entity,
                   const DomainFunction &arg, RangeFunction &dest ) const
  {
    assert( rowMapper_.maxNumDofs() == 1 );

    //typedef typename ElNumberMapType :: iterator mapiterator ;
    //const mapiterator mapend = elNumbers_.end();
    //for( mapiterator rowit = elNumbers_.begin(); rowit != mapend; ++rowit )
    {
      // get current row
      //const int row = (*rowit).first ;
      //ElNumberType& colNumbers = (*rowit).second ;

      ElNumberType& colNumbers = elNumbers_ ;
      const int row = rowMapper_.mapToGlobal( entity , 0 );
      //std::cout << "Dof mult for en = " << row << std::endl;

      typedef typename DomainFunction :: ConstDofBlockPtrType DomainBlockPtrType;
      typedef typename RangeFunction  :: DofBlockPtrType RangeBlockPtrType;
      typedef typename RangeFunction  :: RangeFieldType RangeFieldType ;

      RangeBlockPtrType rangeBlock = dest.block( row );

      // multiply with all column entries
      typedef typename ElNumberType :: iterator iterator ;
#ifndef NDEBUG
      int preCol = -1;
#endif
      const iterator end = colNumbers.end();
      for(iterator it = colNumbers.begin (); it != end ; ++it )
      {
        // get column block
        DomainBlockPtrType domBlock = arg.block( (*it).first );
#ifndef NDEBUG
        assert( preCol < (*it).first );
        preCol = (*it).first;
#endif
        // get corresponding matrix
        const LittleMatrixType& matrix = matrices_[ (*it).second ];

        //std::cout << "Have col = " << (*it).first  << " m = " << (*it).second << "\n";
        // do multiplication
        for( int i = 0; i < littleRows; ++i)
        {
          RangeFieldType sum = 0;
          for( int k = 0 ; k < littleCols; ++k )
          {
            sum += matrix[ i ][ k ] * (*domBlock)[ k ] ;
          }
          (*rangeBlock)[ i ] += sum ;
        }
      }
    }
  }

  //! return reference to stability matrix
  MatrixType & matrix()
  {
    return *this;
  }

  //! set all matrix entries to zero
  void clear()
  {
    // clear stored numbers
    elNumbers_.clear();
  }

  //! return true if precoditioning matrix is provided
  bool hasPreconditionMatrix () const { return false; }

  //! return reference to preconditioner (here also systemMatrix)
  //const PreconditionMatrixType& preconditionMatrix () const { return matrix(); }

  //! reserve memory corresponnding to size of spaces
  void reserve(bool verbose = false )
  {
  }

  //! apply matrix to discrete function
  template< class DomainFunction, class RangeFunction >
  void apply ( const DomainFunction &arg, RangeFunction &dest ) const
  {
    abort();
  }

  //! mult method of matrix object used by oem solver
  void multOEM(const double * arg, double * dest) const
  {
    abort();
  }

  //! communicate data
  void communicate(const double * arg) const
  {
    /*
    if( rowSpace_.grid().comm().size() <= 1 ) return ;

    DestinationType tmp("BlockMatrixObject::communicate_tmp",rowSpace_,arg);
    rowSpace_.communicate( tmp );
    */
  }

  //! resort row numbering in matrix to have ascending numbering
  void resort()
  {
  }

  //! print matrix
  void print(std::ostream & s) const
  {
  }

 //! interface method from LocalMatrixFactory
  ObjectType* newObject() const
  {
    return new ObjectType(*this,
                          rowSpace_,
                          colSpace_);
  }

  //! return local matrix
  LocalMatrixType localMatrix(const EntityType& rowEntity,
                              const EntityType& colEntity) const
  {
    return LocalMatrixType(localMatrixStack_,rowEntity,colEntity);
  }

  //! empty method as we use right preconditioning here
  void createPreconditionMatrix() {}

  /** \brief delete all row belonging to a hanging node and rebuild them */
  template <class HangingNodesType>
  void changeHangingNodes(const HangingNodesType& hangingNodes)
  {
    DUNE_THROW(NotImplemented,"RowwiseMatrixObject: changeHangingNodes not implemented for this matrix object");
  }
};




} // end namespace Dune
#endif
