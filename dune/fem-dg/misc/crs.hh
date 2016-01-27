#ifndef DUNE_CRSMATRIX_HH
#define DUNE_CRSMATRIX_HH

#include <cmath>
#include <cstddef>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/precision.hh>

namespace Dune
{

namespace Fem
{


  template<class K, int ROWS, int COLS>
  class BlockCRSMatrix
  {
    std::vector< double > data_;
    std::vector< int > rows_;
    std::vector< int > cols_;
  public:

    //! export size
    enum {
      //! The number of rows.
      rows = ROWS,
      //! The number of columns.
      cols = COLS
    };

    typedef K Field;
    typedef K value_type ;

    typedef std::vector< value_type > row_type ;

    //===== constructors
    /** \brief Default constructor
     */
    BlockCRSMatrix ()
     : data_(), rows_(rows + 1, int(0) ), cols_() {}

    BlockCRSMatrix ( const BlockCRSMatrix& other )
     : data_( other.data_ ), rows_( other.rows_ ), cols_( other.cols_ ) {}

    template <class Matrix>
    void set( const Matrix& matrix )
    {
      data_.clear();
      cols_.clear();

      assert( matrix.size() == rows );

      int count = 0;
      for( int row=0; row<rows; ++row )
        for( int col=0; col<cols; ++col )
          if( std::abs( matrix[ row ][ col ] ) > 0.0 )
            ++ count;

      data_.resize( count, Field( 0 ) );
      cols_.resize( count, int( 0 ) );

      count = 0 ;
      for( int row=0; row<rows; ++row )
      {
        rows_[ row ] = count ;
        for( int col=0; col<cols; ++col )
        {
          const Field val = matrix[ row ][ col ];
          if( std::abs( val ) > 0.0 )
          {
            assert( count < int( data_.size() ) );
            cols_[ count ] = col ;
            data_[ count ] = val ;
            ++ count ;
          }
        }
      }
      rows_[ rows ] = count ;
    }

    int N () const { return rows; }
    int M () const { return cols; }

    struct Add
    {
      template <class T>
      void operator()( T& a, const T& b ) const
      {
        a += b;
      }
    };

    struct Substract
    {
      template <class T>
      void operator()( T& a, const T& b ) const
      {
        a -= b;
      }
    };

    struct Set
    {
      template <class T>
      void operator()( T& a, const T& b ) const
      {
        a = b;
      }
    };

    template <class X, class Y>
    void mv ( const X& x, Y& y ) const
    {
      mult( x, y, Set() );
    }

    template <class X, class Y>
    void mmv ( const X& x, Y& y ) const
    {
      mult( x, y, Substract() );
    }

    template <class X, class Y>
    void umv ( const X& x, Y& y ) const
    {
      mult( x, y, Add() );
    }

    template <class X, class Y, class Op>
    void mult ( const X& x, Y& y, const Op op ) const
    {
      for( int row=0; row<rows; ++row )
      {
        Field sum = 0;
        for( int c = rows_[ row ]; c < rows_[ row+1 ]; ++ c )
        {
          sum += data_[ c ] * x[ cols_[ c ] ];
        }
        op( y[ row ], sum );
      }
    }

    template <class X, class Y>
    void mvb( const int blockSize, const X& x, Y& y ) const
    {
      // clear result
      y.clear();

      for( int row=0; row<rows; ++row )
      {
        for( int c = rows_[ row ]; c < rows_[ row+1 ]; ++ c )
        {
          Field matrix = data_[ c ];
          for( int r=0, ir = row * blockSize, jr = cols_[ c ] * blockSize ;
                r<blockSize; ++ r, ++ ir, ++ jr )
          {
            y[ ir ] += matrix * x[ jr ] ;
          }
        }
      }
    }

    //! Multiplies M from the right to this matrix
    template < class M >
    BlockCRSMatrix& rightmultiply (const M& )
    {
      abort();
      return *this;
    }

    template < class M >
    BlockCRSMatrix& leftmultiply (const M& )
    {
      abort();
      return *this;
    }

    void invert()
    {
      abort();
    }
  };

} // end namespace Fem

} // end namespace
#endif
