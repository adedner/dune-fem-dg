#ifndef DUNE_BASEFUNCTIONSETS_VECTORCODEGEN_HH
#define DUNE_BASEFUNCTIONSETS_VECTORCODEGEN_HH

#define NEWBASEFCT_CACHING

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>

#include <dune/common/exceptions.hh>
#include <dune/fem/io/io.hh>

namespace Dune
{

namespace Fem { 

  /** \brief default code generator methods */
  template < int simdWidth = 2 >
  struct VectorCodeGenerator 
  {
    static const char* restrictKey() 
    {
      return "__restrict__";
    }

    static void writePreCompHeader(std::ostream& out, const int stage ) 
    {
      const char* codegenPreCompVar = "CODEGEN_COMPILE_INNERLOOPS";
      if( stage == -1 ) 
      {
        out << "#if ! " << codegenPreCompVar << std::endl;
      }
      else if( stage == 0 ) 
      {
        out << "#else" << std::endl;
        out << "#if " << codegenPreCompVar << " == 1" << std::endl;
        out << "extern \"C\" {" << std::endl
            << "  extern inline" << std::endl;
        out << "#endif" << std::endl;
      }
      else if( stage == 1 )
      {
        out << "#if  " << codegenPreCompVar << " == 1" << std::endl;
        out << "  ;" << std::endl;
        out << "}" << std::endl;
        out << "#else" << std::endl;
      }
      else 
      {
        out << "#endif" << std::endl;
        out << "#endif" << std::endl;
      }
    }

    // generate inner loop function name 
    static std::string generateFunctionName( const std::string& prefix, 
                                      const int simdW, const int dimRange, 
                                      const size_t numRows, const size_t numCols )
    {
      std::stringstream funcName;
      funcName << prefix << "_" << simdWidth << "_" << dimRange << "_" << numRows << "_" << numCols ;
      return funcName.str();
    }

    static void writeInnerLoopEval(std::ostream& out, const int simdW, const int dimRange, const size_t numRows, const size_t numCols )
    {
      out << "      for(int row = 0; row < " << numRows << " ; ++row )" << std::endl;
      out << "      {" << std::endl;
      if( simdW == 1 ) 
      {
          out << "        const double phi0 = rangeStorage[ rowMap[ row ] * " << numCols << " + col ][ 0 ];" << std::endl;  
      }
      else 
      {
        for( int i = 0 ; i< simdW ; ++ i ) 
          out << "        const double phi" << i << " = base" << i << "[ row ];" << std::endl;
      }
      for(int r = 0; r < dimRange; ++ r ) 
      {
        out << "        result" << r << "[ row ] += phi0 * dof0" << r;
        for( int i=1; i<simdW; ++i )
          out << " + phi" << i << " * dof" << i << r;
        out << " ;" << std::endl;
      }
      out << "      }" << std::endl;
    }

    static void evaluateCodegen(std::ostream& out, 
                                const int dim, 
                                const int dimRange, 
                                const size_t numRows, 
                                const size_t numCols ) 
    {
      const std::string funcName = 
            generateFunctionName( "evalRangesLoop", simdWidth, dimRange, numRows, numCols );

      writePreCompHeader( out, -1 );

      out << "template <class BaseFunctionSet>" << std::endl;
      out << "struct EvaluateRanges<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
      out << "{" << std::endl;
      out << "  template< class QuadratureType,"<< std::endl;
      out << "            class RangeVectorType," << std::endl;
      out << "            class RangeFactorType," << std::endl;
      out << "            class LocalDofVectorType>" << std::endl;
      out << "  static void eval( const QuadratureType& quad," << std::endl;
      out << "                    const RangeVectorType& rangeStorage," << std::endl; 
      out << "                    const LocalDofVectorType& dofs," << std::endl;
      out << "                    RangeFactorType &rangeVector)" << std::endl;
      out << "  {" << std::endl;
      out << "    typedef typename ScalarRangeType :: field_type field_type;" << std::endl;
      //out << "    typedef typename RangeVectorType :: value_type value_type;" << std::endl; 
      out << "    typedef RangeType value_type;" << std::endl; 

      out << "    // only evaluate cachingPoint once" << std::endl;
      out << "    int rowMap[ " << numRows << " ] ;" << std::endl;
      out << "    for( int row = 0; row < " << numRows << "; ++ row )" << std::endl;
      out << "    {" << std::endl;
      out << "      rowMap[ row ] = quad.cachingPoint( row  );" << std::endl; 
      out << "    }" << std::endl;
      out << std::endl;

      // make length simd conform 
      out << "    field_type resultTmp[ " << numRows * dimRange << " ] = { 0 };" << std::endl << std::endl;

      for(int r=0; r<dimRange ; ++r ) 
      {
        out << "    field_type* result" << r << " = resultTmp + " << r * numRows << " ;" << std::endl; 
      }
      out << std::endl;

      const size_t simdCols = simdWidth * ( numCols / simdWidth );
      if( simdCols > 0 ) 
      {
        for( int i=0; i< simdWidth; ++ i ) 
        {
          out << "    field_type base" << i << "[ " << numRows << " ];" << std::endl;
        }
        out << "    for( int col = 0, dof = 0 ; col < "<< simdCols <<" ; col += " << simdWidth << ", dof += " << simdWidth * dimRange<< " )"<<std::endl;
        out << "    {" << std::endl;
        out << "      for( int row = 0 ; row < " << numRows << " ; ++ row )" << std::endl;
        out << "      {" << std::endl;
        //out << "        const value_type& rangeStorageRow = rangeStorage[ rowMap[ row ] ];" << std::endl;
        out << "        int idx = rowMap[ row ] * " << numCols << " + col;" << std::endl;
        for( int i=0; i< simdWidth; ++ i )
        {
          const char* plusplus = (i == simdWidth-1) ? "  " : "++";
          out << "        base" << i << "[ row ] = rangeStorage[ idx" << plusplus << " ][ 0 ];" << std::endl;
        } 
        out << "      }" << std::endl;
        out << "      " << funcName << "(" << std::endl;
        out << "                ";
        for( int i = 0; i< simdWidth * dimRange; ++i ) 
          out << " dofs[ dof + " << i << " ],";
        out << std::endl << "                 ";
        for( int i=0; i< simdWidth; ++i )
          out << "base" << i << ", ";
        out << std::endl << "                 ";
        for( int r = 0; r < dimRange-1; ++r) 
          out << "result" << r << ", ";
        out << "result" << dimRange-1 << " );" << std::endl;
        out << "    }"<< std::endl;  
        out << std::endl;
      }

      if( numCols > simdCols ) 
      {
        out << "    // remainder iteration" << std::endl;
        out << "    for( int col = " << simdCols << ", dof = " << simdCols * dimRange << " ; col < " << numCols << " ; ++col )" << std::endl;
        out << "    {" << std::endl;
        for( int r=0; r<dimRange; ++r ) 
          out << "      const double dof0" << r << " = dofs[ dof++ ];" << std::endl;
        writeInnerLoopEval( out, 1, dimRange, numRows, numCols ); 
        out << "    }" << std::endl;
        out << std::endl;
      }

      out << "    // store result" << std::endl;   
      out << "    for(int row = 0; row < " << numRows << " ; ++row )" << std::endl;
      out << "    {" << std::endl;
      out << "      RangeType& result = rangeVector[ row ];" << std::endl;
      for( int r = 0 ; r < dimRange; ++ r ) 
      {
        out << "      result[ " << r << " ] = result" << r << "[ row ];" << std::endl;
      }
      out << "    }" << std::endl;
      out << "  }" << std::endl << std::endl;
      out << "};" << std::endl;

      writePreCompHeader( out, 0 );
      out << "  void " << funcName << "(" << std::endl;
      for( int i=0; i<simdWidth; ++i ) 
      {
        out << "        ";
        for( int r=0; r<dimRange; ++ r ) 
          out << "const double dof"<< i << r << ", ";
        out << std::endl;
      }
      for( int i=0; i<simdWidth; ++ i ) 
        out << "        const double* " << restrictKey() << " base" << i << "," << std::endl; 
      for( int r=0; r<dimRange; ++ r ) 
      {
        out << "        double* "<< restrictKey() << " result" << r;
        if( r == dimRange-1 ) out << " )" << std::endl;
        else out << "," << std::endl; 
      }

      writePreCompHeader( out, 1 );
      out << "  {" << std::endl;
      writeInnerLoopEval( out, simdWidth, dimRange, numRows, numCols ); 
      out << "  }" << std::endl;
      writePreCompHeader( out, 2 );
    }

    static void writeInnerLoop(std::ostream& out, const int simdW, const int dimRange, const size_t numCols )
    {
      for( int i=0; i< simdW; ++i ) 
      {
        for( int r=0; r< dimRange; ++r ) 
        {
          out << "      const double fac" << i << r << " = rangeFactor" << i << "[ " << r << " ];" << std::endl;
        }
      }
      if( simdW == 1 )
        out << "      int rowCol = quad.cachingPoint( row ) * " << numCols <<";" << std::endl;
      out << "      for(int col = 0 ; col < " << numCols << " ; ++ col";
      if( simdW == 1 )
        out << ", ++rowCol";
      out << " )" << std::endl;
      out << "      {" << std::endl;
      for( int i = 0 ; i< simdW ; ++ i ) 
      {
        if( simdW == 1 )
          out << "        const double phi" << i << " = rangeStorage[ rowCol ][ 0 ];" << std::endl;  
        else 
          out << "        const double phi" << i << " = base" << i  << " [ col ];" << std::endl;  
      }
      for(int r = 0; r < dimRange; ++ r ) 
      {
        out << "        dofs" << r << "[ col ] += phi0 * fac0" << r;
        for( int i=1; i<simdW; ++i )
        {
          out << " + phi" << i << " * fac" << i << r ; 
        }
        out << " ;" << std::endl;
      }
      out << "      }" << std::endl;
    }

    static void axpyCodegen(std::ostream& out, 
            const int dim, const int dimRange, const size_t numRows, const size_t numCols ) 
    {
      const std::string funcName = 
        generateFunctionName( "axpyRangesLoop", simdWidth, dimRange, numRows, numCols );

      writePreCompHeader( out, -1 );

      out << "template <class BaseFunctionSet>" << std::endl;
      out << "struct AxpyRanges<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
      out << "{" << std::endl;

      out << std::endl;
      out << "  template< class QuadratureType,"<< std::endl;
      out << "            class RangeVectorType," << std::endl;
      out << "            class RangeFactorType," << std::endl;
      out << "            class LocalDofVectorType>" << std::endl;
      out << "  static void axpy( const QuadratureType& quad," << std::endl;
      out << "                    const RangeVectorType& rangeStorage," << std::endl; 
      out << "                    const RangeFactorType& rangeFactors," << std::endl;
      out << "                    LocalDofVectorType& dofs)" << std::endl;
      out << "  {" << std::endl;

      ////////////////////////////////////////////////////
      // axpy 
      ////////////////////////////////////////////////////

      //out << "    typedef typename RangeVectorType :: value_type value_type;" << std::endl; 
      out << "    typedef RangeType value_type;" << std::endl; 
      out << "    typedef typename ScalarRangeType :: field_type field_type;" << std::endl;
      out << std::endl;

      out << "    double dofResult[ " << numCols * dimRange << " ] = { 0 };" << std::endl << std::endl;
      const size_t simdRows  = simdWidth * (numRows / simdWidth) ;

      if( simdRows > 0 ) 
      {
        out << "    for( int row = 0; row < "<< simdRows << " ; row += " << int(simdWidth) << " )" << std::endl;
        out << "    {" << std::endl;
        for( int i=0; i<simdWidth; ++ i ) 
          out << "      const double* rangeFactor" << i << " = &rangeFactors[ row + " << i << " ][ 0 ];" << std::endl;
        out << "      " << funcName << "("; 
        for( int i = 0; i < simdWidth; ++i ) 
          out << " &rangeStorage[ quad.cachingPoint( row + " << i << " ) * " << numCols << " ][ 0 ],";
        out << std::endl; 
        out << "                 rangeFactor0, ";
        for( int i=1; i<simdWidth; ++ i ) 
          out << "rangeFactor" << i << ",";
        out << std::endl;
        out << "                ";
        for( int r = 0; r < dimRange; ++ r ) 
        {
          out << " dofResult + " << r * numCols;
          if( r == dimRange-1 ) 
            out << " );" << std::endl;
          else  
            out << ",";
        }
        out << std::endl;
        out << "    }" << std::endl;
        out << std::endl;
      }

      out << "    double* dofs0 = dofResult;" << std::endl;
      for( int r = 1; r < dimRange; ++ r ) 
        out << "    double* dofs" << r << " = dofResult + " << r * numCols << ";" << std::endl;
      out << std::endl;

      if( numRows > simdRows )
      {
        out << "    // remainder iteration" << std::endl;
        out << "    for( int row = " << simdRows << " ; row < " << numRows << " ; ++row )" << std::endl;
        out << "    {" << std::endl;
        out << "      const double* rangeFactor0 = &rangeFactors[ row ][ 0 ];" << std::endl;
        writeInnerLoop( out, 1, dimRange, numCols ); 
        out << "    }" << std::endl;
        out << std::endl;
      }

      out << "    // sum up results (transform from variable based to point based layout)" << std::endl;
      out << "    for( int col = 0, dof = 0 ; col < "<< numCols << " ; ++col )" << std::endl;
      out << "    {" << std::endl;
      for( int r = 0 ; r < dimRange; ++ r ) 
        out << "      dofs[ dof++ ] += dofs" << r << "[ col ];" << std::endl;
      out << "    }" << std::endl;

      out << "  }" << std::endl << std::endl;
      out << "};" << std::endl;

      ///////////////////////////////////
      //  inner loop 
      ///////////////////////////////////
      writePreCompHeader( out, 0 );
      out << "  void " << funcName << "(" << std::endl;
      out << "       const double* " << restrictKey() << " base0," << std::endl; 
      for( int i=1; i<simdWidth; ++ i ) 
        out << "       const double* " << restrictKey() << " base" << i << "," << std::endl; 
      for( int i=0; i<simdWidth; ++ i ) 
        out << "       const double* " << restrictKey() << " rangeFactor" << i << "," << std::endl; 
      for( int r = 0; r < dimRange; ++r ) 
      {
        out << "       double* " << restrictKey() << " dofs" << r;
        if( r == dimRange-1 ) 
          out << " )" << std::endl;
        else 
          out << "," << std::endl;
      }
      writePreCompHeader( out, 1 );
      out << "  {" << std::endl;
      writeInnerLoop( out, simdWidth, dimRange, numCols ); 
      out << "  }" << std::endl;
      writePreCompHeader( out, 2 );
    }

    static void writeInnerJacEvalLoop(std::ostream& out, const int simdW, const int dim,
                                      const int dimRange, const size_t numRows, const size_t numCols  ) 
    { 
      out << "      for(int row = 0; row < " << numRows << " ; ++row )" << std::endl;
      out << "      {" << std::endl;
      if( simdW == 1 ) 
      {
        out << "        int idx = rowMap[ row ] * " << numCols << " + col ;" << std::endl;
        out << "        const GeometryJacobianType& gjit = geometry.jacobianInverseTransposed( quad.point( row ) );" << std::endl << std::endl;
        for( int i = 0 ; i< simdW ; ++ i ) 
        {
          const char* plusplus = (i == simdW-1) ? "  " : "++";
          out << "        gjit.mv( jacStorage[ idx" << plusplus << " ][ 0 ], gradPhi" << i << " );" << std::endl;
          for( int d = 0 ; d < dim; ++ d ) 
            out << "        const double phi" << i << d << " = gradPhi" << i << "[ " << d << " ];" << std::endl;
        }
        out << std::endl;
      }
      else 
      {
        for( int d = 0; d < dim ; ++ d ) 
        {
          for( int i = 0 ; i< simdW ; ++ i ) 
            out << "        const double phi" << i << d << " = base" << i << d << "[ row ];" << std::endl;
        }
      }
      for( int d = 0; d < dim ; ++ d ) 
      {
        for(int r = 0; r < dimRange; ++ r ) 
        {
          out << "        result" << r << d << "[ row ] += phi0" << d << " * dof0" << r;
          for( int i=1; i<simdW; ++i )
            out << " + phi" << i << d << " * dof" << i << r;
          out << " ;" << std::endl;
        }
      }
      out << "      }" << std::endl;
    }

    static void evaluateJacobiansCodegen(std::ostream& out, 
            const int dim, const int dimRange, const size_t numRows, const size_t numCols ) 
    {
      const std::string funcName = 
            generateFunctionName( "evalJacobianLoop", simdWidth, dimRange, numRows, numCols );

      writePreCompHeader( out, -1 );

      out << "template <class BaseFunctionSet>" << std::endl;
      out << "struct EvaluateJacobians<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
      out << "{" << std::endl;
      out << "  template< class QuadratureType,"<< std::endl;
      out << "            class JacobianRangeVectorType," << std::endl;
      out << "            class LocalDofVectorType," << std::endl;
      out << "            class JacobianRangeFactorType>" << std::endl;
      out << "  static void eval( const QuadratureType&," << std::endl;
      out << "                    const Fem :: EmptyGeometry&," << std::endl; 
      out << "                    const JacobianRangeVectorType&," << std::endl; 
      out << "                    const LocalDofVectorType&," << std::endl;
      out << "                    JacobianRangeFactorType &)" << std::endl;
      out << "  {" << std::endl;
      out << "    std::cerr << \"ERROR: wrong code generated for VectorialBaseFunctionSet::axpyJacobians\" << std::endl;" << std::endl;
      out << "    abort();" << std::endl;
      out << "  }" << std::endl;
      out << "};" << std::endl << std::endl;
      out << "template <class BaseFunctionSet, class Geometry>" << std::endl;
      out << "struct EvaluateJacobians<BaseFunctionSet, Geometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
      out << "{" << std::endl;
      out << "  template< class QuadratureType,"<< std::endl;
      out << "            class JacobianRangeVectorType," << std::endl;
      out << "            class LocalDofVectorType," << std::endl;
      out << "            class JacobianRangeFactorType>" << std::endl;
      out << "  static void eval( const QuadratureType& quad," << std::endl;
      out << "                    const Geometry& geometry," << std::endl; 
      out << "                    const JacobianRangeVectorType& jacStorage," << std::endl; 
      out << "                    const LocalDofVectorType& dofs," << std::endl;
      out << "                    JacobianRangeFactorType& jacFactors)" << std::endl;
      out << "  {" << std::endl;
      out << "    evalJac( quad, geometry, jacStorage, dofs, jacFactors, jacFactors[ 0 ] );" << std::endl;
      out << "  }" << std::endl;
      out << "private:" << std::endl;
      out << "  template< class QuadratureType,"<< std::endl;
      out << "            class JacobianRangeVectorType," << std::endl;
      out << "            class LocalDofVectorType," << std::endl;
      out << "            class JacobianRangeFactorType," << std::endl;
      out << "            class GlobalJacobianRangeType>" << std::endl;
      out << "  static void evalJac( const QuadratureType& quad," << std::endl;
      out << "                       const Geometry& geometry," << std::endl; 
      out << "                       const JacobianRangeVectorType& jacStorage," << std::endl; 
      out << "                       const LocalDofVectorType& dofs," << std::endl;
      out << "                       JacobianRangeFactorType& jacVector," << std::endl;
      out << "                       const GlobalJacobianRangeType& )" << std::endl;
      out << "  {" << std::endl;
      //out << "    typedef typename JacobianRangeVectorType :: value_type  value_type;" << std::endl; 
      out << "    typedef JacobianRangeType value_type;" << std::endl; 
      out << "    typedef typename JacobianRangeType :: field_type field_type;" << std::endl;
      out << "    typedef typename Geometry::Jacobian GeometryJacobianType;" << std::endl;
      out << "    // only evaluate cachingPoint once" << std::endl;
      out << "    int rowMap[ " << numRows << " ] ;" << std::endl;
      out << "    for( int row = 0; row < " << numRows << "; ++ row )" << std::endl;
      out << "    {" << std::endl;
      out << "      rowMap[ row ] = quad.cachingPoint( row  );" << std::endl; 
      out << "    }" << std::endl;
      out << std::endl;

      for( int d = 0; d < dim ; ++ d ) 
      {
        // make length simd conform 
        out << "    field_type resultTmp" << d << "[ " << numRows * dimRange << " ] = { 0 };" << std::endl;
      }
      out << std::endl;

      for( int d = 0; d < dim ; ++ d ) 
      {
        for(int r=0; r<dimRange ; ++r ) 
        {
          out << "    field_type* result" << r << d <<" = resultTmp" << d << " + " << r * numRows << " ;" << std::endl; 
        }
      }
      out << std::endl;

      const size_t simdNumCols = simdWidth * ( numCols / simdWidth );
      if( simdNumCols > 0 ) 
      {
        for( int d = 0; d < dim; ++d )
        {
          for( int i=0; i< simdWidth; ++ i ) 
          {
            out << "    field_type base" << i << d << "[ " << numRows << " ];" << std::endl;
          }
        }

        out << "    typedef typename GlobalJacobianRangeType :: row_type JacobianRangeType;" << std::endl;
        for( int i=0; i< simdWidth; ++ i ) 
        out << "    JacobianRangeType gradPhi" << i << ";" << std::endl;
        out << "    for( int col = 0, dof = 0 ; col < "<< simdNumCols <<" ; col += " << simdWidth << ", dof += " << simdWidth * dimRange<< " )"<<std::endl;
        out << "    {" << std::endl;
        out << "      for( int row = 0; row < " << numRows << " ; ++ row )" << std::endl;
        out << "      {" << std::endl;
        out << "        int idx = rowMap[ row ] * " << numCols << " + col ;" << std::endl;
        out << "        // use reference to GeometryJacobianType to make code compile with SPGrid Geometry" << std::endl;
        out << "        const GeometryJacobianType& gjit = geometry.jacobianInverseTransposed( quad.point( row ) );" << std::endl << std::endl;
        for( int i=0; i< simdWidth; ++ i ) 
        {
          const char* plusplus = (i == simdWidth-1) ? "  " : "++";
          out << "        gjit.mv( jacStorage[ idx" << plusplus << " ][ 0 ], gradPhi" << i << " );" << std::endl;
        }
        out << std::endl;
        for( int d = 0; d < dim; ++ d ) 
        {
          for( int i=0; i< simdWidth; ++ i ) 
          {
            out << "        base" << i << d << "[ row ] = gradPhi"<< i << "[ " << d << " ];" << std::endl;
          }
        }
        out << "      }" << std::endl << std::endl;

        out << "      " << funcName << "(";
        for( int i = 0; i< simdWidth * dimRange; ++i ) 
          out << " dofs[ dof + " << i << " ],";
        out << std::endl << "                 ";
        for( int d = 0; d < dim; ++ d )
        {
          for( int i=0; i< simdWidth; ++i )
            out << "base" << i << d << ", ";
        }
        out << std::endl << "                 ";
        for( int d = 0; d < dim; ++ d )
        {
          for( int r = 0; r < dimRange; ++r) 
          {
            out << "result" << r << d;
            if( r * d == (dim-1)*(dimRange-1) ) out << std::endl;
            else out << ", ";
          }
        }
        out << "      );" << std::endl;
        out << "    }"<< std::endl;  
        out << std::endl;
      }

      // remainder iteration
      if( numCols > simdNumCols ) 
      {
        out << "    // remainder iteration" << std::endl;
        out << "    for( int col = " << simdNumCols << ", dof = " << simdNumCols * dimRange << " ; col < " << numCols << " ; ++col )" << std::endl;
        out << "    {" << std::endl;
        for( int r=0; r<dimRange; ++r ) 
          out << "      const double dof0" << r << " = dofs[ dof++ ];" << std::endl;
        writeInnerJacEvalLoop( out, 1, dim, dimRange, numRows, numCols ); 
        out << "    }" << std::endl;
        out << std::endl;
      }

      out << "    // store result" << std::endl;   
      out << "    for(int row = 0; row < " << numRows << " ; ++row )" << std::endl;
      out << "    {" << std::endl;
      out << "      GlobalJacobianRangeType& result = jacVector[ row ];" << std::endl;
      for( int d = 0 ; d < dim; ++ d ) 
      {
        for( int r = 0 ; r < dimRange; ++ r ) 
        {
          out << "      result[ " << r << " ][ " << d <<" ] = result" << r << d << "[ row ];" << std::endl;
        }
      }
      out << "    }" << std::endl;
      out << "  }" << std::endl << std::endl;
      out << "};" << std::endl;

      writePreCompHeader( out, 0 );
      out << "  void " << funcName << "(" << std::endl;
      for( int i=0; i<simdWidth; ++i ) 
      {
        out << "                        ";
        for( int r=0; r<dimRange; ++ r ) 
          out << " const double dof"<< i << r << ",";
        out << std::endl;
      }
      for( int d=0; d<dim; ++ d ) 
      {
        for( int i=0; i<simdWidth; ++ i ) 
          out << "                         const double* " << restrictKey() << " base" << i << d << "," << std::endl; 
      }
      for( int d=0; d<dim; ++ d ) 
      {
        for( int r=0; r<dimRange; ++ r ) 
        {
          out << "                         double* "<< restrictKey() << " result" << r << d;
          if( d * r == (dim-1)*(dimRange-1) ) out << " )" << std::endl;
          else out << "," << std::endl; 
        }
      }

      writePreCompHeader( out, 1 );
      out << "  {" << std::endl;
      writeInnerJacEvalLoop( out, simdWidth, dim, dimRange, numRows, numCols ); 
      out << "  }" << std::endl;
      writePreCompHeader( out, 2 );
    }

    static void writeInnerLoopAxpyJac(std::ostream& out, const int dim, const int dimRange, const size_t numCols )
    {
      out << "      for( int col = 0; col < " << numCols << " ; ++col )" << std::endl;
      out << "      {" << std::endl;
      for( int d =0; d < dim; ++d )
        out << "        const double phi" << d << " = base" << d << "[ col ];" << std::endl;

      for( int r = 0; r < dimRange; ++r )
      {
        out << "        result" << r << "[ col ]  +=  phi0 * jacFactorInv0" << r;
        for( int d=1; d < dim; ++d )
          out << "  +  phi" << d << " * jacFactorInv" << d << r;
        out << ";" << std::endl;
      }
      out << "      }" << std::endl;
    }

    static void axpyJacobianCodegen(std::ostream& out, 
            const int dim, const int dimRange, const size_t numRows, const size_t numCols ) 
    {
      const std::string funcName = 
            generateFunctionName( "axpyJacobianLoop", simdWidth, dimRange, numRows, numCols );

      writePreCompHeader( out, -1 );

      out << "template <class BaseFunctionSet>" << std::endl;
      out << "struct AxpyJacobians<BaseFunctionSet, Fem :: EmptyGeometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
      out << "{" << std::endl;
      out << "  template< class QuadratureType,"<< std::endl;
      out << "            class JacobianRangeVectorType," << std::endl;
      out << "            class JacobianRangeFactorType," << std::endl;
      out << "            class LocalDofVectorType>" << std::endl;
      out << "  static void axpy( const QuadratureType&," << std::endl;
      out << "                    const Fem :: EmptyGeometry&," << std::endl; 
      out << "                    const JacobianRangeVectorType&," << std::endl; 
      out << "                    const JacobianRangeFactorType&," << std::endl;
      out << "                    LocalDofVectorType&)" << std::endl;
      out << "  {" << std::endl;
      out << "    std::cerr << \"ERROR: wrong code generated for VectorialBaseFunctionSet::axpyJacobians\" << std::endl;" << std::endl;
      out << "    abort();" << std::endl;
      out << "  }" << std::endl;
      out << "};" << std::endl << std::endl;
      out << "template <class BaseFunctionSet, class Geometry>" << std::endl;
      out << "struct AxpyJacobians<BaseFunctionSet, Geometry, " << dimRange << ", " << numRows << ", " << numCols << ">" << std::endl;
      out << "{" << std::endl;

      out << "  template< class QuadratureType,"<< std::endl;
      out << "            class JacobianRangeVectorType," << std::endl;
      out << "            class JacobianRangeFactorType," << std::endl;
      out << "            class LocalDofVectorType>" << std::endl;
      out << "  static void axpy( const QuadratureType& quad," << std::endl;
      out << "                    const Geometry& geometry," << std::endl; 
      out << "                    const JacobianRangeVectorType& jacStorage," << std::endl; 
      out << "                    const JacobianRangeFactorType& jacFactors," << std::endl;
      out << "                    LocalDofVectorType& dofs)" << std::endl;
      out << "  {" << std::endl;
      out << "    typedef JacobianRangeType value_type;" << std::endl; 
      out << "    typedef typename JacobianRangeType :: field_type field_type;" << std::endl;
      const size_t dofs = dimRange * numCols ;
      out << "    field_type result [ " << dofs << " ] = {"; 
      for( size_t dof = 0 ; dof < dofs-1 ; ++ dof ) out << " 0,";
      out << " 0 };" << std::endl << std::endl;
      for( int r=0; r<dimRange; ++r ) 
        out << "    field_type* result" << r << " = result + " << r * numCols << ";" << std::endl; 
      out << std::endl;

      for( int d =0; d < dim; ++d ) 
      {
        out << "    field_type base"<< d << "[ " << numCols << " ] ;" << std::endl; 
      }

      out << "    for( int row = 0; row < " << numRows << " ; ++ row )" << std::endl;
      out << "    {" << std::endl;
      out << "      typedef typename Geometry::Jacobian GeometryJacobianType;" << std::endl;
      out << "      // use reference to GeometryJacobianType to make code compile with SPGrid Geometry" << std::endl;
      out << "      const GeometryJacobianType& gjit = geometry.jacobianInverseTransposed( quad.point( row ) );" << std::endl << std::endl;
      out << "      JacobianRangeType jacFactorTmp;" << std::endl;
      out << "      for( int r = 0; r < " << dimRange << " ; ++r )" << std::endl;
      out << "      {"<<std::endl; 
      out << "        gjit.mtv( jacFactors[ row ][ r ], jacFactorTmp[ r ] );" << std::endl;
      out << "      }" << std::endl << std::endl;

      out << "      for( int col = 0, rowCol = quad.cachingPoint( row ) * " << numCols << "; col < " << numCols << " ; ++ col, ++rowCol )" << std::endl;
      out << "      {" << std::endl;
      for( int d =0; d < dim; ++d ) 
        out << "        base"<< d << "[ col ] = jacStorage[ rowCol ][ 0 ][ " << d << " ];" << std::endl;

      out << "      }" << std::endl;

      out << "      // calculate updates" << std::endl;
      out << "      " << funcName << "(";
      for( int d =0; d < dim; ++d ) out << "base" << d << ", ";
      out << std::endl;
      for( int i =0; i < dim; ++i )
      {
        out << "                               ";
        for( int r = 0; r < dimRange; ++ r ) 
          out << "jacFactorTmp[ " << r  << " ][ " << i << " ], ";
        out << std::endl;
      }
      out << "                               ";
      for( int r = 0; r < dimRange; ++ r ) 
      {
        out << "result" << r;
        if( r == dimRange -1 ) 
          out << " );" << std::endl ;
        else   
          out << ", "; 
      }
      out << "    }" << std::endl << std::endl;

      out << "    // sum up results (transform from variable based to point based layout)" << std::endl;
      out << "    for( int col = 0, dof = 0 ; col < "<< numCols << " ; ++col )" << std::endl;
      out << "    {" << std::endl;
      for( int r = 0 ; r < dimRange; ++ r ) 
        out << "      dofs[ dof++ ]  +=  result" << r << "[ col ];" << std::endl;
      out << "    }" << std::endl;

      out << "  }" << std::endl << std::endl;
      out << "};" << std::endl;


      ///////////////////////////////////
      //  inner loop
      ///////////////////////////////////
      writePreCompHeader( out, 0 );

      out << "  void " << funcName << "(" << std::endl;
      out << "        const double* " << restrictKey() << " base0," << std::endl; 
      for( int i=1; i<dim; ++ i ) 
        out << "        const double* " << restrictKey() << " base" << i << "," << std::endl; 
      for( int i=0; i<dim; ++i ) 
      {
        out << "        ";
        for( int r=0; r<dimRange; ++ r ) 
          out << "const double jacFactorInv"<< i << r << ", ";
        out << std::endl;
      }
      for( int r = 0; r < dimRange; ++r ) 
      {
        out << "        double* " << restrictKey() << " result" << r;
        if( r == dimRange-1 ) 
          out << " )" << std::endl;
        else 
          out << "," << std::endl;
      }
      writePreCompHeader( out, 1 );
      out << "  {" << std::endl;
      writeInnerLoopAxpyJac( out, dim, dimRange, numCols ); 
      out << "  }" << std::endl;
      writePreCompHeader( out, 2 );
    }
  };

  // if this pre processor variable is defined then 
  // we assume that CODEGENERATOR_REPLACEMENT is CodeGenerator of choice 
#ifndef CODEGEN_SIMD_WIDTH 
#define CODEGEN_SIMD_WIDTH 2 
#endif
#define FEM_CODEGENERATOR_REPLACEMENT  VectorCodeGenerator< CODEGEN_SIMD_WIDTH >
} // end namespace Fem 

} // end namespace Dune

#include <dune/fem/space/basefunctions/codegen.hh>
#endif
