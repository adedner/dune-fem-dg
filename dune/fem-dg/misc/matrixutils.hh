#ifndef MATRIXUTILS_HH
#define MATRIXUTILS_HH

#include <string>
#include <iostream>
#include <cmath>

namespace Dune
{
namespace Fem
{
    // NOTE: this is bad coding style
    template< class MatrixImp >
    void testMatrixSymmetry( const MatrixImp& matrix )
    {
      double maxSymError_diagonal = 0;
      double maxSymError_offdiagonal = 0;
      const auto& gridPart = matrix.systemMatrix().rangeSpace().gridPart();
      for( const auto& entity : elements( gridPart ) )
      {
        auto localOpEn = matrix.localMatrix( entity, entity  );
        const auto& baseSet = localOpEn.domainBasisFunctionSet();
        const unsigned int numBasisFunctions = baseSet.size();

        for( unsigned int r = 0; r < numBasisFunctions; ++r )
        {
          for( unsigned int c = 0; c < numBasisFunctions; ++c )
          {
            double v1 = localOpEn.get(r,c);
            double v2 = localOpEn.get(c,r);
             if (std::abs(v1-v2)>1e-8)
              std::cout << "symmetry error on diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
            maxSymError_diagonal = std::max(maxSymError_diagonal,std::abs(v1-v2));
          }
        }
        for (const auto& intersection : intersections( gridPart, entity) )
        {
          if ( !intersection.neighbor() )
            continue;

          const auto& neighbor = intersection.outside();

          auto localOpNb1 = matrix.localMatrix( entity, neighbor );
          auto localOpNb2 = matrix.localMatrix( neighbor, entity );
          for( unsigned int r = 0; r < numBasisFunctions; ++r )
          {
            for( unsigned int c = 0; c < numBasisFunctions; ++c )
            {
              double v1 = localOpNb1.get(r,c);
              double v2 = localOpNb2.get(c,r);
               if (std::abs(v1-v2)>1e-8)
                  std::cout << "symmetry error on off diagonal: " << v1 << " - " << v2 << " = " << v1-v2 << std::endl;
              maxSymError_offdiagonal = std::max(maxSymError_offdiagonal,std::abs(v1-v2));
            }
          }
        }
      }
      if (maxSymError_diagonal > 1e-8)
        std::cout << "non symmetric diagonal: " << maxSymError_diagonal << std::endl;
      if (maxSymError_offdiagonal > 1e-8)
        std::cout << "non symmetric off diagonal: " << maxSymError_offdiagonal << std::endl;
    }

}
}

#endif
