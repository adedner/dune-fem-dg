#ifndef CONS2PRIM_HH
#define CONS2PRIM_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>

#include <dune/fem/quadrature/cachingquadrature.hh>

namespace Dune
{
namespace Fem
{
  /** \brief converts a discrete function of conservative variables to
   *    a discrete function of primitive variables for a visualization purpose only
   *
   *  \param[in] consDF The discrete function of conservative variables
   *  \param[in] model The analytical model
   *  \param[out] primDF The discrete function of primitive variables
   */
  template< class TimeProvider,
            class ConsDiscreteFunctionType,
            class ModelType,
            class PrimDiscreteFunctionType >
  void setupAdditionalVariables( const TimeProvider& tp,
                                 const ConsDiscreteFunctionType& consDF,
                                 const ModelType& model,
                                 PrimDiscreteFunctionType& primDF )
  {
    typedef typename ConsDiscreteFunctionType::Traits::DiscreteFunctionSpaceType
      ConsDiscreteFunctionSpaceType;
    typedef typename PrimDiscreteFunctionType::Traits::DiscreteFunctionSpaceType
      PrimDiscreteFunctionSpaceType;
    typedef typename ConsDiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename ConsDiscreteFunctionSpaceType::IteratorType Iterator;
    typedef typename Iterator :: Entity Entity;
    typedef typename Entity :: Geometry Geometry;
    typedef typename ConsDiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename ConsDiscreteFunctionSpaceType::RangeType ConsRangeType;
    typedef typename PrimDiscreteFunctionSpaceType::RangeType PrimRangeType;

    const ConsDiscreteFunctionSpaceType& space =  consDF.space();

    primDF.clear();

    typedef typename ConsDiscreteFunctionType::LocalFunctionType ConsLocalFuncType;
    typedef typename PrimDiscreteFunctionType::LocalFunctionType PrimLocalFuncType;

    ConsRangeType cons(0.0);
    ConsRangeType cons_bg(0.0);
    PrimRangeType prim(0.0);
    PrimRangeType prim_bg(0.0);

    const Iterator endit = space.end();
    for( Iterator it = space.begin(); it != endit ; ++it)
    {
      // get entity
      const Entity& entity = *it ;
      const Geometry& geo = entity.geometry();

      // get quadrature rule for L2 projection
      Dune::Fem::CachingQuadrature< GridPartType, 0 > quad( entity, 2*space.order()+1 );

      ConsLocalFuncType consLF = consDF.localFunction( entity );
      PrimLocalFuncType primLF = primDF.localFunction( entity );

      const int quadNop = quad.nop();
      for(int qP = 0; qP < quadNop; ++qP)
      {
        const DomainType& xgl = geo.global( quad.point(qP) );
        consLF.evaluate( quad[qP], cons );

        // it is useful to visualize better suited quantities
        bool forVisual = true;
        model.conservativeToPrimitive( tp.time(), xgl, cons, prim, forVisual );

        prim *=  quad.weight(qP);
        primLF.axpy( quad[qP] , prim );
      }
    }
  }

}
}

#endif
