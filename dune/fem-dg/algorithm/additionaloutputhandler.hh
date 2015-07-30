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

  template< class ConsDiscreteFunctionImp, class PrimDiscreteFunctionImp >
  class Cons2PrimCalculator
  {
    public:

    typedef ConsDiscreteFunctionImp                                       InDiscreteFunctionType;
    typedef typename InDiscreteFunctionType::DiscreteFunctionSpaceType    InDiscreteFunctionSpaceType;
    typedef PrimDiscreteFunctionImp                                       OutDiscreteFunctionType;
    typedef typename OutDiscreteFunctionType::DiscreteFunctionSpaceType   OutDiscreteFunctionSpaceType;

    typedef typename InDiscreteFunctionSpaceType::GridPartType  GridPartType;
    typedef typename InDiscreteFunctionSpaceType::IteratorType  Iterator;
    typedef typename Iterator::Entity                           Entity;
    typedef typename Entity::Geometry                           Geometry;
    typedef typename InDiscreteFunctionSpaceType::DomainType    DomainType;
    typedef typename InDiscreteFunctionSpaceType::RangeType     InRangeType;
    typedef typename OutDiscreteFunctionSpaceType::RangeType    OutRangeType;

    typedef typename InDiscreteFunctionType::LocalFunctionType  InLocalFuncType;
    typedef typename OutDiscreteFunctionType::LocalFunctionType OutLocalFuncType;

    typedef std::tuple< OutDiscreteFunctionType* >              TupleType;

    Cons2PrimCalculator( const OutDiscreteFunctionSpaceType& space )
      : outDf_( "_prim", space )
    {}

    /** \brief converts a discrete function of conservative variables to
     *    a discrete function of primitive variables for a visualization purpose only
     *
     *  \param[in] consDF The discrete function of conservative variables
     *  \param[in] model The analytical model
     *  \param[out] primDF The discrete function of primitive variables
     */
    template< class TimeProvider,
              class ModelType >
    void setup( const InDiscreteFunctionType& consDF,
                const TimeProvider& tp,
                const ModelType& model )
    {
      const InDiscreteFunctionSpaceType& space =  consDF.space();
      outDf_.clear();

      InRangeType cons(0.0);
      InRangeType cons_bg(0.0);
      OutRangeType prim(0.0);
      OutRangeType prim_bg(0.0);

      const Iterator endit = space.end();
      for( Iterator it = space.begin(); it != endit ; ++it)
      {
        // get entity
        const Entity& entity = *it ;
        const Geometry& geo = entity.geometry();

        // get quadrature rule for L2 projection
        Dune::Fem::CachingQuadrature< GridPartType, 0 > quad( entity, 2*space.order()+1 );

        InLocalFuncType consLF = consDF.localFunction( entity );
        OutLocalFuncType primLF = outDf_.localFunction( entity );

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

    OutDiscreteFunctionType& result()
    {
      return std::make_tuple( &outDf_ );
    }

    private:
    OutDiscreteFunctionType     outDf_;

  };

  class NoAdditionalOutputHandler
  {
  public:
    typedef std::tuple<>                   TupleType;

    template< class ... Args >
    NoAdditionalOutputHandler( Args& ... )
    {}

    template< class ... Args >
    void setup( Args&... )
    {}

    template< class ... Args >
    TupleType& result( Args&... )
    {
      return tuple_;
    }

  private:
    TupleType tuple_;

  };

}
}

#endif
