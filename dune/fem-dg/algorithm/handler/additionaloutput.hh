#ifndef CONS2PRIM_HH
#define CONS2PRIM_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem-dg/misc/optional.hh>



namespace Dune
{

namespace Fem
{

  template< class DiscreteFunctionImp >
  class Cons2PrimOutputHandler
  {
    public:

    Cons2PrimOutputHandler( const DiscreteFunctionImp& df )
      : solution_( df )
    {}

    /** \brief converts a discrete function of conservative variables to
     *    a discrete function of primitive variables for a visualization purpose only
     */
    template< class TimeProviderImp, class SubStepperImp >
    void step( const TimeProviderImp& tp, const SubStepperImp& stepper )
    {
      typedef typename SubStepperImp::DiscreteFunctionType                  InDiscreteFunctionType;
      typedef typename InDiscreteFunctionType::DiscreteFunctionSpaceType    InDiscreteFunctionSpaceType;
      typedef DiscreteFunctionImp                                           OutDiscreteFunctionType;
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

      const InDiscreteFunctionSpaceType& space =  stepper.solution().space();
      solution_.clear();

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

        InLocalFuncType consLF = stepper.solution().localFunction( entity );
        OutLocalFuncType primLF = solution_.localFunction( entity );

        const int quadNop = quad.nop();
        for(int qP = 0; qP < quadNop; ++qP)
        {
          const DomainType& xgl = geo.global( quad.point(qP) );
          consLF.evaluate( quad[qP], cons );

          // it is useful to visualize better suited quantities
          bool forVisual = true;
          stepper.model().conservativeToPrimitive( tp.time(), xgl, cons, prim, forVisual );

          prim *=  quad.weight(qP);
          primLF.axpy( quad[qP] , prim );
        }
      }
    }
  private:
    DiscreteFunctionImp&     solution_;
  };


  template< class DiscreteFunctionImp >
  class ExactSolutionOutputHandler
  {
    public:
    typedef DiscreteFunctionImp                                                      DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::GridPartType   GridPartType;

    ExactSolutionOutputHandler( DiscreteFunctionType& df, const std::string keyPrefix = "" )
      : solution_( df )
    {}

    template< class TimeProviderImp, class SubStepperImp >
    void step( TimeProviderImp& tp, SubStepperImp& stepper )
    {
      typedef typename SubStepperImp::ProblemType::InstationaryFunctionType      ExactSolutionType;
      typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridFunctionAdapterType;
      auto ftf = stepper.problem().fixedTimeFunction( tp.time() );
      GridFunctionAdapterType adapter( "temporary adapter", ftf , stepper.gridPart(), solution_.space().order()+2 );
      interpolate( adapter, solution_ );
    }

    private:
    DiscreteFunctionType& solution_;
  };

  class NoOutputHandler
  {
  public:
    template< class... Args >
    NoOutputHandler( Args&& ...){}

    template< class... Args >
    void step( Args&& ...  ){}
  };

  template< class Obj >
  class AdditionalOutputHandlerOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    AdditionalOutputHandlerOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

  template<>
  class AdditionalOutputHandlerOptional< void >
    : public OptionalNullPtr< NoOutputHandler >
  {
    typedef OptionalNullPtr< NoOutputHandler >    BaseType;
  public:
    template< class... Args >
    AdditionalOutputHandlerOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

}
}

#endif
