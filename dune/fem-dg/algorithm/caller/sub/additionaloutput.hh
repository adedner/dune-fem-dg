#ifndef FEMDG_SUBNOADDITIONALOUTPUT_HH
#define FEMDG_SUBNOADDITIONALOUTPUT_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>


namespace Dune
{
namespace Fem
{

  class NoOutput
  {
  public:
    template< class... Args >
    NoOutput( Args&& ...){}

    template< class... Args >
    void step( Args&& ...  ){}
  };



  template< class DiscreteFunctionImp >
  class Cons2PrimOutput
  {
    public:

    Cons2PrimOutput( const DiscreteFunctionImp& df )
      : solution_( df )
    {}

    /** \brief converts a discrete function of conservative variables to
     *    a discrete function of primitive variables for a visualization purpose only
     */
    template< class TimeProviderImp, class SubAlgorithmImp >
    void step( const TimeProviderImp& tp, const SubAlgorithmImp& alg )
    {
      typedef typename SubAlgorithmImp::DiscreteFunctionType                InDiscreteFunctionType;
      typedef typename InDiscreteFunctionType::DiscreteFunctionSpaceType    InDiscreteFunctionSpaceType;
      typedef DiscreteFunctionImp                                           OutDiscreteFunctionType;
      typedef typename OutDiscreteFunctionType::DiscreteFunctionSpaceType   OutDiscreteFunctionSpaceType;

      typedef typename InDiscreteFunctionSpaceType::GridPartType  GridPartType;
      typedef typename InDiscreteFunctionSpaceType::RangeType     InRangeType;
      typedef typename OutDiscreteFunctionSpaceType::RangeType    OutRangeType;

      const auto& space =  alg.solution().space();
      solution_.clear();

      InRangeType cons(0.0);
      InRangeType cons_bg(0.0);
      OutRangeType prim(0.0);
      OutRangeType prim_bg(0.0);

      for( const auto& entity : elements( space.gridPart() ) )
      {
        const auto& geo = entity.geometry();

        // get quadrature rule for L2 projection
        Dune::Fem::CachingQuadrature< GridPartType, 0 > quad( entity, 2*space.order()+1 );

        const auto& consLF = alg.solution().localFunction( entity );
        const auto& primLF = solution_.localFunction( entity );

        for( const auto qp : quad )
        {
          const auto& xgl = geo.global( qp.position() );
          consLF.evaluate( qp, cons );

          // it is useful to visualize better suited quantities
          bool forVisual = true;
          alg.model().conservativeToPrimitive( tp.time(), xgl, cons, prim, forVisual );

          prim *=  qp.weight();
          primLF.axpy( qp, prim );
        }
      }
    }
  private:
    DiscreteFunctionImp&     solution_;
  };


  template< class DiscreteFunctionImp >
  class ExactSolutionOutput
  {
    public:
    typedef DiscreteFunctionImp                                                      DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::GridPartType   GridPartType;

    ExactSolutionOutput( DiscreteFunctionType& df, const std::string keyPrefix = "" )
      : solution_( &df )
    {}

    template< class... Args >
    ExactSolutionOutput( Args&&... )
      : solution_( nullptr )
    {}

    template< class TimeProviderImp, class SubAlgorithmImp >
    void step( TimeProviderImp& tp, SubAlgorithmImp& alg )
    {
      if( solution_ )
      {
        typedef typename SubAlgorithmImp::ProblemType::InstationaryFunctionType      ExactSolutionType;
        typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType >  GridFunctionAdapterType;
        auto ftf = alg.problem().fixedTimeFunction( tp.time() );
        GridFunctionAdapterType adapter( "temporary adapter", ftf , solution_->space().gridPart(), solution_->space().order()+2 );
        interpolate( adapter, *solution_ );
      }
    }

    private:
    DiscreteFunctionType* solution_;
  };



  template< class Obj >
  class AdditionalOutputOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    AdditionalOutputOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

  template<>
  class AdditionalOutputOptional< void >
    : public OptionalNullPtr< NoOutput >
  {
    typedef OptionalNullPtr< NoOutput >    BaseType;
  public:
    template< class... Args >
    AdditionalOutputOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

}
}

#endif
