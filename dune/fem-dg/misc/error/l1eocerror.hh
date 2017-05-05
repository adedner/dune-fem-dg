#ifndef FEMDG_L1EOCERROR_HH
#define FEMDG_L1EOCERROR_HH

#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/l1norm.hh>
#include <dune/fem/misc/femeoc.hh>

namespace Dune
{
namespace Fem
{

  class L1EOCError
  {
  public:
    template< class DiscreteFunctionImp >
    L1EOCError( const DiscreteFunctionImp& u )
    : name_( "$L^1$-error("+u.name()+")" )
    {}

    template< class Model, class Solution >
    static void add ( Model& model, Solution &u, int id )
    {
      Dune::Fem::L1Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
      const double error = norm.distance( model.problem().exactSolution( model.time() ), u );
      Dune::Fem::FemEoc::setErrors( id, error );
    }

    const std::string name() const
    {
      return name_;
    }

  private:
    const std::string name_;
  };

}
}

#endif
