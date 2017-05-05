#ifndef FEMDG_L2EOCERROR_HH
#define FEMDG_L2EOCERROR_HH

#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/femeoc.hh>

namespace Dune
{
namespace Fem
{

  class L2EOCError
  {
  public:
    template< class DiscreteFunctionImp >
    L2EOCError( const DiscreteFunctionImp& u )
    : name_( "$L^2$-error("+u.name()+")" )
    {}

    template< class Model, class Solution >
    static void add ( Model& model, Solution &u, int id )
    {
      Dune::Fem::L2Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
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
