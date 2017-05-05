#ifndef FEMDG_DGEOCERROR_HH
#define FEMDG_DGEOCERROR_HH

#include <dune/fem/misc/femeoc.hh>
#include <dune/fem-dg/misc/dgnorm.hh>
#include <dune/fem/misc/femeoc.hh>
namespace Dune
{
namespace Fem
{

  class DGEOCError
  {
  public:
    template< class DiscreteFunctionImp >
    DGEOCError( const DiscreteFunctionImp& u )
    : name_( "DG-error("+u.name()+")" )
    {}

    template< class Model, class Solution >
    static void add ( Model& model, Solution &u, int id )
    {
      Dune::Fem::DGNorm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
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
