#ifndef FEMDG_DGEOCERROR_HH
#define FEMDG_DGEOCERROR_HH

#include <dune/fem/misc/femeoc.hh>
#include <dune/fem-dg/misc/dgnorm.hh>


class DGEOCError
{
public:
  DGEOCError( const std::string name = "DG-error" )
  : id_( Dune::Fem::FemEoc::addEntry( name ) )
  {}

  template< class Solution, class Model, class ExactFunction, class TimeProvider >
  void add ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
  {
    Dune::Fem::DGNorm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
    const double error = norm.distance( model.problem().fixedTimeFunction( tp.time() ), u );
    Dune::Fem::FemEoc::setErrors( id_, error );
  }

  template< class Solution, class ExactFunction >
  void add ( Solution &u, const ExactFunction &f )
  {
    Dune::Fem::DGNorm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
    const double error = norm.distance( f, u );
    Dune::Fem::FemEoc::setErrors( id_, error );
  }

private:
  int id_;
};

#endif
