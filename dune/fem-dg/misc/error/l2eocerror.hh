#ifndef FEMDG_L2EOCERROR_HH
#define FEMDG_L2EOCERROR_HH

#include <dune/fem/misc/l2norm.hh>


class L2EOCError
{
public:
  L2EOCError( const std::string name = "$L^2$-error" )
  : id_( Dune::Fem::FemEoc::addEntry( name ) )
  {}

  template< class Solution, class Model, class ExactFunction, class TimeProvider >
  void add ( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
  {
    Dune::Fem::L2Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
    const double error = norm.distance( model.problem().fixedTimeFunction( tp.time() ), u );
    Dune::Fem::FemEoc::setErrors( id_, error );
  }

  template< class Solution, class ExactFunction >
  void add ( Solution &u, const ExactFunction &f )
  {
    Dune::Fem::L2Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
    const double error = norm.distance( f, u );
    Dune::Fem::FemEoc::setErrors( id_, error );
  }

private:
  int id_;
};

#endif
