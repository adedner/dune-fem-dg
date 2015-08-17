#ifndef FEMDG_L1EOCERROR_HH
#define FEMDG_L1EOCERROR_HH

#include <dune/fem/misc/l1norm.hh>


class L1EOCError
{
public:
  L1EOCError( const std::string name = "$L^1$-error" )
  : id_( Dune::Fem::FemEoc::addEntry( name ) )
  {}

  template< class Solution, class Model, class ExactFunction, class TimeProvider >
  void add( TimeProvider& tp, Solution &u, Model &model, ExactFunction &f )
  {
    Dune::Fem::L1Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
    const double error = norm.distance( model.problem().fixedTimeFunction( tp.time() ), u );
    Dune::Fem::FemEoc::setErrors( id_, error );
  }

  template< class Solution, class ExactFunction >
  void add ( Solution &u, const ExactFunction &f )
  {
    Dune::Fem::H1Norm< typename Solution::DiscreteFunctionSpaceType::GridPartType > norm( u.space().gridPart() );
    const double error = norm.distance( f, u );
    Dune::Fem::FemEoc::setErrors( id_, error );
  }

private:
  int id_;
};

#endif
