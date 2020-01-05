#ifndef DUNE_FEMDG_SOD_HH
#define DUNE_FEMDG_SOD_HH

#include <dune/common/fvector.hh>
#include <dune/fem-dg/examples/euler/problems/chorjo.hh>

template <class DomainType>
Dune::FieldVector<double, DomainType::dimension+2> sod( const double t, const double x0, const DomainType& x )
{
  Dune::FieldVector<double, DomainType::dimension+2> res( 0 );
  return chorin(t,x[0]-x0,res[0],res[1],res[DomainType::dimension+1]);
}

#endif
