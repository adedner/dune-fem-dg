#pragma once

#include <dune/fem/function/localfunction/const.hh>
namespace Dune
{
  namespace Fem
  {
    template< class DiscreteFunction >
    struct TroubledCellIndicatorBase
    {
      virtual ~TroubledCellIndicatorBase() {}

      typedef Dune::Fem::ConstLocalFunction< DiscreteFunction >  LocalFunctionType;
      virtual double operator()( const DiscreteFunction& U, const LocalFunctionType& uEn ) const = 0;
    };
  }
}
