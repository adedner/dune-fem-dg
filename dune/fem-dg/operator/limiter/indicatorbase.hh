#pragma once
namespace Dune
{
  namespace Fem
  {
    template< class DiscreteFunction >
    struct TroubledCellIndicatorBase
    {
      virtual ~TroubledCellIndicatorBase() {}

      typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;
      virtual double operator()( const DiscreteFunction& U, const LocalFunctionType& uEn ) const = 0;
    };
  }
}
