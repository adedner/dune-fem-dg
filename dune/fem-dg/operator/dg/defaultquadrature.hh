#ifndef DUNE_FEMDG_DEFAULTQUADRATURE_HH
#define DUNE_FEMDG_DEFAULTQUADRATURE_HH

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/quadrature/defaultquadratures.hh>
#include <dune/fem/quadrature/interpolationquadrature.hh>

namespace Dune
{
  namespace Fem
  {

    template < class Space >
    struct DefaultQuadrature
    {
      template <class F, int d>
      using DefaultQuadratureTraits = Dune::Fem::DefaultQuadratureTraits< F, d >;

      static int volumeOrder( const int polOrder ) {  return 2 * polOrder; }
      static int faceOrder( const int polOrder )   {  return 2 * polOrder + 1; }
    };

#if HAVE_DUNE_LOCALFUNCTIONS
    template < class FunctionSpace, class GridPart, unsigned int order, template< class > class Storage >
    struct DefaultQuadrature< Dune::Fem::FixedOrderDGLagrangeSpace< FunctionSpace, GridPart, order, Dune::GaussLobattoPointSet, Storage > >
    {
      template <class F, int d>
      using DefaultQuadratureTraits = Dune::Fem::GaussLobattoQuadratureTraits< F, d >;

      static int volumeOrder( const int polOrder ) {  return 2 * polOrder - 1; }
      static int faceOrder( const int polOrder )   {  return 2 * polOrder - 1; }
    };

    template < class FunctionSpace, class GridPart, unsigned int order, template< class > class Storage >
    struct DefaultQuadrature< Dune::Fem::FixedOrderDGLagrangeSpace< FunctionSpace, GridPart, order, Dune::GaussLegendrePointSet, Storage > >
    {
      template <class F, int d>
      using DefaultQuadratureTraits = Dune::Fem::GaussLegendreQuadratureTraits< F, d >;

      static int volumeOrder( const int polOrder ) {  return 2 * polOrder + 1; }
      static int faceOrder( const int polOrder )   {  return 2 * polOrder + 1; }
    };
#endif
  }
}

#endif
