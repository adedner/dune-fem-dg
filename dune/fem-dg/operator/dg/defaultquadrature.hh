#ifndef DUNE_FEMDG_DEFAULTQUADRATURE_HH
#define DUNE_FEMDG_DEFAULTQUADRATURE_HH

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/quadrature/defaultquadratures.hh>
#include <dune/fem/quadrature/interpolationquadrature.hh>

namespace Dune
{
  namespace Fem
  {

    struct DefaultQuadratureGauss
    {
      template <class F, int d>
      using DefaultQuadratureTraits = Dune::Fem::DefaultQuadratureTraits< F, d >;

      static int volumeOrder( const int polOrder ) {  return 2 * polOrder; }
      static int faceOrder( const int polOrder )   {  return 2 * polOrder + 1; }
    };

    template < class Space >
    struct DefaultQuadrature : public DefaultQuadratureGauss
    {
    };

#if HAVE_DUNE_LOCALFUNCTIONS
    struct DefaultQuadratureGaussLobatto
    {
      template <class F, int d>
      using DefaultQuadratureTraits = Dune::Fem::GaussLobattoQuadratureTraits< F, d >;

      static int volumeOrder( const int polOrder ) {  return (polOrder > 0) ? (2 * polOrder - 1) : 0; }
      static int faceOrder( const int polOrder )   {  return (polOrder > 0) ? (2 * polOrder - 1) : 0; }
    };

    template < class FunctionSpace, class GridPart, unsigned int order, class Storage >
    struct DefaultQuadrature< Dune::Fem::FixedOrderDGLagrangeSpace< FunctionSpace, GridPart, order, Dune::GaussLobattoPointSet, Storage > >
    : public DefaultQuadratureGaussLobatto
    {
    };

    template< class LFEMap >
    struct DefaultQuadratureSpec : DefaultQuadratureGauss
    {};

    template< class FunctionSpace, class GridPart, unsigned int order >
    struct DefaultQuadratureSpec< FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, Dune::GaussLobattoPointSet > >
    : public DefaultQuadratureGaussLobatto
    {
    };

    struct DefaultQuadratureGaussLegendre
    {
      template <class F, int d>
      using DefaultQuadratureTraits = Dune::Fem::GaussLegendreQuadratureTraits< F, d >;

      static int volumeOrder( const int polOrder ) {  return 2 * polOrder + 1; }
      static int faceOrder( const int polOrder )   {  return 2 * polOrder + 1; }
    };

    template < class FunctionSpace, class GridPart, unsigned int order, class Storage >
    struct DefaultQuadrature< Dune::Fem::FixedOrderDGLagrangeSpace< FunctionSpace, GridPart, order, Dune::GaussLegendrePointSet, Storage > >
      : public DefaultQuadratureGaussLegendre
    {};

    template< class FunctionSpace, class GridPart, unsigned int order >
    struct DefaultQuadratureSpec< FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, Dune::GaussLegendrePointSet > >
    : public DefaultQuadratureGaussLegendre
    {
    };

    template< class LFEMap, class FunctionSpace, class Storage >
    struct DefaultQuadrature< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      : public DefaultQuadratureSpec< LFEMap >
    {
    };

#endif
  }
}

#endif
