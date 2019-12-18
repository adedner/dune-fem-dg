#ifndef DUNE_FEMDG_PYTHON_OPERATOR_HH
#define DUNE_FEMDG_PYTHON_OPERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/fempy/py/space.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/pybind11/pybind11.hh>

#include <dune/fem-dg/solver/dg.hh>

namespace Dune
{

  namespace FemPy
  {
    template< class DF, class MA, class MD, class Add, class... options >
    inline static void registerOperator ( pybind11::module module,
        pybind11::class_< Fem::DGOperator<DF,MA,MD,Add>, options... > cls )
    {
      using pybind11::operator""_a;
      typedef Fem::DGOperator<DF,MA,MD,Add> Operator;
      typedef typename DF::DiscreteFunctionSpaceType DFSpace;
      typedef typename DF::GridPartType GridPartType;
      typedef typename Fem::SpaceOperatorInterface<DF> Base;
      Dune::FemPy::detail::registerOperator< Operator >( module, cls );
      cls.def( pybind11::init( [] ( const DFSpace &space,
               const MA &advectionModel,
               const MD &diffusionModel,
               const pybind11::dict &parameters )
      {
        return new Operator(space, advectionModel, diffusionModel, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );
      } ), "space"_a, "advectionModel"_a, "diffusionModel"_a, "parameters"_a,
           pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def( "applyLimiter", []( Operator &self, DF &u) { self.limit(u); } );
      cls.def_property_readonly( "fullOperator", [](Operator &self) -> const Base&
          { return self.fullOperator(); } );
      cls.def_property_readonly( "explicitOperator", [](Operator &self) -> const Base&
          { return self.explicitOperator(); } );
      cls.def_property_readonly( "implicitOperator", [](Operator &self) -> const Base&
          { return self.implicitOperator(); } );
      // cls.def( "setTimeStepSize", &Operator::setTimeStepSize );
      // cls.def( "deltaT", &Operator::deltaT );
      // cls.def( "step", &Operator::solve, "target"_a );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMDG_PYTHON_OPERATOR_HH
