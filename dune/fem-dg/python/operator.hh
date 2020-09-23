#ifndef DUNE_FEMDG_PYTHON_OPERATOR_HH
#define DUNE_FEMDG_PYTHON_OPERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/fempy/py/space.hh>
#include <dune/fempy/py/operator.hh>
#include <dune/fempy/pybind11/pybind11.hh>

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
      typedef typename Operator::AdvectionFluxType AdvectionFluxType;
      typedef typename DF::DiscreteFunctionSpaceType DFSpace;
      typedef typename Fem::SpaceOperatorInterface<DF> Base;
      typedef Base FullType;
      typedef Base ExplType;
      typedef Base ImplType;
      Dune::FemPy::detail::registerOperator< Operator >( module, cls );
      cls.def( pybind11::init( [] ( const DFSpace &space,
               const MA &advectionModel,
               const MD &diffusionModel,
               const pybind11::dict &parameters )
      {
        return new Operator(space, advectionModel, diffusionModel, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );
      } ), "space"_a, "advectionModel"_a, "diffusionModel"_a, "parameters"_a,
           pybind11::keep_alive< 1, 2 >(),
           pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >(),
           pybind11::keep_alive< 1, 5 >()
           );
      cls.def( pybind11::init( [] ( const DFSpace &space,
               const MA &advectionModel,
               const MD &diffusionModel )
      {
        return new Operator(space, advectionModel, diffusionModel);
      } ), "space"_a, "advectionModel"_a, "diffusionModel"_a,
           pybind11::keep_alive< 1, 2 >(),
           pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >()
           );
      cls.def( pybind11::init( [] ( const DFSpace &space,
               const MA &advectionModel,
               const MD &diffusionModel,
               const AdvectionFluxType &advectionFlux,
               const pybind11::dict &parameters )
      {
        return new Operator(space, advectionModel, diffusionModel, advectionFlux, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );
      } ), "space"_a, "advectionModel"_a, "diffusionModel"_a, "advectionFlux"_a, "parameters"_a,
           pybind11::keep_alive< 1, 2 >(),
           pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >(),
           pybind11::keep_alive< 1, 5 >(),
           pybind11::keep_alive< 1, 6 >()
           );
      cls.def( pybind11::init( [] ( const DFSpace &space,
               const MA &advectionModel,
               const MD &diffusionModel,
               const AdvectionFluxType &advectionFlux )
      {
        return new Operator(space, advectionModel, diffusionModel, advectionFlux);
      } ), "space"_a, "advectionModel"_a, "diffusionModel"_a, "advectionFlux"_a,
           pybind11::keep_alive< 1, 2 >(),
           pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >(),
           pybind11::keep_alive< 1, 5 >()
           );
      cls.def( "applyLimiter", []( Operator &self, DF &u) { self.applyLimiter(u); } );
      // cls.def( "setTime", &Operator::setTime);
      cls.def( "_setTime", &Operator::setTime);
      cls.def_property_readonly("localTimeStepEstimate",
          [](const Operator &self) -> std::tuple<double,double,double>
          { return {self.timeStepEstimate(),
                    self.explicitOperator().timeStepEstimate(),
                    self.implicitOperator().timeStepEstimate()};
          });
      Dune::Python::insertClass<ExplType>(cls,"ExplType",
          Dune::Python::GenerateTypeName("NotAvailable"),
          Dune::Python::IncludeFiles{});
      Dune::Python::insertClass<ImplType>(cls,"ImplType",
          Dune::Python::GenerateTypeName("NotAvailable"),
          Dune::Python::IncludeFiles{});
      cls.def_property_readonly("fullOperator", [](const Operator &self) -> const FullType&
           { return self.fullOperator(); } );
      cls.def_property_readonly("explicitOperator", [](const Operator &self) -> const ExplType&
           { return self.explicitOperator(); } );
      cls.def_property_readonly("implicitOperator", [](const Operator &self) -> const ImplType&
           { return self.implicitOperator(); } );
      cls.def("counter", [](const Operator &self) -> std::tuple<int,int,int>
           { return self.counter(); } );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMDG_PYTHON_OPERATOR_HH
