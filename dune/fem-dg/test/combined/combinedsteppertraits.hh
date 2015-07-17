#ifndef DUNE_COMBINEDSTEPPERTRAITS_HH
#define DUNE_COMBINEDSTEPPERTRAITS_HH

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/solver/odesolver.hh>
#include <dune/fem/solver/pardginverseoperators.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>

#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>

template <typename A, typename B>
auto func(A a, B b) -> decltype(std::tuple_cat(a,b))
{ return std::tuple_cat(a,b); }

template <typename A, typename B>
struct tuple_appender {
    typedef decltype(func(std::declval<A>(), std::declval<B>())) type;
};

template <class Stepper1,
          class Stepper2>
struct CombinedStepperTraits
{
  typedef typename Stepper1::Traits Traits1;
  typedef typename Stepper2::Traits Traits2;

  //todo: check wether types are the same for both steppers
  typedef typename Traits1::GridType GridType;
  typedef typename Traits1::GridPartType GridPartType;

  typedef Dune::Fem::RestrictProlongTuple< typename Traits1::RestrictionProlongationType,
                                           typename Traits2::RestrictionProlongationType > RestrictionProlongationType;

  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType ;

  typedef typename tuple_appender<typename Traits1::IOTupleType, typename Traits2::IOTupleType>::type IOTupleType;

};

#endif
