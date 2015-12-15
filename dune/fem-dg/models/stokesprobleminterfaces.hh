#ifndef DUNE_FEM_DG_STOKES_PROBLEMINTERFACE_HH
#define DUNE_FEM_DG_STOKES_PROBLEMINTERFACE_HH

#include <dune/common/version.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/fem/misc/gridsolution.hh>

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief describes the interface for a stokes problem
   *
   * \tparam FunctionSpaceImp type of the discrete function space describing the velocity
   * \tparam PressureSpaceImp type of the (scalar) discrete function space describing the pressure
   *
   * \ingroup Problems
   */
  template <class PoissonProblemImp, class StokesProblemImp>
  class StokesProblemInterface
  {
  public:
    typedef PoissonProblemImp                            PoissonProblemType;
    typedef StokesProblemImp                             StokesProblemType;

    typedef std::tuple< PoissonProblemType, StokesProblemType >         ProblemTupleType;

    template< class ProblemTupleImp >
    StokesProblemInterface( ProblemTupleImp problems )
      : problems_( problems )
    {}

    template< int i >
    const typename std::tuple_element<i,ProblemTupleType>::type& get() const
    {
      return std::get<i>( problems_);
    }

  private:
    ProblemTupleType   problems_;
  };
}
}
#endif  /*DUNE_PROBLEMINTERFACE_HH*/
