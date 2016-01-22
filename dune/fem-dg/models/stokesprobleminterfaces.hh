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
namespace Stokes
{

  template< class GridImp>
  class ProblemInterfaceBase
  {
  public:
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, GridImp::dimension > FunctionSpaceType;
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, 1 > PressureFunctionSpaceType;

    typedef ProblemInterface< FunctionSpaceType >         PoissonProblemType;
    typedef ProblemInterface< PressureFunctionSpaceType > StokesProblemType;
  };


  template< class GridImp >
  class ProblemInterface
  {
    typedef ProblemInterfaceBase< GridImp >        BaseType;
  public:

    typedef typename BaseType::FunctionSpaceType         FunctionSpaceType;
    typedef typename BaseType::PressureFunctionSpaceType PressureFunctionSpaceType;

    typedef typename BaseType::PoissonProblemType        PoissonProblemType;
    typedef typename BaseType::StokesProblemType         StokesProblemType;

    typedef std::tuple< PoissonProblemType*, StokesProblemType* >         ProblemTupleType;

    /**
     *  \brief constructor constructing a combined problem of the interface sub problems,
     *  i.e. the poisson and the stokes problem.
     *
     *  \note Use the StokesProblem class to create derived objects.
     */
    ProblemInterface()
      : problems_( std::make_tuple( new PoissonProblemType(), new StokesProblemType() ) )
    {}

    template< int i >
    const typename std::remove_pointer< typename std::tuple_element<i,ProblemTupleType>::type >::type& get() const
    {
      return *(std::get<i>( problems_) );
    }

    template< int i >
    typename std::remove_pointer< typename std::tuple_element<i,ProblemTupleType>::type >::type& get()
    {
      return *(std::get<i>( problems_) );
    }

    // return prefix for data loops
    virtual std::string dataPrefix() const
    {
      return get<0>().dataPrefix();
    }

  protected:
    template< class PoissonProblemPtrImp, class StokesProblemPtrImp >
    void create( const PoissonProblemPtrImp& poisson, const StokesProblemPtrImp& stokes )
    {
      std::get<0>( problems_ ) = poisson;
      std::get<1>( problems_ ) = stokes;
    }

  private:
    mutable ProblemTupleType   problems_;
  };



  /**
   * \brief helper class which helps for the correct (virtual) construction
   * of the problem tuple.
   *
   * \tparam GridImp type of the unterlying grid
   * \tparam StokesProblemImp type of the stokes problem
   *
   * \ingroup Problems
   */
  template< class GridImp,
            template<class> class StokesProblemImp >
  class Problem
    : public ProblemInterface< GridImp >
  {
    typedef ProblemInterface< GridImp >                                      BaseType;

    typedef typename BaseType::PoissonProblemType                            PoissonProblemBaseType;
    typedef typename BaseType::StokesProblemType                             StokesProblemBaseType;

  public:
    typedef typename StokesProblemImp<GridImp>::PoissonProblemType           PoissonProblemType;
    typedef typename StokesProblemImp<GridImp>::StokesProblemType            StokesProblemType;

    Problem()
      : BaseType(),
        poisson_( new PoissonProblemType() ),
        stokes_( new StokesProblemType() )
    {
      BaseType::create( poisson_, stokes_ );
    }

  private:
    mutable PoissonProblemBaseType* poisson_;
    mutable StokesProblemBaseType*  stokes_;

  };


}
}
}
#endif  /*DUNE_PROBLEMINTERFACE_HH*/
