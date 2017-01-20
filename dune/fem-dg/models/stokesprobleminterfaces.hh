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
  class ProblemInterfaceTraits
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
    typedef ProblemInterfaceTraits< GridImp >            Traits;
  public:

    typedef typename Traits::FunctionSpaceType           FunctionSpaceType;
    typedef typename Traits::PressureFunctionSpaceType   PressureFunctionSpaceType;

    typedef typename Traits::PoissonProblemType          PoissonProblemType;
    typedef typename Traits::StokesProblemType           StokesProblemType;

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

    ProblemInterface( PoissonProblemType* poisson, StokesProblemType* stokes )
      : problems_( std::make_tuple( poisson, stokes ) )
    {}

    template< int i >
    const std::remove_pointer_t< std::tuple_element_t<i,ProblemTupleType> >& get() const
    {
      return *(std::get<i>( problems_) );
    }

    template< int i >
    std::remove_pointer_t< std::tuple_element_t<i,ProblemTupleType> >& get()
    {
      return *(std::get<i>( problems_) );
    }

    virtual ~ProblemInterface()
    {
      delete std::get<0>( problems_ );
      delete std::get<1>( problems_ );
    }

    const typename StokesProblemType::ExactSolutionType& exactSolution( const double time=0.0 ) const
    {
      return std::get<1>( problems_ )->exactSolution();
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
      : BaseType( new PoissonProblemType(), new StokesProblemType() )
    {}

  };


}
}
}
#endif  /*DUNE_PROBLEMINTERFACE_HH*/
