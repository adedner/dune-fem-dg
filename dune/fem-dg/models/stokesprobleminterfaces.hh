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

  template< class GridImp>
  class StokesProblemInterfaceBase
  {
  public:
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, GridImp::dimension > FunctionSpaceType;
    typedef Dune::Fem::FunctionSpace< double, double, GridImp::dimension, 1 > PressureFunctionSpaceType;

    typedef ProblemInterface< FunctionSpaceType >         PoissonProblemType;
    typedef ProblemInterface< PressureFunctionSpaceType > StokesProblemType;
  };


  template< class GridImp >
  class StokesProblemInterface
  {
  public:
    typedef StokesProblemInterfaceBase< GridImp >        BaseType;

    typedef typename BaseType::FunctionSpaceType         FunctionSpaceType;
    typedef typename BaseType::PressureFunctionSpaceType PressureFunctionSpaceType;

    typedef typename BaseType::PoissonProblemType        PoissonProblemType;
    typedef typename BaseType::StokesProblemType         StokesProblemType;

    typedef std::tuple< PoissonProblemType*, StokesProblemType* >         ProblemTupleType;

    StokesProblemInterface()
      : problems_( std::make_tuple( new PoissonProblemType(), new StokesProblemType() ) )
    {}

    template< class PoissonProblemImp, class StokesProblemImp >
    void create( const PoissonProblemImp& poisson, const StokesProblemImp& stokes )
    {
      std::get<0>( problems_ ) = poisson;
      std::get<1>( problems_ ) = stokes;
    }


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
  class VirtualStokesProblemCreator
    : public StokesProblemInterface< GridImp >
  {
    typedef StokesProblemInterface< GridImp > BaseType;
  public:
    typedef ProblemInterface< Dune::Fem::FunctionSpace< double, double, GridImp::dimension, GridImp::dimension > > PoissonProblemBaseType;
    typedef ProblemInterface< Dune::Fem::FunctionSpace< double, double, GridImp::dimension, 1  > >                 StokesProblemBaseType;

    typedef typename StokesProblemImp<GridImp>::PoissonProblemType                PoissonProblemType;
    typedef typename StokesProblemImp<GridImp>::StokesProblemType                 StokesProblemType;


    typedef std::tuple< PoissonProblemBaseType*, StokesProblemBaseType* >         ProblemTupleType;

    VirtualStokesProblemCreator()
      : BaseType(),
        problems_( std::make_tuple( new PoissonProblemType(), new StokesProblemType() ) )
    {
      BaseType::create( std::get<0>(problems_), std::get<1>(problems_) );
    }

  private:
    mutable ProblemTupleType   problems_;

  };



}
}
#endif  /*DUNE_PROBLEMINTERFACE_HH*/
