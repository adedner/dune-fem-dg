#ifndef DUNE_FEM_DG_DATAIOMODEL_HH
#define DUNE_FEM_DG_DATAIOMODEL_HH

#include <dune/fem-dg/models/defaultmodel.hh>

namespace Dune
{
namespace Fem
{
  template <class GridImp, class ProblemImp >
  class DefaultNoModelTraits
    : public ProblemImp::FunctionSpaceType
  {
    typedef typename ProblemImp::FunctionSpaceType                        BaseType;

    typedef GridImp                                                       GridType;
  public:
    static const int modelParameterSize = 0;
    typedef typename BaseType::RangeFieldType                             RangeFieldType;
    typedef typename BaseType::DomainFieldType                            DomainFieldType;
    enum { dimRange  = BaseType::dimRange };
    enum { dimDomain = BaseType::dimDomain };

    typedef Dune::FieldVector< DomainFieldType, dimDomain-1 >             FaceDomainType;
    typedef Dune::FieldVector< RangeFieldType, dimDomain * dimRange >     GradientType;
    typedef typename BaseType::JacobianRangeType                          JacobianRangeType;
    typedef typename BaseType::JacobianRangeType                          FluxRangeType;

  };


  template< class GridImp, class ProblemImp >
  class NoModel
    : public DefaultModel< DefaultNoModelTraits< GridImp, ProblemImp > >
  {
    typedef DefaultModel< DefaultNoModelTraits< GridImp, ProblemImp > > BaseType;

  public:
    typedef ProblemImp ProblemType;

    NoModel( const ProblemImp& problem )
      : problem_( problem )
    {}

    inline const ProblemType& problem() const { return problem_; }

  private:
    const ProblemType& problem_;
  };

}
}

#endif
