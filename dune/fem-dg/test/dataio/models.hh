#ifndef DUNE_FEM_DG_DATAIOMODEL_HH
#define DUNE_FEM_DG_DATAIOMODEL_HH

#include <dune/fem-dg/models/defaultmodel.hh>

namespace Dune
{
namespace Fem
{
  template <class GridPart, class ProblemImp >
  class DefaultNoModelTraits
    : public ProblemImp::FunctionSpaceType
  {
    typedef typename ProblemImp::FunctionSpaceType                        BaseType;

    typedef GridPart                                                      GridPartType;
    typedef typename GridPartType::GridType                               GridType;
  public:
    typedef typename GridType::template Codim< 0 >::Entity                EntityType;
    typedef typename GridPartType::IntersectionIteratorType::Intersection IntersectionType;

    typedef typename BaseType::RangeFieldType                             RangeFieldType;
    typedef typename BaseType::DomainFieldType                            DomainFieldType;
    enum { dimRange  = BaseType::dimRange };
    enum { dimDomain = BaseType::dimDomain };

    typedef Dune::FieldVector< DomainFieldType, dimDomain-1 >             FaceDomainType;
    typedef Dune::FieldVector< RangeFieldType, dimDomain * dimRange >     GradientType;
    typedef typename BaseType::JacobianRangeType                          JacobianRangeType;
    typedef typename BaseType::JacobianRangeType                          FluxRangeType;

  };


  template< class GridPartImp, class ProblemImp >
  class NoModel
    : public DefaultModel< DefaultNoModelTraits< GridPartImp, ProblemImp > >
  {

    public:

    NoModel( const ProblemImp& problem )
    {}

  };

}
}

#endif
