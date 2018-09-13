#ifndef PYTHON_MODEL_WRAPPER_HH
#define PYTHON_MODEL_WRAPPER_HH

#include <cmath>
#include <dune/common/fvector.hh>

// Additional
// ----------

namespace Dune {
namespace Fem {

namespace DGTesting {

template< class FunctionSpace >
struct Additional
{
  typedef typename FunctionSpace::DomainType DomainType;
  typedef typename FunctionSpace::RangeType RangeType;
  typedef typename FunctionSpace::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpace::HessianRangeType HessianRangeType;

  template< class Entity, class Point >
  static void advection ( const double &t, const Entity &entity, const Point &x, const RangeType &u, JacobianRangeType &result )
  {
    (result[ 0 ])[ 0 ] = u[ 0 ] * (u[ 1 ] / u[ 0 ]);
    (result[ 0 ])[ 1 ] = u[ 0 ] * (u[ 2 ] / u[ 0 ]);
    (result[ 1 ])[ 0 ] = 0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))) + (u[ 0 ] * (u[ 1 ] / u[ 0 ])) * (u[ 1 ] / u[ 0 ]);
    (result[ 1 ])[ 1 ] = (u[ 0 ] * (u[ 1 ] / u[ 0 ])) * (u[ 2 ] / u[ 0 ]);
    (result[ 2 ])[ 0 ] = (u[ 0 ] * (u[ 1 ] / u[ 0 ])) * (u[ 2 ] / u[ 0 ]);
    (result[ 2 ])[ 1 ] = 0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))) + (u[ 0 ] * (u[ 2 ] / u[ 0 ])) * (u[ 2 ] / u[ 0 ]);
    (result[ 3 ])[ 0 ] = (u[ 3 ] + 0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ])))))) * (u[ 1 ] / u[ 0 ]);
    (result[ 3 ])[ 1 ] = (u[ 3 ] + 0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ])))))) * (u[ 2 ] / u[ 0 ]);
  }

  template< class Entity, class Point >
  static double maxSpeed ( const double &t, const Entity &entity, const Point &x, const DomainType &normal, const RangeType &u )
  {
    double result;
    using std::abs;
    using std::sqrt;
    result = std::abs( normal[ 0 ] * (u[ 1 ] / u[ 0 ]) + normal[ 1 ] * (u[ 2 ] / u[ 0 ]) ) + std::sqrt( (1.4 * (0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))))) / u[ 0 ] );
    return result;
  }

  template< class Entity, class Point >
  static DomainType velocity ( const double &t, const Entity &entity, const Point &x, const RangeType &u )
  {
    DomainType result;
    result[ 0 ] = u[ 1 ] / u[ 0 ];
    result[ 1 ] = u[ 2 ] / u[ 0 ];
    return result;
  }

  template< class Entity, class Point >
  static double physical ( const Entity &entity, const Point &x, const RangeType &u )
  {
    double result;
    result = (u[ 0 ] > 1e-08 ? (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ])))) > 1e-08 ? 1 : 0) : 0);
    return result;
  }

  template< class Intersection, class Point >
  static double jump ( const Intersection& it, const Point &x, const RangeType &u, const RangeType &w )
  {
    double result;
    result = (0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))) + -1 * (0.3999999999999999 * (w[ 3 ] + -1 * (w[ 0 ] * (0.5 * ((w[ 1 ] / w[ 0 ]) * (w[ 1 ] / w[ 0 ]) + (w[ 2 ] / w[ 0 ]) * (w[ 2 ] / w[ 0 ]))))))) / (0.5 * (0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))) + 0.3999999999999999 * (w[ 3 ] + -1 * (w[ 0 ] * (0.5 * ((w[ 1 ] / w[ 0 ]) * (w[ 1 ] / w[ 0 ]) + (w[ 2 ] / w[ 0 ]) * (w[ 2 ] / w[ 0 ])))))));
    return result;
  }

  template< class Entity, class Point >
  static bool boundaryFlux ( const int bndId, const double &t, const Entity& entity, const Point &x, const DomainType &normal, const RangeType &u, RangeType &result )
  {
    switch( bndId )
    {
    default:
      {
        return false;
      }
    }
  }

  template< class Entity, class Point >
  static bool diffusionBoundaryFlux ( const int bndId, const double &t, const Entity& entity, const Point &x, const DomainType &normal, const RangeType &u, const JacobianRangeType &jac, RangeType &result )
  {
    switch( bndId )
    {
    default:
      {
        return false;
      }
    }
  }

  template< class Entity, class Point >
  static bool boundaryValue ( const int bndId, const double &t, const Entity& entity, const Point &x, const RangeType &u, RangeType &result )
  {
    switch( bndId )
    {
    case 1:
      {
        result[ 0 ] = u[ 0 ];
        result[ 1 ] = u[ 1 ];
        result[ 2 ] = u[ 2 ];
        result[ 3 ] = u[ 3 ];
        return true;
      }
      break;
    case 2:
      {
        result[ 0 ] = u[ 0 ];
        result[ 1 ] = u[ 1 ];
        result[ 2 ] = u[ 2 ];
        result[ 3 ] = u[ 3 ];
        return true;
      }
      break;
    case 3:
      {
        result[ 0 ] = u[ 0 ];
        result[ 1 ] = u[ 1 ];
        result[ 2 ] = u[ 2 ];
        result[ 3 ] = u[ 3 ];
        return true;
      }
      break;
    case 4:
      {
        result[ 0 ] = u[ 0 ];
        result[ 1 ] = u[ 1 ];
        result[ 2 ] = u[ 2 ];
        result[ 3 ] = u[ 3 ];
        return true;
      }
      break;
    default:
      {
        return false;
      }
    }
  }

  template< class LimitedRange >
  static void limitedRange ( LimitedRange& limRange )
  {
    {}
  }
  template< class Entity, class Point >
  static void adjustAverageValue ( const Entity& entity, const Point &x, RangeType &u )
  {
    {}
  }
  static const int limitedDimRange = FunctionSpace :: dimRange;
  static const bool hasAdvection = true;
  static const bool hasDiffusion = false;
  static const bool hasStiffSource = false;
  static const bool hasNonStiffSource = false;
  static const bool hasFlux = true;
  static const bool threading = true;
  static const Dune::Fem::Solver::Enum solverId = Dune::Fem::Solver::Enum::fem;
  static const Dune::Fem::Formulation::Enum formId = Dune::Fem::Formulation::Enum::primal;
  static const Dune::Fem::AdvectionLimiter::Enum limiterId = Dune::Fem::AdvectionLimiter::Enum::limited;
  static const Dune::Fem::AdvectionFlux::Enum advFluxId = Dune::Fem::AdvectionFlux::Enum::llf;
  static const Dune::Fem::DiffusionFlux::Enum diffFluxId = Dune::Fem::DiffusionFlux::Enum::none;
  static const Dune::Fem::AdvectionLimiterFunction::Enum limiterFunctionId = Dune::Fem::AdvectionLimiterFunction::Enum::minmod;
};

// Model
// -----

template< class GridPart >
struct PythonModel
{
  typedef GridPart GridPartType;
  typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPart::IntersectionType IntersectionType;
  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 4 > DFunctionSpaceType;
  typedef typename DFunctionSpaceType::DomainFieldType DDomainFieldType;
  typedef typename DFunctionSpaceType::RangeFieldType DRangeFieldType;
  typedef typename DFunctionSpaceType::DomainType DDomainType;
  typedef typename DFunctionSpaceType::RangeType DRangeType;
  typedef typename DFunctionSpaceType::JacobianRangeType DJacobianRangeType;
  typedef typename DFunctionSpaceType::HessianRangeType DHessianRangeType;
  static const int dimDomain = GridPartType::dimensionworld;
  static const int dimD = 4;
  typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 4 > RFunctionSpaceType;
  typedef typename RFunctionSpaceType::DomainFieldType RDomainFieldType;
  typedef typename RFunctionSpaceType::RangeFieldType RRangeFieldType;
  typedef typename RFunctionSpaceType::DomainType RDomainType;
  typedef typename RFunctionSpaceType::RangeType RRangeType;
  typedef typename RFunctionSpaceType::JacobianRangeType RJacobianRangeType;
  typedef typename RFunctionSpaceType::HessianRangeType RHessianRangeType;
  static const int dimR = 4;
  static const int dimLocal = GridPartType::dimension;

  PythonModel ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  {}

  bool init ( const EntityType &entity ) const
  {
    {
      entity_ = &entity;
    }
    return true;
  }

  const EntityType &entity () const
  {
    return *entity_;
  }

  std::string name () const
  {
    return "Model";
  }
  typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
  static const bool symmetric = false;

  template< class Point >
  void source ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
  {
    result[ 0 ] = 0;
    result[ 1 ] = 0;
    result[ 2 ] = 0;
    result[ 3 ] = 0;
  }

  template< class Point >
  void linSource ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
  {
    result[ 0 ] = 0;
    result[ 1 ] = 0;
    result[ 2 ] = 0;
    result[ 3 ] = 0;
  }

  template< class Point >
  void diffusiveFlux ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
  {
    const auto tmp0 = u[ 1 ] / u[ 0 ];
    const auto tmp1 = u[ 0 ] * tmp0;
    const auto tmp2 = u[ 2 ] / u[ 0 ];
    const auto tmp3 = u[ 0 ] * tmp2;
    const auto tmp4 = tmp0 * tmp0;
    const auto tmp5 = tmp2 * tmp2;
    const auto tmp6 = tmp4 + tmp5;
    const auto tmp7 = 0.5 * tmp6;
    const auto tmp8 = u[ 0 ] * tmp7;
    const auto tmp9 = -1 * tmp8;
    const auto tmp10 = u[ 3 ] + tmp9;
    const auto tmp11 = 0.3999999999999999 * tmp10;
    const auto tmp12 = tmp1 * tmp0;
    const auto tmp13 = tmp11 + tmp12;
    const auto tmp14 = tmp1 * tmp2;
    const auto tmp15 = tmp3 * tmp2;
    const auto tmp16 = tmp11 + tmp15;
    const auto tmp17 = u[ 3 ] + tmp11;
    const auto tmp18 = tmp17 * tmp0;
    const auto tmp19 = tmp17 * tmp2;
    (result[ 0 ])[ 0 ] = tmp1;
    (result[ 0 ])[ 1 ] = tmp3;
    (result[ 1 ])[ 0 ] = tmp13;
    (result[ 1 ])[ 1 ] = tmp14;
    (result[ 2 ])[ 0 ] = tmp14;
    (result[ 2 ])[ 1 ] = tmp16;
    (result[ 3 ])[ 0 ] = tmp18;
    (result[ 3 ])[ 1 ] = tmp19;
  }

  template< class Point >
  void linDiffusiveFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
  {
    const auto tmp0 = ubar[ 1 ] / ubar[ 0 ];
    const auto tmp1 = u[ 0 ] * tmp0;
    const auto tmp2 = -1 * tmp1;
    const auto tmp3 = u[ 1 ] + tmp2;
    const auto tmp4 = tmp3 / ubar[ 0 ];
    const auto tmp5 = ubar[ 0 ] * tmp4;
    const auto tmp6 = tmp1 + tmp5;
    const auto tmp7 = ubar[ 2 ] / ubar[ 0 ];
    const auto tmp8 = u[ 0 ] * tmp7;
    const auto tmp9 = -1 * tmp8;
    const auto tmp10 = u[ 2 ] + tmp9;
    const auto tmp11 = tmp10 / ubar[ 0 ];
    const auto tmp12 = ubar[ 0 ] * tmp11;
    const auto tmp13 = tmp8 + tmp12;
    const auto tmp14 = tmp6 * tmp0;
    const auto tmp15 = ubar[ 0 ] * tmp0;
    const auto tmp16 = tmp15 * tmp4;
    const auto tmp17 = tmp14 + tmp16;
    const auto tmp18 = tmp0 * tmp4;
    const auto tmp19 = tmp18 + tmp18;
    const auto tmp20 = tmp7 * tmp11;
    const auto tmp21 = tmp20 + tmp20;
    const auto tmp22 = tmp19 + tmp21;
    const auto tmp23 = 0.5 * tmp22;
    const auto tmp24 = ubar[ 0 ] * tmp23;
    const auto tmp25 = tmp0 * tmp0;
    const auto tmp26 = tmp7 * tmp7;
    const auto tmp27 = tmp25 + tmp26;
    const auto tmp28 = 0.5 * tmp27;
    const auto tmp29 = u[ 0 ] * tmp28;
    const auto tmp30 = tmp24 + tmp29;
    const auto tmp31 = -1 * tmp30;
    const auto tmp32 = u[ 3 ] + tmp31;
    const auto tmp33 = 0.3999999999999999 * tmp32;
    const auto tmp34 = tmp17 + tmp33;
    const auto tmp35 = tmp6 * tmp7;
    const auto tmp36 = tmp15 * tmp11;
    const auto tmp37 = tmp35 + tmp36;
    const auto tmp38 = tmp13 * tmp7;
    const auto tmp39 = ubar[ 0 ] * tmp7;
    const auto tmp40 = tmp39 * tmp11;
    const auto tmp41 = tmp38 + tmp40;
    const auto tmp42 = tmp41 + tmp33;
    const auto tmp43 = u[ 3 ] + tmp33;
    const auto tmp44 = tmp43 * tmp0;
    const auto tmp45 = ubar[ 0 ] * tmp28;
    const auto tmp46 = -1 * tmp45;
    const auto tmp47 = ubar[ 3 ] + tmp46;
    const auto tmp48 = 0.3999999999999999 * tmp47;
    const auto tmp49 = ubar[ 3 ] + tmp48;
    const auto tmp50 = tmp49 * tmp4;
    const auto tmp51 = tmp44 + tmp50;
    const auto tmp52 = tmp43 * tmp7;
    const auto tmp53 = tmp49 * tmp11;
    const auto tmp54 = tmp52 + tmp53;
    (result[ 0 ])[ 0 ] = tmp6;
    (result[ 0 ])[ 1 ] = tmp13;
    (result[ 1 ])[ 0 ] = tmp34;
    (result[ 1 ])[ 1 ] = tmp37;
    (result[ 2 ])[ 0 ] = tmp37;
    (result[ 2 ])[ 1 ] = tmp42;
    (result[ 3 ])[ 0 ] = tmp51;
    (result[ 3 ])[ 1 ] = tmp54;
  }

  template< class Point >
  void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
  {
    const auto tmp0 = u[ 1 ] / u[ 0 ];
    const auto tmp1 = (du[ 0 ])[ 0 ] * tmp0;
    const auto tmp2 = -1 * tmp1;
    const auto tmp3 = (du[ 1 ])[ 0 ] + tmp2;
    const auto tmp4 = tmp3 / u[ 0 ];
    const auto tmp5 = u[ 0 ] * tmp4;
    const auto tmp6 = tmp1 + tmp5;
    const auto tmp7 = u[ 2 ] / u[ 0 ];
    const auto tmp8 = (du[ 0 ])[ 1 ] * tmp7;
    const auto tmp9 = -1 * tmp8;
    const auto tmp10 = (du[ 2 ])[ 1 ] + tmp9;
    const auto tmp11 = tmp10 / u[ 0 ];
    const auto tmp12 = u[ 0 ] * tmp11;
    const auto tmp13 = tmp8 + tmp12;
    const auto tmp14 = tmp6 + tmp13;
    const auto tmp15 = -1 * tmp14;
    const auto tmp16 = tmp6 * tmp0;
    const auto tmp17 = u[ 0 ] * tmp0;
    const auto tmp18 = tmp17 * tmp4;
    const auto tmp19 = tmp16 + tmp18;
    const auto tmp20 = tmp0 * tmp4;
    const auto tmp21 = tmp20 + tmp20;
    const auto tmp22 = (du[ 0 ])[ 0 ] * tmp7;
    const auto tmp23 = -1 * tmp22;
    const auto tmp24 = (du[ 2 ])[ 0 ] + tmp23;
    const auto tmp25 = tmp24 / u[ 0 ];
    const auto tmp26 = tmp7 * tmp25;
    const auto tmp27 = tmp26 + tmp26;
    const auto tmp28 = tmp21 + tmp27;
    const auto tmp29 = 0.5 * tmp28;
    const auto tmp30 = u[ 0 ] * tmp29;
    const auto tmp31 = tmp0 * tmp0;
    const auto tmp32 = tmp7 * tmp7;
    const auto tmp33 = tmp31 + tmp32;
    const auto tmp34 = 0.5 * tmp33;
    const auto tmp35 = (du[ 0 ])[ 0 ] * tmp34;
    const auto tmp36 = tmp30 + tmp35;
    const auto tmp37 = -1 * tmp36;
    const auto tmp38 = (du[ 3 ])[ 0 ] + tmp37;
    const auto tmp39 = 0.3999999999999999 * tmp38;
    const auto tmp40 = tmp19 + tmp39;
    const auto tmp41 = (du[ 0 ])[ 1 ] * tmp0;
    const auto tmp42 = -1 * tmp41;
    const auto tmp43 = (du[ 1 ])[ 1 ] + tmp42;
    const auto tmp44 = tmp43 / u[ 0 ];
    const auto tmp45 = u[ 0 ] * tmp44;
    const auto tmp46 = tmp41 + tmp45;
    const auto tmp47 = tmp46 * tmp7;
    const auto tmp48 = tmp17 * tmp11;
    const auto tmp49 = tmp47 + tmp48;
    const auto tmp50 = tmp40 + tmp49;
    const auto tmp51 = -1 * tmp50;
    const auto tmp52 = tmp13 * tmp7;
    const auto tmp53 = u[ 0 ] * tmp7;
    const auto tmp54 = tmp53 * tmp11;
    const auto tmp55 = tmp52 + tmp54;
    const auto tmp56 = tmp0 * tmp44;
    const auto tmp57 = tmp56 + tmp56;
    const auto tmp58 = tmp7 * tmp11;
    const auto tmp59 = tmp58 + tmp58;
    const auto tmp60 = tmp57 + tmp59;
    const auto tmp61 = 0.5 * tmp60;
    const auto tmp62 = u[ 0 ] * tmp61;
    const auto tmp63 = (du[ 0 ])[ 1 ] * tmp34;
    const auto tmp64 = tmp62 + tmp63;
    const auto tmp65 = -1 * tmp64;
    const auto tmp66 = (du[ 3 ])[ 1 ] + tmp65;
    const auto tmp67 = 0.3999999999999999 * tmp66;
    const auto tmp68 = tmp55 + tmp67;
    const auto tmp69 = tmp6 * tmp7;
    const auto tmp70 = tmp17 * tmp25;
    const auto tmp71 = tmp69 + tmp70;
    const auto tmp72 = tmp68 + tmp71;
    const auto tmp73 = -1 * tmp72;
    const auto tmp74 = (du[ 3 ])[ 0 ] + tmp39;
    const auto tmp75 = tmp74 * tmp0;
    const auto tmp76 = u[ 0 ] * tmp34;
    const auto tmp77 = -1 * tmp76;
    const auto tmp78 = u[ 3 ] + tmp77;
    const auto tmp79 = 0.3999999999999999 * tmp78;
    const auto tmp80 = u[ 3 ] + tmp79;
    const auto tmp81 = tmp80 * tmp4;
    const auto tmp82 = tmp75 + tmp81;
    const auto tmp83 = (du[ 3 ])[ 1 ] + tmp67;
    const auto tmp84 = tmp83 * tmp7;
    const auto tmp85 = tmp80 * tmp11;
    const auto tmp86 = tmp84 + tmp85;
    const auto tmp87 = tmp82 + tmp86;
    const auto tmp88 = -1 * tmp87;
    result[ 0 ] = tmp15;
    result[ 1 ] = tmp51;
    result[ 2 ] = tmp73;
    result[ 3 ] = tmp88;
  }

  template< class Point >
  void alpha ( const Point &x, const DRangeType &u, RRangeType &result ) const
  {
    result[ 0 ] = 0;
    result[ 1 ] = 0;
    result[ 2 ] = 0;
    result[ 3 ] = 0;
  }

  template< class Point >
  void linAlpha ( const DRangeType &ubar, const Point &x, const DRangeType &u, RRangeType &result ) const
  {
    result[ 0 ] = 0;
    result[ 1 ] = 0;
    result[ 2 ] = 0;
    result[ 3 ] = 0;
  }

  bool hasNeumanBoundary () const
  {
    return false;
  }

  bool hasDirichletBoundary () const
  {
    return false;
  }

  bool isDirichletIntersection ( const IntersectionType &intersection, Dune::FieldVector< int, dimR > &dirichletComponent ) const
  {
    return false;
  }

  template< class Point >
  void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
  {
    result = RRangeType( 0 );
  }

private:
  mutable const EntityType *entity_ = nullptr;
};

}}} // namespaces

#endif
