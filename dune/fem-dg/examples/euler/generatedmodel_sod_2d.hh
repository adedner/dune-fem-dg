#ifndef GuardModelImpl_nvc600ee5a9901394cd3682562ac74fdf3_ffdcb84ddeea836e1032a1b0b6306fec
#define GuardModelImpl_nvc600ee5a9901394cd3682562ac74fdf3_ffdcb84ddeea836e1032a1b0b6306fec

#include <dune/fem/misc/boundaryidprovider.hh>
#include <dune/fem-dg/misc/algorithmcreatorenums.hh>

namespace ModelImpl
{

template< class FunctionSpace >
struct Additional
{
  typedef typename FunctionSpace::DomainType DomainType;
  typedef typename FunctionSpace::RangeType RangeType;
  typedef typename FunctionSpace::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpace::HessianRangeType HessianRangeType;
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
  static const Dune::Fem::AdvectionLimiterFunction::Enum limiterFunctionId = Dune::Fem::AdvectionLimiterFunction::Enum::minmod;
  static const Dune::Fem::AdvectionFlux::Enum advFluxId = Dune::Fem::AdvectionFlux::Enum::euler_general;
  static const Dune::Fem::DiffusionFlux::Enum diffFluxId = Dune::Fem::DiffusionFlux::Enum::none;
  static const bool defaultQuadrature = true;
};



  // Model
  // -----

  template< class GridPart >
  struct Model
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

    Model ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
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
      const auto tmp12 = u[ 0 ] * tmp4;
      const auto tmp13 = tmp11 + tmp12;
      const auto tmp14 = tmp0 * tmp2;
      const auto tmp15 = u[ 0 ] * tmp14;
      const auto tmp16 = u[ 0 ] * tmp5;
      const auto tmp17 = tmp11 + tmp16;
      const auto tmp18 = u[ 3 ] + tmp11;
      const auto tmp19 = tmp18 * tmp0;
      const auto tmp20 = tmp18 * tmp2;
      (result[ 0 ])[ 0 ] = tmp1;
      (result[ 0 ])[ 1 ] = tmp3;
      (result[ 1 ])[ 0 ] = tmp13;
      (result[ 1 ])[ 1 ] = tmp15;
      (result[ 2 ])[ 0 ] = tmp15;
      (result[ 2 ])[ 1 ] = tmp17;
      (result[ 3 ])[ 0 ] = tmp19;
      (result[ 3 ])[ 1 ] = tmp20;
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
      const auto tmp14 = tmp0 * tmp4;
      const auto tmp15 = tmp14 + tmp14;
      const auto tmp16 = ubar[ 0 ] * tmp15;
      const auto tmp17 = tmp0 * tmp0;
      const auto tmp18 = u[ 0 ] * tmp17;
      const auto tmp19 = tmp16 + tmp18;
      const auto tmp20 = tmp7 * tmp11;
      const auto tmp21 = tmp20 + tmp20;
      const auto tmp22 = tmp15 + tmp21;
      const auto tmp23 = 0.5 * tmp22;
      const auto tmp24 = ubar[ 0 ] * tmp23;
      const auto tmp25 = tmp7 * tmp7;
      const auto tmp26 = tmp17 + tmp25;
      const auto tmp27 = 0.5 * tmp26;
      const auto tmp28 = u[ 0 ] * tmp27;
      const auto tmp29 = tmp24 + tmp28;
      const auto tmp30 = -1 * tmp29;
      const auto tmp31 = u[ 3 ] + tmp30;
      const auto tmp32 = 0.3999999999999999 * tmp31;
      const auto tmp33 = tmp19 + tmp32;
      const auto tmp34 = tmp7 * tmp4;
      const auto tmp35 = tmp0 * tmp11;
      const auto tmp36 = tmp34 + tmp35;
      const auto tmp37 = ubar[ 0 ] * tmp36;
      const auto tmp38 = tmp0 * tmp7;
      const auto tmp39 = u[ 0 ] * tmp38;
      const auto tmp40 = tmp37 + tmp39;
      const auto tmp41 = ubar[ 0 ] * tmp21;
      const auto tmp42 = u[ 0 ] * tmp25;
      const auto tmp43 = tmp41 + tmp42;
      const auto tmp44 = tmp43 + tmp32;
      const auto tmp45 = u[ 3 ] + tmp32;
      const auto tmp46 = tmp45 * tmp0;
      const auto tmp47 = ubar[ 0 ] * tmp27;
      const auto tmp48 = -1 * tmp47;
      const auto tmp49 = ubar[ 3 ] + tmp48;
      const auto tmp50 = 0.3999999999999999 * tmp49;
      const auto tmp51 = ubar[ 3 ] + tmp50;
      const auto tmp52 = tmp51 * tmp4;
      const auto tmp53 = tmp46 + tmp52;
      const auto tmp54 = tmp45 * tmp7;
      const auto tmp55 = tmp51 * tmp11;
      const auto tmp56 = tmp54 + tmp55;
      (result[ 0 ])[ 0 ] = tmp6;
      (result[ 0 ])[ 1 ] = tmp13;
      (result[ 1 ])[ 0 ] = tmp33;
      (result[ 1 ])[ 1 ] = tmp40;
      (result[ 2 ])[ 0 ] = tmp40;
      (result[ 2 ])[ 1 ] = tmp44;
      (result[ 3 ])[ 0 ] = tmp53;
      (result[ 3 ])[ 1 ] = tmp56;
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
      const auto tmp16 = tmp0 * tmp4;
      const auto tmp17 = tmp16 + tmp16;
      const auto tmp18 = u[ 0 ] * tmp17;
      const auto tmp19 = tmp0 * tmp0;
      const auto tmp20 = (du[ 0 ])[ 0 ] * tmp19;
      const auto tmp21 = tmp18 + tmp20;
      const auto tmp22 = (du[ 0 ])[ 0 ] * tmp7;
      const auto tmp23 = -1 * tmp22;
      const auto tmp24 = (du[ 2 ])[ 0 ] + tmp23;
      const auto tmp25 = tmp24 / u[ 0 ];
      const auto tmp26 = tmp7 * tmp25;
      const auto tmp27 = tmp26 + tmp26;
      const auto tmp28 = tmp17 + tmp27;
      const auto tmp29 = 0.5 * tmp28;
      const auto tmp30 = u[ 0 ] * tmp29;
      const auto tmp31 = tmp7 * tmp7;
      const auto tmp32 = tmp19 + tmp31;
      const auto tmp33 = 0.5 * tmp32;
      const auto tmp34 = (du[ 0 ])[ 0 ] * tmp33;
      const auto tmp35 = tmp30 + tmp34;
      const auto tmp36 = -1 * tmp35;
      const auto tmp37 = (du[ 3 ])[ 0 ] + tmp36;
      const auto tmp38 = 0.3999999999999999 * tmp37;
      const auto tmp39 = tmp21 + tmp38;
      const auto tmp40 = (du[ 0 ])[ 1 ] * tmp0;
      const auto tmp41 = -1 * tmp40;
      const auto tmp42 = (du[ 1 ])[ 1 ] + tmp41;
      const auto tmp43 = tmp42 / u[ 0 ];
      const auto tmp44 = tmp7 * tmp43;
      const auto tmp45 = tmp0 * tmp11;
      const auto tmp46 = tmp44 + tmp45;
      const auto tmp47 = u[ 0 ] * tmp46;
      const auto tmp48 = tmp0 * tmp7;
      const auto tmp49 = (du[ 0 ])[ 1 ] * tmp48;
      const auto tmp50 = tmp47 + tmp49;
      const auto tmp51 = tmp39 + tmp50;
      const auto tmp52 = -1 * tmp51;
      const auto tmp53 = tmp7 * tmp11;
      const auto tmp54 = tmp53 + tmp53;
      const auto tmp55 = u[ 0 ] * tmp54;
      const auto tmp56 = (du[ 0 ])[ 1 ] * tmp31;
      const auto tmp57 = tmp55 + tmp56;
      const auto tmp58 = tmp0 * tmp43;
      const auto tmp59 = tmp58 + tmp58;
      const auto tmp60 = tmp59 + tmp54;
      const auto tmp61 = 0.5 * tmp60;
      const auto tmp62 = u[ 0 ] * tmp61;
      const auto tmp63 = (du[ 0 ])[ 1 ] * tmp33;
      const auto tmp64 = tmp62 + tmp63;
      const auto tmp65 = -1 * tmp64;
      const auto tmp66 = (du[ 3 ])[ 1 ] + tmp65;
      const auto tmp67 = 0.3999999999999999 * tmp66;
      const auto tmp68 = tmp57 + tmp67;
      const auto tmp69 = tmp7 * tmp4;
      const auto tmp70 = tmp0 * tmp25;
      const auto tmp71 = tmp69 + tmp70;
      const auto tmp72 = u[ 0 ] * tmp71;
      const auto tmp73 = (du[ 0 ])[ 0 ] * tmp48;
      const auto tmp74 = tmp72 + tmp73;
      const auto tmp75 = tmp68 + tmp74;
      const auto tmp76 = -1 * tmp75;
      const auto tmp77 = (du[ 3 ])[ 0 ] + tmp38;
      const auto tmp78 = tmp77 * tmp0;
      const auto tmp79 = u[ 0 ] * tmp33;
      const auto tmp80 = -1 * tmp79;
      const auto tmp81 = u[ 3 ] + tmp80;
      const auto tmp82 = 0.3999999999999999 * tmp81;
      const auto tmp83 = u[ 3 ] + tmp82;
      const auto tmp84 = tmp83 * tmp4;
      const auto tmp85 = tmp78 + tmp84;
      const auto tmp86 = (du[ 3 ])[ 1 ] + tmp67;
      const auto tmp87 = tmp86 * tmp7;
      const auto tmp88 = tmp83 * tmp11;
      const auto tmp89 = tmp87 + tmp88;
      const auto tmp90 = tmp85 + tmp89;
      const auto tmp91 = -1 * tmp90;
      result[ 0 ] = tmp15;
      result[ 1 ] = tmp52;
      result[ 2 ] = tmp76;
      result[ 3 ] = tmp91;
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
    typedef Dune::FieldVector<int,4> DirichletComponentType;

    bool hasDirichletBoundary () const
    {
      return false;
    }

    bool isDirichletIntersection ( const IntersectionType &intersection, DirichletComponentType &dirichletComponent ) const
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

  public:
    static constexpr bool hasGamma = true;
    static constexpr double gamma = 1.4;

    template< class Entity, class Point >
    double maxWaveSpeed ( const double &t, const Entity &entity, const Point &x, const DDomainType &normal, const DRangeType &u ) const
    {
      double result;
      using std::abs;
      using std::sqrt;
      result = std::abs( normal[ 0 ] * (u[ 1 ] / u[ 0 ]) + normal[ 1 ] * (u[ 2 ] / u[ 0 ]) ) + std::sqrt( (1.4 * (0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))))) / u[ 0 ] );
      return result;
    }

    template< class Entity, class Point >
    DDomainType velocity ( const double &t, const Entity &entity, const Point &x, const DRangeType &u ) const
    {
      DDomainType result;
      result[ 0 ] = u[ 1 ] / u[ 0 ];
      result[ 1 ] = u[ 2 ] / u[ 0 ];
      return result;
    }

    template< class Entity, class Point>
    double diffusionTimeStep ( const Entity& entity, const Point &x, const DRangeType &u ) const
    {
      double result;
      result = double( (0) );
      return result;
    }
    static constexpr bool hasPhysical = true;

    template< class Entity, class Point >
    double physical ( const Entity &entity, const Point &x, const DRangeType &u ) const
    {
      double result;
      result = (u[ 0 ] < 1e-08 ? 0 : (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ])))) > 1e-08 ? 1 : 0));
      return result;
    }

    template< class Intersection, class Point >
    double jump ( const Intersection& it, const Point &x, const DRangeType &u, const DRangeType &w ) const
    {
      double result;
      result = (0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))) + -1 * (0.3999999999999999 * (w[ 3 ] + -1 * (w[ 0 ] * (0.5 * ((w[ 1 ] / w[ 0 ]) * (w[ 1 ] / w[ 0 ]) + (w[ 2 ] / w[ 0 ]) * (w[ 2 ] / w[ 0 ]))))))) / (0.5 * (0.3999999999999999 * (u[ 3 ] + -1 * (u[ 0 ] * (0.5 * ((u[ 1 ] / u[ 0 ]) * (u[ 1 ] / u[ 0 ]) + (u[ 2 ] / u[ 0 ]) * (u[ 2 ] / u[ 0 ]))))) + 0.3999999999999999 * (w[ 3 ] + -1 * (w[ 0 ] * (0.5 * ((w[ 1 ] / w[ 0 ]) * (w[ 1 ] / w[ 0 ]) + (w[ 2 ] / w[ 0 ]) * (w[ 2 ] / w[ 0 ])))))));
      return result;
    }

    template< class Entity, class Point >
    void adjustAverageValue ( const Entity& entity, const Point &x, DRangeType &u ) const
    {
      {}
    }

    template< class Entity, class Point >
    bool boundaryFlux ( const int bndId, const double &t, const Entity& entity, const Point &x, const DDomainType &normal, const DRangeType &u, RRangeType &result ) const
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
    bool diffusionBoundaryFlux ( const int bndId, const double &t, const Entity& entity, const Point &x, const DDomainType &normal, const DRangeType &u, const DJacobianRangeType &jac, RRangeType &result ) const
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
    bool hasBoundaryValue ( const int bndId, const double &t, const Entity& entity, const Point &x, const DRangeType &u, RRangeType &result ) const
    {
      switch( bndId )
      {
      case 1:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      case 2:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      case 3:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      case 4:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      default:
        {
          return false;
        }
      }
    }

    template< class Entity, class Point >
    bool boundaryValue ( const int bndId, const double &t, const Entity& entity, const Point &x, const DDomainType &normal, const DRangeType &u, RRangeType &result ) const
    {
      switch( bndId )
      {
      case 1:
        {
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = normal[ 0 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ])) + -1 * (normal[ 1 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])));
          result[ 2 ] = normal[ 0 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])) + normal[ 1 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ]));
          result[ 3 ] = u[ 3 ];
          return true;
        }
        break;
      case 2:
        {
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = normal[ 0 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ])) + -1 * (normal[ 1 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])));
          result[ 2 ] = normal[ 0 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])) + normal[ 1 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ]));
          result[ 3 ] = u[ 3 ];
          return true;
        }
        break;
      case 3:
        {
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = normal[ 0 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ])) + -1 * (normal[ 1 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])));
          result[ 2 ] = normal[ 0 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])) + normal[ 1 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ]));
          result[ 3 ] = u[ 3 ];
          return true;
        }
        break;
      case 4:
        {
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = normal[ 0 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ])) + -1 * (normal[ 1 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])));
          result[ 2 ] = normal[ 0 ] * (normal[ 0 ] * u[ 2 ] + u[ 1 ] * (-1 * normal[ 1 ])) + normal[ 1 ] * (-1 * (normal[ 0 ] * u[ 1 ] + normal[ 1 ] * u[ 2 ]));
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
    void limitedRange ( LimitedRange& limRange ) const
    {
      {}
    }

    template< class LimitedRange >
    void obtainBounds ( LimitedRange& globalMin, LimitedRange& globalMax ) const
    {
      globalMin = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};
      globalMax = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    }
  };

} // namespace ModelImpl

#endif // GuardModelImpl_nvc600ee5a9901394cd3682562ac74fdf3_ffdcb84ddeea836e1032a1b0b6306fec
