namespace euler
{




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
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, std::tuple< std::shared_ptr< double > > >::element_type;
    static const std::size_t numConstants = 1;

    Model ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {
      std::get< 0 >( constants_ ) = std::make_shared< double >();
    }

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
      using std::pow;
      double tmp0 = constant< 0 >();
      const auto tmp1 = -1 * tmp0;
      const auto tmp2 = u[ 1 ] / u[ 0 ];
      const auto tmp3 = u[ 0 ] * tmp2;
      const auto tmp4 = tmp1 * tmp3;
      const auto tmp5 = u[ 2 ] / u[ 0 ];
      const auto tmp6 = u[ 0 ] * tmp5;
      const auto tmp7 = tmp1 * tmp6;
      const auto tmp8 = std::pow( u[ 1 ], 2 );
      const auto tmp9 = std::pow( u[ 2 ], 2 );
      const auto tmp10 = tmp8 + tmp9;
      const auto tmp11 = tmp10 / 2;
      const auto tmp12 = -1 * tmp11;
      const auto tmp13 = u[ 3 ] + tmp12;
      const auto tmp14 = 0.3999999999999999 * tmp13;
      const auto tmp15 = std::pow( tmp2, 2 );
      const auto tmp16 = u[ 0 ] * tmp15;
      const auto tmp17 = tmp14 + tmp16;
      const auto tmp18 = tmp17 * tmp1;
      const auto tmp19 = tmp3 * tmp5;
      const auto tmp20 = tmp1 * tmp19;
      const auto tmp21 = std::pow( tmp5, 2 );
      const auto tmp22 = u[ 0 ] * tmp21;
      const auto tmp23 = tmp14 + tmp22;
      const auto tmp24 = tmp23 * tmp1;
      const auto tmp25 = u[ -1 ] + tmp14;
      const auto tmp26 = tmp25 * tmp2;
      const auto tmp27 = tmp1 * tmp26;
      const auto tmp28 = tmp25 * tmp5;
      const auto tmp29 = tmp1 * tmp28;
      (result[ 0 ])[ 0 ] = tmp4;
      (result[ 0 ])[ 1 ] = tmp7;
      (result[ 1 ])[ 0 ] = tmp18;
      (result[ 1 ])[ 1 ] = tmp20;
      (result[ 2 ])[ 0 ] = tmp20;
      (result[ 2 ])[ 1 ] = tmp24;
      (result[ 3 ])[ 0 ] = tmp27;
      (result[ 3 ])[ 1 ] = tmp29;
    }

    template< class Point >
    void linDiffusiveFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      using std::pow;
      const auto tmp0 = ubar[ 1 ] / ubar[ 0 ];
      const auto tmp1 = u[ 0 ] * tmp0;
      const auto tmp2 = -1 * tmp1;
      const auto tmp3 = u[ 1 ] + tmp2;
      const auto tmp4 = tmp3 / ubar[ 0 ];
      const auto tmp5 = ubar[ 0 ] * tmp4;
      const auto tmp6 = tmp1 + tmp5;
      double tmp7 = constant< 0 >();
      const auto tmp8 = -1 * tmp7;
      const auto tmp9 = tmp6 * tmp8;
      const auto tmp10 = ubar[ 2 ] / ubar[ 0 ];
      const auto tmp11 = u[ 0 ] * tmp10;
      const auto tmp12 = -1 * tmp11;
      const auto tmp13 = u[ 2 ] + tmp12;
      const auto tmp14 = tmp13 / ubar[ 0 ];
      const auto tmp15 = ubar[ 0 ] * tmp14;
      const auto tmp16 = tmp11 + tmp15;
      const auto tmp17 = tmp16 * tmp8;
      const auto tmp18 = 2 * tmp4;
      const auto tmp19 = tmp18 * tmp0;
      const auto tmp20 = ubar[ 0 ] * tmp19;
      const auto tmp21 = std::pow( tmp0, 2 );
      const auto tmp22 = u[ 0 ] * tmp21;
      const auto tmp23 = tmp20 + tmp22;
      const auto tmp24 = 2 * u[ 1 ];
      const auto tmp25 = ubar[ 1 ] * tmp24;
      const auto tmp26 = 2 * u[ 2 ];
      const auto tmp27 = ubar[ 2 ] * tmp26;
      const auto tmp28 = tmp25 + tmp27;
      const auto tmp29 = tmp28 / 2;
      const auto tmp30 = -1 * tmp29;
      const auto tmp31 = u[ 3 ] + tmp30;
      const auto tmp32 = 0.3999999999999999 * tmp31;
      const auto tmp33 = tmp23 + tmp32;
      const auto tmp34 = tmp33 * tmp8;
      const auto tmp35 = tmp6 * tmp10;
      const auto tmp36 = ubar[ 0 ] * tmp0;
      const auto tmp37 = tmp36 * tmp14;
      const auto tmp38 = tmp35 + tmp37;
      const auto tmp39 = tmp38 * tmp8;
      const auto tmp40 = 2 * tmp14;
      const auto tmp41 = tmp40 * tmp10;
      const auto tmp42 = ubar[ 0 ] * tmp41;
      const auto tmp43 = std::pow( tmp10, 2 );
      const auto tmp44 = u[ 0 ] * tmp43;
      const auto tmp45 = tmp42 + tmp44;
      const auto tmp46 = tmp45 + tmp32;
      const auto tmp47 = tmp46 * tmp8;
      const auto tmp48 = u[ -1 ] + tmp32;
      const auto tmp49 = tmp48 * tmp0;
      const auto tmp50 = std::pow( ubar[ 1 ], 2 );
      const auto tmp51 = std::pow( ubar[ 2 ], 2 );
      const auto tmp52 = tmp50 + tmp51;
      const auto tmp53 = tmp52 / 2;
      const auto tmp54 = -1 * tmp53;
      const auto tmp55 = ubar[ 3 ] + tmp54;
      const auto tmp56 = 0.3999999999999999 * tmp55;
      const auto tmp57 = ubar[ -1 ] + tmp56;
      const auto tmp58 = tmp57 * tmp4;
      const auto tmp59 = tmp49 + tmp58;
      const auto tmp60 = tmp59 * tmp8;
      const auto tmp61 = tmp48 * tmp10;
      const auto tmp62 = tmp57 * tmp14;
      const auto tmp63 = tmp61 + tmp62;
      const auto tmp64 = tmp63 * tmp8;
      (result[ 0 ])[ 0 ] = tmp9;
      (result[ 0 ])[ 1 ] = tmp17;
      (result[ 1 ])[ 0 ] = tmp34;
      (result[ 1 ])[ 1 ] = tmp39;
      (result[ 2 ])[ 0 ] = tmp39;
      (result[ 2 ])[ 1 ] = tmp47;
      (result[ 3 ])[ 0 ] = tmp60;
      (result[ 3 ])[ 1 ] = tmp64;
    }

    template< class Point >
    void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
    {
      using std::pow;
      const auto tmp0 = u[ 1 ] / u[ 0 ];
      const auto tmp1 = (du[ 0 ])[ 0 ] * tmp0;
      const auto tmp2 = -1 * tmp1;
      const auto tmp3 = (du[ 1 ])[ 0 ] + tmp2;
      const auto tmp4 = tmp3 / u[ 0 ];
      const auto tmp5 = u[ 0 ] * tmp4;
      const auto tmp6 = tmp1 + tmp5;
      double tmp7 = constant< 0 >();
      const auto tmp8 = -1 * tmp7;
      const auto tmp9 = tmp6 * tmp8;
      const auto tmp10 = u[ 2 ] / u[ 0 ];
      const auto tmp11 = (du[ 0 ])[ 1 ] * tmp10;
      const auto tmp12 = -1 * tmp11;
      const auto tmp13 = (du[ 2 ])[ 1 ] + tmp12;
      const auto tmp14 = tmp13 / u[ 0 ];
      const auto tmp15 = u[ 0 ] * tmp14;
      const auto tmp16 = tmp11 + tmp15;
      const auto tmp17 = tmp16 * tmp8;
      const auto tmp18 = tmp9 + tmp17;
      const auto tmp19 = -1 * tmp18;
      const auto tmp20 = 2 * tmp4;
      const auto tmp21 = tmp20 * tmp0;
      const auto tmp22 = u[ 0 ] * tmp21;
      const auto tmp23 = std::pow( tmp0, 2 );
      const auto tmp24 = (du[ 0 ])[ 0 ] * tmp23;
      const auto tmp25 = tmp22 + tmp24;
      const auto tmp26 = 2 * (du[ 1 ])[ 0 ];
      const auto tmp27 = u[ 1 ] * tmp26;
      const auto tmp28 = 2 * (du[ 2 ])[ 0 ];
      const auto tmp29 = u[ 2 ] * tmp28;
      const auto tmp30 = tmp27 + tmp29;
      const auto tmp31 = tmp30 / 2;
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = (du[ 3 ])[ 0 ] + tmp32;
      const auto tmp34 = 0.3999999999999999 * tmp33;
      const auto tmp35 = tmp25 + tmp34;
      const auto tmp36 = tmp35 * tmp8;
      const auto tmp37 = (du[ 0 ])[ 1 ] * tmp0;
      const auto tmp38 = -1 * tmp37;
      const auto tmp39 = (du[ 1 ])[ 1 ] + tmp38;
      const auto tmp40 = tmp39 / u[ 0 ];
      const auto tmp41 = u[ 0 ] * tmp40;
      const auto tmp42 = tmp37 + tmp41;
      const auto tmp43 = tmp42 * tmp10;
      const auto tmp44 = u[ 0 ] * tmp0;
      const auto tmp45 = tmp44 * tmp14;
      const auto tmp46 = tmp43 + tmp45;
      const auto tmp47 = tmp46 * tmp8;
      const auto tmp48 = tmp36 + tmp47;
      const auto tmp49 = -1 * tmp48;
      const auto tmp50 = 2 * tmp14;
      const auto tmp51 = tmp50 * tmp10;
      const auto tmp52 = u[ 0 ] * tmp51;
      const auto tmp53 = std::pow( tmp10, 2 );
      const auto tmp54 = (du[ 0 ])[ 1 ] * tmp53;
      const auto tmp55 = tmp52 + tmp54;
      const auto tmp56 = 2 * (du[ 1 ])[ 1 ];
      const auto tmp57 = u[ 1 ] * tmp56;
      const auto tmp58 = 2 * (du[ 2 ])[ 1 ];
      const auto tmp59 = u[ 2 ] * tmp58;
      const auto tmp60 = tmp57 + tmp59;
      const auto tmp61 = tmp60 / 2;
      const auto tmp62 = -1 * tmp61;
      const auto tmp63 = (du[ 3 ])[ 1 ] + tmp62;
      const auto tmp64 = 0.3999999999999999 * tmp63;
      const auto tmp65 = tmp55 + tmp64;
      const auto tmp66 = tmp65 * tmp8;
      const auto tmp67 = tmp6 * tmp10;
      const auto tmp68 = (du[ 0 ])[ 0 ] * tmp10;
      const auto tmp69 = -1 * tmp68;
      const auto tmp70 = (du[ 2 ])[ 0 ] + tmp69;
      const auto tmp71 = tmp70 / u[ 0 ];
      const auto tmp72 = tmp44 * tmp71;
      const auto tmp73 = tmp67 + tmp72;
      const auto tmp74 = tmp73 * tmp8;
      const auto tmp75 = tmp66 + tmp74;
      const auto tmp76 = -1 * tmp75;
      const auto tmp77 = (du[ -1 ])[ 0 ] + tmp34;
      const auto tmp78 = tmp77 * tmp0;
      const auto tmp79 = std::pow( u[ 1 ], 2 );
      const auto tmp80 = std::pow( u[ 2 ], 2 );
      const auto tmp81 = tmp79 + tmp80;
      const auto tmp82 = tmp81 / 2;
      const auto tmp83 = -1 * tmp82;
      const auto tmp84 = u[ 3 ] + tmp83;
      const auto tmp85 = 0.3999999999999999 * tmp84;
      const auto tmp86 = u[ -1 ] + tmp85;
      const auto tmp87 = tmp86 * tmp4;
      const auto tmp88 = tmp78 + tmp87;
      const auto tmp89 = tmp88 * tmp8;
      const auto tmp90 = (du[ -1 ])[ 1 ] + tmp64;
      const auto tmp91 = tmp90 * tmp10;
      const auto tmp92 = tmp86 * tmp14;
      const auto tmp93 = tmp91 + tmp92;
      const auto tmp94 = tmp93 * tmp8;
      const auto tmp95 = tmp89 + tmp94;
      const auto tmp96 = -1 * tmp95;
      result[ 0 ] = tmp19;
      result[ 1 ] = tmp49;
      result[ 2 ] = tmp76;
      result[ 3 ] = tmp96;
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

    template< std::size_t i >
    const ConstantType< i > &constant () const
    {
      return *std::get< i >( constants_ );
    }

    template< std::size_t i >
    ConstantType< i > &constant ()
    {
      return *std::get< i >( constants_ );
    }

    const double &dt () const
    {
      return *std::get< 0 >( constants_ );
    }

    double &dt ()
    {
      return *std::get< 0 >( constants_ );
    }

  private:
    mutable const EntityType *entity_ = nullptr;
    mutable std::tuple< std::shared_ptr< double > > constants_;
  };

} // namespace euler
