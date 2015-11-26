#ifndef DUNE_FEM_DG_DIFFUSIONFLUXES_HH
#define DUNE_FEM_DG_DIFFUSIONFLUXES_HH

#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>

namespace Dune {


  // DGFluxTupleToVectorConverter
  //-----------------------------

  template <class ArgumentVectorTuple, int passUId >
  class DGFluxTupleToVectorConverter
  {
    integral_constant< int, passUId > uVar;
    const ArgumentVectorTuple& vec_;

    //DGFluxTupleToVectorConverter(const DGFluxTupleToVectorConverter&);
  public:
    typedef typename ArgumentVectorTuple :: value_type TupleType;
    typedef typename TupleType::template Value< integral_constant< int, passUId > >::Type ValueType;

    DGFluxTupleToVectorConverter(const ArgumentVectorTuple& vec)
      : vec_( vec )
    {}

    const TupleType& tuple( const size_t i ) const
    {
      assert( i < vec_.size() );
      return vec_[ i ];
    }

    const ValueType& operator [] ( const size_t i ) const
    {
      return tuple( i )[ uVar ];
    }
  };

  ///////////////////////////////////////////////////////////
  //
  //  Identifier for Diffusion Fluxes for Primal methods
  //
  //////////////////////////////////////////////////////////
  enum DGDiffusionFluxIdentifier {
    method_cdg2    = 0,  // CDG 2 (Compact Discontinuous Galerkin 2)
    method_cdg     = 1,  // CDG (Compact Discontinuous Galerkin)
    method_br2     = 2,  // BR2 (Bassi-Rebay 2)
    method_ip      = 3,  // IP (Interior Penalty)
    method_nipg    = 4,  // NIPG (Non-symmetric Interior  Penalty)
    method_bo      = 5,  // BO (Baumann-Oden)
    method_general = 6,   // general means all methods chosen via parameter file
    method_none    = 7   // no diffusion (advection only)
  };
  enum DGLiftingFluxIdentifier {
    lifting_id_id    = 0,  // int_Omega r([u]).tau  = -int_e [u].{tau}
    lifting_id_A     = 1,  // int_Omega r([u]).tau  = -int_e [u].{Atau}
    lifting_A_A      = 2   // int_Omega r([u]).Atau = -int_e [u].{Atau}
  };

  class DGPrimalFormulationParameters
    : public Dune::Fem::LocalParameter< DGPrimalFormulationParameters, DGPrimalFormulationParameters >
  {
    const std::string keyPrefix_;
  public:
    typedef DGDiffusionFluxIdentifier MethodType;
    typedef DGLiftingFluxIdentifier   LiftingType;

    explicit DGPrimalFormulationParameters( const std::string keyPrefix = "dgdiffusionflux." )
      : keyPrefix_( keyPrefix )
    {}

    static std::string methodNames( const MethodType mthd )
    {
      const std::string method []
        = { "CDG2", "CDG" , "BR2", "IP" , "NIPG", "BO" };
      assert( mthd >= method_cdg2 && mthd < method_general );
      return method[ mthd ];
    }

    virtual MethodType getMethod() const
    {
      const std::string method []
        = { methodNames( method_cdg2 ),
            methodNames( method_cdg ),
            methodNames( method_br2 ),
            methodNames( method_ip ),
            methodNames( method_nipg ),
            methodNames( method_bo )
          };
      return (MethodType) Fem::Parameter::getEnum( keyPrefix_ + "method", method );
    }

    static std::string liftingNames( const LiftingType mthd )
    {
      const std::string method []
        = { "id_id", "id_A" , "A_A" };
      return method[ mthd ];
    }

    virtual LiftingType getLifting() const
    {
      const std::string method []
        = { liftingNames( lifting_id_id ),
            liftingNames( lifting_id_A ),
            liftingNames( lifting_A_A )
          };
      return (LiftingType) Fem::Parameter::getEnum( keyPrefix_ + "lifting", method, 0 );
    }

    virtual double penalty() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "penalty" );
    }

    virtual double liftfactor() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "liftfactor" );
    }

    virtual double theoryparameters() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "theoryparameters", 0. );
    }

    template <class DomainType>
    void upwind( DomainType& upwd ) const
    {
      Fem::Parameter::get(keyPrefix_ + "upwind", upwd, upwd);
    }
  };


  // DGDiffusionFluxBase
  //--------------------

  template <class DiscreteFunctionSpaceImp,
            class Model >
  class DGDiffusionFluxBase
  {
  public:
    enum { evaluateJacobian = false };
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    typedef DGPrimalFormulationParameters ParameterType;
    typedef typename ParameterType :: MethodType   MethodType;
    typedef typename ParameterType :: LiftingType  LiftingType;

    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    enum { dimRange  = DiscreteFunctionSpaceType :: dimRange };

    typedef typename DiscreteFunctionSpaceType :: DomainType           DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType       RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeType            RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType    JacobianRangeType;

    typedef FieldVector< DomainFieldType, dimDomain-1 > FaceDomainType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType         GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename GridPartType :: GridType                          GridType;
    typedef typename DiscreteFunctionSpaceType :: EntityType           EntityType;
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType :: polynomialOrder };

    // type of gradient space
    typedef typename DiscreteFunctionSpaceType ::
        template ToNewDimRange< dimGradRange > :: Type   DiscreteGradientSpaceType;

    typedef typename DiscreteGradientSpaceType :: RangeType GradientType;

    typedef Fem::TemporaryLocalFunction< DiscreteGradientSpaceType > LiftingFunctionType;

    typedef Fem::CachingQuadrature< GridPartType, 0> VolumeQuadratureType ;

    DomainType upwindDefault() const
    {
      DomainType upwind ( M_PI );
      // set upwind to some strange numbers
      if( dimDomain > 1 ) upwind[1] = M_LN2 ;
      if( dimDomain > 2 ) upwind[2] = M_E ;
      return upwind ;
    }
  public:
    /**
     * @brief constructor
     */
    DGDiffusionFluxBase( const Model& mod,
                         const bool initUpwind,
                         const ParameterType& parameters )
    : model_(mod),
      param_( parameters ),
      upwind_( upwindDefault() ),
      // Set CFL number for penalty term (compare diffusion in first pass)
      cflDiffinv_(8. * (polOrd+1.) ),
      dimensionFactor_( 2.0 * ( dimDomain ) ),
      nonconformingFactor_( 2.0 * ( dimDomain - 1 ) )
    {
      if( initUpwind )
      {
        parameter().upwind( upwind_ );
        if( Fem::Parameter :: verbose() )
          std::cout << "Using upwind = " << upwind_ << std::endl;
      }
    }

    //! copy constructor
    DGDiffusionFluxBase( const DGDiffusionFluxBase& other )
    : model_( other.model_ ),
      param_( other.parameter() ),
      upwind_( other.upwind_ ),
      cflDiffinv_( other.cflDiffinv_ ),
      dimensionFactor_( other.dimensionFactor_ ),
      nonconformingFactor_( other.nonconformingFactor_ )
    {
    }

    //! returns true if lifting has to be calculated
    //const bool hasLifting () const { return false; }

  protected:
    bool determineDirection( const bool areaSwitch, const double enVolume, const double nbVolume,
                             const Intersection& intersection ) const
    {
      if (areaSwitch && std::abs( enVolume - nbVolume ) > 1e-8)
      {
        // we need en=K^- and nb=K^+ such that c = (|K^-|/|K^+|) <= 1
        // outside is K^- if outside volume is smaller
        return ( enVolume > nbVolume );
      }
      else
      {
        Fem::Quadrature< DomainFieldType, dimDomain - 1 > quad( intersection.type() , 0 );
        const DomainType normal = intersection.outerNormal( quad.point( 0 ) );
        return determineDirection( normal );
      }
    }

    bool determineDirection( const DomainType& normal ) const
    {
       return ( normal * upwind_ ) < 0 ;
    }

  public:
    void switchUpwind()
    {
      upwind_ *= -1.0;
    }

    bool hasLifting() const { return false; }

    template <class LocalEvaluation, class ArgumentTupleVector >
    void initializeIntersection(const LocalEvaluation& left,
                                const LocalEvaluation& right,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec)
    {
    }

    template <class LocalEvaluation, class ArgumentTupleVector>
    void initializeBoundary(const LocalEvaluation& local,
                            const ArgumentTupleVector& uLeftVec,
                            const std::vector< RangeType >& uRight)
    {
    }

    template <class QuadratureImp, class ArgumentTupleVector >
    void initializeIntersection(const Intersection& intersection,
                                const EntityType& inside,
                                const EntityType& outside,
                                const double time,
                                const QuadratureImp& quadInner,
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec)
    {
    }

    template <class QuadratureImp, class ArgumentTupleVector>
    void initializeBoundary(const Intersection& intersection,
                            const EntityType& entity,
                            const double time,
                            const QuadratureImp& quadInner,
                            const ArgumentTupleVector& uLeftVec,
                            const std::vector< RangeType >& uRight)
    {
    }


  public:
    /**
     * @brief flux function on interfaces between cells for advection and diffusion
     *
     * @param[in] it intersection
     * @param[in] time current time given by TimeProvider
     * @param[in] x coordinate of required evaluation local to \c it
     * @param[in] uLeft DOF evaluation on this side of \c it
     * @param[in] uRight DOF evaluation on the other side of \c it
     * @param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * @param[out] gDiffRight advection flux projected on normal for the other side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * @return wave speed estimate (multiplied with the integration element of the intersection),
     *              to estimate the time step |T|/wave.
     */
    template <class LocalEvaluation>
    double numericalFlux(const LocalEvaluation& left,
                         const LocalEvaluation& right,
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const JacobianRangeType& jacLeft,
                         const JacobianRangeType& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight )
    {
      assert( false );
      abort();
      return 0.0;
    }


    template <class LocalEvaluation>
    double boundaryFlux(const LocalEvaluation& left,
                        const RangeType& uLeft,
                        const RangeType& uRight,
                        const JacobianRangeType& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft )   /*@LST0E@*/
    {

      assert( false );
      abort();
      return 0.0;
    }

    const Model &model () const { return model_; }
    const ParameterType& parameter() const { return param_; }

  protected:
    const Model   &model_;
    const ParameterType& param_;
    DomainType upwind_;
    const double cflDiffinv_;
    const double dimensionFactor_;
    const double nonconformingFactor_;
  }; // end DGPrimalDiffusionFlux

} // end namespace
#endif
