#ifndef DUNE_FEMDG_POISSON_ERROR_ESTIMATOR_HH
#define DUNE_FEMDG_POISSON_ERROR_ESTIMATOR_HH

//- Dune includes
// #include <dune/grid/common/referenceelements.hh>

#include <cmath>

#include <dune/geometry/type.hh>
#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

//- Dune-fem includes
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/nonconformitylevel.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/space/padaptivespace.hh>

#include <dune/fem-dg/pass/context.hh>
#include <dune/fem-dg/operator/adaptation/utility.hh>

namespace Dune
{
namespace Fem
{

  template <class OtherSpace>
  struct HierarchicSimplexDGSpace
  {
    typedef Dune::Fem::DiscontinuousGalerkinSpace< typename OtherSpace :: FunctionSpaceType,
                                                   typename OtherSpace :: GridPartType,
                                                   OtherSpace :: polynomialOrder > Type;

    enum { dimension = OtherSpace :: GridPartType :: dimension };

    static int numDGBaseFcts( int polOrder )
    {
      if( dimension == 2 )
        return  (polOrder + 2) * (polOrder + 1) / 2;
      else if ( dimension == 3 )
        return ((polOrder+1)*(polOrder+2)*(2*polOrder+3)/6  + (polOrder+1)*(polOrder+2)/2)/2;
      else
        abort();
    }

    static int size(int order)
    {
      assert( order > 0 );
      int result = numDGBaseFcts( order );
#ifndef NDEBUG
      if( dimension == 2 )
      {
        assert( order < 12 );
        static int bstart[12] = {1,3,6,10,15,21,28,36,45,55,66,78};
        if( result != bstart[ order ] )
        {
          std::cout << result << " r | b " <<  bstart[ order ]  << std::endl;
          abort();
        }
      }
#endif
      return result;
    }

  };

  template <class OtherSpace>
  struct HierarchicCubeDGSpace
  {
    typedef Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace< typename OtherSpace :: FunctionSpaceType,
                                                                     typename OtherSpace :: GridPartType,
                                                                     OtherSpace :: polynomialOrder > Type;
    typedef typename OtherSpace :: GridType GridType;

    static int size( const int order )
    {
      enum { dimension = GridType :: dimension };
      int result = 1;
      for( int i=0; i<dimension; ++i )
        result *= (order+1);
      return result ;
    }
  };


  template< class DiscreteFunctionImp, template<class> class SigmaDiscreteFunctionChooserImp >
  class SigmaDiscreteFunctionSelector
  {
  public:
    typedef DiscreteFunctionImp                                      DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::GridPartType              GridPartType;
    typedef typename GridPartType::GridType                          GridType;
    static const int polOrder = DiscreteFunctionSpaceType::polynomialOrder;


    typedef typename DiscreteFunctionSpaceType::
      template ToNewDimRange< GridType::dimension * DiscreteFunctionSpaceType::FunctionSpaceType::dimRange >::NewFunctionSpaceType
                                                                     SigmaFunctionSpaceType;

    //- hybrid spaces use PAdaptiveDGSpace
    template <class Grid, int topoId>
    struct SigmaSpaceChooser
    {
      typedef Fem::PAdaptiveDGSpace< SigmaFunctionSpaceType, GridPartType, polOrder > type;
    };

    //- cubes use LegendreDGSpace
    template <class Grid>
    struct SigmaSpaceChooser< Grid, 1 >
    {
      typedef Dune::Fem::LegendreDiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polOrder > type;
    };

    //- cubes use OrthonormalDGSpace
    template <class Grid>
    struct SigmaSpaceChooser< Grid, 0 >
    {
      typedef Dune::Fem::DiscontinuousGalerkinSpace< SigmaFunctionSpaceType, GridPartType, polOrder > type;
    };

    // work arround internal compiler error
    enum { simplexTopoId =  GeometryTypes::cube( GridType::dimension ).id(), //Dune::Impl::SimplexTopology< GridType::dimension >::type::id,
           cubeTopoId    =  GeometryTypes::simplex( GridType::dimension ).id(), //Dune::Impl::CubeTopology< GridType::dimension >::type::id,
           myTopo        =  Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
         };

    enum { topoId = (simplexTopoId == myTopo) ? 0 :  // 0 = simplex, 1 = cube, -1 = hybrid
                      (myTopo == cubeTopoId) ? 1 : -1 };

    typedef typename SigmaSpaceChooser< GridType, topoId >::type         SigmaDiscreteFunctionSpaceType;

  public:
    typedef typename SigmaDiscreteFunctionChooserImp<SigmaDiscreteFunctionSpaceType>::type
                                                                         type;

  };

  // Estimator
  // ---------
  // Template Arguments:
  // UFunction: Type for the function u which is only needed to compute the flux hatK
  //            must implement both evaluate and jacobian
  // SigmaFunction: Main function used in the estimate = grad u + liftings
  //            must implement evaluate
  // DGOperator: Operator providing the flux, needs method
  //             model() to access the diffusion method
  //             numericalFlux(inside, outside,
  //                  uValuesEn, duValuesEn, uValuesNb, duValuesNb,
  //                  fluxEn, dfluxEn, fluxNb, dfluxNb):
  //              method to compute -hatK, only fluxEn and fluxNb is used

  template< class DiscreteFunctionImp, class SigmaDiscreteFunctionImp, class DGOperatorImp >
  class ErrorEstimator
  {
    typedef ErrorEstimator< DiscreteFunctionImp, SigmaDiscreteFunctionImp, DGOperatorImp > ThisType;

  public:
    typedef DiscreteFunctionImp                                      DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::GridPartType              GridPartType;
    typedef typename GridPartType::GridType                          GridType;
    typedef DGOperatorImp                                            DGOperatorType;

    typedef SigmaDiscreteFunctionImp                                 SigmaDiscreteFunctionType;
    typedef typename SigmaDiscreteFunctionType::DiscreteFunctionSpaceType SigmaDiscreteFunctionSpaceType;

    typedef ConstLocalFunction< DiscreteFunctionType >               ConstLocalFunctionType;
    typedef ConstLocalFunctionType                                   LocalFunctionType;

    typedef ConstLocalFunction< SigmaDiscreteFunctionType >          ConstSigmaLocalFunctionType;
    typedef ConstSigmaLocalFunctionType                              SigmaLocalFunctionType;
    typedef typename SigmaLocalFunctionType::RangeType               GradientRangeType;

    typedef typename DiscreteFunctionSpaceType::DomainFieldType      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType       RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::DomainType           DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType            RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType    JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::IteratorType         IteratorType;

    typedef Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType >
                                                                     SigmaConverterType;

    typedef typename GridPartType::IndexSetType                      IndexSetType;
    typedef typename GridPartType::IntersectionIteratorType          IntersectionIteratorType;

    typedef typename GridPartType::IntersectionType                  IntersectionType;

    typedef typename GridType::template Codim< 0 >::Entity           ElementType;
    typedef typename ElementType::Geometry                           GeometryType;

    static const int dimension = GridType::dimension;

    // CACHING
    typedef typename DGOperatorType::FaceQuadratureType              FaceQuadratureType ;
    typedef typename DGOperatorType::VolumeQuadratureType            VolumeQuadratureType ;


    typedef std::vector< double >                                    ErrorIndicatorType;
  protected:
    typedef HierarchicSimplexDGSpace<DiscreteFunctionSpaceType>      HierarchicSimplexDGSpaceType;
    typedef HierarchicCubeDGSpace<DiscreteFunctionSpaceType>         HierarchicCubeDGSpaceType;

    typedef typename HierarchicSimplexDGSpaceType::Type              SimplexDGSpaceType;
    typedef typename HierarchicCubeDGSpaceType::Type                 CubeDGSpaceType ;

    enum PAdaptiveMethodIdentifier {
      none = 0,
      coeffroot = 1,
      prior2p = 2
    };


    const DGOperatorType&            oper_;
    const DiscreteFunctionType&      uh_;
    const SigmaDiscreteFunctionType& sigma_;
    const double                     beta_;
    const DiscreteFunctionSpaceType& dfSpace_;
    GridPartType&                    gridPart_;
    const IndexSetType&              indexSet_;
    GridType&                        grid_;

    std::unique_ptr< SimplexDGSpaceType > dgSimplexSpace_;
    std::unique_ptr< CubeDGSpaceType    > dgCubeSpace_;

    ErrorIndicatorType indicator_;
    ErrorIndicatorType R2_,R1_,Rorth_,smoothness_;
    double             totalIndicator2_,maxIndicator_;
    const double       theta_;
    const bool         maximumStrategy_;
    const int          nonConformityDifference_;
    int                eocId_;
    PAdaptiveMethodIdentifier padaptiveMethod_;


  public:
    // no copy or assignment
    ErrorEstimator( const ErrorEstimator& ) = delete;
    ErrorEstimator& operator= ( const ErrorEstimator& ) = delete;


    ErrorEstimator (DiscreteFunctionType& uh,
                    const DGOperatorType &oper,
                    const SigmaDiscreteFunctionType &sigma,
                    const AdaptationParameters& param = AdaptationParameters() )
     : oper_(oper),
       uh_( uh ),
       sigma_( sigma ),
       beta_(1.),
       dfSpace_( uh_.space() ),
       gridPart_( dfSpace_.gridPart() ),
       indexSet_( gridPart_.indexSet() ),
       grid_( gridPart_.grid() ),
       dgSimplexSpace_(),
       dgCubeSpace_(),
       indicator_( indexSet_.size( 0 )),
       R2_( indexSet_.size( 0 )),
       R1_( indexSet_.size( 0 )),
       Rorth_( indexSet_.size( 0 )),
       totalIndicator2_(0),
       maxIndicator_(0),
       theta_( param.theta() ),
       maximumStrategy_( param.maximumStrategy() ),
       nonConformityDifference_( Dune::Fem::Parameter::getValue("nonconformitydifference",1) ),
       eocId_( -1 ),
       padaptiveMethod_( (PAdaptiveMethodIdentifier) param.padaptiveMethod() )
    {
      clear();
    }

    SimplexDGSpaceType& dgSimplexSpace()
    {
      if( ! dgSimplexSpace_ )
        dgSimplexSpace_.reset( new SimplexDGSpaceType( gridPart_ ) );
      return *dgSimplexSpace_ ;
    }

    CubeDGSpaceType& dgCubeSpace()
    {
      if( ! dgCubeSpace_  )
        dgCubeSpace_.reset( new CubeDGSpaceType( gridPart_ ) );
      return *dgCubeSpace_ ;
    }


    int minimalOrder() const { return 2; }

    // make this class behave as required for a LocalFunctionAdapter
    typedef Dune::Fem::FunctionSpace< double, double, dimension, 7 > FunctionSpaceType;
    void init(const ElementType &en)
    {
      enIndex_ = indexSet_.index(en);
      entity_ = &en;
    }
    template< class PointType >
    void evaluate(const PointType& x,FieldVector<double,7> &ret)
    {
      ret[0] = indicator_[enIndex_];
      ret[1] = R2_[enIndex_];
      ret[2] = R1_[enIndex_];
      ret[3] = Rorth_[enIndex_];
      ret[4] = smoothness_[enIndex_];
      ret[5] = dfSpace_.order(*entity_);
      ret[6] = entity_->geometry().volume();
    }
  private:
    const ElementType *entity_;
    int enIndex_;

  public:
    bool isPadaptive() const
    {
      return (padaptiveMethod_ != none);
    }
    void clear ()
    {
      indicator_.resize( indexSet_.size( 0 ));
      R2_.resize( indexSet_.size( 0 ));
      R1_.resize( indexSet_.size( 0 ));
      Rorth_.resize( indexSet_.size( 0 ));
      smoothness_.resize( indexSet_.size( 0 ));

      std::fill( indicator_.begin(), indicator_.end(), 0);
      std::fill( R2_.begin(), R2_.end(), 0);
      std::fill( R1_.begin(), R1_.end(), 0);
      std::fill( Rorth_.begin(), Rorth_.end(), 0);
      std::fill( smoothness_.begin(), smoothness_.end(), 0);
    }

    template< class RHSFunctionType >
    double estimate ( const RHSFunctionType &rhs )
    {
      if (eocId_ == -1)
      {
        const std::string eocDescription[] = { "R2", "R1", "Rorth" };
        eocId_ = Dune::Fem::FemEoc::addEntry( eocDescription,3);
      }
      clear();
      henMin = 1e10;

      ConstLocalFunctionType uLocal( uh_ );
      ConstLocalFunctionType uOutside( uh_ );

      SigmaLocalFunctionType sigmaLocal( sigma_ );
      SigmaLocalFunctionType sigmaOutside( sigma_ );

      for( const auto& entity : elements( dfSpace_.gridPart() ) )
      {
        auto uGuard = Dune::Fem::bindGuard( uLocal, entity );
        auto sGuard = Dune::Fem::bindGuard( sigmaLocal, entity );

        estimateLocal( rhs, entity, uLocal, sigmaLocal );

        for (const auto& intersection : intersections(gridPart_, entity) )
        {
          if( intersection.neighbor() )
            estimateIntersection( intersection, entity, uLocal, uOutside, sigmaLocal, sigmaOutside );
          else
            estimateBoundary( intersection,entity,uLocal, sigmaLocal);
        }

        calcSmoothness( dfSpace_.order(entity), entity, entity.geometry(), uLocal, smoothness_[indexSet_.index(entity)] );

      }
      std::cout << "*************** hen = " << henMin << std::endl;
      eocIndicator();
      return computeIndicator();
    }

    double computeIndicator()
    {
      totalIndicator2_ = 0.0;
      maxIndicator_ = 0.0;

      if (theta_ == 0)
      {
        for (unsigned int i=0;i<indicator_.size();++i)
        {
          totalIndicator2_ += indicator_[ i ];
          maxIndicator_ = std::max(maxIndicator_,indicator_[i]);
        }
      }
      else if (maximumStrategy_)
      {
        const unsigned int size = indicator_.size();
        for (unsigned int i=0;i<size;++i)
        {
          totalIndicator2_ += indicator_[ i ];
          maxIndicator_ = std::max(maxIndicator_,indicator_[i]);
        }
        // comput global max indicator
        maxIndicator_ = grid_.comm().max( maxIndicator_ );
      }
      else
      {
        typedef std :: multimap<double,int> SortIndicatorType;
        SortIndicatorType sortIndicator_;

        const unsigned int size = indicator_.size();
        for (unsigned int i=0;i<size;++i)
        {
          // if (indicator_[i]<1e-10) continue;
          totalIndicator2_ += indicator_[ i ];
          maxIndicator_ = std::max(maxIndicator_,indicator_[i]);
          sortIndicator_.insert(std::pair<double,int>(indicator_[ i ],i));
        }
        typedef typename SortIndicatorType :: const_reverse_iterator SortIteratorType;
        const SortIteratorType end = sortIndicator_.rend();
        double errorMarked = 0;
        int marked = 0;
        const unsigned int mapSize = sortIndicator_.size();
        for( SortIteratorType it = sortIndicator_.rbegin(); it != end; ++it )
        {
          errorMarked += (*it).first;
          indicator_[ (*it).second ] = -1; // theta>0 and indicator<0: mark
          ++marked;
          if ( marked > theta_*mapSize )
            break;
        }
      }

      return sqrt( totalIndicator2_ );
    }

    void eocIndicator()
    {
      totalIndicator2_ = 0.0;
      maxIndicator_ = 0.0;
      FieldVector<double,3> errorParts(0);

      if (theta_ == 0)
      {
        for (unsigned int i=0;i<indicator_.size();++i)
        {
          errorParts[0] += R2_[i];
          errorParts[1] += R1_[i];
          errorParts[2] += Rorth_[i];
        }
      }
      else
      {
        const unsigned int size = indicator_.size();
        for (unsigned int i=0;i<size;++i)
        {
          if (indicator_[i]<1e-10) continue;
          errorParts[0] += R2_[i];
          errorParts[1] += R1_[i];
          errorParts[2] += Rorth_[i];
          assert( std::abs(indicator_[i] - (R2_[i] + R1_[i] + Rorth_[i]) ) < 1e-8*indicator_[i] );
        }

      }
      errorParts[0] = sqrt(errorParts[0]);
      errorParts[1] = sqrt(errorParts[1]);
      errorParts[2] = sqrt(errorParts[2]);
      Dune::Fem::FemEoc :: setErrors(eocId_, errorParts );

    }

    //! mark all elements due to given tolerance
    bool mark ( const double tolerance ) const
    {
      int marked = 0;
      int hmarked = 0;

      if (tolerance < 0)
      {
        for( const auto& entity : elements( dfSpace_.gridPart() ) )
        {
          grid_.mark( 1, entity );
          ++marked;
        }
      }
      else
      {
        // get local tolerance (note: remove sqrt from error)
        const double localTol2 = ( maximumStrategy_ ) ?
                ( theta_*theta_ * maxIndicator_ ) :
                tolerance*tolerance / (double)indexSet_.size( 0 );

        // loop over all elements
        for( const auto& entity : elements( dfSpace_.gridPart() ) )
        {
          // check local error indicator
          if( (maximumStrategy_ && indicator_[ indexSet_.index( entity ) ] > localTol2)
              ||
              indicator_[ indexSet_.index( entity ) ] == -1 )
          {
            ++marked;
            if ( !isPadaptive() ||
                ( dfSpace_.order(entity) > smoothness_[ indexSet_.index( entity ) ]
                  || dfSpace_.order(entity) == dfSpace_.order() ) )
            {
              ++hmarked;
              grid_.mark( 1, entity );
            }
          }
        }
      }
      if (nonConformityDifference_ >= 0)
        makeNonConformity(dfSpace_.gridPart(),nonConformityDifference_);

      std::cout << "MARKED: " << marked << " from " << indexSet_.size(0) << " - of those " << hmarked << " for h refinement" << std::endl;

      return (marked > 0);
    }
    void closure()
    {
      for( const auto& entity : elements( dfSpace_.gridPart() ) )
      {
        int conf = 0, nconf = 0;
        int lev = entity.level();
        for (const auto& intersection : intersections(gridPart_, entity) )
        {
          if (intersection.neighbor())
          {
            if (intersection.outside().level() > lev)
              ++nconf;
            else
              ++conf;
          }
          else
            ++conf;
        }
        if (conf == 0)
          grid_.mark( 1, entity );
      }
    }
    int newOrder( double tolerance, const ElementType &entity)
    {
      assert( indexSet_.index( entity ) < (int)indicator_.size() );

      const double localTol2 = ( maximumStrategy_ ) ?
              ( theta_*theta_ * maxIndicator_ ) :
              tolerance*tolerance / (double)indexSet_.size( 0 );
      if( indicator_[ indexSet_.index( entity ) ] > localTol2 )
      {
        if (dfSpace_.order(entity) <= smoothness_[ indexSet_.index( entity ) ] )
        {
          return std::min( int(std::ceil(smoothness_[ indexSet_.index( entity ) ])) ,dfSpace_.order() );
          // increase order since smoothness is high
          return std::min( dfSpace_.order(entity)+1,dfSpace_.order() );
        }
        else
        {
          // h ref. was used
          // return std::max( dfSpace_.order(entity)-1,2 );
          return dfSpace_.order(entity);
        }
      }
      // if (indicator_[ indexSet_.index( entity ) ] < maxIndicator_* 0.01 )
      if (smoothness_[ indexSet_.index( entity ) ] < 0 )
      {
        return std::max( dfSpace_.order(entity)-1, minimalOrder() );
      }
      // std::cout << std::endl;

      return dfSpace_.order(entity);
    }

  protected:
    //! caclulate error on element
    template< class RHSFunctionType >
    void estimateLocal ( const RHSFunctionType &rhs, const ElementType &entity,
                         const LocalFunctionType &uLocal, const SigmaLocalFunctionType &sigmaLocal )
    {
      const typename ElementType :: Geometry &geometry = entity.geometry();

      const double volume = geometry.volume();
      const double h2 = ( (dimension == 2 ? volume : std :: pow( volume, 2.0 / (double)dimension )) ) /
                        ( uLocal.order()*uLocal.order() );
      const int index = indexSet_.index( entity );

      const int quadOrder = 2*(dfSpace_.order( entity ) )+2;
      VolumeQuadratureType quad( entity, quadOrder );

      for( const auto qp : quad )
      {
        const auto x = qp.position();
        // R_2 = h^2 |f+laplace u|^2

        RangeType y(0.),tmp(0.);

        const DomainType global = geometry.global( x );
        rhs.f(global,y);

        divergence( entity, geometry, quad, qp.index(), h2, uLocal, sigmaLocal, tmp);

        y+=tmp;

        const double weight = qp.weight() * geometry.integrationElement( x );
        indicator_[ index ] +=h2*weight * (y * y);

        R2_[ index ] += h2*weight * (y * y);
      }

    }

    void calcSmoothness( int pOrder, const ElementType &entity,
                         const typename ElementType::Geometry &geometry,
                         const LocalFunctionType &uLocal, double &smoothness )
    {
      if (geometry.type().isSimplex())
        calcSmoothness<HierarchicSimplexDGSpaceType>(dgSimplexSpace(),pOrder,entity,geometry,uLocal,smoothness);
      else
        calcSmoothness<HierarchicCubeDGSpaceType>(dgCubeSpace(),pOrder,entity,geometry,uLocal,smoothness);
    }
    template <class HDGSpaceType>
    void calcSmoothness( const typename HDGSpaceType::Type &dgSpace,
                         int pOrder, const ElementType &entity,
                         const typename ElementType::Geometry &geometry,
                         const LocalFunctionType &uLocal, double &smoothness )
    {
      typedef typename HDGSpaceType::Type DGSpaceType;
      typedef Dune::Fem::LocalMassMatrix < DGSpaceType, VolumeQuadratureType >  LocalMassMatrixType;
      typedef Dune::Fem::TemporaryLocalFunction< DGSpaceType > DGLFType;

      LocalMassMatrixType localMassMatrix( dgSpace, 2*dgSpace.order(entity) ) ;
      DGLFType udg( dgSpace );
      VolumeQuadratureType quad( entity, 2*(dgSpace.order(entity) )+1 );

      udg.init( entity );
      udg.clear();

      for( const auto qp : quad )
      {
        const auto x = qp.position();
        RangeType uh;
        uLocal.evaluate( qp, uh );
        uh *= qp.weight() * geometry.integrationElement( x );
        udg.axpy( qp, uh );
      }

      localMassMatrix.applyInverse( udg );

      assert(RangeType::dimension == 1);
      // only dofs up to bstart[pOrder] should be used
      // bool correctProjection = true;
      for (int i=HDGSpaceType::size(pOrder); i<udg.numDofs(); ++i)
      {
        if ( std::abs(udg[i]) > 1e-10 )
        {
          std::cout << "for i=" << HDGSpaceType::size(pOrder) << " to " << udg.numDofs() ;
          std::cout << " projection wrong for polorder = " << pOrder
                    << " dof[" << i << "] = " << udg[i] << std::endl;
          //correctProjection = false;
        }
      }
      // assert( correctProjection );

      if (padaptiveMethod_ == coeffroot)
      {
        double a = 0;
        double vol = geometry.volume();
        for (int i=HDGSpaceType::size(pOrder-1) ; i<udg.numDofs(); ++i)
        {
          assert( udg.numDofs() > i );
          a += std::abs(udg[i]);
        }
        a *= sqrt(vol);
        if (std::abs(a)<1e-8 || pOrder==1)
          smoothness = -1; // pOrder+1;
        else
          smoothness = log( (2.*pOrder+1)/(2.*a*a) ) / (2.*log(pOrder))-1;
      }
      else if (padaptiveMethod_ == prior2p)
      {
        if (pOrder<3)
          smoothness = pOrder+1;
        else
        {
          typedef Dune::Fem::H1Norm<GridPartType> H1NormType ;
          typedef typename H1NormType :: IntegratorType  IntegratorType ;

          IntegratorType integrator(2*pOrder);
          typedef typename H1NormType :: template FunctionJacobianSquare< DGLFType >  FctJacType;

          FctJacType udg2( udg );
          typename FctJacType::RangeType eta1(0),eta2(0);
          for (int i=0 ; i<HDGSpaceType::size(pOrder-2); ++i)
            udg[i] = 0;
          integrator.integrate( entity, udg2, eta2 ); // square of H1 norm
          for (int i=HDGSpaceType::size(pOrder-2) ; i<HDGSpaceType::size(pOrder-1); ++i)
            udg[i] = 0;
          integrator.integrate( entity, udg2, eta1 ); // square of H1 norm
          if (std::abs(eta1) < 1e-12 || std::abs(eta2) < 1e-12)
          {
            if (std::abs(eta1) < 1e-12 && std::abs(eta2) < 1e-12)
              smoothness = -1; // the solution is apparently smooth but could be
                                // approximated with a much lower polynomial (2 order
                                // are small) - so coarsen
            else
              smoothness = pOrder+1;
          }
          else
          {
            smoothness = -0.5*std::log(eta1[0]/eta2[0]) / std::log( double(pOrder-1)/double(pOrder-2) );
            assert( smoothness == smoothness && smoothness < 100.);
          }
        }
      }
    }

    void estimateBoundary( const IntersectionType &intersection,
                           const ElementType &inside,
                           const LocalFunctionType &uInside, const SigmaLocalFunctionType &sigmaInside )
    {
      double volume = inside.geometry().volume() ;
      const double h = ( (dimension == 1 ? volume : std::pow(0.5* volume, 1.0 / (double)dimension )) ) /
                       ( uInside.order() );
      const int quadOrder = 2*(dfSpace_.order()+1);
      const int insideIndex = indexSet_.index( inside );
      const  FaceQuadratureType quadInside( gridPart_, intersection, quadOrder, FaceQuadratureType::INSIDE );
      const int numQuadraturePoints = quadInside.nop();

      // obtain all required function values on intersection
      uValuesEn.resize( numQuadraturePoints );
      uValuesNb.resize( numQuadraturePoints );
      uInside.evaluateQuadrature( quadInside, uValuesEn );

      typedef RangeType           RangeTuple;
      typedef SigmaConverterType JacobianTuple;

      typedef QuadratureContext< ElementType, IntersectionType, FaceQuadratureType > ContextType;
      typedef LocalEvaluation< ContextType, RangeTuple, JacobianTuple > LocalEvaluationType;

      ContextType cLocal( inside, intersection, quadInside, volume );

      oper_.boundaryValues(cLocal, uValuesNb);
      duValuesEn.resize( numQuadraturePoints );
      uInside.evaluateQuadrature( quadInside, duValuesEn );
      fluxEn.resize( numQuadraturePoints );
      dfluxEn.resize( numQuadraturePoints );

      oper_.boundaryFlux( cLocal,
                          uValuesEn, duValuesEn, uValuesNb,
                          fluxEn, dfluxEn );

      sigmaValuesEn.resize( numQuadraturePoints );
      sigmaInside.evaluateQuadrature( quadInside, sigmaValuesEn );

      double errorInside = 0.0;
      double faceVol= 0.;
      faceVol = intersection.geometry().volume();

      volume/=3*faceVol;

      for( const auto qp : quadInside )
      {
        const auto idx = qp.index();
        const auto& x = qp.localPosition();

        DomainType unitNormal = intersection.integrationOuterNormal( x );
        const double integrationElement = unitNormal.two_norm();

        unitNormal/=integrationElement;

        // R_1 = h| (d u_l * n_l) + (d u_r * n_r) |^2 = h| (d u_l - d u_r) * n_l |^2
        JacobianRangeType AJacEn;
        fluxEn[idx] /= integrationElement;

        const SigmaConverterType sigmaValueEn( sigmaValuesEn[idx] );
        LocalEvaluationType local( cLocal, uValuesEn[idx], sigmaValueEn );

        oper_.model().diffusion(local[idx], uValuesEn[idx], sigmaValueEn, AJacEn);

        // note that flux=-hatA therefore we compute -hatA+Agrad u
        AJacEn.umv( unitNormal, fluxEn[idx]);

        // R_orth = h^{-1} |u_l-u_r|
        RangeType jump;
        jump=uValuesEn[idx];
        jump-=uValuesNb[idx];

        const auto& weight = qp.weight();

        errorInside  += weight *h* (fluxEn[idx] * fluxEn[idx]) *integrationElement;
        errorInside  += weight *1./h* (jump * jump) *integrationElement;

        R1_[ insideIndex ] += weight *h* (fluxEn[idx]*fluxEn[idx]) *integrationElement;
        Rorth_[ insideIndex ] += weight *1./h* (jump * jump) *integrationElement;
      }
      if( errorInside > 0.0 )
        indicator_[ insideIndex ] +=  errorInside;
    }
    //! caclulate error on boundary intersections
    void estimateIntersection ( const IntersectionType &intersection,
                                const ElementType &inside,
                                const LocalFunctionType &uInside,
                                LocalFunctionType& uOutside,
                                const SigmaLocalFunctionType &sigmaInside,
                                SigmaLocalFunctionType& sigmaOutside )
    {
      const ElementType& outside = intersection.outside();

      const int quadOrder = 2*(dfSpace_.order()+1);

      const int insideIndex = indexSet_.index( inside );
      const int outsideIndex = indexSet_.index( outside );

      const bool isOutsideInterior = (outside.partitionType() == InteriorEntity);
      if( !isOutsideInterior || (insideIndex < outsideIndex) )
      {
        auto uGuard = Dune::Fem::bindGuard( uOutside, outside );
        auto sGuard = Dune::Fem::bindGuard( sigmaOutside, outside );

        // const double volume = std::max( inside.geometry().volume() , outside.geometry().volume() );
        // const double h = 2.*volume / intersection.geometry().volume();
        const double volume = ( inside.geometry().volume() + outside.geometry().volume() );
        const double h = (dimension == 1 ? volume : std::pow(0.5* volume, 1.0 / (double)dimension )) /
                         ( uInside.order() );

        double errorInside, errorOutside;
        if( ! intersection.conforming() )
        {
          estimateIntersection< false > ( intersection, inside, outside, quadOrder, insideIndex, outsideIndex,
                                          uInside, sigmaInside, uOutside, sigmaOutside,
                                          h,volume,
                                          errorInside, errorOutside);
        }
        else
        {
          estimateIntersection< true > ( intersection, inside, outside, quadOrder, insideIndex, outsideIndex,
                                         uInside, sigmaInside, uOutside, sigmaOutside,
                                         h,volume,
                                         errorInside, errorOutside);
        }
        if( errorInside > 0.0 )
        {
          indicator_[ insideIndex ] +=  errorInside;
          if( isOutsideInterior )
            indicator_[ outsideIndex ] +=  errorOutside;
        }
      }
    }

    //! caclulate error on element intersections
    template< bool conforming >
    void estimateIntersection ( const IntersectionType &intersection,
                                const ElementType &inside, const ElementType &outside,
                                const int quadOrder,
                                const int insideIndex,
                                const int outsideIndex,
                                const LocalFunctionType &uInside,
                                const SigmaLocalFunctionType &sigmaInside,
                                const LocalFunctionType &uOutside,
                                const SigmaLocalFunctionType &sigmaOutside,
                                const double h, double volume,
                                double &errorInside, double &errorOutside)
    {
      // use IntersectionQuadrature to create appropriate face quadratures
      typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
      typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;

      // create intersection quadrature
      IntersectionQuadratureType interQuad( gridPart_, intersection, quadOrder );

      // get appropriate references
      const QuadratureImp &quadInside  = interQuad.inside();
      const QuadratureImp &quadOutside = interQuad.outside();
      const int numQuadraturePoints = quadInside.nop();

      // obtain all required function values on intersection
      uValuesEn.resize( numQuadraturePoints );
      uValuesNb.resize( numQuadraturePoints );
      uInside.evaluateQuadrature( quadInside, uValuesEn );
      uOutside.evaluateQuadrature( quadOutside, uValuesNb );
      duValuesEn.resize( numQuadraturePoints );
      duValuesNb.resize( numQuadraturePoints );
      uInside.evaluateQuadrature( quadInside, duValuesEn );
      uOutside.evaluateQuadrature( quadOutside, duValuesNb );
      fluxEn.resize( numQuadraturePoints );
      dfluxEn.resize( numQuadraturePoints );
      fluxNb.resize( numQuadraturePoints );
      dfluxNb.resize( numQuadraturePoints );

      typedef RangeType           RangeTuple;
      typedef SigmaConverterType  JacobianTuple;
      typedef QuadratureContext< ElementType, IntersectionType, QuadratureImp > ContextType;
      typedef LocalEvaluation< ContextType, RangeTuple, JacobianTuple > LocalEvaluationType;

      ContextType cLeft( inside, intersection, quadInside, volume );
      ContextType cRight( outside, intersection, quadOutside, volume );

      oper_.numericalFlux(cLeft, cRight,
                          uValuesEn, duValuesEn, uValuesNb, duValuesNb,
                          fluxEn, dfluxEn, fluxNb, dfluxNb
                          );

      sigmaValuesEn.resize( numQuadraturePoints );
      sigmaValuesNb.resize( numQuadraturePoints );
      sigmaInside.evaluateQuadrature( quadInside, sigmaValuesEn );
      sigmaOutside.evaluateQuadrature( quadOutside, sigmaValuesNb );

      errorInside = 0.0;
      errorOutside = 0.0;
      double faceVol= 0.;
      faceVol = intersection.geometry().volume();

      volume/=3*faceVol;
      for( const auto qp : quadInside )
      {
        const auto idx = qp.index();
        const auto& x = qp.localPosition();

        DomainType unitNormal = intersection.integrationOuterNormal( x );
        const double integrationElement = unitNormal.two_norm();

        unitNormal/=integrationElement;

        // R_1 = h| (d u_l * n_l) + (d u_r * n_r) |^2 = h| (d u_l - d u_r) * n_l |^2
        JacobianRangeType AJacEn,AJacNb;
        fluxEn[idx] /= integrationElement;
        fluxNb[idx] /= integrationElement;

        const SigmaConverterType sigmaValueEn( sigmaValuesEn[idx] );
        const LocalEvaluationType left( cLeft, uValuesEn[idx], sigmaValueEn );
        const SigmaConverterType sigmaValueNb( sigmaValuesNb[idx] );
        const LocalEvaluationType right( cRight, uValuesNb[idx], sigmaValueNb );

        oper_.model().diffusion(left[idx],  uValuesEn[idx], sigmaValueEn, AJacEn);
        oper_.model().diffusion(right[idx], uValuesNb[idx], sigmaValueNb, AJacNb);

        // note that flux=-hatA therefore we compute -hatA+Agrad u
        AJacEn.umv( unitNormal, fluxEn[idx]);
        AJacNb.umv( unitNormal, fluxNb[idx]);

        // R_orth = h^{-1} |u_l-u_r|
        RangeType jump;
        jump=uValuesEn[idx];
        jump-=uValuesNb[idx];

        const auto& weightIn = qp.weight();
        const auto& weightOut = quadOutside.weight( idx );

        errorInside  += weightIn *h* (fluxEn[idx] * fluxEn[idx]) *integrationElement;
        errorInside  += weightIn *1./h* (jump * jump) *integrationElement;
        errorOutside += weightOut *h* (fluxNb[idx] * fluxNb[idx]) *integrationElement;
        errorOutside += weightOut *1./h* (jump * jump) *integrationElement;

        R1_[ insideIndex ] += weightIn *h* (fluxEn[idx]*fluxEn[idx]) *integrationElement;
        Rorth_[ insideIndex ] += weightIn *1./h* (jump * jump) *integrationElement;
        R1_[ outsideIndex ] += weightOut *h* (fluxNb[idx]*fluxNb[idx]) *integrationElement;
        Rorth_[ outsideIndex ] += weightOut *1./h* (jump * jump) *integrationElement;
      }
    }

    mutable double henMin;
    template<class Quadrature>
    void divergence(const ElementType &entity,
                    const GeometryType & geo,
                    const Quadrature &quad,
                    const int qp,
                    double h2,
                    const LocalFunctionType &u_h,
                    const SigmaLocalFunctionType &sigma_h,
                    RangeType& result
                    ) const
    {
#if 0
      // This version would only work for affine geometry and piecewise constant K
      // to do this we need to implement the jacobian on local function adapted
      // describing sigma
      if( geo.affine() )
      {
        result = 0;
        typename SigmaLocalFunctionType::JacobianRangeType hessian;
        sigma_h.jacobian( quad[qp], hessian );
        /*
        typename LocalFunctionType::HessianRangeType hessian;
        for( int r = 0; r < RangeType::dimension; ++r )
        {
          for( int i = 0; i < dimension; ++i )
            result[ r ] += hessian[ r ][ i ][ i ];
        }
        */
        assert( RangeType::dimension == 1 );
        for( int i = 0; i < dimension; ++i )
          result[ 0 ] += hessian[ i ][ i ];
      }
      else
#endif
      typedef PointContext< ElementType > PointContextType;
      const double volume = geo.volume();
      {
        // finite difference approximation
        RangeType ux0,ux1;

        GradientRangeType sigmax0,sigmax1;
        JacobianRangeType Asigmax0,Asigmax1;

        DomainType xGlobal = geo.global( quad.point(qp ) );
        double hen = std::max(1e-12,h2*h2);
        result = 0;
        for( int i = 0; i < dimension; ++i )
        {
          DomainType xgl0(xGlobal);
          DomainType xgl1(xGlobal);
          xgl0[i] -= hen;
          xgl1[i] += hen;
          DomainType x0 = geo.local( xgl0 );
          DomainType x1 = geo.local( xgl1 );

          u_h.evaluate(x0,ux0);
          u_h.evaluate(x1,ux1);
          sigma_h.evaluate(x0,sigmax0);
          sigma_h.evaluate(x1,sigmax1);

          {
            const PointContextType local( entity, xgl0, x0, volume );
            const SigmaConverterType jac0( sigmax0 );
            oper_.model().diffusion(local, ux0, jac0, Asigmax0);
          }
          {
            const PointContextType local( entity, xgl1, x1, volume );
            const SigmaConverterType jac1( sigmax1 );
            oper_.model().diffusion(local, ux1, jac1, Asigmax1);
          }
          for( int r = 0; r < RangeType::dimension; ++r )
          {
            result[r] += (Asigmax1[r][i] - Asigmax0[r][i])/ (2.*hen);
          }
          henMin = std::min(henMin,hen);
        }
      }
    }
    mutable std::vector< RangeType > uValuesEn;
    mutable std::vector< RangeType > uValuesNb;
    mutable std::vector< JacobianRangeType > duValuesEn;
    mutable std::vector< JacobianRangeType > duValuesNb;
    mutable std::vector< RangeType > fluxEn;
    mutable std::vector< JacobianRangeType > dfluxEn;
    mutable std::vector< RangeType > fluxNb;
    mutable std::vector< JacobianRangeType > dfluxNb;
    mutable std::vector< GradientRangeType > sigmaValuesEn;
    mutable std::vector< GradientRangeType > sigmaValuesNb;
  };

}
}

#endif
