#ifndef ERROR_STOKESESTIMATOR_HH
#define ERROR_STOKESESTIMATOR_HH

//- Dune includes
// #include <dune/grid/common/referenceelements.hh>

#include "poissonestimator.hh"

namespace Dune
{
namespace Fem
{

  // Estimator
  // ---------
  // Template Arguments:
  // UFunction: Type for the function u which is only needed to compute the flux hatK
  //            must implement both evaluate and jacobian
  // SigmaFunction: Main function used in the estimate = grad u + liftings
  //            must implement evaluate
  // DGOperator: Operator providing the flux, needs method
  //             model() to access the diffusion method
  //             numericalFlux(dfSpace_.gridPart(),
  //                  intersection, inside, outside, 0, quadInside, quadOutside,
  //                  uValuesEn, duValuesEn, uValuesNb, duValuesNb,
  //                  fluxEn, dfluxEn, fluxNb, dfluxNb):
  //              method to compute -hatK, only fluxEn and fluxNb is used

  template< class DiscreteFunctionImp, class SigmaDiscreteFunctionImp, class DGOperatorImp >
  class StokesErrorEstimator
    : public ErrorEstimator< DiscreteFunctionImp, SigmaDiscreteFunctionImp, DGOperatorImp >
  {
    typedef ErrorEstimator< DiscreteFunctionImp, SigmaDiscreteFunctionImp, DGOperatorImp >
                                                                     BaseType;

  public:
    typedef typename BaseType::DGOperatorType                        DGOperatorType;
    typedef typename BaseType::DiscreteFunctionType                  DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::LocalFunctionType         LocalFunctionType;
    typedef typename BaseType::SigmaDiscreteFunctionType             SigmaDiscreteFunctionType;
    typedef typename BaseType::LocalFunctionType                     SigmaLocalFunctionType;

    typedef typename DiscreteFunctionSpaceType::DomainFieldType      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType       RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::DomainType           DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType            RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType    JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType::GridPartType         GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType         IteratorType;

    typedef typename GridPartType::GridType                          GridType;
    typedef typename GridPartType::IndexSetType                      IndexSetType;
    typedef typename GridPartType::IntersectionIteratorType          IntersectionIteratorType;

    typedef typename GridPartType::IntersectionType                  IntersectionType;

    typedef typename GridType::template Codim< 0 >::Entity           ElementType;
#if not DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    typedef typename GridType::template Codim< 0 >::EntityPointer    ElementPointerType;
#endif
    typedef typename ElementType::Geometry                           GeometryType;

    static const int dimension = GridType::dimension;
    typedef FieldMatrix<double,dimension,dimension>                  JacobianInverseType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >          ElementQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 1 >          FaceQuadratureType;
    typedef std::vector< double >                                    ErrorIndicatorType;

  private:
    using BaseType::uh_;
    using BaseType::sigma_;
    using BaseType::eocId_;
    using BaseType::dfSpace_;
    using BaseType::gridPart_;
    using BaseType::R2_;
    using BaseType::R1_;
    using BaseType::Rorth_;
    using BaseType::indicator_;
    using BaseType::theta_;
    typename BaseType::ErrorIndicatorType Rdiv_;
  public:
    StokesErrorEstimator (DiscreteFunctionType& df,
                          const DGOperatorType &oper,
                          const SigmaDiscreteFunctionType &sigma,
                          const AdaptationParameters& param = AdaptationParameters() )
    : BaseType(df,oper,sigma,param),
      Rdiv_( this->indexSet_.size( 0 ))
    {
    }

    template< class RHSFunctionType >
    double estimate ( const RHSFunctionType &rhs )
    {
      if (eocId_ == -1)
      {
        const std::string eocDescription[] = { "R2", "R1", "Rorth", "Rdiv" };
        eocId_ = Dune::Fem::FemEoc::addEntry( eocDescription,4);
      }

      BaseType::clear();
      Rdiv_.resize( this->indexSet_.size( 0 ));
      std::fill( Rdiv_.begin(), Rdiv_.end(), 0);

      for( const auto& entity : elements( dfSpace_.gridPart() ) )
      {
        const LocalFunctionType uLocal = uh_.localFunction( entity );
        const SigmaLocalFunctionType sigmaLocal = sigma_.localFunction( entity );

        BaseType::estimateLocal( rhs, entity, uLocal, sigmaLocal );
        estimateLocal( rhs, entity, uLocal, sigmaLocal );

        for (const auto& intersection : intersections(gridPart_, entity) )
        {
          if( intersection.neighbor() )
            BaseType::estimateIntersection( intersection, entity, uLocal, sigmaLocal );
           else
             BaseType::estimateBoundary( intersection,entity,uLocal,sigmaLocal);
        }
      }
      indicatorEoc();

      return BaseType::computeIndicator();
    }

    void indicatorEoc()
    {
      FieldVector<double,4> errorParts(0);
      if (theta_ == 0)
      {
        for (unsigned int i=0;i<indicator_.size();++i)
        {
          errorParts[0] += R2_[i];
          errorParts[1] += R1_[i];
          errorParts[2] += Rorth_[i];
          errorParts[3] += Rdiv_[i];
        }
      }
      else
      {
        for (unsigned int i=0;i<indicator_.size();++i)
        {
          errorParts[0] += R2_[i];
          errorParts[1] += R1_[i];
          errorParts[2] += Rorth_[i];
          errorParts[3] += Rdiv_[i];
        }
      }

      errorParts[0] = sqrt(errorParts[0]);
      errorParts[1] = sqrt(errorParts[1]);
      errorParts[2] = sqrt(errorParts[2]);
      errorParts[3] = sqrt(errorParts[3]);
      Dune::Fem::FemEoc :: setErrors(eocId_, errorParts );
    }



  private:
    //! caclulate error on element
    template< class RHSFunctionType >
    void estimateLocal ( const RHSFunctionType &rhs, const ElementType &entity,
                         const LocalFunctionType &uLocal, const SigmaLocalFunctionType &sigmaLocal )
    {
      const typename ElementType :: Geometry &geometry = entity.geometry();

      const double volume = geometry.volume();
      const double h2 = (dimension == 2 ? volume : std :: pow( volume, 2.0 / (double)dimension ))/
                        ( uLocal.order()*uLocal.order() );

      const int index = this->indexSet_.index( entity );

      typename SigmaLocalFunctionType::RangeType sigmaVal;
      ElementQuadratureType quad( entity, 2*(dfSpace_.order() + 2) );
      JacobianInverseType inv;

      for( const auto qp : quad )
      {
        const auto& x = qp.position();
        inv = geometry.jacobianInverseTransposed(x);
        JacobianRangeType uJac(0.);
        uLocal.jacobian( qp,uJac );

        double divergence = 0;
        for (int i=0;i<uJac.rows;++i)
        {
          divergence += uJac[i][i];
        }
        const double weight = qp.weight() * geometry.integrationElement( x );
        indicator_[ index ] += h2*weight * (divergence * divergence);

        Rdiv_[ index ] += h2*weight * (divergence * divergence);
      }

    }

  };

}
}

#endif
