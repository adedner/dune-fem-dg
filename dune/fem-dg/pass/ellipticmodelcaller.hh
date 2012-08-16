#ifndef DUNE_ELLIPTICDISCRETEMODELCALLER_HH
#define DUNE_ELLIPTICDISCRETEMODELCALLER_HH

#include <utility>
#include <memory>

#include <dune/common/fvector.hh> 

#include <dune/fem/pass/callerutility.hh>
#include <dune/fem/pass/ellipticdiscretemodel.hh>
#include <dune/fem/pass/modelcallerdefault.hh>

#include <dune/fem/misc/boundaryidentifier.hh>

namespace Dune
{

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class EllipticDiscreteModelCaller
  : public Fem::DiscreteModelCallerDefault< DiscreteModelImp, ArgumentImp, SelectorImp >
  {
    typedef EllipticDiscreteModelCaller< DiscreteModelImp, ArgumentImp, SelectorImp > ThisType;
    typedef Fem::DiscreteModelCallerDefault< DiscreteModelImp, ArgumentImp, SelectorImp > BaseType; 

  public:
    typedef DiscreteModelImp DiscreteModelType;
    typedef ArgumentImp TotalArgumentType;
    typedef SelectorImp SelectorType;

    typedef typename DiscreteModelType::Traits Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection Intersection;
    typedef typename GridType::template Codim<0>::Entity Entity;

    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;

    typedef Fem::Filter<TotalArgumentType, SelectorType> FilterType;
    typedef typename FilterType::ResultType DiscreteFunctionTupleType;
    typedef Fem::LocalFunctionCreator<DiscreteFunctionTupleType> LFCreator;
    typedef typename LFCreator::ResultType LocalFunctionTupleType;
    typedef Fem::Creator<
      Fem::RangeTypeEvaluator, LocalFunctionTupleType> RangeCreator;
    typedef typename RangeCreator::ResultType RangeTupleType;
    typedef Fem::Creator<
      Fem::JacobianRangeTypeEvaluator, LocalFunctionTupleType> JacobianCreator;
    typedef typename JacobianCreator::ResultType JacobianRangeTupleType;

    typedef Fem::BoundaryIdentifier BoundaryIdentifierType;

    EllipticDiscreteModelCaller( DiscreteModelType &problem )
    : problem_( problem ),
      valuesEn_( RangeCreator::apply() ),
      valuesNeigh_( RangeCreator::apply() ),
      jacobians_( JacobianCreator::apply() )
    {}

    //! return true when a mass matrix has to be build  
    bool hasMass () const  { return false; }

    //! evaluate mass matrix factor 
    template <class MassFactorType>
    void mass(const Entity& en,
              const VolumeQuadratureType& quad,
              const int quadPoint,
              MassFactorType& m)
    {
      m = 0;
      abort();
    }

    void setEntity ( const Entity &entity ) 
    {
      BaseType::setEntity( entity );
      problem_.setEntity( entity );
    }

    void setNeighbor( const Entity &neighbor )
    {
      BaseType::setNeighbor( neighbor );
      problem_.setNeighbor( neighbor );
    }

    // Ensure: entities set correctly before call
    template <class QuadratureType, class CoefficientType>
    void evaluateCoefficientFace(const Intersection& intersection,
                                 const QuadratureType& quadInner, 
                                 const QuadratureType& quadOuter, 
                                 const int quadPoint,
                                 CoefficientType& coeffLeft, 
                                 CoefficientType& coeffRight) 
    {
      // evaluate data functions 
      this->evaluateQuad( quadInner, quadPoint,
                         this->data_->localFunctionsSelf(), this->valuesEn_);
      this->evaluateQuad( quadOuter, quadPoint,
                         this->data_->localFunctionsNeigh(), this->valuesNeigh_);

      problem_.coefficientFace(intersection, this->time_, 
                               quadInner, quadOuter, quadPoint,     
                               this->valuesEn_, 
                               this->valuesNeigh_,
                               coeffLeft,coeffRight);
    }

    // Ensure: entities set correctly before call
    template <class QuadratureType, class CoefficientType>
    void evaluateCoefficientBoundary(const Intersection& intersection,
                                     const QuadratureType& quadInner, 
                                     const int quadPoint,
                                     CoefficientType& coeff) 
    {
      this->evaluateQuad( quadInner, quadPoint,
                          this->data_->localFunctionsSelf(), this->valuesEn_);

      problem_.coefficient(this->data_->self(), 
                           this->time_, 
                           quadInner, quadPoint, 
                           this->valuesEn_, coeff);
    }
      
    template< class CoefficientType >
    void evaluateCoefficient( const Entity &entity,
                              const VolumeQuadratureType &quad, int quadPoint,
                              CoefficientType &coeff )
    {
      evaluateQuad( quad, quadPoint,
                    this->data_->localFunctionsSelf(), this->valuesEn_ );
      problem_.coefficient( entity, this->time_, quad, quadPoint,
                            this->valuesEn_, coeff );
    }

    template <class QuadratureType> 
    BoundaryIdentifierType
    boundaryValue(const Intersection& intersection,
                  const QuadratureType& quad,
                  const int quadPoint,
                  RangeType& boundaryValue) 
    {
      typedef typename Intersection::LocalGeometry Geometry;

      this->evaluateQuad( quad, quadPoint,
                          this->data_->localFunctionsSelf(), this->valuesEn_);
      return problem_.boundaryValue(intersection, this->time_, quad, quadPoint,
                                    this->valuesEn_, boundaryValue);
    }

    void rightHandSide(const Entity& en, 
                       const VolumeQuadratureType& quad, 
                       const int quadPoint, 
                       RangeType& res) 
    {
      this->evaluateQuad( quad, quadPoint,
                         this->data_->localFunctionsSelf(), this->valuesEn_);
      this->evaluateJacobianQuad( quad, quadPoint, 
                                  this->data_->localFunctionsSelf(),
                                  this->jacobians_ );
      problem_.rightHandSide(en, this->time_, quad, quadPoint, 
                             this->valuesEn_,
                             this->jacobians_,
                             res);
    }

  private:
    // our problem 
    DiscreteModelType& problem_;

  protected:
    RangeTupleType valuesEn_;
    RangeTupleType valuesNeigh_;
    JacobianRangeTupleType jacobians_;
  }; // end EllipticDiscreteModelCaller 

} // end namespace Dune 
#endif
