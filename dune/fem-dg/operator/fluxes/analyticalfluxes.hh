#ifndef DUNE_FEM_DG_ANALYTICALFLUXES_HH
#define DUNE_FEM_DG_ANALYTICALFLUXES_HH

/**
 * \brief Defines an interface for analytical fluxes
 *
 * \todo Improve this class. At the moment it is only for documentation purposes
 * since EulerAnalyticalFluxes is the only implementation...
 *
 * \ingroup AnalyticalFluxes
 */
template< class Model >
class AnalyticalFluxBase
{
  public:

  typedef Model                         ModelType;
  typedef typename ModelType::FieldType FieldType;
  typedef typename ModelType::RangeType RangeType;

  AnalyticalFluxBase( const ModelType& mod )
    : model_( mod )
  {}


  template <class FluxRangeType>
  void analyticalFlux( const FieldType gamma,
                       const RangeType& u,
                       FluxRangeType& f) const
  {
  }

  template <class FluxRangeType>
  void jacobian( const FieldType gamma,
                 const RangeType& u,
                 const FluxRangeType& du,
                 RangeType& A) const
  {
  }

  protected:
  const ModelType& model_;
};

#endif
