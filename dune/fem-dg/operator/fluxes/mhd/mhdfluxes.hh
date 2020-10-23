#ifndef DUNE_FEM_DG_MHDFLUXES_HH
#define DUNE_FEM_DG_MHDFLUXES_HH
#warning "Using Mhd NumFluxes"
#define USE_MHDFLUXES_INLINE

#include <cmath>

#ifdef COUNT_FLOPS
#include <double.h>
#endif

#include <dune/common/dynvector.hh>

#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/operator/fluxes/mhd/mhd_eqns.hh>
#include <dune/fem-dg/operator/fluxes/rotator.hh>

namespace Dune
{
namespace Fem
{

  // Dai-Woodward
  template <class Model>
  class DWNumFlux;

  // HLLEM
  template <class Model>
  class HLLEMNumFlux;

  // ************************************************
  template <int dimDomain>
  class ConsVec : public Dune :: FieldVector< double, dimDomain+2>
  {
  public:
    explicit ConsVec (const double& t) : Dune :: FieldVector<double,dimDomain+2>(t) {}
    ConsVec () : Dune :: FieldVector<double,dimDomain+2>(0) {}
  };

  namespace Mhd {
    typedef enum { DW, HLLEM } MhdFluxType;
  }

  // ***********************
  template < class Model, Mhd :: MhdFluxType fluxtype >
  class MHDNumFluxBase
  {
  public:
    typedef Model                                       ModelType;
    static const int dimDomain = Model::dimDomain;
    static const int dimRange = Model::dimRange;
    typedef typename Model::Traits                      Traits;
    typedef typename Traits::GridType                   GridType;
    typedef typename GridType::ctype                    ctype;

    typedef typename Traits::DomainType                 DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::FluxRangeType              FluxRangeType;

    typedef ::Mhd::MhdSolver MhdSolverType;

    typedef double value_t[ 9 ];

  protected:
    MHDNumFluxBase(const Model& mod)
     : model_(mod),
       eos( MhdSolverType::Eosmode::me_ideal ),
       numFlux_(eos,  1.4, 1.0 ), // gamma=mod.c_p() * mod.c_v_inv(), R=1.0 ),
       rotU_(1),
       rotB_(4)
    {
      if( fluxtype == Mhd :: HLLEM )
      {
        if( Dune :: Fem :: Parameter :: verbose () )
          std::cout << "Choosing HLLEM Flux " << std::endl;

        numFlux_.init(::Mhd :: MhdSolver :: mf_rghllem );
      }
    }

  public:
    // Return value: maximum wavespeed*length of integrationOuterNormal
    // gLeft,gRight are fluxed * length of integrationOuterNormal
    double
    numericalFlux( const DomainType& normal,
                   const RangeType& uLeft,
                   const RangeType& uRight,
                   RangeType& gLeft) const
    {
      RangeType ul(uLeft);
      RangeType ur(uRight);

      rotU_.rotateForth(ul, normal);
      rotB_.rotateForth(ul, normal);
      rotU_.rotateForth(ur, normal);
      rotB_.rotateForth(ur, normal);

      value_t res;
      const double dummy[ 3 ] = { 0, 0, 0 };

      const double ldt = numFlux_(ul, ur, dummy, gLeft);

      // rotate flux
      rotU_.rotateBack( gLeft, normal );
      rotB_.rotateBack( gLeft, normal );

      return ldt;
    }

    const Model& model() const { return model_; }
  protected:
    const Model& model_;
    const typename MhdSolverType::Eosmode::meos_t eos;
    mutable MhdSolverType numFlux_;
    Dune::Fem::FieldRotator<DomainType, RangeType> rotU_, rotB_;
    //mutable MhdSolverType::Vec9 ulmhd_, urmhd_, retmhd_;
  };



  //////////////////////////////////////////////////////////
  //
  //  Flux Implementations
  //
  //////////////////////////////////////////////////////////

  template < class Model >
  class DWNumFlux : public MHDNumFluxBase< Model, Mhd::DW >
  {
    typedef MHDNumFluxBase< Model, Mhd::DW > BaseType ;
  public:
    DWNumFlux( const Model& model )
      : BaseType( model )
    {}
    static std::string name () { return "DW (Mhd)"; }
  };

  template < class Model >
  class HLLEMNumFlux : public MHDNumFluxBase< Model, Mhd::HLLEM >
  {
    typedef MHDNumFluxBase< Model, Mhd::HLLEM > BaseType ;
  public:
    HLLEMNumFlux( const Model& model )
      : BaseType( model )
    {}
    static std::string name () { return "HLLEM (Mhd)"; }
  };


  template< class ModelImp >
  class DGAdvectionFlux< ModelImp, AdvectionFlux::Enum::mhd_general >
    : public DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >
  {
    typedef DGAdvectionFluxBase< ModelImp, AdvectionFluxParameters >   BaseType;

    static const int dimRange  = ModelImp::dimRange;
    static const int dimDomain = ModelImp::dimDomain;
    typedef typename ModelImp::DomainType         DomainType;
    typedef typename ModelImp::RangeType          RangeType;
    typedef typename ModelImp::JacobianRangeType  JacobianRangeType;
    typedef typename ModelImp::FluxRangeType      FluxRangeType;
    typedef typename ModelImp::FaceDomainType     FaceDomainType;

  public:
    typedef AdvectionFlux::Enum                   IdEnum;
    typedef typename BaseType::ModelType          ModelType;
    typedef typename BaseType::ParameterType      ParameterType;

    /**
     * \copydoc DGAdvectionFluxBase::DGAdvectionFluxBase()
     */
    template< class ... Args>
    DGAdvectionFlux(  Args&&... args )
      : BaseType( std::forward<Args>(args)... ),
        method_( this->parameter().getMethod() ),
        flux_dw_(this->model() ),
        flux_hllem_( this->model() )
    {}

    /**
     * \copydoc DGAdvectionFluxBase::name()
     */
    static std::string name () { return "AdvectionFlux - MHD (via parameter file)"; }

    /**
     * \copydoc DGAdvectionFluxBase::numericalFlux()
     */
    template< class LocalEvaluation >
    inline double
    numericalFlux( const LocalEvaluation& left,
                   const LocalEvaluation& right,
                   const RangeType& uLeft,
                   const RangeType& uRight,
                   const JacobianRangeType& jacLeft,
                   const JacobianRangeType& jacRight,
                   RangeType& gLeft,
                   RangeType& gRight) const
    {
      Dune::FieldVector< double, dimDomain > normal = left.intersection().integrationOuterNormal( left.localPosition() );
      const double len = normal.two_norm();
      normal *= 1./len;
      double ldt;
      if ( IdEnum::mhd_dw == method_ )
      {
        ldt = flux_dw_.numericalFlux( normal, uLeft, uRight, gLeft );
      }
      else if ( IdEnum::mhd_hllem == method_ )
      {
        ldt = flux_hllem_.numericalFlux( normal, uLeft, uRight, gLeft );
      }
      else
      {
        std::cerr << "Error: Advection flux not chosen via parameter file" << std::endl;
        assert( false );
        std::abort();
      }
      // conservation
      gLeft *= len;
      gRight = gLeft;
      return ldt*len;
    }

  private:
    const IdEnum             method_;
    DWNumFlux< ModelImp >    flux_dw_;
    HLLEMNumFlux< ModelImp > flux_hllem_;
  };

}}
#endif // DUNE_MHDFLUXES_HH
