#ifndef DUNE_FEM_DG_MHDFLUXES_HH
#define DUNE_FEM_DG_MHDFLUXES_HH
#warning "Using Mhd NumFluxes"

#include <cmath>
#include <dune/fem/misc/field.hh>

#ifdef COUNT_FLOPS
#include <double.h>
#endif

#include <dune/fem/storage/vector.hh>
#include <dune/fem/io/parameter.hh>

#include "mhd_eqns.hh"
#include "rotator.hh"

// Dai-Woodward 
template <class Model>
class DWNumFlux;

// HLLEM
template <class Model>
class HLLEMNumFlux;

// ************************************************
template <int dimDomain>
class ConsVec : public FieldVector< double, dimDomain+2> 
{
public:
  explicit ConsVec (const double& t) : FieldVector<double,dimDomain+2>(t) {}
  ConsVec () : FieldVector<double,dimDomain+2>(0) {}
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
  enum { dimDomain = Model::dimDomain };
  enum { dimRange = Model::dimRange };
  typedef typename Model::Traits                      Traits;
  typedef typename Traits::GridType                   GridType;
  typedef typename GridType::ctype                    ctype;
  typedef typename Traits::EntityType                 EntityType;
  typedef typename Traits::EntityPointerType          EntityPointerType;

  typedef typename Traits::DomainType                 DomainType;
  typedef typename Traits::FaceDomainType             FaceDomainType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::FluxRangeType              FluxRangeType;

  typedef Mhd::MhdSolver MhdSolverType;

  typedef double value_t[ 9 ];

protected:  
  MHDNumFluxBase(const Model& mod) 
   : model_(mod),
     eos( MhdSolverType::Eosmode::me_ideal ),
     numFlux_(eos,  mod.c_p() * mod.c_v_inv(), 1.0 ),
     rot_(1) 
  {
    if( fluxtype == Mhd :: HLLEM ) 
    {
      if( Parameter :: verbose () )
        std::cout << "Choosing HLLEM Flux " << std::endl;

      numFlux_.init(Mhd :: MhdSolver :: mf_rghllem ); 
    }
  }
  
public:
  // Return value: maximum wavespeed*length of integrationOuterNormal
  // gLeft,gRight are fluxed * length of integrationOuterNormal
  template< class Intersection, class QuadratureImp >
  inline double
  numericalFlux( const Intersection& intersection,
                 const EntityType& inside,
                 const EntityType& outside,
                 const double time,
                 const QuadratureImp& faceQuadInner,
                 const QuadratureImp& faceQuadOuter,
                 const int quadPoint,
                 const RangeType& uLeft,
                 const RangeType& uRight,
                 RangeType& gLeft,
                 RangeType& gRight) const;

  const Model& model() const { return model_; }
protected:
  const Model& model_;
  const typename MhdSolverType::Eosmode::meos_t eos;
  mutable MhdSolverType numFlux_;
  Adi::FieldRotator<Model> rot_;
  //mutable MhdSolverType::Vec9 ulmhd_, urmhd_, retmhd_;
};


template < class Model, Mhd :: MhdFluxType fluxtype >
template< class Intersection, class QuadratureImp >
inline double MHDNumFluxBase< Model, fluxtype > :: 
numericalFlux( const Intersection& intersection,
               const EntityType& inside,
               const EntityType& outside,
               const double time,
               const QuadratureImp& faceQuadInner,
               const QuadratureImp& faceQuadOuter,
               const int quadPoint,
               const RangeType& uLeft,
               const RangeType& uRight,
               RangeType& gLeft,
               RangeType& gRight) const
{
  DomainType normal = intersection.integrationOuterNormal( faceQuadInner.localPoint( quadPoint ) );
  // double len = normal.two_norm();
  const double len = normal.two_norm();
  normal *= 1./len;

  RangeType ul(uLeft);
  RangeType ur(uRight);

  rot_.rotateForth(ul, normal);
  rot_.rotateForth(ur, normal);

  enum { e = dimDomain + 1 };

  value_t res;
  const double dummy[ 3 ] = { 0, 0, 0 };

  value_t entity = { 0,0,0,0,0,0,0,0 };
  value_t neigh  = { 0,0,0,0,0,0,0,0 };
  for(int i=0; i<e; ++i) 
  {
    entity[ i ] = ul[ i ];
    neigh [ i ] = ur[ i ];
  }

  entity[ 7 ] = ul[ e ];
  neigh [ 7 ] = ur[ e ];

  const double ldt = numFlux_(entity, neigh, dummy, res);

  // copy first components 
  for(int i=0; i<e; ++i) 
    gLeft[ i ] = res[ i ];

  // copy energy 
  gLeft[ e ] = res[ 7 ];

  // rotate flux 
  rot_.rotateBack( gLeft, normal );

  // conservation
  gLeft *= len;
  gRight = gLeft;
  return ldt*len;
}

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

#endif // DUNE_MHDFLUXES_HH
