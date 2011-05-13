#ifndef DUNE_FEM_DG_ROTATOR_HH
#define DUNE_FEM_DG_ROTATOR_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes

namespace EulerFluxes {

  template <class FunctionSpaceT>
  class FieldRotator {
  public:
    //- Global typedefs
    typedef typename FunctionSpaceT::DomainType NormalType;
    typedef typename FunctionSpaceT::RangeType ValueType;

    template <int N> 
    struct FRInt2Type{
      enum { value = N };
    };

    //! Constructor
    //! The vector components to be rotated must be consecutive in the
    //! vector of unknows.
    //! \param startIdx Specifies first component of vector component
    FieldRotator(const int startIdx) : 
      idx_(startIdx) {
      assert(startIdx < ValueType::size-NormalType::size+1);
    }

    //! Rotate data from basic coordinate system into normal coordinate system
    void rotateForth(ValueType& res, 
                     const NormalType& n) const {
      rotateForth(res, n, FRInt2Type<NormalType::size>());
    }

    //! Rotate data from normal coordinate system into basic coordinate system
    void rotateBack(ValueType& res, 
                    const NormalType& n) const {
      rotateBack(res, n, FRInt2Type<NormalType::size>());
    }

  private:
    // Local methods
    void rotateForth(ValueType& res, 
                     const NormalType& n,
                     FRInt2Type<1>) const;
    void rotateForth(ValueType& res, 
                     const NormalType& n,
                     FRInt2Type<2>) const;
    void rotateForth(ValueType& res, 
                     const NormalType& n,
                     FRInt2Type<3>) const;
    void rotateBack(ValueType& res, 
                    const NormalType& n,
                    FRInt2Type<1>) const;
    void rotateBack(ValueType& res, 
                    const NormalType& n,
                    FRInt2Type<2>) const;
    void rotateBack(ValueType& res, 
                    const NormalType& n,
                    FRInt2Type<3>) const;

    
    const int idx_;
    const static double eps_;
  };

  template <class FunctionSpaceT>
  const double FieldRotator<FunctionSpaceT>::eps_ = 1.0e-14;

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(ValueType& res, 
              const NormalType& n,
              FRInt2Type<1>) const {
    // res = arg ;
    res[idx_] = res[idx_] * n[0];
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(ValueType& res, 
              const NormalType& n,
              FRInt2Type<2>) const {
    // res = arg;
    const double a[2] = { res[idx_], res[idx_+1] };
    res[ idx_   ] =  n[0]*a[0] + n[1]*a[1];
    res[ idx_+1 ] = -n[1]*a[0] + n[0]*a[1];
  }
  
  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateForth(ValueType& res, 
              const NormalType& n,
              FRInt2Type<3>) const {
    // res = arg;

    const double a[3] = { res[idx_], res[idx_+1], res[idx_+2] };
 
    const double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8) 
    {
      double d_1 = 1.0/d;
      res[idx_]   =   n[0] * a[0]
                    + n[1] * a[1]
                    + n[2] * a[2];
      res[idx_+1] = - n[1] * d_1 * a[0]
                    + n[0] * d_1 * a[1];
      res[idx_+2] = - n[0] * n[2] * d_1 * a[0]
                    - n[1] * n[2] * d_1 * a[1]
                    + d                 * a[2];
    } 
    else 
    {
      res[idx_]   =   n[2] * a[2];
      res[idx_+1] =          a[1];
      res[idx_+2] = - n[2] * a[0];
    }
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(ValueType& res, 
	           const NormalType& n,
	           FRInt2Type<1>) const 
  {
    // res = arg;
    res[idx_] = res[idx_] * n[0]; 
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(ValueType& res, 
             const NormalType& n,
             FRInt2Type<2>) const 
  {
    // res = arg;

    const double a[2] = { res[idx_], res[idx_+1] };
    res[idx_  ] = n[0]*a[0] - n[1]*a[1];
    res[idx_+1] = n[1]*a[0] + n[0]*a[1];
  }

  template <class FunctionSpaceT>
  inline void FieldRotator<FunctionSpaceT>::
  rotateBack(ValueType& res, 
             const NormalType& n,
             FRInt2Type<3>) const 
  {
    // res = arg;
    double a[3]={res[idx_],res[idx_+1],res[idx_+2]};

    double d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8) 
    {
      double d_1 = 1.0/d;
      res[idx_]   =   n[0]              * a[0]
                    - n[1] * d_1        * a[1]
                    - n[0] * n[2] * d_1 * a[2];
      res[idx_+1] =   n[1]              * a[0]
                    + n[0] * d_1        * a[1]
                    - n[1] * n[2] * d_1 * a[2];
      res[idx_+2] =   n[2]              * a[0]
                    + d                 * a[2];
    } 
    else 
    {
      res[idx_]   = - n[2] * a[2];
      res[idx_+1] =          a[1];
      res[idx_+2] =   n[2] * a[0];
    }
    
  }
} // end namespace EulerFluxes

// depreacted namespace 
namespace Adi {
  using EulerFluxes :: FieldRotator ;
}
#endif
