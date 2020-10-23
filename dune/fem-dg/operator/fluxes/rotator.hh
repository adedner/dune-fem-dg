#ifndef DUNE_FEM_DG_ROTATOR_HH
#define DUNE_FEM_DG_ROTATOR_HH

//- system includes
#include <cmath>


namespace Dune
{
namespace Fem
{

  template <class Domain, class Range>
  class FieldRotator
  {
    template <int N>
    struct RotInt2Type{
      static const int value = N;
    };

  public:
    //- Global typedefs
    typedef Domain NormalType;
    typedef Range  ValueType;
    typedef typename ValueType :: field_type    FieldType;

    //! Constructor
    //! The vector components to be rotated must be consecutive in the
    //! vector of unknows.
    //! \param startIdx Specifies first component of vector component
    FieldRotator(const int startIdx) :
      idx_(startIdx)
    {
      assert(startIdx < (int(ValueType::dimension) - int(NormalType::dimension)+1) );
    }

    //! Rotate data from basic coordinate system into normal coordinate system
    void rotateForth(ValueType& res,
                     const NormalType& n) const
    {
      rotateForth(res, n, RotInt2Type<NormalType::dimension>());
    }

    //! Rotate data from normal coordinate system into basic coordinate system
    void rotateBack(ValueType& res,
                    const NormalType& n) const
    {
      rotateBack(res, n, RotInt2Type<NormalType::dimension>());
    }

  private:
    // Local methods
    void rotateForth(ValueType& res,
                     const NormalType& n,
                     RotInt2Type<1>) const;
    void rotateForth(ValueType& res,
                     const NormalType& n,
                     RotInt2Type<2>) const;
    void rotateForth(ValueType& res,
                     const NormalType& n,
                     RotInt2Type<3>) const;
    void rotateBack(ValueType& res,
                    const NormalType& n,
                    RotInt2Type<1>) const;
    void rotateBack(ValueType& res,
                    const NormalType& n,
                    RotInt2Type<2>) const;
    void rotateBack(ValueType& res,
                    const NormalType& n,
                    RotInt2Type<3>) const;


    const int idx_;
    const static double eps_;
  };

  template <class Domain, class Range>
  const double FieldRotator<Domain, Range>::eps_ = 1.0e-14;

  template <class Domain, class Range>
  inline void FieldRotator<Domain, Range>::
  rotateForth(ValueType& res,
              const NormalType& n,
              RotInt2Type<1>) const {
    // res = arg ;
    res[idx_] = res[idx_] * n[0];
  }

  template <class Domain, class Range>
  inline void FieldRotator<Domain, Range>::
  rotateForth(ValueType& res,
              const NormalType& n,
              RotInt2Type<2>) const {
    // res = arg;
    const FieldType a[2] = { res[idx_], res[idx_+1] };
    res[ idx_   ] =  n[0]*a[0] + n[1]*a[1];
    res[ idx_+1 ] = -n[1]*a[0] + n[0]*a[1];
  }

  template <class Domain, class Range>
  inline void FieldRotator<Domain, Range>::
  rotateForth(ValueType& res,
              const NormalType& n,
              RotInt2Type<3>) const
  {
    const FieldType a[3] = { res[idx_], res[idx_+1], res[idx_+2] };

    const FieldType d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8)
    {
      const FieldType d_1 = 1.0/d;
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

  template <class Domain, class Range>
  inline void FieldRotator<Domain, Range>::
  rotateBack(ValueType& res,
	           const NormalType& n,
	           RotInt2Type<1>) const
  {
    res[idx_] = res[idx_] * n[0];
  }

  template <class Domain, class Range>
  inline void FieldRotator<Domain, Range>::
  rotateBack(ValueType& res,
             const NormalType& n,
             RotInt2Type<2>) const
  {
    const FieldType a[2] = { res[idx_], res[idx_ + 1] };
    res[idx_  ] = n[0]*a[0] - n[1]*a[1];
    res[idx_+1] = n[1]*a[0] + n[0]*a[1];
  }

  template <class Domain, class Range>
  inline void FieldRotator<Domain, Range>::
  rotateBack(ValueType& res,
             const NormalType& n,
             RotInt2Type<3>) const
  {
    // res = arg;
    const FieldType a[3]={res[idx_],res[idx_+1],res[idx_+2]};

    const FieldType d = std::sqrt(n[0]*n[0]+n[1]*n[1]);

    if (d > 1.0e-8)
    {
      const FieldType d_1 = 1.0/d;
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

} // end namespace
} // end namespace

#endif
