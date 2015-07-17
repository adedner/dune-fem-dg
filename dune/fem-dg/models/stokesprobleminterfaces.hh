#ifndef DUNE_FEM_DG_STOKES_PROBLEMINTERFACE_HH
#define DUNE_FEM_DG_STOKES_PROBLEMINTERFACE_HH

#include <dune/common/version.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/fem/misc/gridsolution.hh>

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>

namespace Dune {

template <class FunctionSpaceImp, class PressureSpaceImp>
class StokesProblemInterface : public ProblemInterface < FunctionSpaceImp >
{
public:
  typedef FunctionSpaceImp                                           FunctionSpaceType;
  typedef PressureSpaceImp                                           PressureSpaceType;

  typedef StokesProblemInterface< FunctionSpaceType, PressureSpaceType>  ThisType;

  enum { dimDomain = FunctionSpaceType :: dimDomain };

  typedef typename FunctionSpaceType :: DomainType                   DomainType;
  typedef typename FunctionSpaceType :: RangeType                    RangeType;
  typedef typename FunctionSpaceType :: JacobianRangeType            JacobianRangeType;
  typedef typename FunctionSpaceType :: DomainFieldType              DomainFieldType;
  typedef typename FunctionSpaceType :: RangeFieldType               RangeFieldType;

  typedef typename PressureSpaceType :: RangeType                    PressureRangeType;
  typedef typename PressureSpaceType :: JacobianRangeType            PressureJacobianRangeType;

  typedef FieldMatrix< RangeFieldType, dimDomain, dimDomain >        DiffusionMatrixType;

public:

  //! destructor
  virtual ~StokesProblemInterface() {}

  //! the exact pressure
  virtual void p(const DomainType& x, PressureRangeType& ret) const {}//= 0;

  //! the pressure boundary data


  virtual void gp(const DomainType& x, PressureRangeType& ret) const
  {
    ret = 0;
  }


  //! the gradient of the exact solution
  //virtual void gradient(const DomainType& x,
  //                      JacobianRangeType& grad) const = 0;

protected:

  StokesProblemInterface() : exactPressure_(*this){}


  //! the exact pressure to the problem for EOC calculation
  class ExactPressure
    : public Fem::Function< PressureSpaceType, ExactPressure >
  {
  private:
    typedef Fem::Function< PressureSpaceType, ExactPressure >      BaseType;

    typedef StokesProblemInterface< FunctionSpaceType,PressureSpaceType>   DataType;
  protected:
    PressureSpaceType  functionSpace_;
    const DataType &data_;

  public:
    inline ExactPressure (  const ThisType& data )
      : BaseType( ),
        functionSpace_(),
        data_(data)
    {
    }

    inline void evaluate ( const DomainType &x, PressureRangeType &ret ) const
    {
       data_.p( x, ret );
    }

    inline void jacobian ( const DomainType &x,PressureJacobianRangeType &ret ) const
    {
    //   data_.gradient( x, ret );
    }

    inline void evaluate (const DomainType &x,
                          const double time, PressureRangeType &phi ) const
    {
			evaluate( x, phi );
    }
  }; // end class ExactSolutionX

public:
  //! type of function converter for exact solution and gradient
  typedef ExactPressure ExactPressureType;
//protected:
  ExactPressureType exactPressure_;
public:
  const ExactPressureType& exactPressure() const { return exactPressure_; }

};
}
#endif  /*DUNE_PROBLEMINTERFACE_HH*/
