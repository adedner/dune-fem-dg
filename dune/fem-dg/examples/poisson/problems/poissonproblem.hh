#ifndef DUNE_FEM_DG_POISSONPROBLEM_HH
#define DUNE_FEM_DG_POISSONPROBLEM_HH

#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/fem-dg/models/defaultprobleminterfaces.hh>
#include "benchmarkproblems.hh"

namespace Dune
{
namespace Fem
{
namespace Poisson
{

  /**
   * \brief Interface class for the Poisson problem.
   *
   * \ingroup PoissonProblems.
   */
  template <class GridType, int dimRange>
  struct Problem : public ProblemInterface<
       Dune :: Fem :: FunctionSpace< double, double, GridType::dimension, dimRange> >
  {
  public:
    typedef ProblemInterface<
         Dune :: Fem::FunctionSpace< double, double,
         GridType::dimension, dimRange >
            >                    BaseType;

    enum{ dimDomain = BaseType :: dimDomain };
    enum{ dim = dimDomain };
    typedef typename BaseType :: DomainType                DomainType;
    typedef typename BaseType :: RangeType                 RangeType;
    typedef typename BaseType :: JacobianRangeType         JacobianRangeType;
    typedef typename BaseType :: DiffusionMatrixType   DiffusionMatrixType;
    typedef typename GridType :: ctype   FieldType ;

    typedef DataFunctionIF< dimDomain, FieldType, FieldType > DataFunctionType;

    // not copy or assignment
    Problem( const Problem& ) = delete;
    Problem& operator= ( const Problem& ) = delete;

    /**
     * \brief define problem parameters
     */
    Problem(const int problemNumber) :
      BaseType (),
      data_(0),
      myName_("")
    {
      FieldType shift = Dune :: Fem::Parameter :: getValue< double > ("globalshift", 0);
      FieldType factor = Dune :: Fem::Parameter :: getValue< double > ("factor", 1);
      FieldType advection = Dune :: Fem::Parameter :: getValue< double > ("advection", 1);
      if( problemNumber == 0 )
      {
        data_ = new BenchMark_1<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 1 )
      {
        data_ = new BenchMark_1_2<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 2 )
      {
        data_ = new BenchMark_2<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 3 )
      {
        data_ = new BenchMark_3<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 4 )
      {
        data_ = new BenchMark_4<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 5 )
      {
        data_ = new BenchMark_5<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 6 )
      {
        data_ = new BenchMark_6<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 7 )
      {
        data_ = new BenchMark_7<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 8 )
      {
        data_ = new CurvedRidges< dim, FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 9 )
      {
        data_ = new BenchMark_9<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 10 )
      {
        data_ = new SinSin<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 11 )
      {
        data_ = new Excercise2_3< dim, FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 12 )
      {
        data_ = new CastilloProblem<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 13 )
      {
        data_ = new ReentrantCorner<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 14 )
      {
        data_ = new RiviereProblem<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 15 )
      {
        data_ = new HeatProblem<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 16 )
      {
        data_ = new AlbertaProblem<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 17 )
      {
        data_ = new SinSin<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 18 )
      {
        data_ = new CosCos<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 19 )
      {
        data_ = new SinSinSin<dim,FieldType,FieldType> (shift,factor,advection);
      }
      if( problemNumber == 20 )
      {
        data_ = new Hump<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 21 )
      {
        data_ = new BenchMark3d_1<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 23 )
      {
        data_ = new BenchMark3d_3<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 24 )
      {
        data_ = new BenchMark3d_4<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 25 )
      {
        data_ = new BenchMark3d_5<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 31 )
      {
        data_ = new BoundaryLayerProblem<dim,FieldType,FieldType> (shift,factor);
      }
      if( problemNumber == 32 )
      {
        data_ = new FicheraCorner<dim,FieldType,FieldType> (shift,factor);
      }

      if( data_ == 0 )
      {
        std::cerr << "ERROR: wrong problem number " << std::endl;
        abort();
      }

      myName_ = "Poisson Eqn.";
    }

    const DataFunctionType& data() const { assert( data_ ); return *data_; }

    //! this problem has no source term
    bool hasStiffSource() const { return true; }
    bool hasNonStiffSource() const { return false; }

    void f(const DomainType& arg,
           RangeType& res) const
    {
      res = data().rhs( &arg[ 0 ] );
    }

    double stiffSource(const DomainType& arg,
                       const double t,
                       const RangeType& u,
                       RangeType& res) const
    {
      res = data().rhs( &arg[ 0 ] );
      return 0;
    }

    double nonStiffSource(const DomainType& arg,
                  const double t,
                  const RangeType& u,
                  RangeType& res) const
    {
      abort();
      return 0.0;
    }

    //! return start time
    double startTime() const { return 0; }

    double constantAdvection() const { return data().advection(); }

    bool constantK() const { return data().constantLocalK();  }

    void K(const DomainType& x, DiffusionMatrixType& m) const
    {
      double k[ dimDomain ][ dimDomain ];
      data().K( &x[0], k );

      for(int i=0; i<dimDomain; ++i)
        for(int j=0; j<dimDomain; ++j)
           m[i][j] = k[i][j];
    }

    /**
     * \brief getter for the velocity
     */
    void velocity(const DomainType& x, DomainType& v) const
    {
      double vv[ dimDomain ];
      data().velocity( &x[0], vv );
      for (int i=0;i<dimDomain; ++i)
        v[i] = vv[i];
    }

    void u(const DomainType& arg, RangeType& res) const
    {
      evaluate(arg, res );
    }

    /**
     * \brief evaluates \f$ u_0(x) \f$
     */
    void evaluate(const DomainType& arg, RangeType& res) const
    {
      evaluate(arg, 0, res);
    }

    /**
     * \brief evaluate exact solution
     */
    void evaluate(const DomainType& arg, const double t, RangeType& res) const
    {
      res = data().exact( &arg[ 0 ] );
    }

    void gradient(const DomainType& x,
                  JacobianRangeType& grad) const
    {
      for (int i=0;i<RangeType::dimension;++i)
        data().gradExact( &x[ 0 ], &grad[ i ][ 0 ] );
    }

    /**
     * \brief latex output for EocOutput
     */
    std::string description() const
    {
      std::ostringstream ofs;

      ofs << "Problem: " << myName_ ;
      return ofs.str();
    }

  protected:
    DataFunctionType* data_;
    std::string myName_;
  };

}
}
}
#endif  /*DUNE_PROBLEM_HH__*/

