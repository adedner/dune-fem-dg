#ifndef DUNE_FEM_DG_POISSON_MODEL_HH
#define DUNE_FEM_DG_POISSON_MODEL_HH

#include <dune/fem/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem-dg/models/defaultmodel.hh>

namespace Dune
{
  template <class ModelType>
  class UpwindFlux;
}

/**********************************************
 * Analytical model                           *
 *********************************************/
/**
 * @brief Traits class for PoissonModel
 */
template <class GridPart,
          class ProblemImp>
class PoissonModelTraits
{
public:
  typedef ProblemImp  ProblemType;
  typedef GridPart                                                   GridPartType;
  typedef typename GridPartType :: GridType                          GridType;
  static const int dimDomain = GridType::dimensionworld;
  static const int dimRange = ProblemType :: dimRange;
  static const int dimGradRange = dimRange * dimDomain;

  typedef double RangeFieldType;
  typedef double DomainFieldType;
  // Definition of domain and range types
  typedef Dune::FieldVector< double, dimDomain >                     DomainType;
  typedef Dune::FieldVector< double, dimDomain-1 >                   FaceDomainType;
  typedef Dune::FieldVector< double, dimRange >                      RangeType;
  typedef Dune::FieldVector< double, dimGradRange >                  GradientType;
  // ATTENTION: These are matrices (c.f. PoissonModel)
  typedef Dune::FieldMatrix< double, dimRange, dimDomain >           FluxRangeType;
  typedef Dune::FieldMatrix< double, dimRange, dimDomain >           JacobianRangeType;
  typedef Dune::FieldMatrix< double, dimGradRange, dimDomain >       DiffusionRangeType;
  typedef typename GridType :: template Codim< 0 > :: Entity         EntityType;
  typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
  typedef typename IntersectionIterator :: Intersection              IntersectionType;

};

/**
 * @brief describes the analytical model
 *
 * This is an description class for the problem
 * \f{eqnarray*}{ V + \nabla a(U)      & = & 0 \\
 * \partial_t U + \nabla (F(U)+A(U,V)) & = & 0 \\
 *                          U          & = & g_D \\
 *                   \nabla U \cdot n  & = & g_N \f}
 *
 * where each class methods describes an analytical function.
 * <ul>
 * <li> \f$F\f$:   advection() </li>
 * <li> \f$a\f$:   diffusion1() </li>
 * <li> \f$A\f$:   diffusion2() </li>
 * <li> \f$g_D\f$  boundaryValue() </li>
 * <li> \f$g_N\f$  boundaryFlux1(), boundaryFlux2() </li>
 * </ul>
 *
 * \attention \f$F(U)\f$ and \f$A(U,V)\f$ are matrix valued, and therefore the
 * divergence is defined as
 *
 * \f[ \Delta M = \nabla \cdot (\nabla \cdot (M_{i\cdot})^t)_{i\in 1\dots n} \f]
 *
 * for a matrix \f$M\in \mathbf{M}^{n\times m}\f$.
 *
 * @param GridPart GridPart for extraction of dimension
 * @param ProblemType Class describing the initial(t=0) and exact solution
 */


////////////////////////////////////////////////////////
//
//  Analytical model for the Heat Equation
//      dx(u) + div(uV) - epsilon*lap(u)) = 0
//  where V is constant vector
//
////////////////////////////////////////////////////////
template <class GridPartType, class ProblemImp>
class PoissonModel : public DefaultModel< PoissonModelTraits< GridPartType, ProblemImp > >
{
public:
  typedef ProblemImp  ProblemType ;

  typedef typename GridPartType :: GridType                          GridType;
  static const int dimDomain = GridType::dimensionworld;
  static const int dimRange = ProblemType :: dimRange;
  typedef PoissonModelTraits< GridPartType, ProblemType >             Traits;
  typedef typename Traits :: DomainFieldType                         DomainFieldType;
  typedef typename Traits :: RangeFieldType                          RangeFieldType;
  typedef typename Traits :: DomainType                              DomainType;
  typedef typename Traits :: RangeType                               RangeType;
  typedef typename Traits :: GradientType                            GradientType;
  typedef typename Traits :: FluxRangeType                           FluxRangeType;
  typedef typename Traits :: DiffusionRangeType                      DiffusionRangeType;
  typedef typename Traits :: FaceDomainType                          FaceDomainType;
  typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

  typedef typename ProblemType :: DiffusionMatrixType   DiffusionMatrixType ;

  typedef typename Traits :: EntityType                       EntityType;
  typedef typename Traits :: IntersectionType                 IntersectionType;

public:
  static const int ConstantVelocity = false;
  /**
   * @brief Constructor
   *
   * initializes model parameter
   *
   * @param problem Class describing the initial(t=0) and exact solution
   */
  PoissonModel(const ProblemType& problem) : problem_(problem)
  {
  }

  inline bool hasFlux() const { return true ; }

  inline bool hasStiffSource() const { return true ; }
  inline bool hasNonStiffSource() const { return false ; }

  template <class LocalEvaluation>
  inline double stiffSource( const LocalEvaluation& local,
                        const RangeType& u,
                        const JacobianRangeType& du,
                        RangeType & s) const
  {
    DomainType xgl = local.entity().geometry().global( local.point() );
    problem_.f( xgl, s );
    return 0.0;
  }

  template <class LocalEvaluation>
  inline double nonStiffSource( const LocalEvaluation& local,
                        const RangeType& u,
                        const JacobianRangeType& du,
                        RangeType & s) const
  {
    DomainType xgl = local.entity().geometry().global( local.point() );
    problem_.f( xgl, s );
    return 0.0;
  }


  /**
   * @brief advection term \f$F\f$
   *
   * @param en entity on which to evaluate the advection term
   * @param time current time of TimeProvider
   * @param x coordinate local to entity
   * @param u \f$U\f$
   * @param f \f$f(U)\f$
   */
  template <class LocalEvaluation>
  inline  void advection(const LocalEvaluation& local,
                         const RangeType& u,
                         const JacobianRangeType& jacu,
                         FluxRangeType & f) const
  {
    // evaluate velocity V
    DomainType v;
    velocity( local.entity(), local.time(), local.point(), v );

    //f = uV;
    for( int r=0; r<dimRange; ++r )
      for( int d=0; d<dimDomain; ++d )
        f[r][d] = v[ d ] * u[ r ];
  }

  bool hasDirichletBoundary () const
  {
    return true ;
  }

  bool isDirichletPoint( const DomainType& global ) const
  {
    return true ;
  }

protected:

  /**
   * @brief velocity calculation, is called by advection()
   */
  inline  void velocity(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        DomainType& v) const
  {
    problem_.velocity(en.geometry().global(x), v);
  }

public:

  /**
   * @brief diffusion term \f$a\f$
   */
  template <class LocalEvaluation>
  inline void jacobian(const LocalEvaluation& local,
                        const RangeType& u,
                        DiffusionRangeType& a) const
  {
    a = 0;

    assert( a.rows == dimRange * dimDomain );
    assert( a.cols == dimDomain );

    for (int r=0;r<dimRange;r++)
      for (int d=0;d<dimDomain;d++)
        a[dimDomain*r+d][d] = u[r];
  }

  template <class T>
  T SQR( const T& a ) const
  {
    return (a * a);
  }

  template <class LocalEvaluation>
  inline void eigenValues(const LocalEvaluation& local,
                          const RangeType& u,
                          RangeType& maxValue) const
  {
    DiffusionMatrixType K ;
    DomainType xgl = local.entity().geometry().global( local.point() );
    problem_.K( xgl, K );

    DomainType values ;
    // calculate eigenvalues
    Dune::FMatrixHelp :: eigenValues( K, values );

    maxValue = SQR(values[ dimDomain -1 ]) /values[0];
    return ;

    // take max eigenvalue
    maxValue = values.infinity_norm();
  }

  inline double lambdaK( const DiffusionMatrixType& K ) const
  {
    DomainType values ;
    // calculate eigenvalues
    Dune::FMatrixHelp :: eigenValues( K, values );

    // value[ 0 ] is smallest ev
    return SQR(values[ dimDomain -1 ]) / values[ 0 ];
  }

  template <class LocalEvaluation>
  inline double penaltyFactor( const LocalEvaluation& left,
                               const LocalEvaluation& right,
                               const RangeType& uLeft,
                               const RangeType& uRight ) const
  {
    DiffusionMatrixType K ;
    double betaK, betaInside, betaOutside ;
    if( problem_.constantK() )
    {
      DiffusionMatrixType Kinside ;
      DiffusionMatrixType Koutside;

      const DomainType xglIn = left.entity().geometry().center();
      problem_.K( xglIn , Kinside );
      const DomainType xglOut = right.entity().geometry().center();
      problem_.K( xglOut , Koutside );

      K = Kinside ;
      K += Koutside ;
      K *= 0.5 ;

      betaInside  = lambdaK( Kinside );
      betaOutside = lambdaK( Koutside );

      betaK = lambdaK( K );
    }
    else
    {
      const DomainType xgl = left.entity().geometry().global( left.point() );
      problem_.K( xgl , K );

      betaK = lambdaK( K );
      betaInside = betaOutside = betaK;
    }

    const double jump = std::tanh( std::abs( betaInside - betaOutside ) );

    // only for small values of betS apply betS in smooth cases
    const double betaN = std :: min( betaK , 1.0 );

    // betS becomes 1 if the eigen values of both matrices are the same
    betaK = betaK * jump + (1.0 - jump) * betaN;

    return betaK ;
  }

  inline double penaltyBoundary( const EntityType& inside,
                                 const double time,
                                 const DomainType& xInside,
                                 const RangeType& uLeft ) const
  {
    return penaltyFactor( inside, inside, time, xInside, uLeft, uLeft );
  }

  /**
   * @brief diffusion term \f$A\f$
   */
  template <class LocalEvaluation>
  inline void diffusion(const LocalEvaluation& local,
                        const RangeType& u,
                        const JacobianRangeType& jac,
                        FluxRangeType& A) const
  {
    // for constant K evalute at center (see Problem 4)
    const DomainType xgl = ( problem_.constantK() ) ?
      local.entity().geometry().center () : local.entity().geometry().global(local.point())  ;

    DiffusionMatrixType K ;

    // fill diffusion matrix
    problem_.K( xgl, K );

    // apply diffusion
    for( int r =0; r<dimRange; ++r )
      K.mv( jac[ r ] , A[ r ] );
  }

  template <class LocalEvaluation>
  inline double diffusionTimeStep( const LocalEvaluation& local,
                                   const double circumEstimate,
                                   const RangeType& u ) const
  {
    return 0;
  }

public:
  /**
   * @brief checks for existence of dirichlet boundary values
   */
  template <class LocalEvaluation>
  inline bool hasBoundaryValue(const LocalEvaluation& local ) const
  {
    return true;
  }

  /**
   * @brief dirichlet boundary values
   */
  template <class LocalEvaluation>
  inline  void boundaryValue(const LocalEvaluation& local,
                             const RangeType& uLeft,
                             RangeType& uRight) const
  {
    if ( local.intersection().boundaryId() == 99) // Dirichlet zero boundary conditions
    {
      uRight = 0;
    }
    else
    {
      DomainType xgl = local.entity().geometry().global( local.point() );
      problem_.g(xgl, uRight);
    }
  }

  /**
   * @brief diffusion boundary flux
   */
  template <class LocalEvaluation>
  inline double diffusionBoundaryFlux( const LocalEvaluation& local,
                                       const RangeType& uLeft,
                                       const JacobianRangeType& gradLeft,
                                       RangeType& gLeft ) const
  {
  }

  const ProblemType& problem () const { return problem_; }

 protected:
  const ProblemType& problem_;
  friend class Dune::UpwindFlux<PoissonModel>;
};

#endif
