#ifndef DUNE_FEM_DG_DEFAULTMODEL_HH
#define DUNE_FEM_DG_DEFAULTMODEL_HH

#include<limits>
#include<tuple>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem-dg/misc/integral_constant.hh>
#include <dune/fem-dg/misc/error/eocerrorlist.hh>
#include <dune/fem-dg/misc/error/l2eocerror.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{
namespace Fem
{


  template< int i, class Arg >
  class PassIdGenerator;

  template< int i, class... bools >
  class PassIdGenerator< i, std::tuple< bools... > >
  {
    typedef std::tuple< bools... > T;

    template<int ii>
    static constexpr typename std::enable_if< (ii<std::tuple_size<T>::value), int>::type getId()
    { return std::tuple_element_t<ii,T>::value; }

    template<int ii>
    static constexpr typename std::enable_if<!(ii<std::tuple_size<T>::value), int>::type getId()
    { return 0; }

    static const int thisId = getId<i>();
  public:
    static const int size = PassIdGenerator<i-1,T>::size + thisId;
    static const int id = (thisId <= 0) ? -1 : size-1;
  };

  template< bool... b >
  class PassIdGenerator< -1, std::tuple< _bool<b>... > >
  {
    typedef std::tuple< _bool<b>... > T;
  public:
    //start id
    static const int size = 0;
    static const int id = -1;
  };


  /**
   *  \brief Default traits class for models
   */
  template <class GridImp, class ProblemImp,class ActivationImp=std::tuple<>, int maxModelParameterSize = 0 >
  struct DefaultModelTraits
  {
    typedef ProblemImp                                                    ProblemType;
    typedef typename ProblemType::FunctionSpaceType                       FunctionSpaceType;
    typedef GridImp                                                       GridType;
    static constexpr int dimDomain = FunctionSpaceType::dimDomain;
    static constexpr int dimRange  = FunctionSpaceType::dimRange;
    static constexpr int dimGradRange = dimRange * dimDomain;

    typedef typename FunctionSpaceType::DomainFieldType                   DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType                    RangeFieldType;

    typedef typename FunctionSpaceType::DomainType                        DomainType;
    typedef typename FunctionSpaceType::RangeType                         RangeType;
    typedef Dune::FieldVector< DomainFieldType, dimDomain-1 >             FaceDomainType;

    typedef typename FunctionSpaceType::JacobianRangeType                 JacobianRangeType;
    typedef JacobianRangeType                                             FluxRangeType;
    typedef Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain >     DiffusionMatrixType;
    typedef Dune::FieldMatrix< RangeFieldType, dimGradRange, dimDomain >  DiffusionRangeType;

    //id generator for pass ids
    template< int i >
    using IdGenerator = PassIdGenerator<i,ActivationImp>;

    //function space short cut
    template< int... dims >
    using ParameterSpaces = std::tuple< Dune::Fem::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dims >... >;

    //default
    typedef ParameterSpaces<>                                             ParameterSpacesType;

    static const int modelParameterSize = IdGenerator<maxModelParameterSize>::size;
  };

  struct IsFemDGModel
  {};

  /**
   * \brief Describes a general analytical model
   *
   * \ingroup AnalyticalModels
   *
   * The main goal of this class is provide all the analytical data of
   * the partial differential equation.
   *
   * This model class describes the data for the following partial differential equation
   *
   * Let \f$ Omega\subset\mathbb{R}^d \f$ be a domain and \f$ \partial\Omega = \Gamma_D \cup \Gamma_N \f$
   * be the boundary of the domain which can be composed into a
   * Dirichlet and a Neumann boundary part.
   *
   * Then we are searching a solution \f$ (u,v) \f$ such that
   *
   * \f{eqnarray*}{ v + \nabla a(u)                 & = & 0             & \text{on }\Omega \\
   * \partial_t u + \nabla \cdot (F(u)+A(u,v)) +S_1(u) + S_2(u) & = & 0 & \text{on }\Omega \\
   *                                u                & = & g_D          & \text{on }\Gamma_D\\
   *                                   F(u) \cdot n  & = & g_{N,1}      & \text{on }\Gamma_N\\
   *                                  A(u,v) \cdot n & = & g_{N,2}      & \text{on }\Gamma_N
   *                               \f}
   *
   * where each class methods describes an analytical function.
   *
   * | Method                     | formular                                  |
   * | -------------------------- | ----------------------------------------- |
   * | mass()                     | \f$ R   \f$                               |
   * | stiffSource()              | \f$ S_1 \f$                               |
   * | nonStiffSource()           | \f$ S_2 \f$                               |
   * | boundaryValue()            | \f$ g_D\f$                                |
   * | hasBoundaryValue()         | true if \f$ x \in \Gamma_D \f$            |
   * | hasRobinBoundaryValue()    | true if \f$ x \in \Gamma_R \f$            |
   * | diffusion()                | \f$ A \f$                                 |
   * | advection()                | \f$ F \f$                                 |
   * | jacobian()                 | \f$ a \f$                                 |
   * | boundaryFlux()             | \f$ g_{N,1} \f$                           |
   * | diffusionBoundaryFlux()    | \f$ g_{N,2} \f$                           |
   * | hasFlux()                  | true if \f$ F\neq 0\f$, false otherwise   |
   * | hasStiffSource()           | true if \f$ S_1\neq 0\f$, false otherwise |
   * | hasNonStiffSource()        | true if \f$ S_2\neq 0\f$, false otherwise |
   *
   * \attention \f$F(u)\f$ and \f$A(u,v)\f$ are matrix valued, and therefore the
   * divergence is defined as
   * \f[ \Delta M = \nabla \cdot (\nabla \cdot (M_{i\cdot})^t)_{i\in 1\dots n} \f]
   * for a matrix \f$M\in \mathbf{M}^{n\times m}\f$.
   *
   * \tparam Traits traits class
   */
  template < class Traits >
  class DefaultModel : public IsFemDGModel
  {
  public:
    static const int modelParameterSize = Traits::modelParameterSize;
    static const bool scalingLimiter = false ;

    static const int dimDomain = Traits::dimDomain;
    static const int dimRange  = Traits::dimRange;

    typedef typename Traits::FunctionSpaceType           FunctionSpaceType;
    typedef FunctionSpaceType                            DFunctionSpaceType;
    typedef FunctionSpaceType                            RFunctionSpaceType;
    typedef typename Traits::DomainType                  DomainType;
    typedef typename Traits::RangeType                   RangeType;
    typedef typename Traits::GradientType                GradientType;
    typedef typename Traits::FluxRangeType               FluxRangeType;
    typedef typename Traits::FaceDomainType              FaceDomainType;
    typedef typename Traits::JacobianRangeType           JacobianRangeType;

    // default is 0,...,dimRange-1
    typedef Dune::FieldVector< int, dimRange > LimitedRangeType;

    explicit DefaultModel( double time = 0 )
      : time_( time )
    {
      for( int d=0; d<dimRange; ++d )
        defaultLimitedRange_[ d ] = d;
    }

    /**
     * \brief return set of components to be modified by the limiter */
    const LimitedRangeType& limitedRange() const
    {
      return defaultLimitedRange_;
    }

    /**
     * \brief returns whether the mass term is not the identity
     */
    inline bool hasMass () const { return false ; }

    /**
     * \brief returns whether the advection term is zero
     * or not, i.e. \f$ F\neq 0 \f$
     */
    inline bool hasFlux () const { return false ; }

    /**
     *  \brief returns whether the stiff source term is zero
     *  or not, i.e.\ returns whether \f$ S_1\neq 0 \f$
     */
    inline bool hasStiffSource () const { return false ; }

    /**
     *  \brief returns whether the non stiff source term is zero
     *  or not, i.e.\f$ S_2\neq 0 \f$
     */
    inline bool hasNonStiffSource () const { return false ; }

    /**
     *  \brief returns whether the shock indicator needs to be computed when the
     *  LimitPass is used.
     */
    bool calculateIndicator () const { return true; }

    /**
     *  \brief returns whether the shock indicator needs to be computed when the
     *  LimitPass is used.
     */
    void obtainBounds( RangeType& globalMin, RangeType& globalMax ) const
    {
      globalMin =  1e308;
      globalMax = -1e308;
    }

    /**
     *  \brief compute pressure and temperature for compressible flow
     */
    template <class State>
    void pressureTemperatur( const State& u, double& p, double& T ) const
    {
      p = 0.0;
      T = 0.0;
    }

    /**
     * \brief returns the mass factor term \f$ R \f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] diagonal of the mass term, i.e. the result \f$ R(u) \f$
     *
     */
    template <class LocalEvaluation>
    inline void mass (const LocalEvaluation& local,
                      const RangeType& u,
                      RangeType& diag) const
    {
      diag = 1 ;
    }

    /**
     * \brief returns the stiff source term \f$ S_1 \f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  du evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] s the result \f$ S_1(u) \f$
     *
     * \returns The time step restriction which is given
     * by the stiff source.
     */
    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               const GradientType& du,
                               RangeType & s) const
    {
      return stiffSource( local, u, s );
    }

    /**
     * \brief returns the stiff source term \f$ S_1 \f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  du evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] s the result \f$ S_1(u) \f$
     *
     * \returns The time step restriction which is given
     * by the stiff source.
     */
    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& jac,
                               RangeType & s) const
    {
      return stiffSource( local, u, s );
    }

    /**
     * \brief returns the stiff source term \f$ S_1 \f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] s the result \f$ S_1(u) \f$
     *
     * \returns The time step restriction which is given
     * by the stiff source.
     */
    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

    /**
     * \brief returns the non stiff source term \f$ S_2 \f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  du evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] s the result \f$ S_2(u) \f$
     *
     * \returns The time step restriction which is given
     * by the non stiff source.
     */
    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  const GradientType& du,
                                  RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

    /**
     * \brief returns the non stiff source term \f$ S_2 \f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  du evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] s the result \f$ S_2(u) \f$
     *
     * \returns The time step restriction which is given
     * by the non stiff source.
     */
    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& jac,
                                  RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

    /**
     * \brief returns the non stiff source term \f$ S_2 \f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] s the result \f$ S_2(u) \f$
     *
     * \returns The time step restriction which is given
     * by the non stiff source.
     */
    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

   /**
     * \brief advection term \f$F\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] f the result \f$F(u)\f$
     */
    template <class LocalEvaluation>
    inline void advection (const LocalEvaluation& local,
                           const RangeType& u,
                           FluxRangeType & f) const
    {
      f = 0;
    }

    /**
     * \brief velocity calculation
     *
     * \note This method is internally called by by advection() and
     * externally called by some numerical fluxes (i.e. UpwindFlux etc.).
     *
     * \param[in] local local evaluation
     */
    template <class LocalEvaluation>
    DomainType velocity (const LocalEvaluation& local,
                         const RangeType& u ) const
    {
      return (DomainType)0;
    }

    /**
     * \brief diffusion term \f$a\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] a the result \f$ a(u) \f$
     */
    template <class LocalEvaluation >
    inline void jacobian (const LocalEvaluation& local,
                          const RangeType& u,
                          JacobianRangeType& a) const
    {
      a = 0;

      assert( a.rows == dimRange * dimDomain );
      assert( a.cols == dimDomain );

      for (int r=0;r<dimRange;r++)
        for (int d=0;d<dimDomain;d++)
          a[dimDomain*r+d][d] = u[r];
    }



   /**
     * \brief diffusion term \f$A\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  jac evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] A The result \f$ A(u) \f$
     */
    template <class LocalEvaluation, class JacobianType>
    inline void diffusion (const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianType& jac,
                           FluxRangeType& A) const
    {}

    /**
     * \brief diffusion term \f$A\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[in]  vecJac evaluation of the gradient of the local function, i.e. \f$\nabla u_E( \hat{x} )\f$
     * \param[out] A The result \f$ A(u) \f$
     */
    template <class LocalEvaluation>
    inline void diffusion (const LocalEvaluation& local,
                           const RangeType& u,
                           const GradientType& vecJac,
                           FluxRangeType& A) const
    {
      Dune::Fem::FieldMatrixConverter< GradientType, FluxRangeType> jac( vecJac );
      diffusion( local, u, jac, A );
    }


    /**
     * \brief return the maximal wave speed due to the advection part
     *
     * \todo Add example for advspeed != totalspeed or explain the difference.
     *
     * \param[in]  local local evaluation
     * \param[in]  normal the normal of the intersection
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] advspeed the maximal wave speed, i.e. the largest eigenvalue of the flux jacobian.
     *                      Needed for viscosity terms, for example local Lax-Friedrich flux.
     * \param[out] totalspeed the maximal wave speed, i.e. the largest eigenvalue of the flux jacobian.
     *                        Needed for the time step estimation.
     */
    template <class LocalEvaluation>
    inline void maxWaveSpeed (const LocalEvaluation&,
                              const DomainType& normal,
                              const RangeType& u,
                              double& advspeed,
                              double& totalspeed ) const
    {
      advspeed = totalspeed = 0;
    }

    template <class LocalEvaluation>
    [[deprecated("Use maxWaveSpeed instead!")]]
    inline void maxSpeed (const LocalEvaluation&,
                          const DomainType& normal,
                          const RangeType& u,
                          double& advspeed,
                          double& totalspeed ) const
    {
      DUNE_THROW(NotImplemented,"Deprecated: use maxWaveSpeed instead!");
    }

    /**
     * \brief returns the maximal time step size which is given through
     * the diffusion term.
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     */
    template <class LocalEvaluation>
    inline double diffusionTimeStep( const LocalEvaluation& local,
                                     const RangeType& u ) const
    {
      return 0;
    }

    /**
     * \brief checks for existence of Robin boundary values
     *
     * \param[in] local local evaluation
     */
    template< class LocalEvaluation >
    inline bool hasRobinBoundaryValue(const LocalEvaluation& local) const
    {
      return false;
    }

    /**
     * \brief checks for existence of Dirichlet boundary values
     *
     * \param[in] local local evaluation
     */
    template< class LocalEvaluation >
    inline bool hasBoundaryValue(const LocalEvaluation& local) const
    {
      return true;
    }

    /**
     * \brief Neuman boundary values \f$g_{N,1}\f$ for pass2
     *
     * \param[in]  local local evaluation
     * \param[in]  uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[in]  jac evaluation of the gradient of the local function, i.e. \f$\nabla u_{E^+}( \hat{x} )\f$
     * \param[out] gLeft the Neumann boundary value \f$g_{N,1}\f$.
     */
    template <class LocalEvalution>
    inline double boundaryFlux (const LocalEvalution& local,
                                const RangeType& uLeft,
                                const JacobianRangeType& jac,
                                RangeType& gLeft) const
    {
      gLeft = 0.;
      return 0.;
    }

    /**
     * \brief Neuman boundary values \f$g_{N,1}\f$ for pass2
     *
     * \param[in]  local local evaluation
     * \param[in]  uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[out] gLeft the Neumann boundary value \f$g_{N,1}\f$.
     */
    template <class LocalEvalution>
    inline double boundaryFlux (const LocalEvalution& local,
                                const RangeType& uLeft,
                                RangeType& gLeft) const
    {
      gLeft = 0.;
      return 0.;
    }

    /**
     * \brief returns the diffusion boundary flux \f$g_{N,2}\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[in]  gradLeft evaluation of the gradient of the local function, i.e. \f$\nabla u_{E^+}( \hat{x} )\f$
     * \param[out] gLeft the diffusion boundary flux \f$g_{N,2}\f$
     */
    template <class LocalEvaluation>
    inline double diffusionBoundaryFlux (const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const GradientType& gradLeft,
                                         RangeType& gLeft) const
    {
      Dune::Fem::FieldMatrixConverter< GradientType, JacobianRangeType> jacLeft( gradLeft );
      return diffusionBoundaryFlux( local, uLeft, jacLeft, gLeft );
    }

    /**
     * \brief returns the diffusion boundary flux \f$g_{N,2}\f$
     *
     * \param[in]  local local evaluation
     * \param[in]  uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[in]  gradLeft evaluation of the gradient of the local function, i.e. \f$\nabla u_{E^+}( \hat{x} )\f$
     * \param[out] gLeft the diffusion boundary flux \f$g_{N,2}\f$
     */
    template <class LocalEvaluation, class JacobianRangeImp>
    inline double diffusionBoundaryFlux (const LocalEvaluation& local,
                                         const RangeType& uLeft,
                                         const JacobianRangeImp& jacLeft,
                                         RangeType& gLeft) const
    {
      std::cerr <<"diffusionBoundaryFlux shouldn't be used in this model" <<std::endl;
      abort();
    }



    /**
     * \brief returns the Dirichlet boundary values
     *
     * \param[in]  local local evaluation
     * \param[in]  uLeft evaluation of the local function, i.e. \f$ u_{E^+}( \hat{x} ) \f$
     * \param[out] uRight the Dirichlet boundary value
     */
    template <class LocalEvaluation>
    inline void boundaryValue (const LocalEvaluation& local,
                               const RangeType& uLeft,
                               RangeType& uRight) const
    {
      std::cerr << "Method DefaultModel::boundaryValue not implemented." << std::endl;
      abort();
    }

    /** \brief check with the problem setup whether or not a cell is allowed to be refined
     *
     *  \param[in] local local evaluation
     *  \return true if the cell can be refined
     */
    template <class LocalEvaluation >
    inline bool allowsRefinement( const LocalEvaluation& local ) const
    {
      return true;
    }

    /**
     * \brief returns the maximum eigen value of \f$ K \f$.
     *
     * \param[in]  local local evaluation
     * \param[in]  u evaluation of the local function, i.e. \f$ u_E( \hat{x} ) \f$
     * \param[out] the maximal eigen value
     */
    template <class LocalEvaluation>
    inline void eigenValues (const LocalEvaluation& local,
                             const RangeType& u,
                             RangeType& maxValue) const
    {
      DUNE_THROW(Dune::NotImplemented,"DefaultModel::eigenValues is not implemented");
    }

    template< class LocalEvaluation >
    inline double penaltyBoundary( const LocalEvaluation& local,
                                   const RangeType& uLeft ) const
    {
      return penaltyFactor( local, local, uLeft, uLeft );
    }

    template <class LocalEvaluation>
    inline double penaltyFactor (const LocalEvaluation& left,
                                 const LocalEvaluation& right,
                                 const RangeType& uLeft,
                                 const RangeType& uRight) const
    {
      DUNE_THROW(Dune::NotImplemented,"DefaultModel::penaltyValues is not implemented");
      return 0.0;
    }

    template <class Entity>
    void setEntity( const Entity& entity ) const
    {
    }

    void setTime (double time)
    {
      time_ = time;
    }

    double time () const
    {
      return time_;
    }

    template< class DiscreteFunction, class... DiscreteFunctions >
    void eocErrors( const DiscreteFunction& df, const DiscreteFunctions&... dfs ) const
    {
      //default l2 error
      EOCErrorList::setErrors<L2EOCError>( *this, df );
      eocErrors( dfs... );
    }

    void eocErrors() const
    {}

  protected:
    double time_;
    LimitedRangeType defaultLimitedRange_;
  };

}
}
#endif
