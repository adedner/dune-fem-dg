#ifndef DUNE_FEM_DG_DEFAULTMODEL_HH
#define DUNE_FEM_DG_DEFAULTMODEL_HH

#include<limits>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/fmatrixconverter.hh>

/**********************************************
 * Default model
 *********************************************/
template < class Traits >
class DefaultModel
{
public:
  typedef std::tuple <>  ModelParameter;

  static const int dimDomain = Traits :: dimDomain ;
  static const int dimRange  = Traits :: dimRange ;

  typedef typename Traits :: DomainType                              DomainType;
  typedef typename Traits :: RangeType                               RangeType;
  typedef typename Traits :: GradientType                            GradientType;
  typedef typename Traits :: FluxRangeType                           FluxRangeType;
  typedef typename Traits :: FaceDomainType                          FaceDomainType;
  typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

  typedef typename Traits :: EntityType                       EntityType;
  typedef typename Traits :: IntersectionType                 IntersectionType;

  void setTime ( double time ) {}

  inline bool hasFlux() const { return false ; }

  inline bool hasStiffSource() const { return false ; }
  inline bool hasNonStiffSource() const { return false ; }


  inline double nonStiffSource( const EntityType& en,
                                const double time,
                                const DomainType& x,
                                const RangeType& u,
                                const GradientType& du,
                                RangeType & s) const
  {
    s = 0 ;
    return std::numeric_limits< double > :: max ();
  }

  inline double nonStiffSource( const EntityType& en,
                                const double time,
                                const DomainType& x,
                                const RangeType& u,
                                const JacobianRangeType& jac,
                                RangeType & s) const
  {
    s = 0 ;
    return std::numeric_limits< double > :: max ();
  }

  inline double nonStiffSource( const EntityType& en,
                                const double time,
                                const DomainType& x,
                                const RangeType& u,
                                RangeType & s) const
  {
    s = 0 ;
    return std::numeric_limits< double > :: max ();
  }

  inline double stiffSource( const EntityType& en,
                             const double time,
                             const DomainType& x,
                             const RangeType& u,
                             const GradientType& du,
                             RangeType & s) const
  {
    return stiffSource( en, time, x, u, s );
  }


  inline double stiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const JacobianRangeType& jac,
                        RangeType & s) const
  {
    return stiffSource( en, time, x, u, s );
  }

  inline double stiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        RangeType & s) const
  {
    s = 0 ;
    return std::numeric_limits< double > :: max ();
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
  inline  void advection(const EntityType& en,
                         const double time,
                         const DomainType& x,
                         const RangeType& u,
                         FluxRangeType & f) const
  {
    f = 0;
  }

  /**
   * @brief velocity calculation, is called by advection()
   */
  inline  void velocity(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        DomainType& v) const
  {
    v = 0 ;
  }

  /**
   * @brief diffusion term \f$a\f$
   */
  inline void jacobian(const EntityType& en,
                       const double time,
                       const DomainType& x,
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

  template <class LocalEvaluation>
  inline void eigenValues(const LocalEvaluation& local,
                          const RangeType& u,
                          RangeType& maxValue) const
  {
    DUNE_THROW(Dune::NotImplemented,"DefaultModel::eigenValues is not implemented");
  }

  inline double penaltyFactor( const double time,
                               const DomainType& xInside,
                               const EntityType& inside,
                               const RangeType& uLeft ) const
                               //const EntityType& outside = inside,
                               //const RangeType& uRight = uLeft ) const
  {
    DUNE_THROW(Dune::NotImplemented,"DefaultModel::penaltyValues is not implemented");
    return 0.0;
  }

  /**
   * @brief diffusion term \f$A\f$
   */
  template <class JacobianType>
  inline void diffusion(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const JacobianType& jac,
                        FluxRangeType& A) const
  {
  }

  inline void diffusion(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const GradientType& vecJac,
                        FluxRangeType& A) const
  {
    Dune::Fem::FieldMatrixConverter< GradientType, FluxRangeType> jac( vecJac );
    diffusion( en, time, x, u, jac, A );
  }

  //! maximal wave speed due to advection part
  template <class LocalEvaluation>
  inline void maxSpeed( const LocalEvaluation&,
                        const DomainType& normal,
                        const RangeType& u,
                        double& advspeed,
                        double& totalspeed ) const
  {
    advspeed = totalspeed = 0;
  }

  //! maximal wave speed due to diffusion part
  inline double diffusionTimeStep( const IntersectionType &it,
                                   const double enVolume,
                                   const double circumEstimate,
                                   const double time,
                                   const FaceDomainType& x,
                                   const RangeType& u ) const
  {
    return 0;
  }

public:
  /**
   * @brief checks for existence of dirichlet boundary values
   */
  inline bool hasBoundaryValue(const IntersectionType& it,
                               const double time,
                               const FaceDomainType& x) const
  {
    return true;
  }

  /**
   * @brief neuman boundary values \f$g_N\f$ for pass1
   */
  template <class LocalEvalution>
  inline double boundaryFlux(const LocalEvalution&,
                             const RangeType& uLeft,
                             const JacobianRangeType& jac,
                             RangeType& gLeft) const
  {
    gLeft = 0.;
    return 0.;
  }


  /**
   * @brief neuman boundary values \f$g_N\f$ for pass1
   */
  template <class LocalEvalution>
  inline double boundaryFlux(const LocalEvalution&,
                             const RangeType& uLeft,
                             RangeType& gLeft) const
  {
    gLeft = 0.;
    return 0.;
  }

  /**
   * @brief diffusion boundary flux
   */
  template <class LocalEvalution>
  inline double diffusionBoundaryFlux( const LocalEvalution&,
                                       const RangeType&,
                                       const JacobianRangeType&,
                                       RangeType& ) const
  {
    std::cerr <<"Method DefaultModel::diffusionBoundaryFlux not implemented." <<std::endl;
    abort();
  }



  /**
   * @brief dirichlet boundary values
   */
  template <class LocalEvalution>
  inline void boundaryValue(const LocalEvalution&,
                             const RangeType&,
                             RangeType& ) const
  {
    std::cerr << "Method DefaultModel::boundaryValue not implemented." << std::endl;
    abort();
  }

  /** \brief check with the problem setup whether or not a cell is allowed to be refined
   *
   *  \param[in] it Intersection
   *  \param[in] time Current model time
   *  \param[in] x Point w.r.t. intersection \a it
   *  \return true if the cell can be refined
   */
  inline bool allowsRefinement( const IntersectionType& it,
                                const double time,
                                const FaceDomainType& x ) const
  {
    return true;
  }
};
#endif
