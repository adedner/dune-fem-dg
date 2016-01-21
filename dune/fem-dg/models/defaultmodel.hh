#ifndef DUNE_FEM_DG_DEFAULTMODEL_HH
#define DUNE_FEM_DG_DEFAULTMODEL_HH

#include<limits>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/fmatrixconverter.hh>

namespace Dune
{
namespace Fem
{

  /**
   *  \brief Default traits class for models
   */
  template <class GridPartImp, class ProblemImp>
  struct DefaultModelTraits
  {
    typedef ProblemImp                                                     ProblemType;
    typedef GridPartImp                                                    GridPartType;
    typedef typename GridPartType::GridType                                GridType;
    static const int dimDomain = GridType::dimensionworld;
    static const int dimRange = ProblemType::dimRange;
    static const int dimGradRange = dimRange * dimDomain;

    typedef typename GridType::ctype                                       RangeFieldType;
    typedef typename GridType::ctype                                       DomainFieldType;

    typedef Dune::FieldVector< DomainFieldType, dimDomain >                DomainType;
    typedef Dune::FieldVector< DomainFieldType, dimDomain-1 >              FaceDomainType;
    typedef Dune::FieldVector< RangeFieldType, dimRange >                  RangeType;

    typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain >       FluxRangeType;
    typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain >       JacobianRangeType;
    typedef Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain >      DiffusionMatrixType;
    typedef Dune::FieldMatrix< RangeFieldType, dimGradRange, dimDomain >   DiffusionRangeType;

    typedef typename GridType::template Codim< 0 >::Entity                 EntityType;
    typedef typename GridPartType::IntersectionIteratorType                IntersectionIterator;
    typedef typename IntersectionIterator::Intersection                    IntersectionType;
  };

  /**
   * \brief Analytical model interface
   *
   * \ingroup AnalyticalModels
   */
  template < class Traits >
  class DefaultModel
  {
  public:
    typedef std::tuple <>  ModelParameter;

    static const int dimDomain = Traits::dimDomain;
    static const int dimRange  = Traits::dimRange;

    typedef typename Traits::DomainType                  DomainType;
    typedef typename Traits::RangeType                   RangeType;
    typedef typename Traits::GradientType                GradientType;
    typedef typename Traits::FluxRangeType               FluxRangeType;
    typedef typename Traits::FaceDomainType              FaceDomainType;
    typedef typename Traits::JacobianRangeType           JacobianRangeType;

    typedef typename Traits::EntityType                  EntityType;
    typedef typename Traits::IntersectionType            IntersectionType;

    void setTime (double time) {}

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

    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  const GradientType& du,
                                  RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  const JacobianRangeType& jac,
                                  RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

    template <class LocalEvaluation>
    inline double nonStiffSource (const LocalEvaluation& local,
                                  const RangeType& u,
                                  RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               const GradientType& du,
                               RangeType & s) const
    {
      return stiffSource( local, u, s );
    }


    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               const JacobianRangeType& jac,
                               RangeType & s) const
    {
      return stiffSource( local, u, s );
    }

    template <class LocalEvaluation>
    inline double stiffSource (const LocalEvaluation& local,
                               const RangeType& u,
                               RangeType & s) const
    {
      s = 0 ;
      return std::numeric_limits< double > :: max ();
    }

    /**
     * \brief advection term \f$F\f$
     *
     * \param local local evaluation
     * \param u \f$U\f$
     * \param f \f$F(U)\f$
     */
    template <class LocalEvaluation>
    inline void advection (const LocalEvaluation& local,
                           const RangeType& u,
                           FluxRangeType & f) const
    {
      f = 0;
    }

    /**
     * \brief velocity calculation, is called by advection()
     */
    template <class LocalEvaluation>
    DomainType velocity (const LocalEvaluation& local) const
    {
      return (DomainType)0;
    }

    /**
     * \brief diffusion term \f$a\f$
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

    template <class LocalEvaluation>
    inline void eigenValues (const LocalEvaluation& local,
                             const RangeType& u,
                             RangeType& maxValue) const
    {
      DUNE_THROW(Dune::NotImplemented,"DefaultModel::eigenValues is not implemented");
    }

    inline double penaltyFactor( const double time,
                                 const DomainType& xInside,
                                 const EntityType& inside,
                                 const RangeType& uLeft ) const
    {
      DUNE_THROW(Dune::NotImplemented,"DefaultModel::penaltyValues is not implemented");
      return 0.0;
    }

    /**
     * \brief diffusion term \f$A\f$
     */
    template <class LocalEvaluation, class JacobianType>
    inline void diffusion (const LocalEvaluation& local,
                           const RangeType& u,
                           const JacobianType& jacu,
                           FluxRangeType& A) const
    {}

    template <class LocalEvaluation>
    inline void diffusion (const LocalEvaluation& local,
                           const RangeType& u,
                           const GradientType& vecJac,
                           FluxRangeType& A) const
    {
      Dune::Fem::FieldMatrixConverter< GradientType, FluxRangeType> jac( vecJac );
      diffusion( local, u, jac, A );
    }

    //! maximal wave speed due to advection part
    template <class LocalEvaluation>
    inline void maxSpeed (const LocalEvaluation&,
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

    /**
     * \brief checks for existence of dirichlet boundary values
     */
    inline bool hasBoundaryValue(const IntersectionType& it,
                                 const double time,
                                 const FaceDomainType& x) const
    {
      return true;
    }

    /**
     * \brief neuman boundary values \f$g_N\f$ for pass1
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
     * \brief neuman boundary values \f$g_N\f$ for pass1
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
     * \brief diffusion boundary flux
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

    /** \brief boundary flux for the diffusion part
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
     * \brief dirichlet boundary values
     */
    inline void boundaryValue (const IntersectionType& it,
                               const double time,
                               const FaceDomainType& x,
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
  };

}
}
#endif
