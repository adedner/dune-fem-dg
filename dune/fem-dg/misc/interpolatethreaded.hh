#ifndef DUNE_FEM_DG_INTERPOLATE_HH
#define DUNE_FEM_DG_INTERPOLATE_HH

#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/misc/threads/threaditerator.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    // interpolate
    // -----------


    template <class GridFunction, class DiscreteFunction, class Iterators >
    struct RunIterpolateThreaded
    {
      const GridFunction& u_;
      DiscreteFunction& v_;
      Iterators& iterators_;

      RunIterpolateThreaded( const GridFunction& u, DiscreteFunction& v, Iterators& iterators )
        : u_( u ), v_( v ), iterators_( iterators )
      {
        // init caches
        (*this)( true );
      }

      void operator() ( const bool stopAfterFirst = false )
      {
        // reserve memory for local dof vector
        std::vector< typename DiscreteFunction::RangeFieldType > ldv;
        ldv.reserve( v_.space().blockMapper().maxNumDofs() * DiscreteFunction::DiscreteFunctionSpaceType::localBlockSize );

        typename GridFunction::LocalFunctionType uLocal( u_ );
        LocalInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType > interpolation( v_.space() );

        const auto endit = iterators_.end();
        // iterate over selected partition
        for( auto it = iterators_.begin(); it != endit; ++it )
        {
          const auto& entity = *it;

          // obtain local interpolation
          auto guard = bindGuard( interpolation, entity );

          // resize local dof vector
          ldv.resize( v_.space().basisFunctionSet( entity ).size() );

          // interpolate u locally
          uLocal.init( entity );
          interpolation( uLocal, ldv );

          // write local dofs into v
          v_.setLocalDofs( entity, ldv );

          // if true return since only one element is needed for
          // initialization of singleton caches
          if( stopAfterFirst ) return ;
        }
      }
    };

    /**
     * \function interpolate
     * \ingroup  DiscreteFunctionSpace
     * \brief    perform native interpolation of a discrete function space
     *
     * By definition of its degrees of freedom, each discrete function space
     * has a native interpolation, which can be computed very quickly.
     *
     * For example, the native interpolation of a Lagrange discrete function
     * space is the evaluation in its Lagrange points.
     * An orthonormal DG space would instead perform an \f$L^2\f$-Projection.
     *
     * The actual implementation must locally be provided by the discrete
     * function space through the method
     * \code
     * template< class LocalFunction, class LocalDofVector >
     * void interpolate ( const LocalFunction &f, LocalDofVector &dofs ) const;
     * \endcode
     *
     * \param[in]   u  grid function to interpolate
     * \param[out]  v  discrete function to represent the interpolation
     */
    template< class GridFunction, class DiscreteFunction >
    static inline void interpolateThreaded ( const GridFunction &u, DiscreteFunction &v )
    {
      typedef typename DiscreteFunction :: DiscreteFunctionSpaceType ::
        GridPartType GridPartType;
      typedef Dune::Fem::ThreadIterator< GridPartType, Dune::Interior_Partition> ThreadIteratorType;

      // create thread range iterators
      ThreadIteratorType iterators( v.space().gridPart() );
      iterators.update();

      RunIterpolateThreaded< GridFunction, DiscreteFunction, ThreadIteratorType>
        runThread( u, v, iterators );

      ////////////////////////////////////////////////////////////
      // BEGIN PARALLEL REGION, first stage, element integrals
      ////////////////////////////////////////////////////////////
      try
      {
        // see threadpool.hh
        Fem :: MPIManager :: run( runThread );
      }
      catch ( const Dune::Fem::SingleThreadModeError& e )
      {
        Dune::Fem::interpolate( u, v );
      }

      /////////////////////////////////////////////////
      // END PARALLEL REGION
      /////////////////////////////////////////////////
    }

    template< class Function, class DiscreteFunction >
    static inline std::enable_if_t< !std::is_convertible< Function, HasLocalFunction >::value >
    interpolateThreaded ( const Function &u, DiscreteFunction &v )
    {
      typedef typename DiscreteFunction :: DiscreteFunctionSpaceType :: GridPartType  GridPartType;
      typedef GridFunctionAdapter< Function, GridPartType > GridFunctionType;

      GridFunctionType uGrid( "uGrid", u, v.space().gridPart() );

      interpolateThreaded( uGrid, v );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
