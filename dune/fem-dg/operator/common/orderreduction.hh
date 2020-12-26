#ifndef DUNE_FEM_DG_ORDERREDUCTION_HH
#define DUNE_FEM_DG_ORDERREDUCTION_HH

#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/hpdg/legendre.hh>
#include <dune/fem/space/hpdg/orthogonal.hh>

#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/grid/utility/persistentcontainer.hh>

namespace Dune {

  namespace Fem {

    /** \brief OrderReduction reduces the polynomial order of a discrete
     *         function by projecting from P_E -> P-1_E on each element.
     *
     *  \note  This is only implemented for hierarchic spaces currently.
     */
    template <class DF>
    class OrderReduction
    {
    public:
      typedef DF  DomainFunctionType;
      typedef DF  RangeFunctionType;

      typedef typename DomainFunctionType :: DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename DomainSpaceType :: GridType GridType;
      typedef typename DomainSpaceType :: EntityType  EntityType;

      typedef DomainSpaceType RangeSpaceType;

      struct ElementOrder
      {
        int order_;
        ElementOrder() : order_( -1 ) {}
        int order() const { return order_; }
        void set( const int order ) { order_ = order; }
      };

      typedef Dune::PersistentContainer< GridType, ElementOrder >
        PersistentContainerType;

      //static_assert(Dune::Fem::Capabilities::isHierarchic< DomainSpaceType > :: v, "OrderReduction is only implemented for hierarchic discrete function spaces");

      //mutable PersistentContainerType maxRelOrder_;

      OrderReduction( const DomainSpaceType& spc )
        //: maxRelOrder_( spc.gridPart().grid(), 0 )
      {
        //maxRelOrder_.resize();
      }

      OrderReduction(const GridType& grid )
        //: maxRelOrder_( grid, 0 )
      {
        //maxRelOrder_.resize();
      }

      void operator () (const DomainFunctionType& arg, RangeFunctionType& dest )
      {
        const auto& basisSets = arg.space().basisFunctionSets();

        typedef typename RangeFunctionType::RangeFieldType  RangeFieldType;

        Dune::DynamicVector< RangeFieldType > localDofs;
        ConstLocalFunction< DomainFunctionType > hlf( arg );

        for( const auto& entity : arg.space() )
        {
          auto guard = bindGuard( hlf, entity );

          const int numDofs = hlf.numDofs();
          const int order   = hlf.order();

          const int lowerOrder = std::max( order-1, 0 );

          // get basis function set with one order lower
          const auto& lowerBase = basisSets.basisFunctionSet( entity, lowerOrder );

          localDofs.resize( numDofs );
          // copy dofs
          localDofs = hlf.localDofVector();

          // erase higher order moments
          for( int i = lowerBase.size(); i<numDofs; ++i )
          {
            localDofs[ i ] = 0;
          }

          // set local dofs of destination
          dest.setLocalDofs( entity, localDofs );
        }
      }

      int maxRelevantOrder( const EntityType& entity ) const
      {
        return -1;
        //return maxRelOrder_[ entity ].order();
      }
    };

    /** \brief Specialization of OrderReduction for TupleDiscreteFunction.
     *         The implementation default to the above for each tuple member.a
     *
     *  \note  This is only implemented for hierarchic spaces currently.
     */
    template < class... DFs >
    class OrderReduction< TupleDiscreteFunction< DFs... > >
    {
    public:
      typedef TupleDiscreteFunction< DFs... > DomainFunctionType;
      typedef DomainFunctionType  RangeFunctionType;

      typedef typename DomainFunctionType :: DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename DomainSpaceType :: EntityType  EntityType;
      typedef DomainSpaceType RangeSpaceType;

      typedef typename DomainSpaceType :: GridType GridType;
      typedef typename DomainFunctionType :: template SubDiscreteFunction< 0 >::Type SubFunction;
      typedef typename OrderReduction< SubFunction > :: PersistentContainerType  PersistentContainerType;

      //static_assert(Dune::Fem::Capabilities::isHierarchic< DomainSpaceType > :: v, "OrderReduction is only implemented for hierarchic discrete function spaces");

      OrderReduction( const DomainSpaceType& )
      {
      }

      void operator () (const DomainFunctionType& arg, RangeFunctionType& dest )
      {
        Hybrid::forEach( typename DomainFunctionType::Sequence{}, [ & ]( auto i )
            {
              typedef typename DomainFunctionType :: template SubDiscreteFunction< i >::Type DF;
              OrderReduction< DF >( arg.space().gridPart().grid() )
                    ( arg.template subDiscreteFunction< i >(), dest.template subDiscreteFunction< i >() );
            } );
      }

      int maxRelevantOrder( const EntityType& entity ) const
      {
        return -1;
        //return maxRelOrder_[ entity ].order();
      }
    };

  } // end namespace Fem
} // end namespace Dune

#endif
