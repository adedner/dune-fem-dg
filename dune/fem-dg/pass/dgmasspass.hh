#ifndef DUNE_FEM_DG_MASS_INVERSE_PASS_HH
#define DUNE_FEM_DG_MASS_INVERSE_PASS_HH

#include <cassert>
#include <iosfwd>
#include <sstream>
#include <utility>

#include <dune/common/typetraits.hh>

#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/common/tupleutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/capabilities.hh>

#include <dune/fem-dg/pass/dginversemass.hh>
#include <dune/fem-dg/misc/crs.hh>

namespace Dune
{

  namespace Fem
  {

    // DGMassInverseMassPass
    // ---------------------

    /** \brief Pass applying the local inverse mass matrix on each element
     *
     *  \ingroup Pass
     *
     *  \tparam  functionalId   pass id of functional to convert
     *  \tparam  PreviousPass   type of previous pass
     *  \tparam  id             pass id
     */
    template< class ScalarDiscreteFunctionSpace, bool inverse >
    class DGMassInverseMassImplementation
    {
    public:
      static_assert((ScalarDiscreteFunctionSpace::dimRange == 1), "dimRange of ScalarSpace > 1");
      //! \brief discrete function space type
      typedef ScalarDiscreteFunctionSpace ScalarDiscreteFunctionSpaceType ;
      typedef typename ScalarDiscreteFunctionSpaceType::RangeType ScalarRangeType;
      typedef typename ScalarDiscreteFunctionSpaceType::GridPartType GridPartType;

      static const int numScalarBasis = ScalarDiscreteFunctionSpaceType :: localBlockSize ;

      class Key
      {
      protected:
        const GridPartType & gridPart_;
        const int numScalarBasis_;
        const bool inverse_ ;
      public:
        //! constructor taking space
        Key(const GridPartType& gridPart)
          : gridPart_(gridPart),
            numScalarBasis_( numScalarBasis ),
            inverse_( inverse ) {}

        //! copy constructor
        //! returns true if indexSet pointer and numDofs are equal
        bool operator == (const Key& otherKey) const
        {
          // mapper of space is singleton
          return (&(gridPart_.indexSet()) == & (otherKey.gridPart_.indexSet()) )
                 && ( numScalarBasis_ == otherKey.numScalarBasis_ ) && ( inverse_ == otherKey.inverse_ );
        }

        //! return reference to grid part for construction
        const GridPartType& gridPart() const { return gridPart_; }
      };
      typedef Key KeyType;

    private:
      typedef typename ScalarDiscreteFunctionSpaceType :: RangeFieldType  ctype;
      typedef Dune::FieldMatrix< ctype, numScalarBasis, numScalarBasis >         LocalMatrixType ;
      typedef Dune::Fem::BlockCRSMatrix< ctype, numScalarBasis, numScalarBasis > SparseLocalMatrixType ;

      typedef Fem::CachingQuadrature< GridPartType, 0 > VolumeQuadratureType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      typedef typename ScalarDiscreteFunctionSpaceType :: BasisFunctionSetType  BasisFunctionSetType ;

      typedef typename EntityType :: Geometry  Geometry ;
      typedef typename Geometry :: GlobalCoordinate  GlobalCoordinate ;

      typedef LocalMassMatrix< ScalarDiscreteFunctionSpaceType, VolumeQuadratureType > LocalMassMatrixType ;
      typedef LocalMassMatrixImplementationDgOrthoNormal<
            ScalarDiscreteFunctionSpaceType, VolumeQuadratureType >  DgOrthoNormalMassMatrixType ;

      static const bool isOrthoNormal = false ; //Conversion< LocalMassMatrixType, DgOrthoNormalMassMatrixType > :: exists ;

      struct VectorEntry
      {
        SparseLocalMatrixType  matrix_ ;
        GlobalCoordinate center_ ;
        VectorEntry () : matrix_(), center_( std::numeric_limits< double > :: max() )
        {}
      };

      typedef std::vector< VectorEntry* >      MatrixStorageType;

    public:
      template <class Key>
      explicit DGMassInverseMassImplementation ( const Key& key )
      : scalarSpace_( const_cast< GridPartType& > (key.gridPart()) ),
        volumeQuadratureOrder_( 2 * scalarSpace_.order() ),
        matrices_( Fem :: MPIManager :: maxThreads() ),
        localMassMatrix_( scalarSpace_, [this](const int order) { return
            Capabilities::DefaultQuadrature< ScalarDiscreteFunctionSpaceType >::volumeOrder(order); }  ),
        sequence_( -1 )
      {
        assert( Fem::MPIManager::singleThreadMode() );
        setup();
      }

      ~DGMassInverseMassImplementation ()
      {
        const int maxThreads = Fem :: MPIManager :: maxThreads() ;
        for( int thread = 0; thread < maxThreads ; ++ thread )
        {
          for( size_t i=0; i<matrices_[ thread ].size(); ++i )
            delete matrices_[ thread ][ i ];
        }
      }

      //! interface method
      template < class DestinationType >
      void prepare( const DestinationType& argument, DestinationType &destination, const bool doSetup ) const
      {
        if( doSetup )
          setup();
      }

      //! interface method
      template < class DestinationType >
      void finalize( const DestinationType &argument, DestinationType &destination ) const
      {
      }

      //! apply inverse mass matrix locally
      template < class DestinationType >
      void applyLocal( const EntityType& entity,
                       const DestinationType& argument,
                       DestinationType& destination ) const
      {
        typedef typename DestinationType :: LocalFunctionType LocalFunctionType;
        const Geometry& geometry = entity.geometry();
        LocalFunctionType localDest = destination.localFunction( entity );
        const LocalFunctionType localArg = argument.localFunction( entity );
        if( isOrthoNormal && geometry.affine() )
        {
          // get inverse mass factor
          const double massFactorInv = localMassMatrix_.getAffineMassFactor( geometry );
          const double factor = ( ! inverse ) ? massFactorInv : 1.0/massFactorInv ;

          // apply factor * localArg and store in dest
          localDest.axpy( factor, localArg );
        }
        else
        {
          const int idx = scalarSpace_.gridPart().indexSet().index( entity );
          const int thread = setup( idx, entity, geometry );
          enum { dimRange = DestinationType :: DiscreteFunctionSpaceType :: dimRange };
          assert( matrices_[ thread ][ idx ] );
          multiplyBlock( dimRange, matrices_[ thread ][ idx ]->matrix_, localArg, localDest );
        }
      }

      //! apply inverse mass matrix to local function
      template< class Caller, class LocalFunction >
      void applyInverse( Caller& caller, const EntityType& entity, LocalFunction& localDest ) const
      {
        // this pass only can handle no mass
        assert( ! caller.hasMass() );
        applyInverse( entity, localDest );
      }

      //! apply inverse mass matrix to local function
      template< class LocalFunction >
      void applyInverse( LocalFunction& localDest ) const
      {
        applyInverse( localDest.entity(), localDest );
      }

      //! apply inverse mass matrix to local function
      template< class LocalFunction >
      void applyInverse( const EntityType& entity, LocalFunction& localDest ) const
      {
        enum { dimRange = LocalFunction :: dimRange };
        const int idx = scalarSpace_.gridPart().indexSet().index( entity );
        const Geometry& geometry = entity.geometry();

        const int thread = setup( idx, entity, geometry );

        static const int numBasis = dimRange * numScalarBasis ;
        assert( localDest.numDofs() == numBasis );

        // vector holding basis function evaluation
        std::array< ctype, numBasis > lfX;
        for( int i=0; i<numBasis; ++i )
          lfX[ i ] = localDest[ i ];

        assert( matrices_[ thread ][ idx ] );
        multiplyBlock( dimRange, matrices_[ thread ][ idx ]->matrix_, lfX, localDest );
      }

      template < class DestinationType >
      void compute ( const DestinationType& argument, DestinationType &destination ) const
      {
        // make sure this method is not called in multi thread mode
        assert( Fem :: MPIManager :: singleThreadMode() );

        // set pointer
        prepare( argument, destination, true );

        for( const auto& en : elements( scalarSpace_.gridPart(), Dune::Partitions::all ) )
        {
          applyLocal( en, argument, destination );
        }

        // remove pointer
        finalize( argument, destination );
      }

    protected:
      template < class ConstLocalFunction, class LocalFunctionType >
      void multiplyBlock( const int dimRange,
                          const SparseLocalMatrixType& matrix,
                          const ConstLocalFunction& arg,
                          LocalFunctionType& dest ) const
      {
        // clear destination
        // dest.clear();
        matrix.mvb( dimRange, arg, dest );

        /*
        static const int rows = LocalMatrixType::rows;
        static const int cols = LocalMatrixType::cols;

        assert( int(rows * dimRange) == dest.numDofs() );

        // apply matrix to arg and store in dest
        for(int i=0; i<rows; ++i )
        {
          for(int j=0; j<cols; ++j )
          {
            for( int r=0, ir = i * dimRange, jr = j * dimRange ;
                 r<dimRange; ++ r, ++ ir, ++ jr )
            {
              dest[ ir ] += matrix[ i ][ j ] * arg[ jr ] ;
            }
          }
        }
        */
      }

      void setup() const
      {
        assert( Fem :: MPIManager :: singleThreadMode() );
        const int gpSequence = scalarSpace_.gridPart().sequence();
        if( sequence_ != gpSequence )
        {
          const GridPartType &gridPart = scalarSpace_.gridPart();
          const int gpSize = gridPart.indexSet().size( 0 ) ;
          const int maxThreads = Fem :: MPIManager :: maxThreads() ;
          for( int thread = 0; thread < maxThreads ; ++ thread )
          {
            matrices_[ thread ].resize( gpSize, (VectorEntry *) 0 );
          }
          sequence_ = gpSequence ;
        }
      }

      int setup( const int idx,
                 const EntityType& entity,
                 const Geometry& geometry ) const
      {
        const int thread = Fem :: MPIManager :: thread();
        bool computeMatrix = false ;
        if( matrices_[ thread ][ idx ] == 0 )
        {
          VectorEntry* entry = new VectorEntry();
          matrices_[ thread ][ idx ] = entry;
          entry->center_ = geometry.center();
          computeMatrix = true ;

        }
        else
        {
          VectorEntry* entry = matrices_[ thread ][ idx ] ;
          const GlobalCoordinate center = geometry.center();
          GlobalCoordinate diff ( center );
          diff -= entry->center_;
          // for cells on the boundary we have to compute this all the time
          computeMatrix = diff.infinity_norm() > 1e-12 ;
          if( computeMatrix ) entry->center_ = center ;
        }

        if( computeMatrix )
        {
          VectorEntry* entry = matrices_[ thread ][ idx ];
          setup( idx, entity,
                 geometry,
                 scalarSpace_.basisFunctionSet( entity ),
                 entry->matrix_ );
        }
        return thread ;
      }

      void setup( const int idx,
                  const EntityType& entity,
                  const Geometry& geo,
                  const BasisFunctionSetType& set,
                  SparseLocalMatrixType& storedMatrix ) const
      {
        // clear matrix
        LocalMatrixType matrix( 0 );

        // vector holding basis function evaluation
        std::array< ScalarRangeType, numScalarBasis > phi ;

        // make sure that the number of basis functions is correct
        assert( numScalarBasis == int(set.size()) );

        // create appropriate quadrature
        VolumeQuadratureType volQuad( entity, volumeQuadratureOrder_ );

        const int volNop = volQuad.nop();
        for(int qp=0; qp<volNop; ++qp)
        {
          // calculate integration weight
          const double intel = ( volQuad.weight(qp) * geo.integrationElement(volQuad.point(qp)) );

          // eval basis functions
          set.evaluateAll(volQuad[ qp ], phi);

          for(int m=0; m<numScalarBasis; ++m)
          {
            const ScalarRangeType& phi_m = phi[ m ];
            const ctype val = intel * (phi_m * phi_m);
            matrix[ m ][ m ] += val;

            for(int k=m+1; k<numScalarBasis; ++k)
            {
              const ctype val = intel * (phi_m * phi[ k ]);
              matrix[ m ][ k ] += val;
              matrix[ k ][ m ] += val;
            }
          }
        }

        for(int m=0; m<numScalarBasis; ++m)
        {
          for(int k=0; k<numScalarBasis; ++k)
          {
            if( std::abs( matrix[ m ][ k ] ) < 1e-14 )
              matrix[ m ][ k ] = 0;
          }
        }

        // if we are looking for the inverse, then invert matrix
        if( inverse )
          matrix.invert();

        // store matrix
        storedMatrix.set( matrix );
      }

    protected:
      ScalarDiscreteFunctionSpaceType          scalarSpace_;
      const int                                volumeQuadratureOrder_;
      mutable std::vector< MatrixStorageType > matrices_;

      LocalMassMatrixType                      localMassMatrix_;

      mutable int                              sequence_ ;
    };


    /** \brief Pass applying the local inverse mass matrix on each element
     *
     *  \ingroup Pass
     *
     *  \tparam  functionalId   pass id of functional to convert
     *  \tparam  PreviousPass   type of previous pass
     *  \tparam  id             pass id
     */
    template< int functionalId, class PreviousPass, int id, bool inverse >
    class DGMassInverseMassPass
    : public Dune::Fem::LocalPass< DGInverseMassPassDiscreteModel< functionalId, PreviousPass >, PreviousPass, id >
    {
      typedef DGMassInverseMassPass< functionalId, PreviousPass, id, inverse > ThisType;
      typedef Dune::Fem::LocalPass< DGInverseMassPassDiscreteModel< functionalId, PreviousPass >, PreviousPass, id > BaseType;

    public:
      //! type of the discrete model used
      typedef DGInverseMassPassDiscreteModel< functionalId, PreviousPass > DiscreteModelType;

      //! pass ids up to here (tuple of integral constants)
      typedef typename BaseType::PassIds PassIds;

      //! \brief argument type
      typedef typename BaseType::TotalArgumentType TotalArgumentType;
      //! \brief destination type
      typedef typename BaseType::DestinationType DestinationType;

      //! \brief discrete function space type
      typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      // type of communication manager object which does communication
      typedef typename DiscreteFunctionSpaceType::template ToNewDimRange< 1 >::Type ScalarDiscreteFunctionSpaceType;
      typedef DGMassInverseMassImplementation< ScalarDiscreteFunctionSpaceType, inverse > MassInverseMassType ;
      typedef typename MassInverseMassType :: KeyType KeyType;
      typedef Fem::SingletonList< KeyType, MassInverseMassType >  MassInverseMassProviderType;

    private:
      static const std::size_t functionalPosition = DiscreteModelType::functionalPosition;

      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      typedef typename DestinationType :: LocalFunctionType    LocalFunctionType;

    public:
      using BaseType::passNumber;
      using BaseType::space;

      explicit DGMassInverseMassPass ( PreviousPass &previousPass,
                                       const DiscreteFunctionSpaceType &spc )
      : BaseType( previousPass, spc, "DGMassInverseMassPass" ),
        massInverseMass_( MassInverseMassProviderType :: getObject( KeyType( spc.gridPart() ) ) ),
        arg_( 0 ),
        dest_( 0 )
      {
      }

      //! constructor for use with thread pass
      DGMassInverseMassPass ( const DiscreteModelType& discreteModel,
                              PreviousPass &previousPass,
                              const DiscreteFunctionSpaceType &spc,
                              const int volQuadOrd  = -1,
                              const int faceQuadOrd = -1)
      : BaseType( previousPass, spc, "DGMassInverseMassPass" ),
        massInverseMass_( MassInverseMassProviderType :: getObject( KeyType( spc.gridPart() ) ) ),
        arg_( 0 ),
        dest_( 0 )
      {
      }

      ~DGMassInverseMassPass() { MassInverseMassProviderType :: removeObject( massInverseMass_ ); }

      void printTexInfo ( std::ostream &out ) const
      {
        BaseType::printTexInfo( out );
        out << "DGMassInverseMassPassname() :" << "\\\\" << std::endl;
      }

      //! this pass needs no communication
      bool requireCommunication () const { return false ; }

      //! interface method
      void prepare( const TotalArgumentType &argument, DestinationType &destination ) const
      {
        prepare( argument, destination, true );
      }

      //! prepare for ThreadPass
      void prepare( const TotalArgumentType &argument, DestinationType &destination, const bool doSetup ) const
      {
        arg_  = &argument ;
        dest_ = &destination;
        massInverseMass_.prepare( *(std::get< functionalPosition >( argument )), destination, doSetup );
      }

      //! finalize for ThreadPass
      void finalize( const TotalArgumentType &argument, DestinationType &destination, const bool ) const
      {
        finalize( argument, destination );
      }

      //! interface method
      void finalize( const TotalArgumentType &argument, DestinationType &destination ) const
      {
        massInverseMass_.finalize( *(std::get< functionalPosition >( argument )), destination );
        arg_  = 0;
        dest_ = 0;
      }

      //! apply inverse mass matrix locally
      void applyLocal( const EntityType& entity ) const
      {
        assert( arg_  );
        assert( dest_ );
        massInverseMass_.applyLocal( entity, *(std::get< functionalPosition >( *arg_ )), *dest_ );
      }

      //! apply local with neighbor checker (not needed here)
      template <class NBChecker>
      void applyLocal( const EntityType& entity, const NBChecker& ) const
      {
        applyLocal( entity );
      }

      /** \brief  apply local for all elements that do not need information from
       *          other processes (here all elements)
       */
      template <class NBChecker>
      void applyLocalInterior( const EntityType& entity, const NBChecker& ) const
      {
        applyLocal( entity );
      }

      /** \brief  apply local for all elements that need information from
       *          other processes (here no elements)
       */
      template <class NBChecker>
      void applyLocalProcessBoundary( const EntityType& entity, const NBChecker& ) const
      {
        DUNE_THROW(InvalidStateException,"DGInverseMassPass does not need a second phase for ThreadPass");
      }

    protected:
      void compute ( const TotalArgumentType &argument, DestinationType &destination ) const
      {
        massInverseMass_.compute( *(std::get< functionalPosition >( argument )), destination );
      }

    protected:
      MassInverseMassType& massInverseMass_ ;

      mutable const TotalArgumentType* arg_ ;
      mutable DestinationType*   dest_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_APPLYLOCALOPERATOR_HH
