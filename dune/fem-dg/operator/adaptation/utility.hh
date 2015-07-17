#ifndef DUNE_ADAPTATION_UTILITY_HH
#define DUNE_ADAPTATION_UTILITY_HH

// include restricion, prolongation and adaptation operator classes for discrete functions
#include <dune/fem/io/parameter.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/misc/compatibility.hh>

namespace Dune
{
  struct AdaptationParameters
    : public Fem::LocalParameter< AdaptationParameters, AdaptationParameters >
  {
    protected:
    const std::string keyPrefix_;
    int markStrategy_;



    public:
    //! constructor
    AdaptationParameters( const std::string keyPrefix = "fem.adaptation." )
      : keyPrefix_( keyPrefix ),
        markStrategy_( getStrategy() )
    {}

    int getStrategy () const
    {
      const std::string names[] = { "shockind", "apost", "grad" };
      // default value is gradient
      return Fem::Parameter::getEnum( keyPrefix_ + "markingStrategy", names, 2 );
    }

    //! simulation end time
    virtual double endTime () const
    {
      return Fem::Parameter::getValue< double >( "femdg.stepper.endtime" );
    }

    //! retujrn refinement tolerance
    virtual double refinementTolerance () const
    {
      return Fem::Parameter::getValue< double >( keyPrefix_ + "refineTolerance" );
    }

    //! return percentage of refinement tolerance used for coarsening tolerance
    virtual double coarsenPercentage () const
    {
      return Fem::Parameter::getValue< double >( keyPrefix_ + "coarsenPercent", 0.1 );
    }

    //! return product of refinementTolerance and coarsenPercentage
    virtual double coarsenTolerance () const
    {
      return refinementTolerance() * coarsenPercentage();
    }

    //! return maximal level achieved by refinement
    virtual int finestLevel ( const int refineStepsForHalf = 1 ) const
    {
      return refineStepsForHalf *
             Fem::Parameter::getValue< int >( keyPrefix_ + "finestLevel" );
    }

    //! return minimal level achieved by refinement
    virtual int coarsestLevel ( const int refineStepsForHalf = 1 ) const
    {
      return refineStepsForHalf *
             Fem::Parameter::getValue< int >( keyPrefix_ + "coarsestLevel", 0 );
    }

    //! return depth for refining neighbors of a cell marked for refinement
    virtual int neighborRefLevel () const
    {
      return Fem::Parameter::getValue< int >( keyPrefix_ + "grad.neighborRefLevel", 1 );
    }

    //! return true if marking strategy is based on shock indicator
    virtual bool shockIndicator () const
    {
      return markStrategy_ == 0;
    }

    //! return true if marking strategy is based on shock indicator
    virtual bool gradientBasedIndicator () const
    {
      return markStrategy_ == 2;
    }

    //! return true if aposteriori indicator is enabled
    virtual bool aposterioriIndicator () const
    {
      return markStrategy_ == 1;
    }

    virtual int adaptCount() const
    {
      return Fem::Parameter::getValue<int>( keyPrefix_ + "adaptcount", 1 );
    }

    virtual bool method () const
    {
      const std::string names[] = { "none", "generic", "callback" };
      // default value is generic
      return Fem::Parameter::getEnum( keyPrefix_ + "method", names, 1 );
    }

    virtual double theta() const
    {
      return Fem::Parameter::getValue<double>( keyPrefix_ + "theta", 0.5);
    }

    virtual bool maximumStrategy() const
    {
      return Fem::Parameter::getValue<bool>( keyPrefix_ + "maximumstrategy", true );
    }

    virtual bool adaptive () const
    {
      return method() != 0;
    }

    //! return true if verbosity mode is enabled
    virtual bool verbose () const { return Fem::Parameter::getValue< bool >( keyPrefix_ + "verbose", false ); }
  };

  class ComputeMinMaxVolume
  {
  protected:
    typedef std::pair< double, double > VolumePairType;

    template <class Entity>
    double findCoarsestVolume( const Entity& entity, const bool hasGridHierarchy ) const
    {
      // go to father, if possible
      // otherwise write min and max volume on backup/restore
      if( hasGridHierarchy && entity.level() > 0 )
      {
        return findCoarsestVolume( Fem::make_entity( entity.father() ), hasGridHierarchy );
      }
      else  // return entity's volume
        return entity.geometry().volume();
    }

    template <class GridPart>
    VolumePairType computeMinMaxVolume( const GridPart& gridPart,
                                        const int coarsestLevel,
                                        const int finestLevel )
    {
      typedef typename GridPart :: GridType GridType ;
      typedef typename GridPart :: template Codim< 0 > :: IteratorType IteratorType;

      const bool hasGridHierarchy = Fem :: GridPartCapabilities :: hasGrid< GridPart > :: v;

      double weight = Dune::DGFGridInfo<GridType>::refineWeight();
      // if weight is not set, use 1/(2^d)
      if( weight < 0 )
        weight = 1.0/std::pow( 2.0, double( GridType :: dimension) );

      double minVolume[ 1 ] = { std::numeric_limits< double > ::max() };
      double maxVolume[ 1 ] = { std::numeric_limits< double > ::min() };

      // if grid is not empty compute smallest and biggest volume
      const IteratorType end = gridPart.template end< 0 >();
      for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
      {
        const double volume = findCoarsestVolume( *it, hasGridHierarchy );
        minVolume[ 0 ] = std::min( minVolume[ 0 ], volume );
        maxVolume[ 0 ] = std::max( maxVolume[ 0 ], volume );
      }

      // compute global min and max of values
      computeGlobalMinMax( gridPart.grid().comm(), 1, maxVolume, minVolume );

      double finestWeight   = std::pow( weight, double(finestLevel) );
      double coarsestWeight = std::pow( weight, double(coarsestLevel) );
      // set local variables
      return VolumePairType( double(maxVolume[ 0 ] * coarsestWeight),
                             double(minVolume[ 0 ] * finestWeight) );
    }

    template <class CommunicatorType, class Vector>
    void computeGlobalMinMax( const CommunicatorType& comm, const int size, Vector& max, Vector& min ) const
    {
      //std::vector< double > buffer( 2 * size, 0.0 );
      double buffer[ 2 * size ];
      // store max
      for( int i=0; i<size; ++i )
        buffer[ i ] = max[ i ];

      // store 1/min
      const double eps = std::numeric_limits< double > :: epsilon ();
      for( int i=0, ib=size; i<size; ++i, ++ib )
        buffer[ ib ] = (std::abs( min[ i ] ) > eps) ? 1.0/min[ i ] : 0;

      // compute global maximum
      comm.max( &buffer[ 0 ], 2*size );

      // store max again
      for( int i=0; i<size; ++i )
        max[ i ] = buffer[ i ];
      // store 1/min again
      for( int i=0, ib=size; i<size; ++i, ++ib )
        min[ i ] = (std::abs( buffer[ ib ] ) > eps) ? 1.0/buffer[ ib ] : 0;
    }

    //! constructor computing coarsest and finest volume
    template <class GridPart>
    ComputeMinMaxVolume( const GridPart& gridPart,
                         const int coarsestLevel,
                         const int finestLevel )
      : volumes_( computeMinMaxVolume( gridPart, coarsestLevel, finestLevel ) )
    {
    }

    const VolumePairType volumes_;

  public:
    double coarsestVolume() const { return volumes_.first; }
    double finestVolume() const { return volumes_.second; }
  };

} // end namespace Dune
#endif
