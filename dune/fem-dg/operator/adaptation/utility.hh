#ifndef DUNE_ADAPTATION_UTILITY_HH
#define DUNE_ADAPTATION_UTILITY_HH

// include restricion, prolongation and adaptation operator classes for discrete functions
#include <dune/fem/io/parameter.hh>
#include <dune/fem/gridpart/common/capabilities.hh>

namespace Dune
{
namespace Fem
{

  /**
   * \brief Parameter class describing adaptation strategies.
   *
   * \ingroup ParameterClass
   */
  struct AdaptationParameters
    : public Fem::LocalParameter< AdaptationParameters, AdaptationParameters >
  {
    typedef Fem::Parameter ParameterType;

    protected:
    const std::string keyPrefix_;
    int markStrategy_;
    const Dune::Fem::ParameterReader parameter_;

    public:
    //! constructor
    AdaptationParameters( const std::string keyPrefix = "fem.adaptation.",
                          const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : keyPrefix_( keyPrefix ),
        markStrategy_( getStrategy() ),
        parameter_(parameter)
    {}
    AdaptationParameters( const Dune::Fem::ParameterReader &parameter )
      : keyPrefix_( "fem.adaptation." ),
        markStrategy_( getStrategy() ),
        parameter_(parameter)
    {}

    ////TODO:
    //AdaptationParameters( std::initializer_list<std::pair<std::string,std::string> > paramList, const std::string keyPrefix = "fem.adaptation."  )
    //  : keyPrefix_( keyPrefix ),
    //    markStrategy_( getStrategy() )
    //{
    //  for( auto it : paramList )
    //    ParameterType::append( it->first, it->second );
    //}

    int getStrategy () const
    {
      const std::string names[] = { "shockind", "apost", "grad" };
      // default value is gradient
      return ParameterType::getEnum( keyPrefix_ + "markingStrategy", names, 1 );
    }

    //! simulation end time
    virtual double endTime () const
    {
      return ParameterType::getValue< double >( "femdg.stepper.endtime" );
    }

    //! retujrn refinement tolerance
    virtual double refinementTolerance () const
    {
      return ParameterType::getValue< double >( keyPrefix_ + "refineTolerance", 0.05 );
    }

    //! return percentage of refinement tolerance used for coarsening tolerance
    virtual double coarsenPercentage () const
    {
      return ParameterType::getValue< double >( keyPrefix_ + "coarsenPercent", 0.1 );
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
             ParameterType::getValue< int >( keyPrefix_ + "finestLevel", 0 );
    }

    //! return minimal level achieved by refinement
    virtual int coarsestLevel ( const int refineStepsForHalf = 1 ) const
    {
      return refineStepsForHalf *
             ParameterType::getValue< int >( keyPrefix_ + "coarsestLevel", 0 );
    }

    //! return depth for refining neighbors of a cell marked for refinement
    virtual int neighborRefLevel () const
    {
      return ParameterType::getValue< int >( keyPrefix_ + "grad.neighborRefLevel", 1 );
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
      return ParameterType::getValue<int>( keyPrefix_ + "adaptcount", 1 );
    }

    virtual int method () const
    {
      const std::string names[] = { "none", "generic", "callback" };
      // default value is generic
      return ParameterType::getEnum( keyPrefix_ + "method", names, 1 );
    }

    virtual double theta() const
    {
      return ParameterType::getValue<double>( keyPrefix_ + "theta", 0.5);
    }

    virtual bool maximumStrategy() const
    {
      return ParameterType::getValue<bool>( keyPrefix_ + "maximumstrategy", true );
    }

    virtual int padaptiveMethod () const
    {
      const std::string names[] = { "none", "coeffroot", "prior2p" };
      // default value is generic
      return ParameterType::getEnum( keyPrefix_ + "padaptive.method", names, 1 );
    }

    virtual bool adaptive () const
    {
      return (method() != 0 && adaptCount() != 0);
    }

    //! return true if verbosity mode is enabled
    virtual bool verbose () const { return ParameterType::getValue< bool >( keyPrefix_ + "verbose", false ); }
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
        return findCoarsestVolume( entity.father(), hasGridHierarchy );
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

      const bool hasGridHierarchy = Fem :: GridPartCapabilities :: hasGrid< GridPart > :: v;

      double weight = Dune::DGFGridInfo<GridType>::refineWeight();
      // if weight is not set, use 1/(2^d)
      if( weight < 0 )
        weight = 1.0/std::pow( 2.0, double( GridType :: dimension) );

      double minVolume[ 1 ] = { std::numeric_limits< double > ::max() };
      double maxVolume[ 1 ] = { std::numeric_limits< double > ::min() };

      // if grid is not empty compute smallest and biggest volume

      for( const auto& en : elements( gridPart ) )
      {
        const double volume = findCoarsestVolume( en, hasGridHierarchy );
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

} // end namespace
} // end namespace Dune
#endif
