#ifndef FEMDG_ALGORITHM_GRIDINITIALIZER_HH
#define FEMDG_ALGORITHM_GRIDINITIALIZER_HH

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/common/timer.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem-dg/misc/parameterkey.hh>

#include <dune/fem-dg/operator/adaptation/utility.hh>
#include <dune/fem-dg/algorithm/caller/checkpoint.hh>
#include <dune/fem-dg/algorithm/caller/adapt.hh>

namespace Dune
{
namespace Fem
{

  template< class GridImp,
            class CheckPointCallerImp = GridCheckPointCaller< GridImp > >
  class DefaultGridInitializer
  {
    typedef CheckPointCallerImp  CheckPointCallerType;
    public:

    static Dune::GridPtr< GridImp > initialize( const std::string name = "" )
    {
    #ifdef ALUGRID_CONSTRUCTION_WITH_STREAMS
      // precision of out streams (here ALUGrid backup streams)
      const int precision = Dune::Fem::Parameter::getValue< int > ("fem.io.precision", 16);
      ALU3DSPACE ALUGridExternalParameters :: precision () = precision;
    #endif

      // grid pointer object
      Dune::GridPtr< GridImp > gridptr;

      if( CheckPointCallerType::checkPointExists("") )
        gridptr = CheckPointCallerType::restoreGrid("");
      else  // normal new start
      {
        // ----- read in runtime parameters ------
        const std::string filekey = Dune::Fem::IOInterface::defaultGridKey( GridImp::dimension );
        const std::string filename = Dune::Fem::Parameter::commonInputPath() + "/" + Dune::Fem::Parameter::getValue< std::string >(  filekey );

        // initialize grid with given macro file name
        gridptr = Dune::GridPtr< GridImp >( filename );
        Dune::Fem::Parameter::appendDGF( filename );

        // load balance grid in case of parallel runs
        gridptr->loadBalance();

        // and refine the grid globally
        AdaptationParameters adaptParam( ParameterKey::generate( "", "fem.adaptation." ) );
        const int startLevel = adaptParam.coarsestLevel();
        for(int level=0; level < startLevel ; ++level)
          Dune::Fem::GlobalRefine::apply(*gridptr, 1 );
      }

      return gridptr;
    }

  };

}
}

#endif
