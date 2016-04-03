#ifndef FEMDG_ADVDIFFPROBLEMS_HH
#define FEMDG_ADVDIFFPROBLEMS_HH

#include "problems/problem.hh"
#include "problems/problemQuasiHeatEqn.hh"
#include "problems/pulse.hh"
#include "problems/sin.hh"
#include "problems/deformationalflow.hh"

namespace Dune
{
namespace Fem
{
  template< class FunctionSpaceImp, class GridImp >
  struct AnalyticalAdvDiffProblemCreator
  {
    typedef EvolutionProblemInterface< FunctionSpaceImp > ProblemInterfaceType;
    static ProblemInterfaceType* apply()
    {
      // choice of explicit or implicit ode solver
      static const std::string probString[]  = { "heat" ,"sin", "quasi", "pulse" };
      const int probNr = Fem::Parameter::getEnum( "problem", probString, 0 );
      if( probNr == 0 )
        return new U0< GridImp, FunctionSpaceImp::dimRange > ();
      else if ( probNr == 1 )
        return new U0Sin< GridImp, FunctionSpaceImp::dimRange > ();
      else if ( probNr == 2 )
        return new QuasiHeatEqnSolution< GridImp, FunctionSpaceImp::dimRange > ();
      else if ( probNr == 3 )
        return new Pulse< GridImp, FunctionSpaceImp::dimRange > ();
      else
      {
        abort();
        return 0;
      }
    }
  };

}
}
#endif
