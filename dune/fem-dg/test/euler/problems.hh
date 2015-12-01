#ifndef FEMDG_EULERPROBLEMS_HH
#define FEMDG_EULERPROBLEMS_HH

#include "problems/problems.hh"

namespace Dune
{
namespace Fem
{
  template< class GridImp >
  struct AnalyticalEulerProblemCreator
  {
    typedef ProblemBase< GridImp >                      ProblemInterfaceType;

    static ProblemInterfaceType* apply()
    {
      const std::string problemNames []
          = { "sod" , "withman", "withmansmooth", "smooth1d" , "ffs" , "diffraction" , "shockbubble", "p123" };

      const int problemNumber = Fem :: Parameter :: getEnum ( "problem" , problemNames );

      if( problemNames[ problemNumber ] == "sod" )
        return new U0Sod< GridImp > ( );
      else if( problemNames[ problemNumber ] == "smooth1d" )
        return new U0Smooth1D< GridImp > ();
      else if( problemNames[ problemNumber ] == "ffs" )
        return new U0FFS< GridImp > ();
      else if( problemNames[ problemNumber ] == "p123" )
        return new U0P123< GridImp >();
      std::cerr << "Error: Problem " << problemNames[ problemNumber ]
                << " not implemented." << std::endl;

      // choice of explicit or implicit ode solver
      return new U0Smooth1D< GridImp > ();
    }
  };

}
}

#endif
