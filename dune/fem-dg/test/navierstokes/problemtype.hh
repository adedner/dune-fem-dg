#ifndef FEMHOWTO_NSEQ_PROBLEMTYPE_HH
#define FEMHOWTO_NSEQ_PROBLEMTYPE_HH

#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>

///////////////////////////////////////
// AVAILABLE PROBLEMS
///////////////////////////////////////
#if PROBLEM==2
   // a nonstationary exact solution 
   // for Navier-Stokes equations
   #include "nswaves.hh"
   typedef NSWaves< GridSelector :: GridType >  NSProblemType;
   #define PROBLEM_HAS_SOLUTION

#elif PROBLEM==3
   // a nonstationary exact solution 
   // for Navier-Stokes equations
   #include "nssmooth.hh"
   typedef NSSmoothSolution< GridSelector :: GridType >  NSProblemType;
   #define PROBLEM_HAS_SOLUTION

#else
   #error "No valid problem number specified"
#endif

#endif
