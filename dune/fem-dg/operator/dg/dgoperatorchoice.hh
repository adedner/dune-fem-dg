#ifndef NS_DGOPCHOICE_HH
#define NS_DGOPCHOICE_HH

#ifdef HEADERCHECK
#define PRIMALDG
#endif

#ifdef PRIMALDG
#warning "Primal formulation based operators (CDG2,CDG,BR2,IP,NIPG,BO) available"
#include "primaloperator.hh"

#elif (defined DUALDG) || (defined FLUXDG)
#warning "Flux formulation based operator (LDG) available"
#include "fluxoperator.hh"

#else
#error "Please define either PRIMALDG or FLUXDG in AM_CXXFLAGS"
#endif

#endif
