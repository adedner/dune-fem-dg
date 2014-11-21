# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/fem-dg/operator/limiter/limiter.hh])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-fem
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-fem]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
# preset variable to path such that #include <dune/...> works
#AC_SUBST(AM_CPPFLAGS, '\${DUNE_CPPFLAGS} -I$(top_srcdir)')
#AC_SUBST(AM_LDFLAGS, '\${DUNE_LDFLAGS} \${DUNE_LIBS}')
AC_SUBST(AM_CPPFLAGS, '${DUNE_CPPFLAGS} -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '${DUNE_LDFLAGS} ${DUNE_LIBS}')

AC_CONFIG_FILES([
  Makefile
  m4/Makefile
  dune/Makefile
  dune/fem-dg/Makefile
  dune/fem-dg/pass/Makefile
  dune/fem-dg/operator/Makefile
  dune/fem-dg/operator/adaptation/Makefile
  dune/fem-dg/operator/limiter/Makefile
  dune/fem-dg/operator/dg/Makefile
  dune/fem-dg/operator/fluxes/Makefile
  dune/fem-dg/stepper/Makefile
  dune/fem-dg/solver/Makefile
  dune/fem-dg/main/Makefile
  dune/fem-dg/misc/Makefile
  dune/fem-dg/models/Makefile
  dune/fem-dg/test/Makefile
  dune/fem-dg/test/dataio/Makefile
  dune/fem-dg/test/advdiff/Makefile
  dune/fem-dg/test/euler/Makefile
  dune/fem-dg/test/poisson/Makefile
  dune/fem-dg/test/navierstokes/Makefile
  dune/fem-dg/test/poisson/Makefile
  dune/fem-dg/assemble/Makefile
  dune-fem-dg.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
