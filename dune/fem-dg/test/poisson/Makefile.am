# install these headers 
poissondir=$(includedir)/test/poisson
poisson_HEADERS = models.hh  problemcreator.hh  problems.hh

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

BASEDIR=../../main

#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE = YASPGRID
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
#GRIDTYPE=CARTESIANGRID_ALUGRID_CUBE
#GRIDTYPE = SPGRID
GRIDTYPE=ALUGRID_CUBE
POLORDER = 2
GRIDDIM  = 2
FLUX = 1    # LLF 1, HLL 2
WANT_ISTL=1

#USE_OMP=-fopenmp

EXTRA_PROGRAMS = poisson_all # quasiadvdiff quasiadvdiff12 plaplace plaplace12
check_PROGRAMS = poisson # quasiadvdiff quasiadvdiff12 plaplace plaplace12

poisson_all_SOURCES = main.cc
poisson_SOURCES = main.cc

AM_CPPFLAGS = $(USE_OMP)  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DDGSCHEME -DONB \
              -DFLUX=$(FLUX) -DPRIMALDG -DUSE_ASSERT_THROW -DWANT_ISTL=$(WANT_ISTL)
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)

poisson_all_CPPFLAGS = $(AM_CPPFLAGS)
poisson_CPPFLAGS = $(poisson_all_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = parameter

CLEANFILES=manager.*.log eoc_*.tex

PROG=poisson
# codegeneration 
generatecode:
	$(MAKE) -i clean
	$(MAKE) CXXFLAGS="-O0 -Wall -DNDEBUG -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
	./$(PROG) femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none codegenparameter

# compile fast code  
compilecode:
	$(MAKE) clean 
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)

codegen:
	$(MAKE) generatecode
	$(MAKE) compilecode

include $(top_srcdir)/am/global-rules
