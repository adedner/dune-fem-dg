# install these headers 
stokesdir=$(includedir)/test/stokes
stokes_HEADERS = problem.hh corner.hh problemcreator.hh stokesassembler.hh

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

EXTRA_PROGRAMS = stokes_all 
check_PROGRAMS = stokes 

stokes_all_SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_1.cc $(BASEDIR)/main_2.cc \
                  $(BASEDIR)/main_0.cc $(BASEDIR)/main_3.cc
stokes_SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_pol.cc

AM_CPPFLAGS = $(USE_OMP)  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DDGSCHEME -DONB \
              -DFLUX=$(FLUX) -DPRIMALDG -DUSE_ASSERT_THROW -DWANT_ISTL=$(WANT_ISTL) -DSTOKES
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)

stokes_CPPFLAGS = $(AM_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = parameter

CLEANFILES=manager.*.log eoc_*.tex

PROG=stokes
# codegeneration 
generatecode:
	$(MAKE) -i clean
	$(MAKE) CXXFLAGS="-O0 -Wall -DNDEBUG -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
	./$(PROG) femhowto.maximaltimesteps:1 fem.io.outputformat:none codegenparameter

# compile fast code  
compilecode:
	$(MAKE) clean 
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)

codegen:
	$(MAKE) generatecode
	$(MAKE) compilecode

include $(top_srcdir)/am/global-rules
