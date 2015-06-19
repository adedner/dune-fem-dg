# install these headers 
advdiffdir=$(includedir)/test/advdiff
advdiff_HEADERS = models.hh  problemcreator.hh steppertraits.hh

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

BASEDIR=../../main

#GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = YASPGRID
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
#GRIDTYPE=CARTESIANGRID_ALUGRID_CUBE
#GRIDTYPE = SPGRID
POLORDER = 0
GRIDDIM  = 2
DIMRANGE = 4
FLUX = 1    # LLF 1, HLL 2

#USE_OMP=-fopenmp

# INFO SPACE OPERATOR:
#     1. define PRIMALDG for IP, BR2, CDG, CDG2
#     2. define DUALDG for LDG
# INFO TRACK LIFTING:
#     1. define LOCALDEBUG to calculate \sum_e\int_Omega(r_e*l_e) and 
#        \sum_e\int_Omega(r_e*l_e). They will be output to std::cout from the Stepper
# INFO TESTOPERATOR:
#     1. define TESTOPERATOR for linear advdiff equation to output various
#        information on space operator
EXTRA_PROGRAMS = advdiff pulse # quasiadvdiff quasiadvdiff12 plaplace plaplace12
check_PROGRAMS = advdiffonep pulseonep # quasiadvdiff quasiadvdiff12 plaplace plaplace12

advdiff_SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_1.cc $(BASEDIR)/main_2.cc \
                  $(BASEDIR)/main_0.cc $(BASEDIR)/main_3.cc
advdiffonep_SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_pol.cc

pulse_SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_1.cc $(BASEDIR)/main_2.cc \
                $(BASEDIR)/main_0.cc $(BASEDIR)/main_3.cc
pulseonep_SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_pol.cc

AM_CPPFLAGS = $(USE_OMP) -I../../problems/advdiff  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) \
              -DDIMRANGE=$(DIMRANGE) -DFLUX=$(FLUX) -DPRIMALDG -DUSE_ASSERT_THROW
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)

advdiff_CPPFLAGS = $(AM_CPPFLAGS)
advdiffonep_CPPFLAGS = $(advdiff_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)

pulse_CPPFLAGS = $(AM_CPPFLAGS)
pulseonep_CPPFLAGS = $(advdiff_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = parameter

CLEANFILES=manager.*.log eoc_*.tex

PROG=advdiffonep
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
