# install these headers 
advdiffdir=$(includedir)/test/advdiff
advdiff_HEADERS = models.hh  problemcreator.hh

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

#GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = YASPGRID
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
#GRIDTYPE=CARTESIANGRID_ALUGRID_CUBE
#GRIDTYPE = SPGRID
POLORDER = 2
GRIDDIM  = 2
DIMRANGE = 1
FLUX = 1    # LLF 1, HLL 2

#USE_OMP=-fopenmp

# INFO SPACE OPERATOR:
#     1. define PRIMALDG for IP, BR2, CDG, CDG2
#     2. define DUALDG for LDG
# INFO TRACK LIFTING:
#     1. define LOCALDEBUG to calculate \sum_e\int_Omega(r_e*l_e) and 
#        \sum_e\int_Omega(r_e*l_e). They will be output to std::cout from the Algorithm
# INFO TESTOPERATOR:
#     1. define TESTOPERATOR for linear advdiff equation to output various
#        information on space operator
EXTRA_PROGRAMS = advdiff_all
check_PROGRAMS = advdiff

advdiff_SOURCES = main.cc 
advdiff_all_SOURCES = main.cc

AM_CPPFLAGS = $(USE_OMP) -I../../problems/advdiff  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) \
              -DDIMRANGE=$(DIMRANGE) -DFLUX=$(FLUX) -DPRIMALDG
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)

advdiff_all_CPPFLAGS = $(AM_CPPFLAGS)
advdiff_CPPFLAGS = $(advdiff_all_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = parameter

CLEANFILES=manager.*.log eoc_*.tex

PROG=advdiff
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
