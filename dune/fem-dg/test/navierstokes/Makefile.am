# install these headers 
nseqdir=$(includedir)/dune/fem-dg/test/navierstokes
nseq_HEADERS= ns_model.hh problemcreator.hh nswaves.hh \
ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh

#USE_OMP=-fopenmp
#USE_OMP=-DUSE_PTHREADS=1

LDADD=$(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

BASEDIR=../../main

SOURCESONEP = $(BASEDIR)/main.cc $(BASEDIR)/main_pol.cc 
SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_0.cc $(BASEDIR)/main_1.cc \
          $(BASEDIR)/main_2.cc $(BASEDIR)/main_3.cc $(BASEDIR)/main_4.cc
SOURCES12 = $(BASEDIR)/main.cc $(BASEDIR)/main_1.cc $(BASEDIR)/main_2.cc

#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE=PARALLELGRID_ALUGRID_CUBE
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
#GRIDTYPE=SPGRID_COUNT_FLOPS
#GRIDTYPE=SPGRID
#GRIDTYPE=CARTESIANGRID_ALUGRID_CUBE
#GRIDTYPE = ALUGRID_CUBE
GRIDTYPE=YASPGRID
#GRIDTYPE = YASPGRID
GRIDDIM=2
POLORDER=2

PROBLEM = 2 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

DIFFFLUX=PRIMALDG
#DIFFFLUX=FLUXDG

# INFO SPACE OPERATOR:
#     1. define PRIMALDG for various space operators in primal formulation
#     2. define DUALDG for the LDG space operator in dual formulation
# INFO TRACK LIFTING:
#     1. define LOCALDEBUG to calculate \sum_e\int_Omega(r_e*l_e) and 
#        \sum_e\int_Omega(r_e*l_e). They will be output to std::cout from the Stepper
# INFO LIMITER
#     1. define LIMITER to limit the advection part of the solution (not checked)
AM_CPPFLAGS = $(USE_OMP) -I../../problems/nseq  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) \
              -DFLUX=$(FLUX) -D$(DIFFFLUX) -DLEGENDREBASIS # -DCOUNT_FLOPS # -DLEGENDREBASIS # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)

EXTRA_PROGRAMS = nseq nseq12
check_PROGRAMS = nseqonep

nodist_nseq_SOURCES = $(SOURCES) 
nodist_nseqonep_SOURCES = $(SOURCESONEP) 
nodist_nseq12_SOURCES = $(SOURCES12) 

nseq_CPPFLAGS = $(AM_CPPFLAGS)
nseqonep_CPPFLAGS = $(nseq_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)
nseq12_CPPFLAGS = $(nseq_CPPFLAGS) -DONLY_P1_P2

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = parameter

CLEANFILES = manager.*.log 

PROG=nseqonep
# width of SIMD unit (usually 2 or 4)
SIMDWIDTH=4

# codegeneration 
generatecode:
	$(MAKE) -i clean
	$(MAKE) CXXFLAGS="-g -Wall -Wfatal-errors -DBASEFUNCTIONSET_CODEGEN_GENERATE -DCODEGEN_SIMD_WIDTH=$(SIMDWIDTH)" $(PROG)
	$(MAKE) generate

generate:
	./$(PROG) fem.eoc.steps:1 fem.io.outputformat:none fem.ode.order:1 parameter 

compilecode:
	$(MAKE) clean 
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)

codegen:
	$(MAKE) generatecode
	$(MAKE) compilecode

include $(top_srcdir)/am/global-rules
