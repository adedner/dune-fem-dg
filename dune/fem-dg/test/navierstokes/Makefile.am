# install these headers 
nseqdir=$(includedir)/dune/fem-dg/test/navierstokes
nseq_HEADERS= models.hh problemcreator.hh thermodynamics.hh

#USE_OMP=-fopenmp
#USE_OMP=-DUSE_PTHREADS=1

LDADD=$(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

SOURCES = main.cc

#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE=PARALLELGRID_ALUGRID_CUBE
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
# GRIDTYPE=SPGRID
#GRIDTYPE=SPGRID_COUNT_FLOPS
#GRIDTYPE=CARTESIANGRID_ALUGRID_CUBE
# GRIDTYPE = ALUGRID_CUBE
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
#        \sum_e\int_Omega(r_e*l_e). They will be output to std::cout from the Algorithm
# INFO LIMITER
#     1. define LIMITER to limit the advection part of the solution (not checked)
AM_CPPFLAGS = $(USE_OMP) -I../../problems/nseq  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) \
              -DFLUX=$(FLUX) -D$(DIFFFLUX) # -DCOUNT_FLOPS # -DLEGENDREBASIS # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP) 

EXTRA_PROGRAMS = navierstokes_all
check_PROGRAMS = navierstokes

nodist_navierstokes_all_SOURCES = $(SOURCES) 
nodist_navierstokes_SOURCES = $(SOURCES) 

navierstokes_all_CPPFLAGS = $(AM_CPPFLAGS)
navierstokes_CPPFLAGS = $(navierstokes_all_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = parameter

CLEANFILES = manager.*.log 

PROG=navierstokes
# codegeneration 
generatecode:
	$(MAKE) -i clean
	$(MAKE) CODEGEN_OBJFILE=  CXXFLAGS="-g -Wall -Wfatal-errors -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
	$(MAKE) generate

generate:
	./$(PROG) fem.eoc.steps:1 femdg.stepper.maximaltimesteps:1 fem.io.outputformat:none fem.ode.order:1 parameter 

compilecode:
	$(MAKE) clean 
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)

codegen:
	$(MAKE) generatecode
	$(MAKE) compilecode

include $(top_srcdir)/am/global-rules
