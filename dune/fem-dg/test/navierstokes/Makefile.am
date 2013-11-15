# install these headers 
nseqdir=$(includedir)/src/problems/nseq
nseq_HEADERS= ns_model.hh problemcreator.hh nswaves.hh \
ns_model_spec.hh constant.hh problemtype.hh nssmooth.hh thermodynamics.hh passtraits.hh

#USE_OMP=-fopenmp
#USE_OMP=-DUSE_PTHREADS=1

LDADD=$(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)

BASEDIR=$(DUNE_FEM_DG_ROOT)/dune/fem-dg/main

SOURCES = $(BASEDIR)/main.cc $(BASEDIR)/main_0.cc $(BASEDIR)/main_1.cc \
						$(BASEDIR)/main_2.cc $(BASEDIR)/main_3.cc $(BASEDIR)/main_4.cc
SOURCES12 = $(BASEDIR)/main.cc $(BASEDIR)/main_1.cc $(BASEDIR)/main_2.cc
SOURCESONEP = $(BASEDIR)/main.cc $(BASEDIR)/main_pol.cc 

#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE=PARALLELGRID_ALUGRID_CUBE
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
GRIDTYPE=SPGRID
#GRIDTYPE=CARTESIANGRID_ALUGRID_CUBE
#GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE=YASPGRID
#GRIDTYPE = YASPGRID
GRIDDIM=2
POLORDER=2

PROBLEM = 2 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

DIFFFLUX=PRIMALDG
#DIFFFLUX=FLUXDG

# INFO SPACE OPERATOR:
# 		1. define PRIMALDG for various space operators in primal formulation
# 	  2. define DUALDG for the LDG space operator in dual formulation
# INFO TRACK LIFTING:
# 		1. define LOCALDEBUG to calculate \sum_e\int_Omega(r_e*l_e) and 
# 			 \sum_e\int_Omega(r_e*l_e). They will be output to std::cout from the Stepper
# INFO LIMITER
# 		1. define LIMITER to limit the advection part of the solution (not checked)
AM_CPPFLAGS = $(USE_OMP) -I../../problems/nseq  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) \
              -DFLUX=$(FLUX) -D$(DIFFFLUX) -DLAGRANGEBASIS # -DLEGENDREBASIS # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP) 


check_PROGRAMS = constant constantonep constant12 \
								 nseq nseqonep nseq12
nodist_constant_SOURCES = $(SOURCES)
nodist_constantonep_SOURCES = $(SOURCESONEP)
nodist_constant12_SOURCES = $(SOURCES12) 
nodist_nseq_SOURCES = $(SOURCES) 
nodist_nseqonep_SOURCES = $(SOURCESONEP) 
nodist_nseq12_SOURCES = $(SOURCES12) 

constant_CPPFLAGS = $(AM_CPPFLAGS)
constantonep_CPPFLAGS = $(constant_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)
constant12_CPPFLAGS = $(constant_CPPFLAGS) -DONLY_P1_P2
nseq_CPPFLAGS = $(AM_CPPFLAGS)
nseqonep_CPPFLAGS = $(nseq_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)
nseq12_CPPFLAGS = $(nseq_CPPFLAGS) -DONLY_P1_P2

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=nseqonep
# codegeneration 
generatecode:
	$(MAKE) -i clean
	$(MAKE) CODEGEN_OBJFILE=  CXXFLAGS="-g -Wall -Wfatal-errors -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
	$(MAKE) generate

generate:
	./$(PROG) femhowto.eocSteps:1 femhowto.maximaltimesteps:1 fem.io.outputformat:none fem.ode.order:1 paramFiles/paramNSWaves 

compilecode:
	$(MAKE) clean 
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)

codegen:
	$(MAKE) generatecode
	$(MAKE) compilecode

include $(top_srcdir)/am/global-rules
