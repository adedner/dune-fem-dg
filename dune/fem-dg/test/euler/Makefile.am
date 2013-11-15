# install these headers 
eulerdir=$(includedir)/dune/fem-dg/test/euler
euler_HEADERS= eulermodel.hh  problems.hh problemcreator.hh \
chorjo.hh chorjo.cc 							 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS)


BASEDIR=../../main

# helper variables
SOURCEBASE= $(BASEDIR)/main.cc chorjo.cc
SOURCES = $(SOURCEBASE) $(BASEDIR)/main_0.cc $(BASEDIR)/main_1.cc \
					$(BASEDIR)/main_2.cc $(BASEDIR)/main_3.cc 
SOURCES12 = $(SOURCEBASE) $(BASEDIR)/main_1.cc $(BASEDIR)/main_2.cc
SOURCEONEP = $(SOURCEBASE) $(BASEDIR)/main_pol.cc 

# INFO FLUX
#  set FLUX to 1 for LLF
#  set FLUX to 2 for HLL (not tested!)
# INFO TRACK LIFTING:
#  define LOCALDEBUG to calculate \sum_e\int_Omega(r_e*l_e) and 
#   \sum_e\int_Omega(r_e*l_e). They will be output to std::cout from the Stepper
# INFO LIMITER
#  define LIMITER to apply limiting of the numerical solution (not tested!)
# INFO WELL-BALANCING
#   define WELLBALANCE to apply well-balance of the numerical solution (not tested!)
# INFO DIFFFLUXTYPE
#   define PRIMALDG to use the BR2, IP, CDG, CDG2, BO methods
#   define DUALDG to use the LDG method

#GRIDTYPE = ALUGRID_CONFORM
#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE = ALBERTAGRID
#GRIDTYPE = YASPGRID
#GRIDTYPE = SPGRID
GRIDTYPE = ALUGRID_CUBE
GRIDDIM = 2
FLUX=1       # see comment below
DIFFFLUXTYPE=PRIMALDG
POLORDER=2

#USE_OMP=-fopenmp
#USE_OMP=-DUSE_PTHREADS=1

AM_CPPFLAGS = $(USE_OMP) -I../../problems/coldbubble  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) \
              -DFLUX=$(FLUX) -D$(DIFFFLUXTYPE) -DWANT_CACHED_COMM_MANAGER=0
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)

check_PROGRAMS = euler euler12 euleronep 

euler_SOURCES = $(SOURCES)
euler12_SOURCES = $(SOURCES12)
euleronep_SOURCES = $(SOURCEONEP)

euler_CPPFLAGS = $(AM_CPPFLAGS) -DLIMITER
euler12_CPPFLAGS = $(euler_CPPFLAGS) -DONLY_P1_P2
euleronep_CPPFLAGS = $(euler_CPPFLAGS) -DONLY_ONE_P -DPOLORDER=$(POLORDER)

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"

EXTRA_DIST = param* 

CLEANFILES=manager.*.log 

PROG=euler
# codegeneration 
generatecode:
	$(MAKE) -i clean
	$(MAKE) CXXFLAGS="-O0 -Wall -DNDEBUG -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
	./$(PROG) femhowto.maximaltimesteps:1 fem.io.outputformat:none parameter 

compilecode:
	$(MAKE) clean 
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)

codegen:
	$(MAKE) generatecode
	$(MAKE) compilecode

include $(top_srcdir)/am/global-rules
