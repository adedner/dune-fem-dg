#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE = CARTESIANGRID_ALUGRID_CUBE
GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = SPGRID
#GRIDTYPE = YASPGRID
GRIDDIM = 3
POLORDER = 1
DIMRANGE=1
WANT_SIONLIB=0

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DWANT_SIONLIB=$(WANT_SIONLIB) \
              -DPOLORDER=$(POLORDER) \
              -DDIMRANGE=$(DIMRANGE) # -DENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) 

EXTRA_DIST = parameter

# install these headers 
testdir=$(includedir)/dune/fem-dg/test
test_HEADERS= checkpointing.hh

check_PROGRAMS = checkpointing dataconvert 
checkpointing_SOURCES = checkpointing.cc
dataconvert_SOURCES = dataconvert.cc
