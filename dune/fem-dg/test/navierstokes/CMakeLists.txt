# install these headers
set(nseqdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/test/navierstokes)
set(nseq_HEADERS  ns_model.hh problemcreator.hh nswaves.hh
ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh)
#USE_OMP=-fopenmp
#USE_OMP=-DUSE_PTHREADS=1
set(LDADD ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${DUNEMPILDFLAGS} ${DUNEMPILIBS})
set(BASEDIR ../../main)
set(SOURCESONEP  ${BASEDIR}/main.cc ${BASEDIR}/main_pol.cc)
set(SOURCES  ${BASEDIR}/main.cc ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc
          ${BASEDIR}/main_2.cc ${BASEDIR}/main_3.cc ${BASEDIR}/main_4.cc)
set(SOURCES12  ${BASEDIR}/main.cc ${BASEDIR}/main_1.cc ${BASEDIR}/main_2.cc)
#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE=PARALLELGRID_ALUGRID_CUBE
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
# GRIDTYPE=SPGRID
#GRIDTYPE=CARTESIANGRID_ALUGRID_CUBE
set(GRIDTYPE  ALUGRID_CUBE)
#GRIDTYPE=YASPGRID
#GRIDTYPE = YASPGRID
set(GRIDDIM 3)
set(POLORDER 2)
install(FILES ${nseq_HEADERS} DESTINATION ${nseqdir})