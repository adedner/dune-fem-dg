# install these headers
set(eulerdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/test/euler)
set(euler_HEADERS  eulermodel.hh  problems.hh problemcreator.hh
chorjo.hh chorjo.cc passtraits.hh)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${DUNEMPILDFLAGS} ${DUNEMPILIBS})
set(BASEDIR ../../main)
# helper variables
set(SOURCEBASE  ${BASEDIR}/main.cc chorjo.cc)
set(SOURCES  ${SOURCEBASE} ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc
                                        ${BASEDIR}/main_2.cc ${BASEDIR}/main_3.cc)
set(SOURCES12  ${SOURCEBASE} ${BASEDIR}/main_1.cc ${BASEDIR}/main_2.cc)
set(SOURCEONEP  ${SOURCEBASE} ${BASEDIR}/main_pol.cc)
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
set(GRIDTYPE  ALUGRID_CUBE)
set(GRIDDIM  2)
install(FILES ${euler_HEADERS} DESTINATION ${eulerdir})