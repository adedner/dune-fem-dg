# install these headers
set(eulerdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/test/euler)
set(euler_HEADERS  eulermodel.hh  
                   problems.hh 
                   problemcreator.hh
                   chorjo.hh 
                   chorjo.cc 
                   passtraits.hh)

set(BASEDIR ../../main)

# helper variables
set(SOURCEBASE  ${BASEDIR}/main.cc chorjo.cc)
set(SOURCES  ${SOURCEBASE} ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc
                                        ${BASEDIR}/main_2.cc ${BASEDIR}/main_3.cc)
set(SOURCES12  ${SOURCEBASE} ${BASEDIR}/main_1.cc ${BASEDIR}/main_2.cc)
set(SOURCEONEP  ${SOURCEBASE} ${BASEDIR}/main_pol.cc)

set(GRIDTYPE YASPGRID)
set(GRIDDIM  2)
set(FLUX 1)
set(DIFFFLUXTYPE PRIMALDG)

set(TESTS euler
          euler12
          euleronep) 



add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-D${DIFFFLUXTYPE}" )
add_definitions( "-DWANT_CACHED_COMM_MANAGER=0" )


include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/euler")

add_executable( euler ${BASEDIR}/main.cc ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc ${BASEDIR}/main_2.cc ${BASEDIR}/main_3.cc )
dune_target_link_libraries(euler ${DUNE_LIBS})
set_property(TARGET euler APPEND PROPERTY 
  COMPILE_DEFINITIONS "LIMITER")

add_executable( euler12 ${BASEDIR}/main.cc ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc )
dune_target_link_libraries(euler12 ${DUNE_LIBS})
set_property(TARGET euler12 APPEND PROPERTY 
  COMPILE_DEFINITIONS "ONLY_P1_P2")

add_executable( euleronep ${BASEDIR}/main_pol.cc )
dune_target_link_libraries(euleronep ${DUNE_LIBS})
set_property(TARGET euleronep APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")

# We do not want want to build the tests during make all
# but just build them on demand
install(FILES ${euler_HEADERS} DESTINATION ${eulerdir})
foreach(i ${TESTS})
  add_dune_alugrid_flags(${i})
  add_test(${i} ${i})
endforeach(i ${TESTS})
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
