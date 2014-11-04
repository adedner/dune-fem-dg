# install these headers
set(nseqdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/test/navierstokes)
set(nseq_HEADERS  ns_model.hh 
                  problemcreator.hh 
                  nswaves.hh
                  ns_model_spec.hh 
                  nssmooth.hh 
                  thermodynamics.hh 
                  passtraits.hh)

set(BASEDIR ../../main)

# helper variables
set(SOURCEBASE  ${BASEDIR}/main.cc chorjo.cc)
set(SOURCES  ${SOURCEBASE} ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc
                                        ${BASEDIR}/main_2.cc ${BASEDIR}/main_3.cc)
set(SOURCES12  ${SOURCEBASE} ${BASEDIR}/main_1.cc ${BASEDIR}/main_2.cc)
set(SOURCEONEP  ${SOURCEBASE} ${BASEDIR}/main_pol.cc)

set(GRIDTYPE YASPGRID)
set(GRIDDIM  3)
set(POLORDER 2)
set(PROBLEM 2)
set(FLUX 1)
set(DIFFFLUX PRIMALDG)

set(TESTS nseq
          nseq12
          nseqonep) 



add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-D${DIFFFLUX}" )
add_definitions( "-DPROBLEM=${PROBLEM}" )


include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/navierstokes")

add_executable( nseq ${BASEDIR}/main.cc ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc ${BASEDIR}/main_2.cc ${BASEDIR}/main_3.cc )
dune_target_link_libraries(nseq ${DUNE_LIBS})
set_property(TARGET nseq APPEND PROPERTY 
  COMPILE_DEFINITIONS "LIMITER")

add_executable( nseq12 ${BASEDIR}/main.cc ${BASEDIR}/main_0.cc ${BASEDIR}/main_1.cc )
dune_target_link_libraries(nseq12 ${DUNE_LIBS})
set_property(TARGET nseq12 APPEND PROPERTY 
  COMPILE_DEFINITIONS "ONLY_P1_P2")

add_executable( nseqonep ${BASEDIR}/main_pol.cc )
dune_target_link_libraries(nseqonep ${DUNE_LIBS})
set_property(TARGET nseqonep APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")

# We do not want want to build the tests during make all
# but just build them on demand
install(FILES ${nseq_HEADERS} DESTINATION ${nseqdir})
foreach(i ${TESTS})
  add_dune_alugrid_flags(${i})
  add_test(${i} ${i})
endforeach(i ${TESTS})
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
