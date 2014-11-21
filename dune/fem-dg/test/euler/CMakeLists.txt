# install these headers
set(eulerdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/test/euler)
set(euler_HEADERS  eulermodel.hh  
                   problems.hh 
                   problemcreator.hh
                   chorjo.hh 
                   chorjo.cc 
                   passtraits.hh)

set(GRIDTYPE YASPGRID)
set(GRIDDIM  2)
set(POLORDER 2)
set(FLUX 1)
set(DIFFFLUXTYPE PRIMALDG)

set(TESTS euler
          eulerall)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-D${DIFFFLUXTYPE}" )
add_definitions( "-DWANT_CACHED_COMM_MANAGER=0" )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/euler")

add_executable( eulerall ${SOURCEALL} chorjo.cc )
dune_target_link_libraries(eulerall ${DUNE_LIBS})
set_property(TARGET eulerall APPEND PROPERTY 
  COMPILE_DEFINITIONS "LIMITER")

add_executable( euler ${SOURCEONE} chorjo.cc )
dune_target_link_libraries(euler ${DUNE_LIBS})
set_property(TARGET euler APPEND PROPERTY 
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
