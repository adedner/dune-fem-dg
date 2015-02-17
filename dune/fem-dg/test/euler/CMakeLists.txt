dune_install( eulermodel.hh problems.hh problemcreator.hh chorjo.hh chorjo.cc passtraits.hh)

set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/euler/")

set(GRIDTYPE YASPGRID)
set(GRIDDIM  2)
set(POLORDER 2)
set(FLUX 1)
set(DIFFFLUXTYPE PRIMALDG)

configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-D${DIFFFLUXTYPE}" )
add_definitions( "-DWANT_CACHED_COMM_MANAGER=0" )

add_executable( eulerall "${MAININCLUDE_PATH}/main.cc" 
                         "${MAININCLUDE_PATH}/main_0.cc" 
                         "${MAININCLUDE_PATH}/main_1.cc"
                         "${MAININCLUDE_PATH}/main_2.cc" 
                         "${MAININCLUDE_PATH}/main_3.cc" chorjo.cc )
dune_target_link_libraries(eulerall ${DUNE_LIBS})
set_property(TARGET eulerall APPEND PROPERTY 
  COMPILE_DEFINITIONS "LIMITER")

add_executable( euler "${MAININCLUDE_PATH}/main.cc" "${MAININCLUDE_PATH}/main_pol.cc" chorjo.cc )
dune_target_link_libraries(euler ${DUNE_LIBS})
set_property(TARGET euler APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")

dune_add_test( euler eulerall )
