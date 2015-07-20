dune_install( eulermodel.hh problems.hh problemcreator.hh chorjo.hh chorjo.cc passtraits.hh)

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/euler/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

set(POLORDER 2)
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DWANT_CACHED_COMM_MANAGER=0" )

add_executable( eulerall main.cc chorjo.cc )
dune_target_link_libraries(eulerall ${DUNE_LIBS})
set_property(TARGET eulerall APPEND PROPERTY 
             COMPILE_DEFINITIONS "LIMITER")

add_executable( euler main.cc chorjo.cc )
dune_target_link_libraries(euler ${DUNE_LIBS})
set_property(TARGET euler APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")

add_code_generate_targets( euler eulerall )

dune_add_test( euler eulerall )
