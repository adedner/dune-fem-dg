dune_install( models.hh  problemcreator.hh  problems.hh )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/euler/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )

add_executable( eulerall main.cc problems/chorjo.cc )
dune_target_enable_all_packages( eulerall )
set_property(TARGET eulerall APPEND PROPERTY 
             COMPILE_DEFINITIONS "LIMITER")

add_executable( euler main.cc problems/chorjo.cc )
dune_target_enable_all_packages( euler )
set_property(TARGET euler APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=2")

add_code_generate_targets( euler eulerall )

dune_fem_add_test( euler eulerall )
